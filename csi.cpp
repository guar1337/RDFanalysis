#include "csi.hh"

R__LOAD_LIBRARY(libgsl.so);
R__LOAD_LIBRARY(/home/zalewski/aku/ELC/build/libEloss.so);
R__LOAD_LIBRARY(/home/zalewski/aku/TELoss/libTELoss.so);

void loadParameters()
{
	TString fileName;
	std::string dummy;
	Double_t tmpA, tmpB, tmpC;

	for (auto &&iii : vecCalibratedColumns)
	{
		
		fileName = "/home/zalewski/dataTmp/calibrationParameters/geo123/" + iii + ".cal";
		if (cs::runNo==5)
		{
			fileName = "/home/zalewski/dataTmp/calibrationParameters/geo5/" + iii + ".cal";
		}
		
		
		//std::cout<<fileName<<std::endl;

		ifstream instream(fileName.Data());

		if (!instream) 
		{
			printf ("#Cannot open %s coefficient file\n",
			fileName.Data());
		}

		std::getline(instream, dummy);
		int chNo = std::stoi( dummy );
		Acoef.clear();
		Bcoef.clear();

		for (int jjj = 0; jjj < chNo; jjj++)
		{
			instream>>tmpA>>tmpB>>tmpC;
			Acoef.push_back(tmpA);
			Bcoef.push_back(tmpB);
			//std::cout<<tmpA<< "\t" <<tmpB<<std::endl;
		}

		vecAcoef.push_back(Acoef);
		vecBcoef.push_back(Bcoef);

		//cout<<endl;
		instream.close();
	}
}

Int_t getStripNumber(ROOT::VecOps::RVec<Double_t> &inputArray)
{
	return std::distance(inputArray.begin(),std::max_element(inputArray.begin(), inputArray.end()));
}

Double_t calculateToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return ToF;
}

TVector3 getBeamVector(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t m_kinE)
{
	Double_t m_fMWPC_1_X = m_MWPC_1_X + cs::MWPC1_X_displacement;
	Double_t m_fMWPC_1_Y = m_MWPC_1_Y + cs::MWPC1_Y_displacement;
	Double_t m_fMWPC_1_Z = -816.0;
	Double_t m_fMWPC_2_X = m_MWPC_2_X + cs::MWPC2_X_displacement;
	Double_t m_fMWPC_2_Y = m_MWPC_2_Y + cs::MWPC2_Y_displacement;
	Double_t m_fMWPC_2_Z = -270.0;

	TVector3 m_beamVector(m_fMWPC_2_X - m_fMWPC_1_X, m_fMWPC_2_Y - m_fMWPC_1_Y, m_fMWPC_2_Z - m_fMWPC_1_Z);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

Double_t getKineticEnergy(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	return cs::mass6He*(m_gamma-1.0);
}

bool filterToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return (ToF>100 && ToF<200);
}

struct filterMWPC
{
	int MWPC_ID;
	//create object containing ID of the column which is to be calibrated
	filterMWPC(int MWPC_number) 
	{
		MWPC_ID = MWPC_number;
	};

	bool operator()(ROOT::VecOps::RVec<unsigned short> &MWPCHitsArray, unsigned short MWPCHitsCount)
	{		
		if (MWPCHitsCount!=0)
		{
			int sizeof_clust = 1;
			//one wire got signal - simplest solution
			if (MWPCHitsCount==1)
			{
				//printf("simplest solution %i\n", MWPCHitsArray[0]);
				if (MWPC_ID == 0 || MWPC_ID == 2)
				{
					//*MWPC_pos = (15.5 - MWPCHitsArray[0])*1.25;
				}

				else
				{
					//*MWPC_pos = (-15.5 + MWPCHitsArray[0])*1.25;
				}
				return 1;
			}

			//more than one but is it one cluster?
			else
			{	//checking...
				for (int iii = 1; iii < MWPCHitsCount; iii++)
				{
					if ((MWPCHitsArray[iii] - MWPCHitsArray[iii-1])==1)
					{
						sizeof_clust++;
					}
				}
				//
				if (sizeof_clust==MWPCHitsCount)
				{

					if (MWPC_ID == 0 || MWPC_ID == 2)
					{
						//*MWPC_pos = (15.5 - (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25;
					}

					else
					{
						//*MWPC_pos = (-15.5 + (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25;
					}
					return 1;

				}
				else if (sizeof_clust<MWPCHitsCount)
				{
					return 0;
				}
			}
		}
		else
		{
			return 0;
		}
		return 0;
	}

};

struct getMWPCpos
{
	int MWPC_ID;
	//create object containing ID of the column which is to be calibrated
	getMWPCpos(int MWPC_number) 
	{
		MWPC_ID = MWPC_number;
	};

	Double_t operator()(ROOT::VecOps::RVec<unsigned short> &MWPCHitsArray, unsigned short MWPCHitsCount)
	{		
		if (MWPCHitsCount!=0)
		{
			int sizeof_clust = 1;
			//one wire got signal - simplest solution
			if (MWPCHitsCount==1)
			{
				//printf("simplest solution %i\n", MWPCHitsArray[0]);
				if (MWPC_ID == 0 || MWPC_ID == 2)
				{
					return (15.5 - MWPCHitsArray[0])*1.25 + rnd->Uniform(0.0,1.25)-0.625;
				}

				else
				{
					return (-15.5 + MWPCHitsArray[0])*1.25 + rnd->Uniform(0.0,1.25)-0.625;
				}
			}

			//more than one but is it one cluster?
			else
			{	//checking...
				for (int iii = 1; iii < MWPCHitsCount; iii++)
				{
					if ((MWPCHitsArray[iii] - MWPCHitsArray[iii-1])==1)
					{
						sizeof_clust++;
					}
				}
				//
				if (sizeof_clust==MWPCHitsCount)
				{

					if (MWPC_ID == 0 || MWPC_ID == 2)
					{
						return (15.5 - (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25 + rnd->Uniform(0.0,1.25)-0.625;
					}

					else
					{
						return (-15.5 + (MWPCHitsArray[0]+MWPCHitsArray[MWPCHitsCount-1])/2.0)*1.25 + rnd->Uniform(0.0,1.25)-0.625;
					}
				}
				else if (sizeof_clust<MWPCHitsCount)
				{
					std::cout<<"Rogue MWPC event got through. No. of MWPC Hits: "<<MWPCHitsCount<<std::endl;
				}
			}
		}
		else
		{
			std::cout<<"There is 0 MWPC hits"<<std::endl;
		}
		return 0;
	}

};

struct Calibrator
{
	int columnID;
	//create object containing ID of the column which is to be calibrated
	Calibrator(int columnToCalibrate) 
	{
		columnID = columnToCalibrate;
	};

	ROOT::VecOps::RVec<Double_t> operator()(ROOT::VecOps::RVec<unsigned short> &inputArray)
	{		
		auto mySize = inputArray.size();
		ROOT::VecOps::RVec<Double_t> outputArray(mySize);
		for (size_t iii = 0; iii < mySize; iii++)
		{
			//std::cout<<vecAcoef.at(columnID).at(iii)<<" "<<"My ID is: "<<vecCalibratedColumns.at(columnID)<<std::endl;
			outputArray.at(iii) = vecAcoef.at(columnID).at(iii) + vecBcoef.at(columnID).at(iii) * inputArray[iii];
		}

		return outputArray;
	}
};

ROOT::RDF::RNode ApplyDefines(	ROOT::RDF::RNode df, 
								const std::vector<std::string> &colNames,
								unsigned int iii = 0)
{
	//recursive function updating RDataFrame over scope of std::vector with column names for calibration
	if (iii == colNames.size())
	{
		return df.Define("tof", calculateToF, {"tdcF3", "tdcF5"})
				 .Define("aF5", "(F5[0]+F5[1]+F5[2]+F5[3])/4.0")
				 .Define("MWPC_1_X", getMWPCpos(0),{"x1","nx1"})
				 .Define("MWPC_1_Y", getMWPCpos(1),{"y1","ny1"})
				 .Define("MWPC_2_X", getMWPCpos(2),{"x2","nx2"})
				 .Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"})
				 .Define("SQX_L_strip", getStripNumber, {"cal_SQX_L"})
				 .Define("SQX_R_strip", getStripNumber, {"cal_SQX_R"})
				 .Define("CsI_L_strip", getStripNumber, {"cal_CsI_L"})
				 .Define("CsI_R_strip", getStripNumber, {"cal_CsI_R"})
				 .Define("sqlde", [](ROOT::VecOps::RVec<Double_t> &SQX_L, Int_t stripNumber){return SQX_L[stripNumber];}, {"cal_SQX_L", "SQX_L_strip"})
				 .Define("sqrde", [](ROOT::VecOps::RVec<Double_t> &SQX_R, Int_t stripNumber){return SQX_R[stripNumber];}, {"cal_SQX_R", "SQX_R_strip"})
				 .Define("sqletot", [](ROOT::VecOps::RVec<Double_t> &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"cal_CsI_L", "CsI_L_strip"})
				 .Define("sqretot", [](ROOT::VecOps::RVec<Double_t> &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"cal_CsI_R", "CsI_R_strip"})
				 .Define("sqletotr", [](ROOT::VecOps::RVec<UShort_t> &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"CsI_L", "CsI_L_strip"})
				 .Define("sqretotr", [](ROOT::VecOps::RVec<UShort_t> &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"CsI_R", "CsI_R_strip"});
	}

	std::string inputColumn = colNames[iii];
	inputColumn.insert(0, "cal_");
	//std::cout<<colNames[iii]<<"\t"<<inputColumn<<std::endl;
	//calling 'Calibrator' so I can pass other info
	return ApplyDefines(df.Define(inputColumn, Calibrator(iii), {colNames[iii]}), colNames, iii + 1);
}

void calibratorCaller(TString inFileName)
{	
	inFileName.ReplaceAll("raw","cln");
	std::cout<<"Calibrating "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("cleaned", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cln","cal");
	// Print columns' names
	auto dfWithDefines = ApplyDefines(inDF, vecCalibratedColumns);
	auto c = dfWithDefines.Count();
	Int_t customClusterSize = *c / numberOfThreads;
	std::cout << customClusterSize << std::endl;
	ROOT::RDF::RSnapshotOptions myOpts;
	myOpts.fAutoFlush=customClusterSize;
	dfWithDefines.Snapshot("calibrated", outFilename.Data(), dfWithDefines.GetColumnNames(), myOpts);
}


void translator(TString inFileName)
{
	std::cout<<"Translating "<<inFileName<<std::endl;
	TFile inputFile(inFileName.Data(), "READ");
	TTree *inputTree;
	inputTree = (TTree*)inputFile.Get("AnalysisxTree");
	TString outFilename = inFileName.ReplaceAll("raw","simp");

	inputTree->Process(selector, outFilename.Data());
	
	inputFile.Close();	
}

void cleaner(TString inFileName)
{
	inFileName.ReplaceAll("raw","simp");
	std::cout<<"Cleaning "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("simplified", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("simp","cln");

	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
						  .Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
						  .Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
						  .Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
						  .Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
						  .Filter(filterMWPC(0),{"x1","nx1"})
						  .Filter(filterMWPC(1),{"y1","ny1"})
						  .Filter(filterMWPC(2),{"x2","nx2"})
						  .Filter(filterMWPC(3),{"y2","ny2"})
						  .Filter(filterToF,{"tdcF3", "tdcF5"});

	outDF.Snapshot("cleaned", outFilename.Data());
}

void analysis(TString inFileName)
{
	inFileName.ReplaceAll("raw","cal");
	ROOT::RDataFrame inDF("calibrated", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cal","dE");
	std::cout<<"Analysing "<<inFileName<<std::endl;
	inDF.Define("CsI_L_strip", getStripNumber, {"cal_CsI_L"})
		.Define("CsI_R_strip", getStripNumber, {"cal_CsI_R"})
		.Define("sqletot", [](ROOT::VecOps::RVec<Double_t> &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"CsI_L", "CsI_L_strip"})
		.Define("sqretot", [](ROOT::VecOps::RVec<Double_t> &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"CsI_R", "CsI_R_strip"})
		.Define("kinE", [](Double_t ToF){return getKineticEnergy(ToF);},{"tof"});



}

void csi()
{
	ROOT::EnableImplicitMT();
	numberOfThreads = 20;
	TStopwatch *stopwatch = new TStopwatch();
	
	rnd = new TRandom3();
	selector = TSelector::GetSelector("/home/zalewski/aku/analysis/simplifier123.C");
	std::cout<<"..."<<std::endl<<".."<<std::endl<<"."<<std::endl;

	loadParameters();
	TString str_name;
	TSystemDirectory *dir_data = new TSystemDirectory("data",raw_data_path.Data());
	TIter bluster = dir_data->GetListOfFiles();
	while (TObject *obj = bluster())
	{
		str_name = obj->GetName();
		//std::cout<<str_name<<std::endl;
		if (str_name.Contains("raw_csi") &&  //if we want to ommit some phrase (np.li9_3_21_0)
			!str_name.Contains("li9") //  "!"  sign is important
								)
		{
			TString inputFilePath = raw_data_path + str_name;
			//printf("fName:\t%s\n",inputFilePath.Data());
			//translator(inputFilePath);
			//cleaner(inputFilePath);
			calibratorCaller(inputFilePath);
			//analysis(inputFilePath);

		}
	}
	stopwatch->Print();
}