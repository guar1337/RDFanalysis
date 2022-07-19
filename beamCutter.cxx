#include "/home/zalewski/aku/analysis/constants.h"
#include "TVector3.h"
#include "TLorentzVector.h"

std::vector<std::string> columnList{"MWPC_1_X",
									"MWPC_1_Y",
									"MWPC_2_X",
									"MWPC_2_Y",
									"lvBeam"};
std::vector<Double_t> parameters(11);
enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos, sLang1, sLang2, sLang3, sRang, sDistL, sDistR};
enum sCuts {mX1, X1, mY1, Y1, mX2, X2, mY2, Y2, mX3, X3, mY3, Y3};
std::vector<std::string> corrNames = {"-X1", "X1", "-Y1", "Y1", "-X2", "X2", "-Y2", "Y2", "-X3", "X3", "-Y3", "Y3"};
std::vector<Int_t> vCut;
Double_t tarPos[3], tarAngle[3];

Double_t ToFconstant = cs::tof_const;
Double_t tdcBinning = 0.125;
TRandom3 *rnd;
Double_t	m_MWPC_1_Z = -816.0;
Double_t	m_MWPC_2_Z = -270.0;
Bool_t verbosity=false;

TVector3 getBeamVector(std::vector<Double_t> rvecMWPC, Double_t m_kinE)
{
	TVector3 m_beamVector(rvecMWPC[3] - rvecMWPC[0], rvecMWPC[4] - rvecMWPC[1], rvecMWPC[5] - rvecMWPC[2]);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

void loadCutsCorrectionParameters(Int_t m_runNo=1)
{
	vCut.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/dataOut/tarVertexCuts.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("loadCutsCorrectionParameters:\tFailed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	Int_t tmpContainer;
    if (verbosity==true) printf("%d.\t", m_runNo);
	for (int iii = 0; iii < 12; iii++)
	{
		outStreamGenerated>>tmpContainer;
		vCut.push_back(tmpContainer);
		if (verbosity==true) printf("%s = %d\t", corrNames[iii].c_str(), vCut[iii]);
	}
	if (verbosity==true) std::cout<<std::endl;
	
}

std::vector<Double_t> getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y)
{
	std::vector<Double_t> rvecMWPC(6);
	rvecMWPC.at(0) = m_MWPC_1_X + parameters[sMWPC_1_X];
	rvecMWPC.at(1) = m_MWPC_1_Y + parameters[sMWPC_1_Y];
	rvecMWPC.at(2) = -816.0;
	rvecMWPC.at(3) = m_MWPC_2_X + parameters[sMWPC_2_X];
	rvecMWPC.at(4) = m_MWPC_2_Y + parameters[sMWPC_2_Y];
	rvecMWPC.at(5) = -270.0;

	return rvecMWPC;
}

TVector3 getTarVertex(std::vector<Double_t> rvecMWPC, Int_t m_geo)
{
	Double_t m_tarPos = tarPos[m_geo-1];
	Double_t m_tarAngle = tarAngle[m_geo-1];

	Double_t m_dX = rvecMWPC[3] - rvecMWPC[0];
	Double_t m_dY = rvecMWPC[4] - rvecMWPC[1];
	Double_t m_dZ = rvecMWPC[5] - rvecMWPC[2];
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, m_tarPos);
	TVector3 m_beamPoint(rvecMWPC[3], rvecMWPC[4], rvecMWPC[5]);
	TVector3 m_tarPerpendicular(sin(m_tarAngle), 0.0, cos(m_tarAngle));
	Double_t m_dCoeff = ((m_tarPoint-m_beamPoint).Dot(m_tarPerpendicular))/(m_vBeam.Dot(m_tarPerpendicular));
		
	Double_t m_evX = rvecMWPC[3] + m_dX * m_dCoeff;
	Double_t m_evY = rvecMWPC[4] + m_dY * m_dCoeff;
	Double_t m_evZ = rvecMWPC[5] + m_dZ * m_dCoeff;

	return TVector3(m_evX, m_evY, m_evZ);
}

Double_t getToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tdcF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tdcF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tdcF5-tdcF3)*tdcBinning+ToFconstant;
	return ToF;
}

bool filterToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return (ToF>160 && ToF<182);
}

Double_t getKinE(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	return cs::mass6He*(m_gamma-1.0);
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
				//printdcF("simplest solution %i\n", MWPCHitsArray[0]);
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
				//printdcF("simplest solution %i\n", MWPCHitsArray[0]);
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
					std::cout<<"Rogue MWPC event got through"<<std::endl;
				}
			}
		}
		else
		{
			std::cout<<"Rogue MWPC event got through"<<std::endl;
		}
		return 0;
	}

};

TLorentzVector getLVbeam(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t kinE)
{
	TVector3 m_vBeam;
	TLorentzVector m_lvBeam;
	Double_t m_dX = m_MWPC_2_X - m_MWPC_1_X;
	Double_t	m_dY = m_MWPC_2_Y - m_MWPC_1_Y;
	Double_t	m_dZ = m_MWPC_2_Z - m_MWPC_1_Z;

	m_vBeam.SetXYZ(m_dX, m_dY, m_dZ);
	Double_t ene_beam = cs::mass6He + kinE;
	Double_t mom_beam = sqrt(ene_beam*ene_beam - cs::mass6He*cs::mass6He);
						
	m_vBeam.SetMag(mom_beam);
	m_lvBeam.SetVectM(m_vBeam, cs::mass6He);

	return m_lvBeam;
}

int beamCutter()
{
	loadCutsCorrectionParameters();

	parameters[sMWPC_1_X] = -1.0;
	parameters[sMWPC_1_Y] = -2.1375;
	parameters[sMWPC_2_X] = 0.2; 
	parameters[sMWPC_2_Y] = -1.125;
	parameters[sTarPos] = 10.0;
	parameters[sLang1] = 0.0;
	parameters[sLang2] = 0.0;
	parameters[sLang3] = -1.0;
	parameters[sRang] = 0.0;
	parameters[sDistL] = -15.0;
	parameters[sDistR] = -25.0;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos] + 2.0;

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();
	rnd = new TRandom3();
	ROOT::EnableImplicitMT();
	TString inFname = TString::Format("/home/zalewski/dataTmp/simp/geo%d/simp_geo%d.root", cs::runNo, cs::runNo);

	ROOT::RDataFrame inDF("simplified", inFname.Data());
	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
						.Filter("trigger==1")
						.Filter("tdcF3[0]*tdcF3[1]*tdcF3[2]*tdcF3[3]*tdcF5[0]*tdcF5[1]*tdcF5[2]*tdcF5[3]")
						.Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
						.Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
						.Filter(filterToF,{"tdcF3", "tdcF5"})
						.Filter(filterMWPC(0),{"x1","nx1"})
						.Filter(filterMWPC(1),{"y1","ny1"})
						.Filter(filterMWPC(2),{"x2","nx2"})
						.Filter(filterMWPC(3),{"y2","ny2"});

		 
	auto newDF = outDF.Define("tof", getToF, {"tdcF3", "tdcF5"})
							.Define("kinE", getKinE, {"tof"})
							.Define("geo", "cs::runNo")
							.Define("aF5", "(F5[0]+F5[1]+F5[2]+F5[3])/4.0")
							.Define("MWPC_1_X", getMWPCpos(0),{"x1","nx1"})
							.Define("MWPC_1_Y", getMWPCpos(1),{"y1","ny1"})
							.Define("MWPC_2_X", getMWPCpos(2),{"x2","nx2"})
							.Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"})
							.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
							.Define("tarVertex", getTarVertex, {"MWPC", "geo"})
							.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
							.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
							.Define("tarCut1", "tarVertex.X()>vCut[mX1] && tarVertex.X()<vCut[X1] && tarVertex.Y()>vCut[mY1] && tarVertex.Y()<vCut[Y1]")
							.Define("tarCut2", "tarVertex.X()>vCut[mX2] && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]")
							.Define("tarCut3", "tarVertex.X()>vCut[mX3] && tarVertex.X()<vCut[X3] && tarVertex.Y()>vCut[mY3] && tarVertex.Y()<vCut[Y3]");
	auto filteredDF = newDF.Filter("tof>160 && tof<180");

	TString outFname = TString::Format("/home/zalewski/dataTmp/beamSource/beamIntegral%d_v1.root", cs::runNo);
	filteredDF.Snapshot("beamSource", outFname.Data());

	return 0;
}