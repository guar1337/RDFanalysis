#include "/home/zalewski/aku/analysis/constants.h"
#include "TVector3.h"
#include "TLorentzVector.h"

std::vector<std::string> columnList{"MWPC_1_X",
												"MWPC_1_Y",
												"MWPC_2_X",
												"MWPC_2_Y",
												"lvBeam"};

Double_t ToFconstant = cs::tof_const;
Double_t tdcBinning = 0.125;
TRandom3 *rnd;
Double_t	m_MWPC_1_Z = -816.0;
Double_t	m_MWPC_2_Z = -270.0;

Double_t getToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tdcF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tdcF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tdcF5-tdcF3)*tdcBinning+ToFconstant;
	return ToF;
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
	rnd = new TRandom3();
	ROOT::EnableImplicitMT();
	Int_t numberOfThreads = 20;
	ROOT::RDataFrame inDF("simplified", "/home/zalewski/dataTmp/simp/geo3/simp_geo3.root");
	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
						.Filter("trigger==1")
						.Filter("tdcF3[0]*tdcF3[1]*tdcF3[2]*tdcF3[3]*tdcF5[0]*tdcF5[1]*tdcF5[2]*tdcF5[3]")
						.Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
						.Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
						.Filter(filterMWPC(0),{"x1","nx1"})
						.Filter(filterMWPC(1),{"y1","ny1"})
						.Filter(filterMWPC(2),{"x2","nx2"})
						.Filter(filterMWPC(3),{"y2","ny2"});

		 
	auto newDF = outDF.Define("tof", getToF, {"tdcF3", "tdcF5"})
							.Define("kinE", getKinE, {"tof"})
							.Define("aF5", "(F5[0]+F5[1]+F5[2]+F5[3])/4.0")
							.Define("MWPC_1_X", getMWPCpos(0),{"x1","nx1"})
							.Define("MWPC_1_Y", getMWPCpos(1),{"y1","ny1"})
							.Define("MWPC_2_X", getMWPCpos(2),{"x2","nx2"})
							.Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"})
							.Define("geo", "cs::runNo")
							.Define("lvBeam", getLVbeam, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y", "kinE"});
	auto filteredDF = newDF.Filter("tof>165 && tof<182 && aF5>600 && aF5<2800");

	filteredDF.Snapshot("beamSource", "/home/zalewski/dataTmp/beamSource/beamSource3.root", columnList);

	return 0;
}