#include "/home/zalewski/aku/analysis/ui.hh"

R__LOAD_LIBRARY(libgsl.so);
R__LOAD_LIBRARY(/home/zalewski/aku/ELC/build/libEloss.so);
R__LOAD_LIBRARY(/home/zalewski/aku/TELoss/libTELoss.so);

//d,t reaction
Double_t getDTLang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass5He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	m_lv3H.SetTheta(m_thetaCM);
	m_lv3H.Boost(m_boostVect);

	return m_lv3H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getDTRang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm5HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm5Heene = m_cm5HekinE + cs::mass5He;
	Double_t m_cm5Hemom = sqrt(m_cm5Heene*m_cm5Heene - cs::mass5He*cs::mass5He);

	TLorentzVector m_lv5He(0.0,0.0,m_cm5Hemom, m_cm5Heene);

	m_lv5He.SetTheta(TMath::Pi() - m_thetaCM);
	m_lv5He.Boost(m_boostVect);
	return m_lv5He.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

//p,t reaction
Double_t getPTLang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + m_Q;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);
	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	m_lv3H.SetTheta(m_thetaCM);
	m_lv3H.Boost(m_boostVect);

	return m_lv3H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getPTRang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + m_Q;

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	m_lv4He.SetTheta(TMath::Pi() - m_thetaCM);
	m_lv4He.Boost(m_boostVect);
	return m_lv4He.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t altPTLang(Double_t m_thetaCM, TLorentzVector m_lvBeam, Double_t m_Q)
{
	Int_t multi = (rnd->Integer(2)==1) ? 1 : -1;
	TLorentzVector m_lv6He(m_lvBeam);
	// /m_lv6He.SetTheta(0.0);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + m_Q;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	Double_t m_betaCM = m_lvCM.Beta();
	Double_t m_gammaCM = m_lvCM.Gamma();
	Double_t m_beta3HCM = m_lv3H.Beta();
	Double_t m_beta4HeCM = m_lv4He.Beta();
	m_lv3H.SetTheta(m_thetaCM);
	m_lv3H.Boost(m_boostVect);

	Double_t m_betaCM3H = (sqrt(m_cm3HkinE*m_cm3HkinE+2*m_cm3HkinE*cs::mass3H))/(m_cm3HkinE+cs::mass3H);
	Double_t m_betaCM4He = (sqrt(m_cm4HekinE*m_cm4HekinE+2*m_cm4HekinE*cs::mass4He))/(m_cm4HekinE+cs::mass4He);
	Double_t m_beta3HRatio = m_betaCM/m_betaCM3H;
	Double_t m_beta4HeRatio = m_betaCM/m_betaCM4He;
	Double_t m_beta3HRatio2 = m_beta3HRatio*m_beta3HRatio;

	//calculating CM angles
	Double_t m_gammaTan2 = std::pow(m_gammaCM*tan(m_lv3H.Angle(m_lvBeam.Vect())),2);
	Double_t m_cosThetaCM3H = (-m_beta3HRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta3HRatio2)))/(1.0+m_gammaTan2);
	Double_t m_ThetaCM3H = acos(m_cosThetaCM3H);
	return m_ThetaCM3H;
}

//p,p reaction
Double_t getPPLang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0,0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv1H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv1H.SetTheta(m_thetaCM);

	//m_lv6He.Boost(m_boostVect);
	m_lv1H.Boost(m_boostVect);
	Double_t m_sqlangpp = 180.0 * (m_lvBeam.Angle(m_lv1H.Vect()))/double(TMath::Pi());
	//Double_t m_sqrangpp = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqlangpp;
}

Double_t getPPRang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0,0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	//m_lv1H.Boost(-m_boostVect);

	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	//m_lv1H.SetTheta(m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	//m_lv1H.Boost(m_boostVect);
	//Double_t m_sqlangpp = 180.0 * (m_lvBeam.Angle(m_lv1H.Vect()))/double(TMath::Pi());
	Double_t m_sqrangpp = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqrangpp;
}

//d,d reaction
Double_t getDDRang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	//m_lv2H.Boost(-m_boostVect);

	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	//m_lv2H.SetTheta(m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	//m_lv2H.Boost(m_boostVect);
	//Double_t m_sqlangdd = 180.0 * (m_lvBeam.Angle(m_lv2H.Vect()))/double(TMath::Pi());
	Double_t m_sqrangdd = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqrangdd;
}

Double_t getDDLang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv2H.SetTheta(m_thetaCM);

	//m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);
	Double_t m_sqlangdd = 180.0 * (m_lvBeam.Angle(m_lv2H.Vect()))/double(TMath::Pi());
	//Double_t m_sqrangdd = 180.0 * (m_lvBeam.Angle(m_lv6He.Vect()))/double(TMath::Pi());
	return m_sqlangdd;
}

Int_t getStripNumber(ROOT::VecOps::RVec<Double_t> &inputArray)
{
	return std::distance(inputArray.begin(),std::max_element(inputArray.begin(), inputArray.end()));
}

TVector3 getLeftDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	m_xStrip+=rnd->Uniform(0.0,1.0)-0.5;
	m_yStrip+=rnd->Uniform(0.0,1.0)-0.5;

	// coordinates of hit in LAB system	
	Double_t X2hDet = -cs::widthStripX * m_xStrip * cos(sqlAng);
	Double_t Y2hDet = cs::widthStripY * m_yStrip;
	Double_t Z2hDet = cs::widthStripX * m_xStrip * sin(sqlAng);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	m_xStrip+=rnd->Uniform(0.0,1.0)-0.5;
	m_yStrip+=rnd->Uniform(0.0,1.0)-0.5;

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * m_xStrip * cos(sqrAng);
	Double_t Y6HeDet = cs::widthStripY * m_yStrip;
	Double_t Z6HeDet = cs::widthStripX * m_xStrip * sin(sqrAng);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
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

TVector3 getTarVertex(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y)
{
	Double_t m_fMWPC_1_X = m_MWPC_1_X + cs::MWPC1_X_displacement;
	Double_t m_fMWPC_1_Y = m_MWPC_1_Y + cs::MWPC1_Y_displacement;
	Double_t m_fMWPC_1_Z = -816.0;
	Double_t m_fMWPC_2_X = m_MWPC_2_X + cs::MWPC2_X_displacement;
	Double_t m_fMWPC_2_Y = m_MWPC_2_Y + cs::MWPC2_Y_displacement;
	Double_t m_fMWPC_2_Z = -270.0;
	Double_t m_dX = m_fMWPC_2_X - m_fMWPC_1_X;
	Double_t m_dY = m_fMWPC_2_Y - m_fMWPC_1_Y;
	Double_t m_dZ = m_fMWPC_2_Z - m_fMWPC_1_Z;
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, tarPos);
	TVector3 m_beamPoint(m_fMWPC_2_X, m_fMWPC_2_Y, m_fMWPC_2_Z);
	TVector3 m_tarPerpendicular(sin(tarAngle), 0.0, cos(tarAngle));
	Double_t m_dCoeff = ((m_tarPoint-m_beamPoint).Dot(m_tarPerpendicular))/(m_vBeam.Dot(m_tarPerpendicular));
		
	Double_t m_evX = m_fMWPC_2_X + m_dX * m_dCoeff;
	Double_t m_evY = m_fMWPC_2_Y + m_dY * m_dCoeff;
	Double_t m_evZ = m_fMWPC_2_Z + m_dZ * m_dCoeff;

	return TVector3(m_evX, m_evY, m_evZ);
}

Double_t getKineticEnergy(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	return cs::mass6He*(m_gamma-1.0);
}

bool filterLeftDetector(ROOT::VecOps::RVec<Double_t> &leftDetectorArray)
{
	int aboveThresholdEventCounter(0);
	auto leftDetectorArraySize = leftDetectorArray.size();
	for (size_t iii = 0; iii < leftDetectorArraySize; iii++)
	{
		if (leftDetectorArray[iii]>1.2)
		{
			aboveThresholdEventCounter++;
		}
	}
	if (aboveThresholdEventCounter==1)
	{
		return true;
	}

	else
	{
		return false;
	}
	
}

bool filterRightDetector(ROOT::VecOps::RVec<Double_t> &rightDetectorArray)
{
	int aboveThresholdEventCounter(0);
	auto rightDetectorArraySize = rightDetectorArray.size();
	for (size_t iii = 0; iii < rightDetectorArraySize; iii++)
	{
		if (rightDetectorArray[iii]>5.0)
		{
			aboveThresholdEventCounter++;
		}
	}
	if (aboveThresholdEventCounter==1)
	{
		return true;
	}

	else
	{
		return false;
	}
	
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

bool filterToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return (ToF>160 && ToF<182);
}

Double_t calculateToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return ToF;
}

bool filterSQ(ROOT::VecOps::RVec<unsigned short> &SQ, Double_t threshold)
{
	Bool_t sqFlag=false;
	Double_t calibratedSQ;

	for (int iii = 0; iii < 32; iii++)
	{
		//SQX_L filtering
		if (threshold==1.2)
		{
			//vecAcoef.at(0) is array with calibration parameters of SQX_L
			calibratedSQ = vecAcoef.at(0).at(iii) + vecBcoef.at(0).at(iii) * SQ[iii];
		}

		//SQX_R filtering
		else
		{
			//vecAcoef.at(3) is array with calibration parameters of SQX_R
			calibratedSQ = vecAcoef.at(3).at(iii) + vecBcoef.at(3).at(iii) * SQ[iii];
		}
		
		
		if (calibratedSQ>threshold)
		{
			sqFlag=true;
		}
	}
	return sqFlag;
}

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

void cleaner(TString inFileName)
{
	Double_t MWPCposContainer;
	if (cs::runNo==5)
	{
		ToFconstant = cs::tof_const_5;
		tdcBinning = 0.0625;
	}
	
	Double_t sqlThreshold = 0.1;
	Double_t sqrThreshold = 2.0;

	inFileName.ReplaceAll("raw","simp");
	std::cout<<"Cleaning "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("simplified", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("simp","cln");
	//get rid of broken MWPC events
	auto outDF = inDF.Filter("nx1<100 && nx2<100 && ny1<100 && ny2<100")
						.Filter("(tdcF3[0]-tdcF3[1]) > -50.0 && (tdcF3[0]-tdcF3[1]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[2]) > -50.0 && (tdcF3[0]-tdcF3[2]) < 50.0")
						.Filter("(tdcF3[0]-tdcF3[3]) > -50.0 && (tdcF3[0]-tdcF3[3]) < 50.0")
						.Filter("(tdcF5[0]-tdcF5[1]) > -50.0 && (tdcF5[0]-tdcF5[1]) < 50.0")
						.Filter(filterMWPC(0),{"x1","nx1"})
						.Filter(filterMWPC(1),{"y1","ny1"})
						.Filter(filterMWPC(2),{"x2","nx2"})
						.Filter(filterMWPC(3),{"y2","ny2"})
						.Filter(filterToF,{"tdcF3", "tdcF5"})
						.Filter([sqlThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_L){return filterSQ(SQX_L, sqlThreshold);},{"SQX_L"})
						.Filter([sqrThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_R){return filterSQ(SQX_R, sqrThreshold);},{"SQX_R"});

	outDF.Snapshot("cleaned", outFilename.Data());
}

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
				 .Define("MWPC_2_Y", getMWPCpos(3),{"y2","ny2"});
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

void analysis(TString inFileName)
{
	inFileName.ReplaceAll("raw","cal");
	ROOT::RDataFrame inDF("calibrated", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cal","dE");
	std::cout<<"Analysing "<<inFileName<<std::endl;

	//load Graphical cuts
	TFile cutgFile("/home/zalewski/aku/analysis/gcuts.root","READ");
	GCutHe4 = (TCutG*)cutgFile.Get("dehe4");
	GCutHe6 = (TCutG*)cutgFile.Get("dehe6");
	GCut2H = (TCutG*)cutgFile.Get("dd");
	GCut3H = (TCutG*)cutgFile.Get("tt");

	AELC *h1_Si = new ELC(1, 1, si_Nel, 2.35, si_A, si_Z, si_W, 500.,1500);
	AELC *h2_Si = new ELC(2, 1, si_Nel, 2.35, si_A, si_Z, si_W, 500.,1500);
	AELC *h3_Si = new ELC(3, 1, si_Nel, 2.35, si_A, si_Z, si_W, 500.,1500);

	AELC *he4_Si = new ELC(4, 2, si_Nel, 2.35, si_A, si_Z, si_W, 500.,1500);
	AELC *he6_Si = new ELC(6, 2, si_Nel, 2.35, si_A, si_Z, si_W, 500.,1500);

	siEloss1H.SetEL(1, 2.330); // density in g/cm3
	siEloss1H.AddEL(14., 28.086, 1);  //Z, mass
	siEloss1H.SetZP(1., 1.);		//Z, A
	siEloss1H.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	siEloss1H.SetDeltaEtab(100);

	siEloss2H.SetEL(1, 2.330); // density in g/cm3
	siEloss2H.AddEL(14., 28.086, 1);  //Z, mass
	siEloss2H.SetZP(1., 2.);		//Z, A
	siEloss2H.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	siEloss2H.SetDeltaEtab(100);

	siEloss3H.SetEL(1, 2.330); // density in g/cm3
	siEloss3H.AddEL(14., 28.086, 1);  //Z, mass
	siEloss3H.SetZP(1., 3.);		//Z, A
	siEloss3H.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	siEloss3H.SetDeltaEtab(100);

	siEloss4He.SetEL(1, 2.330); // density in g/cm3
	siEloss4He.AddEL(14., 28.086, 1);  //Z, mass
	siEloss4He.SetZP(2., 4.);		//Z, A
	siEloss4He.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	siEloss4He.SetDeltaEtab(100);

	siEloss6He.SetEL(1, 2.330); // density in g/cm3
	siEloss6He.AddEL(14., 28.086, 1);  //Z, mass
	siEloss6He.SetZP(2., 6.);		//Z, A
	siEloss6He.SetEtab(100000, 200.);	// ?, MeV calculate ranges
	siEloss6He.SetDeltaEtab(100);

	Double_t qDT_45 = cs::Qdt;
	Double_t qDT_0 = 0.0;

	Double_t qPT_75 = cs::Qpt;
	Double_t qPT_0 = 0.0;

	//load geometrical parameters
	switch (cs::runNo)
	{
		case 0:
		{
			printf("Wrong geometry\n");
			sqlAng = (0.0) * TMath::DegToRad();
			sqrAng = (0.0) * TMath::DegToRad();
			tarPos = 0.0;
			break;
		}

		case 1:
		{
			tarPos = 10.0;
			tarThickness = 80.0 + cs::tarThicknessShift;
			tarAngle = 45.0 * TMath::DegToRad();

			sqlAng = (65.0 + 0.0) * TMath::DegToRad();
			sqlDist = 170.0;
			
			sqrAng = (15.0 + -0.192) * TMath::DegToRad();
			sqrDist = 250.0;
			break;
		}

		case 2:
		{
			tarPos = 10.0 + 2.0;
			tarThickness = 160 + 2*cs::tarThicknessShift;
			tarAngle = 6.0 * TMath::DegToRad();

			sqlAng = (50.0 + 0.536) * TMath::DegToRad();
			sqlDist = 170.0;

			sqrAng = (15.0 + 0.275) * TMath::DegToRad();
			sqrDist = 250.0;
			break;
		}

		case 3:
		{
			tarPos = 10.0 + 2.0;
			tarThickness = 160 + 2*cs::tarThicknessShift;
			tarAngle = 0.0 * TMath::DegToRad();

			sqlAng = (35.0 + 0.536) * TMath::DegToRad();
			sqlDist = 170.0;

			sqrAng = (15.0 + 0.275) * TMath::DegToRad();
			sqrDist = 250.0;
			break;
		}
	}

	Double_t X2Hlab = sqlDist*sin(sqlAng) + (cs::sqlXzero + cs::leftDetShift) * cos(sqlAng);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(sqlAng) - (cs::sqlXzero + cs::leftDetShift) * sin(sqlAng);
	TVector3 leftDetPosition(X2Hlab, Y2Hlab, Z2Hlab);

	Double_t X6Helab = sqrDist*sin(-sqrAng) - (cs::sqrXzero + cs::rightDetShift) * cos(sqrAng);
	Double_t Y6Helab = cs::sqrYstart + cs::widthStripY;
	Double_t Z6Helab = sqrDist*cos(sqrAng) - (cs::sqrXzero + cs::rightDetShift) * sin(sqrAng);
	TVector3 rightDetPosition(X6Helab, Y6Helab, Z6Helab);

	inDF.Filter(filterLeftDetector, {"cal_SQX_L"})
		.Filter(filterLeftDetector, {"cal_SQY_L"})
		.Filter(filterRightDetector, {"cal_SQX_R"})
		.Filter(filterRightDetector, {"cal_SQY_R"});
	auto outDF = inDF.Define("kinE", [](Double_t ToF){return getKineticEnergy(ToF);},{"tof"})
					 .Define("vBeam", getBeamVector, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y", "kinE"})
					 .Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
					 .Define("tarVertex", getTarVertex, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
					 .Define("geo", "cs::runNo")
					 .Define("SQ300_strip", getStripNumber, {"cal_SQ300"})
					 .Define("SQX_L_strip", getStripNumber, {"cal_SQX_L"})
					 .Define("SQY_L_strip", getStripNumber, {"cal_SQY_L"})
					 .Define("SQX_R_strip", getStripNumber, {"cal_SQX_R"})
					 .Define("SQY_R_strip", getStripNumber, {"cal_SQY_R"})
					 .Define("CsI_L_strip", getStripNumber, {"cal_CsI_L"})
					 .Define("CsI_R_strip", getStripNumber, {"cal_CsI_R"})
					 .Define("sq300", [](ROOT::VecOps::RVec<Double_t> &SQ300, Int_t stripNumber){return SQ300[stripNumber];}, {"cal_SQ300", "SQ300_strip"})
					 .Define("sqlde", [](ROOT::VecOps::RVec<Double_t> &SQX_L, Int_t stripNumber){return SQX_L[stripNumber];}, {"cal_SQX_L", "SQX_L_strip"})
					 .Define("sqrde", [](ROOT::VecOps::RVec<Double_t> &SQX_R, Int_t stripNumber){return SQX_R[stripNumber];}, {"cal_SQX_R", "SQX_R_strip"})
					 .Define("sqletot", [](ROOT::VecOps::RVec<Double_t> &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"cal_CsI_L", "CsI_L_strip"})
					.Define("sqretot", [](ROOT::VecOps::RVec<Double_t> &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"cal_CsI_R", "CsI_R_strip"})

					 .Define("sqletotp", [](Double_t sqlde, Double_t sq300){return siEloss1H.GetE0dE(sqlde+sq300, 1000.0);}, {"sqlde", "sq300"})
					 //.Define("sqletotd", [](Double_t sqlde, Double_t sq300){return siEloss3H.GetE0dE(sqlde+sq300, 1000.0);}, {"sqlde", "sq300"})
					 .Define("sqletott", [](Double_t sqlde, Double_t sq300){return siEloss3H.GetE0dE(sqlde+sq300, 1000.0);}, {"sqlde", "sq300"})
					 .Define("sqretot4", [](Double_t sqrde){return siEloss4He.GetE0dE(sqrde, 1000.0);}, {"sqrde"})
					 //.Define("sqretot6", [](Double_t sqrde){return siEloss6He->GetE(sqrde, 1000.0);}, {"sqrde"})

					 .Define("sqltime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_L, Int_t stripNumber){return tSQX_L[stripNumber];}, {"tSQX_L", "SQX_L_strip"})
					 .Define("sqrtime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_R, Int_t stripNumber){return tSQX_R[stripNumber];}, {"tSQX_R", "SQX_R_strip"})

					 .Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip"})
					 .Define("leftLabVertex", [leftDetPosition](TVector3 leftDetVertex){return leftDetVertex+leftDetPosition;}, {"leftDetVertex"})
					 .Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
					 .Define("rightLabVertex", [rightDetPosition](TVector3 rightDetVertex){return rightDetVertex+rightDetPosition;}, {"rightDetVertex"})
					 .Define("v2H", [](TVector3 leftLabVertex, TVector3 tarVertex){return TVector3(leftLabVertex-tarVertex);}, {"leftLabVertex", "tarVertex"})
					 .Define("v6He", [](TVector3 rightLabVertex, TVector3 tarVertex){return TVector3(rightLabVertex-tarVertex);}, {"rightLabVertex", "tarVertex"})
					 .Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
					 .Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
					 .Define("he4", [](Double_t sqretot, Double_t sqrde){return GCutHe4->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"})
					 .Define("he6", [](Double_t sqretot, Double_t sqrde){return GCutHe6->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"})
					 .Define("dd", [](Double_t sqletot, Double_t sqlde){return GCut2H->IsInside(sqletot,sqlde);}, {"sqletot","sqlde"})
					 .Define("tt", [](Double_t sqletot, Double_t sqlde){return GCut3H->IsInside(sqletot,sqlde);}, {"sqletot","sqlde"})
					 .Define("thetaCM", "rnd->Uniform(0.0,TMath::Pi())")
					 .Define("sqlangpp", getPPLang, {"thetaCM","lvBeam"})
					 .Define("sqrangpp", getPPRang, {"thetaCM","lvBeam"})
					 .Define("sqlangdd", getDDLang, {"thetaCM","lvBeam"})
					 .Define("sqrangdd", getDDRang, {"thetaCM","lvBeam"})

					 .Define("sqlangdt_45", [qDT_45](Double_t thetaCM, TLorentzVector lvBeam){return getDTLang(thetaCM, lvBeam, qDT_45);}, {"thetaCM","lvBeam"})
					 .Define("sqrangdt_45", [qDT_45](Double_t thetaCM, TLorentzVector lvBeam){return getDTRang(thetaCM, lvBeam, qDT_45);}, {"thetaCM","lvBeam"})
					 .Define("sqlangdt_0", [qDT_0](Double_t thetaCM, TLorentzVector lvBeam){return getDTLang(thetaCM, lvBeam, qDT_0);}, {"thetaCM","lvBeam"})
					 .Define("sqrangdt_0", [qDT_0](Double_t thetaCM, TLorentzVector lvBeam){return getDTRang(thetaCM, lvBeam, qDT_0);}, {"thetaCM","lvBeam"})

					 .Define("sqlangpt_75", [qPT_75](Double_t thetaCM, TLorentzVector lvBeam){return getPTLang(thetaCM, lvBeam, qPT_75);}, {"thetaCM","lvBeam"})
					 .Define("sqrangpt_75", [qPT_75](Double_t thetaCM, TLorentzVector lvBeam){return getPTRang(thetaCM, lvBeam, qPT_75);}, {"thetaCM","lvBeam"})
					 .Define("sqlangpt_0", [qPT_0](Double_t thetaCM, TLorentzVector lvBeam){return getPTLang(thetaCM, lvBeam, qPT_0);}, {"thetaCM","lvBeam"})
					 .Define("sqrangpt_0", [qPT_0](Double_t thetaCM, TLorentzVector lvBeam){return getPTRang(thetaCM, lvBeam, qPT_0);}, {"thetaCM","lvBeam"});


	//inDF.Define("vBeamLength", getBeamVector, {"vBeam"});
	//Define("fs_mc_weight", [] (double luminosity, andallotherargs...) { return lumiWeight(sampleTag, luminosity, genWeight, eventWeight, leptonWeight, jvtWeight, bTagWeight, pileupWeight, FFWeight), {"luminosity", "andallothercolumns"})
/*

	AELC *h2_Si = new ELC(2, 1, si_Nel, 2.35, si_A, si_Z, si_W, 100.,1500);
	std::cout<<"Energy "<<h2_Si->GetE(10.0, -3.5)<<std::endl;
*/



outDF.Snapshot("analyzed", outFilename.Data());
}

void ui()
{
	ROOT::EnableImplicitMT();
	numberOfThreads = 20;
	TStopwatch *stopwatch = new TStopwatch();
	
	rnd = new TRandom3();
	selector = TSelector::GetSelector("/home/zalewski/aku/analysis/simplifier123.C");
	std::cout<<"..."<<std::endl<<".."<<std::endl<<"."<<std::endl;

	if (cs::runNo==5)
	{
		selector = TSelector::GetSelector("/home/zalewski/aku/analysis/simplifier5.C");
	}

	loadParameters();
	TString str_name;
	TSystemDirectory *dir_data = new TSystemDirectory("data",raw_data_path.Data());
	TIter bluster = dir_data->GetListOfFiles();
	while (TObject *obj = bluster())
	{
		str_name = obj->GetName();
		//std::cout<<str_name<<std::endl;
		if (str_name.Contains("raw_geo") &&  //if we want to ommit some phrase (np.li9_3_21_0)
			!str_name.Contains("li9") //  "!"  sign is important
								)
		{
			TString inputFilePath = raw_data_path + str_name;
			//printf("fName:\t%s\n",inputFilePath.Data());
			//translator(inputFilePath);
			//cleaner(inputFilePath);
			//calibratorCaller(inputFilePath);
			//analysis(inputFilePath);
		}
	}
	stopwatch->Print();
}