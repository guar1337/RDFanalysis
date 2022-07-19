#include "/home/zalewski/aku/analysis/constants.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TCanvas.h"

Int_t varX, varZ;

std::string smallFile = "/home/zalewski/dataTmp/small/smallE.root";
enum sPar {sLang1, sLang2, sLang3, sRang, sTarPos, sDistL, sDistR};
std::vector<std::string> parNames = {"sLang1", "sLang2", "sLang3", "sRang", "sTarPos", "sDistL", "sDistR"};

Double_t tarPos[3], leftAngle[3], tarAngle[3];
Double_t sqlDist, sqrDist, varLDist, varRDist;
Double_t rightAngle;
Bool_t calculate;

const Int_t numberOfParameters = 7;
Double_t lowerParamBound[numberOfParameters], upperParamBound[numberOfParameters];
Int_t firstRun, lastRun;

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
	Double_t m_sqlangpp = m_lvBeam.Angle(m_lv1H.Vect())*TMath::RadToDeg();
	//Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
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
	//Double_t m_sqlangpp = m_lvBeam.Angle(m_lv1H.Vect())*TMath::RadToDeg();
	Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqrangpp;
}

Double_t getEnePP(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	m_lvBeam.SetE(cs::mass6He+155.0);
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
	Double_t m_sqlEpp = m_lv1H.E()-m_lv1H.M();
	//Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqlEpp;
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
	//Double_t m_sqlangdd = m_lvBeam.Angle(m_lv2H.Vect())*TMath::RadToDeg();
	Double_t m_sqrangdd = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
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
	Double_t m_sqlangdd = m_lvBeam.Angle(m_lv2H.Vect())*TMath::RadToDeg();
	//Double_t m_sqrangdd = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqlangdd;
}

Double_t getEneDD(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	Double_t m_ene = 155.0 + cs::mass6He;
	Double_t m_mom = sqrt(m_ene*m_ene - cs::mass6He*cs::mass6He);
	TLorentzVector m_lv6He(0.0, 0.0, m_mom, m_ene);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv2H.SetTheta(m_thetaCM);

	//m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);
	Double_t m_sqlEdd = m_lv2H.E()-m_lv2H.M();
	//Double_t m_sqrangdd = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqlEdd;
}

std::vector<Double_t> getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_1_Z, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t m_MWPC_2_Z)
{
	std::vector<Double_t> rvecMWPC(6);
	rvecMWPC.at(0) = 0.0 + m_MWPC_1_X - 1.0;
	rvecMWPC.at(1) = 0.0 + m_MWPC_1_Y - 2.1375;
	rvecMWPC.at(2) = 0.0 + m_MWPC_1_Z;
	rvecMWPC.at(3) = 0.0 + m_MWPC_2_X + 0.2; 
	rvecMWPC.at(4) = 0.0 + m_MWPC_2_Y - 1.125;
	rvecMWPC.at(5) = 0.0 + m_MWPC_2_Z;

	return rvecMWPC;
}

TVector3 getBeamVector(std::vector<Double_t> rvecMWPC, Double_t m_kinE)
{
	TVector3 m_beamVector(rvecMWPC[3] - rvecMWPC[0], rvecMWPC[4] - rvecMWPC[1], rvecMWPC[5] - rvecMWPC[2]);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

TVector3 getTarVertex(std::vector<Double_t> rvecMWPC, Int_t m_geo, const Double_t *m_pars)
{
	Double_t m_tarPos = 0.0;
	Double_t m_tarAngle;

	switch (m_geo)
	{
	case 1:
		m_tarAngle = 45.0*TMath::DegToRad();
		m_tarPos = tarPos[0];
		break;

	case 2:
		m_tarAngle = 6.0*TMath::DegToRad();
		m_tarPos = tarPos[1];
		break;

	case 3:
		m_tarAngle = 0.0*TMath::DegToRad();
		m_tarPos = tarPos[2];
		break;
	
	default:
		break;
	}

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

TVector3 getLeftDetVertex(Int_t m_xStrip, Int_t m_yStrip, Int_t m_geo)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system	
	Double_t X2hDet = -cs::widthStripX * xStrip * cos(leftAngle[m_geo-1]);
	Double_t Y2hDet = cs::widthStripY * yStrip;
	Double_t Z2hDet = cs::widthStripX * xStrip * sin(leftAngle[m_geo-1]);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip, Int_t m_geo)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * xStrip * cos(rightAngle);
	Double_t Y6HeDet = cs::widthStripY * yStrip;
	Double_t Z6HeDet = cs::widthStripX * xStrip * sin(rightAngle);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

TVector3 getLeftDetPosition(Int_t m_geo, const Double_t *m_pars)
{

	Double_t X2Hlab = sqlDist*sin(leftAngle[m_geo-1]) + (cs::sqlXzero) * cos(leftAngle[m_geo-1]);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(leftAngle[m_geo-1]) - (cs::sqlXzero) * sin(leftAngle[m_geo-1]);

	return TVector3(X2Hlab, Y2Hlab, Z2Hlab);
}

TVector3 getRightDetPosition(Int_t m_geo, const Double_t *m_pars)
{

	Double_t X6Helab = sqrDist*sin(-rightAngle) - (cs::sqlXzero) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart + cs::widthStripY;
	Double_t Z6Helab = sqrDist*cos(rightAngle) - (cs::sqlXzero) * sin(rightAngle);

	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

Double_t myAngAngFit(Double_t m_leftAngle, TLorentzVector m_lvBeam, Double_t tarMass)
{
	TLorentzVector m_lvBeamCopy(m_lvBeam);
	TLorentzVector m_TarCM{0.0,0.0,0.0,tarMass};
	TLorentzVector m_lvCM = m_lvBeam + m_TarCM;
	Double_t m_thetaCM = m_lvCM.Theta();
	TVector3 boostVect = m_lvCM.BoostVector();
	
	Double_t gammaSquare = m_lvCM.Gamma() *  m_lvCM.Gamma();
	Double_t tanSquare = pow(tan(m_leftAngle*TMath::DegToRad()),2);
	Double_t cosLeftAng = (1.0 - gammaSquare*tanSquare)/(1 + gammaSquare*tanSquare);
	Double_t thetaCM = TMath::Pi() - (acos(cosLeftAng)+m_thetaCM);

	m_lvBeam.Boost(-boostVect);
	//printf("TMath::Pi(): %f\tthetaCM: %f\tlAng: %f\n",TMath::Pi()*TMath::RadToDeg(), thetaCM*TMath::RadToDeg(), m_leftAngle);
	m_lvBeam.SetTheta(thetaCM);
	m_lvBeam.Boost(boostVect);
	/*
	m_TarCM.Boost(-boostVect);
	m_TarCM.SetTheta(acos(cosLeftAng)-m_Theta);
	m_TarCM.Boost(boostVect);
	//printf("sqlangIN: %f\tsqlangOUT: %f\tdiff: %f\tm_Theta: %f\n",m_leftAngle, m_TarCM.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg(), m_leftAngle-m_TarCM.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg(),m_Theta*TMath::RadToDeg());
	*/
	return m_lvBeam.Vect().Angle(m_lvBeamCopy.Vect())*TMath::RadToDeg();
}

Double_t drawPT_3H(Double_t m_thetaCM)
{
	Double_t beamEne = cs::mass6He + 155.0;
	Double_t beamMom = sqrt(beamEne*beamEne - cs::mass6He*cs::mass6He);
	TLorentzVector m_lv6He(0.0, 0.0, beamMom, beamEne);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + cs::Qpt;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2.0*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2.0*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	m_lv3H.SetTheta(m_thetaCM);
	m_lv4He.SetTheta(TMath::Pi()-m_thetaCM);

	m_lv3H.Boost(m_boostVect);
	m_lv4He.Boost(m_boostVect);
	return m_lv3H.Theta()*TMath::RadToDeg();
}

Double_t drawPT_4He(Double_t m_thetaCM)
{
	Double_t beamEne = cs::mass6He + 155.0;
	Double_t beamMom = sqrt(beamEne*beamEne - cs::mass6He*cs::mass6He);
	TLorentzVector m_lv6He(0.0, 0.0, beamMom, beamEne);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + cs::Qpt;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2*cs::mass4He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	Double_t m_cm4HekinE = m_finTcm*(m_finTcm+2*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm4Heene = m_cm4HekinE + cs::mass4He;
	Double_t m_cm4Hemom = sqrt(m_cm4Heene*m_cm4Heene - cs::mass4He*cs::mass4He);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	TLorentzVector m_lv4He(0.0,0.0,m_cm4Hemom, m_cm4Heene);

	m_lv3H.SetTheta(m_thetaCM);
	m_lv4He.SetTheta(TMath::Pi()-m_thetaCM);

	m_lv3H.Boost(m_boostVect);
	m_lv4He.Boost(m_boostVect);
	return m_lv4He.Theta()*TMath::RadToDeg();
}

Double_t recoPT(Double_t m_sqlang, TLorentzVector m_lvBeam, Int_t m_lang)
{
	//generate pt reaction (lv3H actually) with:
	//			theta in range 0:3.14
	//			beam vector
	//			q of the reaction
	Int_t multi = (gRandom->Integer(2)==1) ? 1 : -1;
	multi = 1;
	Double_t m_myAngle;
	TLorentzVector m_lv6He(m_lvBeam);
	// /m_lv6He.SetTheta(0.0);
	TLorentzVector m_lv1H(0.0,0.0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv1H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv1H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass1H) + cs::Qpt;

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
	//m_lv3H.SetTheta(m_thetaCM);
	//m_lv3H.Boost(m_boostVect);

	//reconstruct pt reaction (angle) knowing:
	//			beta and gamma of the CM, beta of 3H in CM, beta of 4He in CM
	//			angle between 3H and beam vectors
	Double_t m_betaCM3H = (sqrt(m_cm3HkinE*m_cm3HkinE+2*m_cm3HkinE*cs::mass3H))/(m_cm3HkinE+cs::mass3H);
	Double_t m_betaCM4He = (sqrt(m_cm4HekinE*m_cm4HekinE+2*m_cm4HekinE*cs::mass4He))/(m_cm4HekinE+cs::mass4He);
	Double_t m_beta3HRatio = m_betaCM/m_betaCM3H;
	Double_t m_beta4HeRatio = m_betaCM/m_betaCM4He;
	Double_t m_beta3HRatio2 = m_beta3HRatio*m_beta3HRatio;
	Double_t m_beta4HeRatio2 = m_beta4HeRatio*m_beta4HeRatio;

	Double_t m_gammaTan2 = std::pow(m_gammaCM*tan(m_sqlang*TMath::DegToRad()),2);

	//calculating CM angles based on left angle (m_lang) or right angle (!m_lang)
	if (m_lang==0)
	{
		multi = 1;
		Double_t m_cosThetaCM4He = (-m_beta4HeRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta4HeRatio2)))/(1.0+m_gammaTan2);
		Double_t m_thetaCM4He = acos(m_cosThetaCM4He);
		Double_t m_thetaCM3H = TMath::Pi() - m_thetaCM4He;
		m_myAngle = atan(sin(m_thetaCM3H)/(m_gammaCM*(cos(m_thetaCM3H)+m_betaCM/m_betaCM3H)));
	}
	
	else if(m_lang==1)
	{
		multi = 1;
		Double_t m_cosThetaCM3H = (-m_beta3HRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta3HRatio2)))/(1.0+m_gammaTan2);
		Double_t m_thetaCM3H = acos(m_cosThetaCM3H);
		Double_t m_thetaCM4He = TMath::Pi() - m_thetaCM3H;
		m_myAngle = atan(sin(m_thetaCM4He)/(m_gammaCM*(cos(m_thetaCM4He)+m_betaCM/m_betaCM4He)));
	}

	return m_myAngle*TMath::RadToDeg();
}

void drawer(int input1 = 0, int input2 = 0)
{

	varX = input1;
	varZ = input2;
	Double_t finalPars[7];
	finalPars[sLang1] = 0.3;
	finalPars[sLang2] = 0.0;
	finalPars[sLang3] = 0.0;
	finalPars[sRang] = -1.0;
	finalPars[sTarPos] = 10.0;
	finalPars[sDistL] = -20.0;
	finalPars[sDistR] = -30.0;

	varLDist = 0.0;
	varRDist = 0.0;

	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;

	tarPos[0] = finalPars[sTarPos];
	tarPos[1] = finalPars[sTarPos];
	tarPos[2] = finalPars[sTarPos];

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + finalPars[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + finalPars[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + finalPars[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + finalPars[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + finalPars[sDistL];
	sqrDist = 250.0 + finalPars[sDistR];
	ROOT::EnableImplicitMT();
	ROOT::RDataFrame smallDF("smallReal", smallFile);

		auto newDF = smallDF.Define("X1",[&finalPars](Double_t MWPC_1_X){return (MWPC_1_X);}, {"MWPC_1_X"})
									.Define("Y1",[&finalPars](Double_t MWPC_1_Y){return (MWPC_1_Y);}, {"MWPC_1_Y"})
									.Define("Z1",[&finalPars](){return -816.0;})
									.Define("X2",[&finalPars](Double_t MWPC_2_X){return (MWPC_2_X );}, {"MWPC_2_X"})
									.Define("Y2",[&finalPars](Double_t MWPC_2_Y){return (MWPC_2_Y);}, {"MWPC_2_Y"})
									.Define("Z2",[&finalPars](){return -270.0;})
									.Define("MWPC", getMWPC, {"X1", "Y1", "Z1", "X2", "Y2", "Z2"})
				.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
				.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
				.Define("tarVertex", [finalPars](std::vector<Double_t> rvecMWPC, Int_t geo){return getTarVertex(rvecMWPC, geo, finalPars);}, {"MWPC", "geo"})

							.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
							.Define("leftLabVertex", [finalPars](Int_t geo){return getLeftDetPosition(geo, finalPars);}, {"geo"})
							.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
							.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip", "geo"})
							.Define("rightLabVertex", [finalPars](Int_t geo){return getRightDetPosition(geo, finalPars);}, {"geo"})
							.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

							.Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
							.Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
							.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
							.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"});


	auto ppDF = newDF.Filter("pp").Cache<Double_t, Double_t, TLorentzVector, Int_t, Double_t>({"sqlang", "sqrang", "lvBeam", "geo", "resqlde1"});
	auto ddDF = newDF.Filter("dd").Cache<Double_t, Double_t, TLorentzVector, Int_t, Double_t>({"sqlang", "sqrang", "lvBeam", "geo", "resqlde2"});
	auto ptDF = newDF.Filter("pt").Cache<Double_t, Double_t, TLorentzVector, Int_t>({"sqlang", "sqrang", "lvBeam", "partPT"});
	auto ptDFLeft = newDF.Filter("pt && partPT==1").Cache<Double_t, Double_t, TLorentzVector, Int_t>({"sqlang", "sqrang", "lvBeam", "partPT"});
	auto ptDFRight = newDF.Filter("pt && partPT==0").Cache<Double_t, Double_t, TLorentzVector, Int_t>({"sqlang", "sqrang", "lvBeam", "partPT"});
	auto vecDF = newDF.Cache<TVector3, TVector3, TVector3, TVector3, TVector3, TVector3, TVector3>({"leftDetVertex", "leftLabVertex", "leftGlobVertex", "rightDetVertex", "rightLabVertex", "rightGlobVertex", "tarVertex"});

	auto nEventsPP = ppDF.Count().GetValue();
	auto sumPP = ppDF.Define("resqrang",[tarMass1H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass1H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");
						
	auto nEventsDD = ddDF.Count().GetValue();
	auto sumDD = ddDF.Define("resqrang",[tarMass2H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass2H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");

	auto nEventsPT = ptDF.Count().GetValue();
	auto sumPT = ptDF.Define("resqrang", recoPT, {"sqlang", "lvBeam", "partPT"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");

	auto nEventsPTLang = ptDFLeft.Count().GetValue()/(1000);
	auto sumPTLang = ptDFLeft.Define("resqrang", recoPT, {"sqlang", "lvBeam", "partPT"})
									 .Define("difSqrang", "pow(resqrang-sqrang,2)");

	auto nEventsPTRang = ptDFRight.Count().GetValue()/(1000);
	auto sumPTRight = ptDFRight.Define("resqlang", recoPT, {"sqrang", "lvBeam", "partPT"})
										.Define("difSqlang", "pow(resqlang-sqrang,2)");

	TGraph ppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ppAngAng_geo1 = sumPP.Filter("geo==1").Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ppAngAng_geo3 = sumPP.Filter("geo==3").Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ppAngE = sumPP.Filter("geo==1").Graph<Double_t, Double_t>("sqlang", "resqlde1").GetValue();

	TGraph ddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ddAngAng_geo1 = sumDD.Filter("geo==1").Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ddAngAng_geo3 = sumDD.Filter("geo==3").Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ddAngE = sumDD.Filter("geo==1").Graph<Double_t, Double_t>("sqlang", "resqlde2").GetValue();

	TGraph ptAngAngLeft = sumPTLang.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ptAngAngRight = sumPTRight.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	
	//TGraph fppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	//TGraph fddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	TGraph fptAngAngLeft = sumPTLang.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	TGraph fptAngAngRight = sumPTRight.Graph<Double_t, Double_t>("resqlang", "sqrang").GetValue();

	auto lDF = vecDF.Define("lX", "leftGlobVertex.X()")
					.Define("lDetX", "leftLabVertex.X()")
					.Define("lZ", "leftGlobVertex.Z()")
					.Define("lDetZ", "leftLabVertex.Z()");
	
	auto rDF = vecDF.Define("rX", "rightGlobVertex.X()")
					.Define("rDetX", "rightLabVertex.X()")
					.Define("rZ", "rightGlobVertex.Z()")
					.Define("rDetZ", "rightLabVertex.Z()");

	auto tarDF = vecDF.Define("tarX", "tarVertex.X()")
					  .Define("tarZ", "tarVertex.Z()");
	

	TGraph setupL = lDF.Graph<Double_t, Double_t>("lX", "lZ").GetValue();
	TGraph detL = lDF.Graph<Double_t, Double_t>("lDetX", "lDetZ").GetValue();
	TGraph setupR = rDF.Graph<Double_t, Double_t>("rX", "rZ").GetValue();
	TGraph detR = rDF.Graph<Double_t, Double_t>("rDetX", "rDetZ").GetValue();
	TGraph setupTar = tarDF.Graph<Double_t, Double_t>("tarX", "tarZ").GetValue();

	TProfile difPPSqrang1 = sumPP.Filter("geo==1").Profile1D<Double_t, Double_t>({"difPPSqrang", "difPPSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();
	TProfile difPPSqrang2 = sumPP.Filter("geo==2").Profile1D<Double_t, Double_t>({"difPPSqrang", "difPPSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();
	TProfile difPPSqrang3 = sumPP.Filter("geo==3").Profile1D<Double_t, Double_t>({"difPPSqrang", "difPPSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();

	TProfile difDDSqrang1 = sumDD.Filter("geo==1").Profile1D<Double_t, Double_t>({"difDDSqrang", "difDDSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();
	TProfile difDDSqrang2 = sumDD.Filter("geo==2").Profile1D<Double_t, Double_t>({"difDDSqrang", "difDDSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();
	TProfile difDDSqrang3 = sumDD.Filter("geo==3").Profile1D<Double_t, Double_t>({"difDDSqrang", "difDDSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();

	TProfile difPTSqrang = sumPT.Profile1D<Double_t, Double_t>({"difPTSqrang", "difPTSqrang", 45, 25, 85}, "sqlang", "difSqrang").GetValue();

	std::vector<Double_t> vAng3H;
	std::vector<Double_t> vAng4He;
	std::vector<Double_t> vAngPP_1H;
	std::vector<Double_t> vAngPP_6He;
	std::vector<Double_t> vEPP_1H;
	std::vector<Double_t> vAngDD_2H;
	std::vector<Double_t> vAngDD_6He;
	std::vector<Double_t> vEDD_2H;

	Double_t beamEneTmp = cs::mass6He + 160.0;
	Double_t beamMomTmp = sqrt(beamEneTmp*beamEneTmp - cs::mass6He*cs::mass6He);
	TLorentzVector lvBeamTmp(0.0, 0.0, beamMomTmp, beamEneTmp);

	for (int iii = 0; iii < 5000; iii++)
	{
		Double_t thetaCM = gRandom->Uniform(0.0,TMath::Pi());
		vAng3H.push_back(drawPT_3H(thetaCM));
		vAng4He.push_back(drawPT_4He(thetaCM));

		vAngPP_1H.push_back(getPPLang(thetaCM,lvBeamTmp));
		vAngPP_6He.push_back(getPPRang(thetaCM,lvBeamTmp));
		vEPP_1H.push_back(getEnePP(thetaCM,lvBeamTmp));

		vAngDD_2H.push_back(getDDLang(thetaCM,lvBeamTmp));
		vAngDD_6He.push_back(getDDRang(thetaCM,lvBeamTmp));
		vEDD_2H.push_back(getEneDD(thetaCM,lvBeamTmp));
	}
	
	TGraph ptAngAng(vAng3H.size(), &vAng3H[0], &vAng4He[0]);
	TGraph fppAngAng(vAngPP_1H.size(), &vAngPP_1H[0], &vAngPP_6He[0]);
	TGraph fddAngAng(vAngDD_2H.size(), &vAngDD_2H[0], &vAngDD_6He[0]);
	TGraph fppAngE(vAngPP_1H.size(), &vAngPP_1H[0], &vEPP_1H[0]);
	TGraph fddAngE(vAngDD_2H.size(), &vAngDD_2H[0], &vEDD_2H[0]);



	TCanvas myCanvas("myCanvas", "Minimize results", 1920, 1080);
	//myCanvas.SetBatch();
	myCanvas.Divide(2,2);

	myCanvas.cd(1);
	ppAngAng.GetXaxis()->SetLimits(25.0, 80.0);
	ppAngAng.GetYaxis()->SetRangeUser(3.0, 22.0);
	ppAngAng.GetXaxis()->SetTitle("Angle of {}^{1}H [LAB deg]");
	ppAngAng.GetYaxis()->SetTitle("Angle of {}^{6}He [LAB deg]");
	ppAngAng.SetTitle("Angle-angle relation for elastic scattering");
	ppAngAng.SetMarkerColor(kRed);
	ppAngAng.Draw("AP");
	ppAngAng_geo1.SetMarkerColor(kBlue);
	ppAngAng_geo1.SetMarkerStyle(7);
	ppAngAng_geo1.Draw("P,same");
	ppAngAng_geo3.SetMarkerColor(kBlue);
	ppAngAng_geo3.Draw("P,same");
	fppAngAng.SetMarkerStyle(7);
	fppAngAng.Draw("P,same");

/*
	TLatex *tex1 = new TLatex(26.05426,10.62081,"Elastic scattering on protons");
   tex1->SetLineWidth(2);
   tex1->Draw("same");
	TLatex *tex2 = new TLatex(39.47308,20.15142,"Elastic scattering on deuterons");
   tex2->SetLineWidth(2);
   tex2->Draw("same");


	myCanvas.cd(2);
	ddAngAng.GetXaxis()->SetLimits(30.0, 85.0);
	ddAngAng.GetYaxis()->SetRangeUser(0.0, 25.0);
	ddAngAng.GetXaxis()->SetTitle("Angle of {}^{1}H [LAB deg]");
	ddAngAng.GetYaxis()->SetTitle("Angle of {}^{6}He [LAB deg]");
	ddAngAng.SetTitle("Angle-angle relation for elastic scattering");*/
	ddAngAng.SetMarkerColor(kRed);
	ddAngAng.Draw("P,same");
	ddAngAng_geo1.SetMarkerStyle(7);
	ddAngAng_geo1.SetMarkerColor(kBlue);
	ddAngAng_geo1.Draw("P,same");
	ddAngAng_geo3.SetMarkerStyle(7);
	ddAngAng_geo3.SetMarkerColor(kBlue);
	ddAngAng_geo3.Draw("P,same");
	//fddAngAng.SetMarkerColor(kRed);
	//fddAngAng.SetMarkerStyle(7);
	fddAngAng.Draw("P,same");

	myCanvas.cd(2);
	ppAngE.GetXaxis()->SetLimits(60.0, 80.0);
	ppAngE.GetYaxis()->SetRangeUser(0.0, 20.0);
	ppAngE.SetMarkerColor(kRed);
	ppAngE.SetMarkerStyle(53);	
	ppAngE.Draw("AP");

	ddAngE.SetMarkerStyle(53);
	ddAngE.SetMarkerColor(kRed);
	ddAngE.Draw("P, same");

	fppAngE.SetMarkerStyle(7);
	fddAngE.SetMarkerStyle(7);
	fppAngE.Draw("P, same");
	fddAngE.Draw("P, same");

	myCanvas.cd(3);
	ptAngAngLeft.GetXaxis()->SetLimits(5.0, 45.0);
	ptAngAngLeft.GetYaxis()->SetRangeUser(0.0, 35.0);
	ptAngAngLeft.Draw("AP");
	ptAngAngRight.Draw("P, same");
	fptAngAngLeft.SetMarkerColor(kRed);
	fptAngAngRight.SetMarkerColor(kRed);
	fptAngAngLeft.Draw("P,same");
	fptAngAngRight.Draw("P,same");
	ptAngAng.SetMarkerColor(kBlue);
	ptAngAng.Draw("P,same");

	myCanvas.cd(4);
	setupL.GetXaxis()->SetLimits(-100.0, 200.0);
	setupL.GetYaxis()->SetRangeUser(-50.0, 300.0);
	setupL.Draw("AP");
	detL.SetMarkerStyle(53);
	detL.SetMarkerSize(2);
	detL.SetMarkerColor(kRed);
	detR.SetMarkerStyle(53);
	detR.SetMarkerSize(2);
	detR.SetMarkerColor(kRed);

	setupR.Draw("P,same");
	setupTar.Draw("P,same");
	detL.Draw("P,same");
	detR.Draw("P,same");




/*
	difPPSqrang1.GetXaxis()->SetLimits(25.0, 85.0);
	difPPSqrang1.GetYaxis()->SetRangeUser(-3.0, 3.0);
	difPPSqrang1.SetMarkerStyle(53);
	difPPSqrang1.SetMarkerSize(2);
	difPPSqrang1.SetMarkerColor(kRed);
	difPPSqrang1.Draw("P");

	difPPSqrang2.GetXaxis()->SetLimits(25.0, 85.0);
	difPPSqrang2.GetYaxis()->SetRangeUser(-3.0, 3.0);
	difPPSqrang2.SetMarkerStyle(53);
	difPPSqrang2.SetMarkerSize(2);
	difPPSqrang2.SetMarkerColor(kGreen);
	difPPSqrang2.Draw("P, same");

	difPPSqrang3.SetMarkerStyle(53);
	difPPSqrang3.SetMarkerSize(2);
	difPPSqrang3.SetMarkerColor(kBlue);
	difPPSqrang3.Draw("P, same");

	difDDSqrang1.SetMarkerStyle(54);
	difDDSqrang1.SetMarkerSize(2);
	difDDSqrang1.SetMarkerColor(kRed);
	difDDSqrang1.Draw("P, same");

	difDDSqrang2.SetMarkerStyle(54);
	difDDSqrang2.SetMarkerSize(2);
	difDDSqrang2.SetMarkerColor(kGreen);
	difDDSqrang2.Draw("P, same");

	difDDSqrang3.SetMarkerStyle(54);
	difDDSqrang3.SetMarkerSize(2);
	difDDSqrang3.SetMarkerColor(kBlue);
	difDDSqrang3.Draw("P, same");

	difPTSqrang.SetMarkerStyle(20);
	difPTSqrang.SetMarkerSize(1);
	difPTSqrang.SetMarkerColor(kBlack);
	difPTSqrang.Draw("P, same");
*/
	//TString outFname = "/home/zalewski/Desktop/6He/analysis/experimental/myData" + myRunNumber + "cvs";
	myCanvas.Print(TString::Format("/home/zalewski/Desktop/6He/analysis/exp/myData_%d_%d.png", varX, varZ));
	myCanvas.SaveAs("/home/zalewski/Desktop/macro.C");


//gPad->Print("/home/zalewski/Desktop/myPP.jpeg");

}