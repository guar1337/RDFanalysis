#include "/home/zalewski/aku_dec/analysis/ui.hh"

R__LOAD_LIBRARY(libgsl.so);
R__LOAD_LIBRARY(/home/zalewski/aku_dec/ELC/build/libEloss.so);
//R__LOAD_LIBRARY(/home/zalewski/aku/TELoss/libTELoss.so);
Bool_t flag=true;

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

void makeSmallFile()
{
	TChain smallChain("small");
	smallChain.Add("/home/zalewski/dataTmp/small/small1.root");
	smallChain.Add("/home/zalewski/dataTmp/small/small2.root");
	smallChain.Add("/home/zalewski/dataTmp/small/small3.root");

	ROOT::RDataFrame smallDF(smallChain);
	auto c = smallDF.Count();
	Int_t customClusterSize = *c / numberOfThreads;
	ROOT::RDF::RSnapshotOptions myOpts;
	myOpts.fAutoFlush=customClusterSize;

	smallDF.Snapshot("small", "/home/zalewski/dataTmp/small/small.root", columnList, myOpts);
}

Double_t getDiff(Double_t m_mc2H, TLorentzVector m_lv2H_CM, TLorentzVector m_flvBeam)
{
	Double_t realAngle = m_lv2H_CM.Vect().Angle(m_flvBeam.Vect())*TMath::RadToDeg();
	return m_mc2H-realAngle;//m_lv2H_CM.Theta()*TMath::RadToDeg();

}

void joinDE()
{
	TChain dEChain("analyzed");
	dEChain.Add("/home/zalewski/dataTmp/dE/geo1/dE_geo1.root");
	dEChain.Add("/home/zalewski/dataTmp/dE/geo2/dE_geo2.root");
	dEChain.Add("/home/zalewski/dataTmp/dE/geo3/dE_geo3.root");

	ROOT::RDataFrame dEDF(dEChain);
	dEDF.Snapshot("analyzed", "/home/zalewski/dataTmp/dE/dE_geo.root");

}

//d,p reaction
Double_t getDPLang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	Double_t m_Q = 0.0;
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm1HkinE = m_finTcm*(m_finTcm+2*cs::mass7He)/(2.0*m_Ecm);
	Double_t m_cm1Hene = m_cm1HkinE + cs::mass1H;
	Double_t m_cm1Hmom = sqrt(m_cm1Hene*m_cm1Hene - cs::mass1H*cs::mass1H);

	TLorentzVector m_lv1H(0.0,0.0,m_cm1Hmom, m_cm1Hene);
	m_lv1H.SetTheta(m_thetaCM);
	m_lv1H.Boost(m_boostVect);

	return m_lv1H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t getDPRang(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	Double_t m_Q = 0.0;
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_finTcm = m_Ecm - (cs::mass6He + cs::mass2H) + m_Q;

	Double_t m_cm7HekinE = m_finTcm*(m_finTcm+2*cs::mass2H)/(2.0*m_Ecm);
	Double_t m_cm7Heene = m_cm7HekinE + cs::mass7He;
	Double_t m_cm7Hemom = sqrt(m_cm7Heene*m_cm7Heene - cs::mass7He*cs::mass7He);

	TLorentzVector m_lv7He(0.0,0.0,m_cm7Hemom, m_cm7Heene);

	m_lv7He.SetTheta(TMath::Pi() - m_thetaCM);
	m_lv7He.Boost(m_boostVect);
	return m_lv7He.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}


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

Double_t reco1HEne(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv1H(0,0,0,cs::mass1H);
	TLorentzVector m_lvCM = m_lv6He+m_lv1H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv1H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv1H.SetTheta(2.0*m_sqlang*TMath::DegToRad());

	//m_lv6He.Boost(m_boostVect);
	m_lv1H.Boost(m_boostVect);
	Double_t m_ene1H = m_lv1H.E() - cs::mass1H;
	//Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_ene1H;
}

Double_t reco2HEne(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;

	TVector3 m_boostVect = m_lvCM.BoostVector();

	//m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);

	//m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);
	m_lv2H.SetTheta(2.0*m_sqlang*TMath::DegToRad());

	//m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);
	Double_t m_ene2H = m_lv2H.E() - cs::mass2H;
	//Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_ene2H;
}

Double_t getAngCorr(TLorentzVector m_lv2H, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvCM = m_lvBeam + m_lv2H;
	return m_lvCM.Theta()*TMath::RadToDeg();
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

Double_t calibrateCsIL_2H(ROOT::VecOps::RVec<unsigned short> rawCsI_L, int id)
{
	Double_t energy = csIcalPars[id][0] + 
					  csIcalPars[id][1]*rawCsI_L[id] + 
					  csIcalPars[id][2]*pow(rawCsI_L[id],2) + 
					  csIcalPars[id][3]*pow(rawCsI_L[id],3) + 
					  csIcalPars[id][4]*pow(rawCsI_L[id],4);
	return energy;
}

Double_t recoDT(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	//generate pt reaction (lv3H actually) with:
	//			theta in range 0:3.14
	//			beam vector
	//			q of the reaction
	Int_t multi = (rnd->Integer(2)==1) ? 1 : -1;
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

	//calculating CM angles
	Double_t m_gammaTan2 = std::pow(m_gammaCM*tan(m_sqlang*TMath::DegToRad()),2);
	Double_t m_cosThetaCM3H = (-m_beta3HRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta3HRatio2)))/(1.0+m_gammaTan2);
	Double_t m_thetaCM3H = acos(m_cosThetaCM3H);
	Double_t m_thetaCM4He = TMath::Pi() - m_thetaCM3H;
	Double_t m_4He_lab = atan(sin(m_thetaCM4He)/(m_gammaCM*(cos(m_thetaCM4He)+m_betaCM/m_betaCM4He)));
	Double_t m_3H_lab = atan(sin(m_thetaCM3H)/(m_gammaCM*(cos(m_thetaCM3H)+m_betaCM/m_betaCM3H)));
	return m_4He_lab*TMath::RadToDeg();
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

Double_t getEnePP(Double_t m_thetaCM, TLorentzVector m_lvBeam)
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
	Double_t m_sqlEpp = m_lv1H.E()-cs::mass1H;
	//Double_t m_sqrangpp = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqlEpp;
}

Double_t getEneDD(Double_t m_thetaCM, TLorentzVector m_lvBeam)
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
	Double_t m_sqlEdd = m_lv2H.E()-cs::mass2H;
	//Double_t m_sqrangdd = m_lvBeam.Angle(m_lv6He.Vect())*TMath::RadToDeg();
	return m_sqlEdd;
}

Double_t getCMAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getCMAngle1H_mod(TLorentzVector m_lv2H, TLorentzVector m_lvBeam, Double_t fTheta)
{
	TVector3 v2H = m_lv2H.Vect();
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-boostVect);
	Double_t m_LAng = v2H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();

	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng))*TMath::RadToDeg();
	
	Double_t beamXLab = m_lvBeam.Z();
	Double_t beamYLab = m_lvBeam.Y();
	Double_t beamZLab = m_lvBeam.Z();

	TRotation myRot;
	TVector3 newX(beamXLab,	0,			0);
	TVector3 newY(0,		beamYLab,	0);
	TVector3 newZ(0,		0,			beamZLab);
	//myRot.RotateAxes(newX,newY,newZ);
	TVector3 vBeam = m_lvBeam.Vect();
	myRot.SetZAxis(vBeam);
	//myRot.Invert();
	v2H.Transform(myRot);
	Double_t m_LAng_rot = v2H.Theta()*TMath::RadToDeg();

	Double_t m_tanSquare_rot = tan(m_LAng_rot*TMath::DegToRad()) * tan(m_LAng_rot*TMath::DegToRad());
	Double_t m_cosLeftAng_rot = (1.0 - m_gammaSquare2H*m_tanSquare_rot)/(1.0 + m_gammaSquare2H*m_tanSquare_rot);
	Double_t m_sqlangCM_rot = (acos(m_cosLeftAng_rot))*TMath::RadToDeg();

	//printf("%f\t%f\t%f\t%f\t%f\t%f\n", m_LAng, m_LAng_rot, m_LAng_rot-m_LAng, 180.0-m_sqlangCM, 180.0-m_sqlangCM_rot, fTheta);
	return m_sqlangCM;
}

Double_t getCMAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getCMAngle2H_mod(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Int_t getStripNumber(ROOT::RVecD &inputArray)
{
	Int_t m_stripNo = std::distance(inputArray.begin(),std::max_element(inputArray.begin(), inputArray.end()));
	return m_stripNo;
}

TLorentzVector getLV1H(Double_t m_sqlde, Double_t m_sqletot, TVector3 m_v1H)
{
	Double_t m_ene1H = m_sqlde + m_sqletot + cs::mass1H;
	Double_t m_mom1H = sqrt(m_ene1H*m_ene1H - cs::mass1H*cs::mass1H);
	m_v1H.SetMag(m_mom1H);
	return TLorentzVector(m_v1H, m_ene1H);
}

TLorentzVector getreLV1H(Double_t m_sqlde, Double_t m_sqletot, TVector3 m_v1H)
{
	Double_t m_reSqlde = h1_Si->GetE(m_sqlde, -(ionDeadLayer+3.0));
	Double_t m_ene1H = m_reSqlde + cs::mass1H;
	Double_t m_mom1H = sqrt(m_ene1H*m_ene1H - cs::mass1H*cs::mass1H);
	m_v1H.SetMag(m_mom1H);
	return TLorentzVector(m_v1H, m_ene1H);
}

TLorentzVector getLV2H(Double_t m_sqlde, Double_t m_sqletot, TVector3 m_v2H)
{
	Double_t m_ene2H = m_sqlde + cs::mass2H;
	Double_t m_mom2H = sqrt(m_ene2H*m_ene2H - cs::mass2H*cs::mass2H);
	m_v2H.SetMag(m_mom2H);
	return TLorentzVector(m_v2H, m_ene2H);
}

TLorentzVector getLV3H(Double_t m_sqlde, Double_t m_sqletot, TVector3 m_v2H)
{
	Double_t m_ene3H = m_sqlde + m_sqletot + cs::mass3H;
	Double_t m_mom3H = sqrt(m_ene3H*m_ene3H - cs::mass3H*cs::mass3H);
	m_v2H.SetMag(m_mom3H);
	return TLorentzVector(m_v2H, m_ene3H);
}

TLorentzVector getLV4He(Double_t m_sqrde, Double_t m_sqretot, TVector3 m_v4He)
{
	Double_t m_ene4He = m_sqrde + m_sqretot + cs::mass4He;
	Double_t m_mom4He = sqrt(m_ene4He*m_ene4He - cs::mass4He*cs::mass4He);
	m_v4He.SetMag(m_mom4He);
	return TLorentzVector(m_v4He, m_ene4He);
}

TLorentzVector getLV5He(Double_t m_sqrde, Double_t m_sqretot, TVector3 m_v5He)
{
	Double_t m_ene5He = m_sqrde + m_sqretot + cs::mass5He;
	Double_t m_mom5He = sqrt(m_ene5He*m_ene5He - cs::mass5He*cs::mass5He);
	m_v5He.SetMag(m_mom5He);
	return TLorentzVector(m_v5He, m_ene5He);
}

Double_t getCMAngleDT(Double_t m_sqrang, TLorentzVector m_lvBeam)
{
	Int_t m_lang = 0;
	Int_t multi = (gRandom->Integer(2)==1) ? 1 : -1;
	multi = 1;
	Double_t m_myAngle;
	TLorentzVector m_lv6He(m_lvBeam);
	// /m_lv6He.SetTheta(0.0);
	TLorentzVector m_lv2H(0.0,0.0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();
	m_lv2H.Boost(-m_boostVect);
	m_lv6He.Boost(-m_boostVect);

	Double_t m_Ecm = m_lv2H.E() + m_lv6He.E();
	Double_t m_Q = (cs::mass6He + cs::mass2H) - (cs::mass3H+cs::mass5He);
	Double_t m_iniTcm = m_Ecm - (cs::mass6He + cs::mass2H);
	Double_t m_finTcm = m_iniTcm + cs::Qdt;

	Double_t m_cm3HkinE = m_finTcm*(m_finTcm+2.0*cs::mass5He)/(2.0*m_Ecm);
	Double_t m_cm3Hene = m_cm3HkinE + cs::mass3H;
	Double_t m_cm3Hmom = sqrt(m_cm3Hene*m_cm3Hene - cs::mass3H*cs::mass3H);

	Double_t m_cm5HekinE = m_finTcm*(m_finTcm+2.0*cs::mass3H)/(2.0*m_Ecm);
	Double_t m_cm5Heene = m_cm5HekinE + cs::mass4He;
	Double_t m_cm5Hemom = sqrt(m_cm5Heene*m_cm5Heene - cs::mass5He*cs::mass5He);

	TLorentzVector m_lv3H(0.0,0.0,m_cm3Hmom, m_cm3Hene);
	TLorentzVector m_lv5He(0.0,0.0,m_cm5Hemom, m_cm5Heene);

	Double_t m_betaCM = m_lvCM.Beta();
	Double_t m_gammaCM = m_lvCM.Gamma();
	Double_t m_beta3HCM = m_lv3H.Beta();
	Double_t m_beta4HeCM = m_lv5He.Beta();
	//m_lv3H.SetTheta(m_thetaCM);
	//m_lv3H.Boost(m_boostVect);

	//reconstruct pt reaction (angle) knowing:
	//			beta and gamma of the CM, beta of 3H in CM, beta of 4He in CM
	//			angle between 3H and beam vectors
	Double_t m_betaCM3H = (sqrt(m_cm3HkinE*m_cm3HkinE+2.0*m_cm3HkinE*cs::mass3H))/(m_cm3HkinE+cs::mass3H);
	Double_t m_betaCM5He = (sqrt(m_cm5HekinE*m_cm5HekinE+2.0*m_cm5HekinE*cs::mass5He))/(m_cm5HekinE+cs::mass5He);
	Double_t m_beta3HRatio = m_betaCM/m_betaCM3H;
	Double_t m_beta5HeRatio = m_betaCM/m_betaCM5He;
	Double_t m_beta3HRatio2 = m_beta3HRatio*m_beta3HRatio;
	Double_t m_beta5HeRatio2 = m_beta5HeRatio*m_beta5HeRatio;

	Double_t m_gammaTan2 = std::pow(m_gammaCM*tan(m_sqrang*TMath::DegToRad()),2);
	Double_t m_thetaCM5He;
	//calculating CM angles based on left angle (m_lang) or right angle (!m_lang)
	if (m_lang==0)
	{
		Double_t m_cosThetaCM5He = (-m_beta5HeRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta5HeRatio2)))/(1.0+m_gammaTan2);
		m_thetaCM5He = acos(m_cosThetaCM5He);
		Double_t m_thetaCM3H = TMath::Pi() - m_thetaCM5He;
		m_myAngle = atan(sin(m_thetaCM3H)/(m_gammaCM*(cos(m_thetaCM3H)+m_betaCM/m_betaCM3H)));
	}
	
	else if(m_lang==1)
	{
		Double_t m_cosThetaCM3H = (-m_beta3HRatio*m_gammaTan2+multi*sqrt(1.0+m_gammaTan2*(1.0-m_beta3HRatio2)))/(1.0+m_gammaTan2);
		Double_t m_thetaCM3H = acos(m_cosThetaCM3H);
		m_thetaCM5He = TMath::Pi() - m_thetaCM3H;
		m_myAngle = atan(sin(m_thetaCM5He)/(m_gammaCM*(cos(m_thetaCM5He)+m_betaCM/m_betaCM5He)));
	}

	return m_thetaCM5He*TMath::RadToDeg();
}

Double_t getCMAngleDTene(Double_t m_sqlde, Double_t m_sqletot, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He_DT(m_lvBeam);
	TLorentzVector m_lv2H_DT(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM_DT = m_lv6He_DT + m_lv2H_DT;

	TVector3 m_boostVect_DT = m_lvCM_DT.BoostVector();
	m_lv6He_DT.Boost(-m_boostVect_DT);
	m_lv2H_DT.Boost(-m_boostVect_DT);
	//Let's go d,t now!
	Double_t Qdt = (cs::mass2H + cs::mass6He) - (cs::mass3H + cs::mass5He);
	Double_t eneCM = m_lv2H_DT.E() + m_lv6He_DT.E();
	Double_t kinECMBefore = m_lv2H_DT.E() - m_lv2H_DT.M() + m_lv6He_DT.E() - m_lv6He_DT.M();
	Double_t kinECMAfter = kinECMBefore + Qdt;
	Double_t kinE3H	= (kinECMAfter/(2.0*eneCM))*(kinECMAfter+2.0*cs::mass5He);
	Double_t m_beta = m_lvCM_DT.Beta();
	Double_t m_gamma = m_lvCM_DT.Gamma();
	Double_t m_cpCM = sqrt(kinE3H*(kinE3H+2.0*cs::mass3H));
	Double_t cosAng = m_sqlde + m_sqletot - (m_gamma-1.0)*cs::mass3H - m_gamma*kinE3H;
	Double_t CMangle = acos(cosAng/(m_gamma*m_beta*m_cpCM)) * TMath::RadToDeg();
	return 180.0 - CMangle;
}

TLorentzVector getreLV2H(Double_t m_sqlde, Double_t m_sqletot, TVector3 m_v2H)
{
	Double_t m_reSqlde = h2_Si->GetE(m_sqlde, -(ionDeadLayer+3.0));
	Double_t m_ene2H = m_reSqlde + cs::mass2H;
	if (m_sqletot>0.0)
	{
		m_ene2H += m_sqletot;
	}
	
	Double_t m_mom2H = sqrt(m_ene2H*m_ene2H - cs::mass2H*cs::mass2H);
	m_v2H.SetMag(m_mom2H);
	return TLorentzVector(m_v2H, m_ene2H);
}

Double_t getMissingMass1H(TLorentzVector m_lv1H, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lv6He = (m_lvBeam + m_lvTar) - m_lv1H;
	return m_lv6He.M() - cs::mass6He;
}

Double_t getMissingMass2H(TLorentzVector m_lv2H, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	if (flag)
	{
		printf("m2H: %f\tm6He: %f\n", m_lv2H.M(), m_lvBeam.M());
		flag=false;
	}
	TLorentzVector m_lv6He = (m_lvBeam + m_lvTar) - m_lv2H;
	return m_lv6He.M() - cs::mass6He;
}

Double_t getMissingMass3H(TLorentzVector m_lv3H, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	if (flag)
	{
		printf("m2H: %f\tm6He: %f\n", m_lv3H.M(), m_lvBeam.M());
		flag=false;
	}
	TLorentzVector m_lv6He = (m_lvBeam + m_lvTar) - m_lv3H;
	return m_lv6He.M() - cs::mass5He;
}

Double_t getMissingMass4He(TLorentzVector m_lv4He, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	if (flag)
	{
		printf("m2H: %f\tm6He: %f\n", m_lv4He.M(), m_lvBeam.M());
		flag=false;
	}
	TLorentzVector m_lv3H = (m_lvBeam + m_lvTar) - m_lv4He;
	return m_lv3H.M() - cs::mass3H;
}

Double_t getMissingMassN(TLorentzVector m_lv4He, TLorentzVector m_lv3H, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	if (flag)
	{
		printf("m2H: %f\tm6He: %f\n", m_lv4He.M(), m_lvBeam.M());
		flag=false;
	}
	TLorentzVector m_neutron = (m_lvBeam + m_lvTar) - (m_lv4He + m_lv3H);
	return m_neutron.M() - cs::massN;
}

Double_t inelSQLang2H_18(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);
	//Let's go inelastic now!
	//kinetic energy in the system
	Double_t eneCM = m_lv2H.E() + m_lv6He.E();
	Double_t kinECMBefore = m_lv2H.E() - m_lv2H.M() + m_lv6He.E() - m_lv6He.M();
	Double_t kinECMAfter = kinECMBefore - 1.8;
	Double_t kinE2H	= (kinECMAfter/(2.0*eneCM))*(kinECMAfter+2.0*(cs::mass6He+1.8));
	Double_t ene2H	= kinE2H + cs::mass2H;
	Double_t kinE6He = (kinECMAfter/(2.0*eneCM))*(kinECMAfter+2.0*cs::mass2H);
	Double_t ene6He = kinE6He + cs::mass6He + 1.8;
	Double_t newMom2H	= sqrt(ene2H*ene2H-cs::mass2H*cs::mass2H);
	Double_t newMom6He = sqrt(ene6He*ene6He-(cs::mass6He+1.8)*(cs::mass6He+1.8));

	TVector3 v2H(m_lv2H.Vect());
	TVector3 v6He(m_lv6He.Vect());
	v2H.SetMag(newMom2H);
	v6He.SetMag(newMom6He);
	m_lv2H.SetVectMag(v2H, cs::mass2H);
	m_lv2H.SetTheta(m_thetaCM);
	m_lv6He.SetVectMag(v6He, cs::mass6He+1.8);
	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);


	//NEUTRON EMISSION - negligible lifetime of excited 6He allows for immidiate decay
	//TLorentzVector lv6He_exc_IN = m_lv6He;
	TVector3 boostVect_6He = m_lv6He.BoostVector();
	//neutrons are sitting and waiting for some action
	Double_t decayE = cs::mass6He + 1.8 - (cs::mass4He + 2.0*cs::massN);
	Double_t momRatio = rnd->Uniform(0.1,1);

	//Setting neutron vectors
	//4He + 2n decay of 6He is classical - differences are negligible
	Double_t neut1Theta_CM = acos(rnd->Uniform(-1.0,1.0));
	Double_t neut1Phi_CM = rnd->Uniform(0.0,2.0*TMath::Pi());
	Double_t neut2Theta_CM = acos(rnd->Uniform(-1.0,1.0));
	Double_t neut2Phi_CM = rnd->Uniform(0.0,2.0*TMath::Pi());

	TVector3 vneu1;
	vneu1.SetPtThetaPhi(1.0,neut1Theta_CM,neut1Phi_CM);
	TVector3 vneu2;
	vneu2.SetPtThetaPhi(1.0,neut2Theta_CM,neut2Phi_CM);


	Double_t cos_angle_between_neutrons = cos(vneu1.Angle(vneu2));
	Double_t L_fac = momRatio * momRatio + 2.0*momRatio * cos_angle_between_neutrons + 1.0;
	Double_t momNeu1=sqrt(2*decayE*cs::mass4He*cs::massN/(cs::mass4He*(1.0+momRatio*momRatio)+cs::massN*L_fac));
	Double_t momNeu2 = momRatio*momNeu1;
	Double_t mom4He	= momNeu1*sqrt(L_fac);
	Double_t eneNeu1 = momNeu1*momNeu1/(2.0*cs::massN);
	Double_t eneNeu2 = momNeu2*momNeu2/(2.0*cs::massN);
	Double_t ene4He	= mom4He*mom4He/(2.0*cs::mass4He);

	vneu1.SetMag(momNeu1);
	vneu2.SetMag(momNeu2);
	TVector3 v4He(-(vneu1+vneu2));
	//printf("Mag from math: %f\tMag from neut: %f\n", mom4He, v4He.Mag());
	
	//TLorentzVector lvNeut1(vneu1,cs::massN+eneNeu1);
	//TLorentzVector lvNeut2(vneu2,cs::massN+eneNeu2);
	TLorentzVector lv4He;
	lv4He.SetE(cs::mass4He+ene4He);

	lv4He.Boost(boostVect_6He);
	//lvNeut1.Boost(boostVect_6He);
	//lvNeut2.Boost(boostVect_6He);
	
	return m_lv2H.Angle(m_lvBeam.Vect())*TMath::RadToDeg();
}

Double_t inelSQRang2H_18(Double_t m_thetaCM, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lv6He(m_lvBeam);
	TLorentzVector m_lv2H(0,0,0,cs::mass2H);
	TLorentzVector m_lvCM = m_lv6He+m_lv2H;
	TVector3 m_boostVect = m_lvCM.BoostVector();

	m_lv6He.Boost(-m_boostVect);
	m_lv2H.Boost(-m_boostVect);
	//Let's go inelastic now!
	//kinetic energy in the system
	Double_t newEcm	= m_lv2H.E() + m_lv6He.E();
	Double_t newE2H	= (1.0/(2.0*newEcm))*(newEcm*newEcm - (cs::mass6He+1.8)*(cs::mass6He+1.8) + cs::mass2H*cs::mass2H);
	Double_t newE6He	= (1.0/(2.0*newEcm))*(newEcm*newEcm + (cs::mass6He+1.8)*(cs::mass6He+1.8) - cs::mass2H*cs::mass2H);
	Double_t newMom2H	= sqrt(newE2H*newE2H-cs::mass2H*cs::mass2H);
	Double_t newMom6He = sqrt(newE6He*newE6He-(cs::mass6He+1.8)*(cs::mass6He+1.8));

	m_lv2H.SetE(newE2H);
	m_lv2H.SetTheta(m_thetaCM);
	m_lv6He.SetE(newE6He);
	m_lv6He.SetTheta(TMath::Pi()-m_thetaCM);

	m_lv6He.Boost(m_boostVect);
	m_lv2H.Boost(m_boostVect);


	//NEUTRON EMISSION - negligible lifetime of excited 6He allows for immidiate decay
	//TLorentzVector lv6He_exc_IN = m_lv6He;
	TVector3 boostVect_6He = m_lv6He.BoostVector();
	//neutrons are sitting and waiting for some action
	Double_t decayE = cs::mass6He + 1.8 - (cs::mass4He + 2.0*cs::massN);
	Double_t momRatio = rnd->Uniform(0.1,1);

	//Setting neutron vectors
	//4He + 2n decay of 6He is classical - differences are negligible
	Double_t neut1Theta_CM = acos(rnd->Uniform(-1.0,1.0));
	Double_t neut1Phi_CM = rnd->Uniform(0.0,2.0*TMath::Pi());
	Double_t neut2Theta_CM = acos(rnd->Uniform(-1.0,1.0));
	Double_t neut2Phi_CM = rnd->Uniform(0.0,2.0*TMath::Pi());

	TVector3 vneu1;
	vneu1.SetPtThetaPhi(1.0,neut1Theta_CM,neut1Phi_CM);
	TVector3 vneu2;
	vneu2.SetPtThetaPhi(1.0,neut2Theta_CM,neut2Phi_CM);


	Double_t cos_angle_between_neutrons = cos(vneu1.Angle(vneu2));
	Double_t L_fac = momRatio * momRatio + 2.0*momRatio * cos_angle_between_neutrons + 1.0;
	Double_t momNeu1=sqrt(2*decayE*cs::mass4He*cs::massN/(cs::mass4He*(1.0+momRatio*momRatio)+cs::massN*L_fac));
	Double_t momNeu2 = momRatio*momNeu1;
	Double_t mom4He	= momNeu1*sqrt(L_fac);
	Double_t eneNeu1 = momNeu1*momNeu1/(2.0*cs::massN);
	Double_t eneNeu2 = momNeu2*momNeu2/(2.0*cs::massN);
	Double_t ene4He	= mom4He*mom4He/(2.0*cs::mass4He);

	vneu1.SetMag(momNeu1);
	vneu2.SetMag(momNeu2);
	TVector3 v4He(-(vneu1+vneu2));
	//printf("Mag from math: %f\tMag from neut: %f\n", mom4He, v4He.Mag());
	
	//TLorentzVector lvNeut1(vneu1,cs::massN+eneNeu1);
	//TLorentzVector lvNeut2(vneu2,cs::massN+eneNeu2);
	TLorentzVector lv4He;
	lv4He.SetE(cs::mass4He+ene4He);

	lv4He.Boost(boostVect_6He);
	//lvNeut1.Boost(boostVect_6He);
	//lvNeut2.Boost(boostVect_6He);
	
	return lv4He.Vect().Angle(m_lvBeam.Vect())*TMath::RadToDeg();
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

TVector3 getLeftDetPosition(Int_t m_geo)
{
	Double_t X2Hlab = sqlDist*sin(leftAngle[m_geo-1]) + (cs::sqlXzero) * cos(leftAngle[m_geo-1]);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripX;
	Double_t Z2Hlab = sqlDist*cos(leftAngle[m_geo-1]) - (cs::sqlXzero) * sin(leftAngle[m_geo-1]);
	return TVector3(X2Hlab, Y2Hlab, Z2Hlab);
}

TVector3 getRightDetPosition(Int_t m_geo)
{
	Double_t X6Helab = sqrDist*sin(-rightAngle) - (cs::sqlXzero) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart;
	Double_t Z6Helab = sqrDist*cos(rightAngle) - (cs::sqlXzero) * sin(rightAngle);
	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

ROOT::RVecD getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_1_Z, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t m_MWPC_2_Z)
{
	ROOT::RVecD rvecMWPC(6);
	rvecMWPC.at(0) = m_MWPC_1_X;
	rvecMWPC.at(1) = m_MWPC_1_Y;
	rvecMWPC.at(2) = m_MWPC_1_Z;
	rvecMWPC.at(3) = m_MWPC_2_X;
	rvecMWPC.at(4) = m_MWPC_2_Y;
	rvecMWPC.at(5) = m_MWPC_2_Z;

	return rvecMWPC;
}

TVector3 getBeamVector(ROOT::RVecD rvecMWPC, Double_t m_kinE)
{
	TVector3 m_beamVector(rvecMWPC[3] - rvecMWPC[0], rvecMWPC[4] - rvecMWPC[1], rvecMWPC[5] - rvecMWPC[2]);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

TVector3 getTarVertex(ROOT::RVecD rvecMWPC, Int_t m_geo)
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

TVector3 getTarOnDet(ROOT::RVecD rvecMWPC)
{
	Double_t m_detPos = detDist;
	Double_t m_detAngle = 0.0;

	Double_t m_dX = rvecMWPC[3] - rvecMWPC[0];
	Double_t m_dY = rvecMWPC[4] - rvecMWPC[1];
	Double_t m_dZ = rvecMWPC[5] - rvecMWPC[2];
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, m_detPos);
	TVector3 m_beamPoint(rvecMWPC[3], rvecMWPC[4], rvecMWPC[5]);
	TVector3 m_tarPerpendicular(sin(m_detAngle), 0.0, cos(m_detAngle));
	Double_t m_dCoeff = ((m_tarPoint-m_beamPoint).Dot(m_tarPerpendicular))/(m_vBeam.Dot(m_tarPerpendicular));
		
	Double_t m_evX = rvecMWPC[3] + m_dX * m_dCoeff;
	Double_t m_evY = rvecMWPC[4] + m_dY * m_dCoeff;
	Double_t m_evZ = rvecMWPC[5] + m_dZ * m_dCoeff;

	return TVector3(m_evX, m_evY, m_evZ);
}

Double_t getSQlde(ROOT::RVecD &leftDetectorArray)
{
	Double_t myEne = 0.0;
	for (Double_t sqlde : leftDetectorArray)
	{
		if (sqlde>1.0) myEne+=sqlde;
	}
	return myEne;

}

TVector3 getZeroDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * xStrip;
	Double_t Y6HeDet = cs::widthStripY * yStrip;
	Double_t Z6HeDet = cs::widthStripX * xStrip;
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

TVector3 getZeroLabVertex()
{
	Double_t X6Helab = -cs::sqlXzero;
	Double_t Y6Helab = cs::sqrYstart;
	Double_t Z6Helab = detDist;
	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

Double_t getKineticEnergy(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);	
	return he6_Si->GetE(cs::mass6He*(m_gamma-1.0), beamDeadLayer-20.0/cos(TMath::PiOver4()));
	//return cs::mass6He*(m_gamma-1.0);
}

Double_t getRawKineticEnergy(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);	
	return he6_Si->GetE(cs::mass6He*(m_gamma-1.0), beamDeadLayer+20.0/cos(TMath::PiOver4()));
	//return cs::mass6He*(m_gamma-1.0);
}

Double_t getEloss(Double_t m_tof)
{
	Double_t m_beta_squared= pow((cs::tofBase/m_tof)/cs::c, 2.0);
	Double_t m_gamma=1.0/sqrt(1.0-m_beta_squared);
	Double_t m_kinE = cs::mass6He*(m_gamma-1.0);
	return cs::mass6He*(m_gamma-1.0) - he6_Si->GetE(m_kinE, beamDeadLayer);
}

bool filterLeftDetector(ROOT::RVecD &leftDetectorArray)
{
	int aboveThresholdEventCounter(0);
	auto leftDetectorArraySize = leftDetectorArray.size();
	for (size_t iii = 0; iii < leftDetectorArraySize; iii++)
	{
		if (leftDetectorArray[iii]>0.1)
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

bool filterRightDetector(ROOT::RVecD &rightDetectorArray)
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
	return (ToF>165 && ToF<182);
}

Double_t calculateToF(ROOT::VecOps::RVec<unsigned short> &timeF3, ROOT::VecOps::RVec<unsigned short> &timeF5)
{
	Double_t tF3 = ((timeF3[0]+timeF3[1]+timeF3[2]+timeF3[3])/4.0);
	Double_t tF5 = ((timeF5[0]+timeF5[1])/2.0);
	Double_t ToF = (tF5-tF3)*tdcBinning+ToFconstant;
	return ToF;
}

ROOT::RDF::RNode ApplyGraphicalCuts(ROOT::RDF::RNode df, Int_t cutsType)
{
	if (cutsType==elastic)
	{
		TFile cutgFile("/home/zalewski/aku_dec/analysis/gcuts.root","READ");
		GCutHe4 = (TCutG*)cutgFile.Get("dehe4");
		GCutHe6 = (TCutG*)cutgFile.Get("dehe6");
		GCutP = (TCutG*)cutgFile.Get("p");
		GCutD = (TCutG*)cutgFile.Get("d");
		GCutT = (TCutG*)cutgFile.Get("t");
		GCutangAngPP = (TCutG*)cutgFile.Get("angAngPP");
		GCutangAngDD = (TCutG*)cutgFile.Get("angAngDD");
		GCutlangEPP = (TCutG*)cutgFile.Get("langEPP");
		GCutlangEDD = (TCutG*)cutgFile.Get("langEDD");
		GCtimeCutR = (TCutG*)cutgFile.Get("timeCutR");
		GCtimeCutL = (TCutG*)cutgFile.Get("timeCutL");

		TFile dtCutFile("/home/zalewski/aku_dec/analysis/dtCutFile.root","READ");
		GCutHe4 = (TCutG*)dtCutFile.Get("dehe4");
		GCutdtLeneLang = (TCutG*)dtCutFile.Get("dtLeneLang");
		GCutdtLeneRang = (TCutG*)dtCutFile.Get("dtLeneRang");		
		GCnoProtDeut = (TCutG*)dtCutFile.Get("noProtDeut");
		GCutdtEneEne = (TCutG*)dtCutFile.Get("dtEneEne");
		GCutdtAngAng = (TCutG*)dtCutFile.Get("dtAngAng");
		GCbHe6 = (TCutG*)dtCutFile.Get("bHe6");
		GCdtTT = (TCutG*)dtCutFile.Get("dtTT");
		GCdtLTl = (TCutG*)dtCutFile.Get("dtLTl");
		GCdtLTr = (TCutG*)dtCutFile.Get("dtLTr");
		GCdtRTl = (TCutG*)dtCutFile.Get("dtRTl");
		GCdtRTr = (TCutG*)dtCutFile.Get("dtRTr");
		GCdts1 = (TCutG*)dtCutFile.Get("dts1");
		GCdts2 = (TCutG*)dtCutFile.Get("dts2");

		return df.Define("he4", [](Double_t sqretot, Double_t sqrde){return GCutHe4->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"})
				 .Define("he6", [](Double_t sqretot, Double_t sqrde){return GCutHe6->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"})
				 .Define("p", [](Double_t sqletot, Double_t sqlde){return GCutP->IsInside(sqletot,sqlde);}, {"sqletot","sqlde"})
				 .Define("d", [](Double_t sqletot, Double_t sqlde){return GCutD->IsInside(sqletot,sqlde);}, {"sqletot","sqlde"})
				 .Define("t", [](Double_t sqletot, Double_t sqlde){return GCutT->IsInside(sqletot,sqlde);}, {"sqletot","sqlde"})
				 .Define("angAngPP", [](Double_t sqrang, Double_t sqlang){return GCutangAngPP->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
				 .Define("angAngDD", [](Double_t sqrang, Double_t sqlang){return GCutangAngDD->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
				 .Define("langEPP", [](Double_t sqletot, Double_t sqlde, Double_t sqlang){return GCutlangEPP->IsInside(sqlang, sqletot+sqlde);}, {"sqletot","sqlde", "sqlang"})
				 .Define("langEDD", [](Double_t sqletot, Double_t sqlde, Double_t sqlang){return GCutlangEDD->IsInside(sqlang, sqletot+sqlde);}, {"sqletot","sqlde", "sqlang"})
				 .Define("timeCutR", [](Int_t sqrtime, ROOT::RVec<unsigned short> &tdcF5){return GCtimeCutR->IsInside(sqrtime, tdcF5[0]-sqrtime);}, {"sqrtime", "tdcF5"})
				 .Define("timeCutL", [](Int_t sqltime, ROOT::RVec<unsigned short> &tdcF5){return GCtimeCutL->IsInside(sqltime, tdcF5[0]-sqltime);}, {"sqltime", "tdcF5"})
				 .Define("pp", "he6 && angAngPP && langEPP && timeCutL && timeCutR")
				 .Define("dd", "he6 && angAngDD && langEDD && timeCutL && timeCutR")
				 .Define("sumCut", "pp||dd");
	}
	
	else if (cutsType==MC)
	{
		TFile MCcutgFile("/home/zalewski/aku_dec/analysis/mcCuts.root","READ");
		GCutmcHe6 = (TCutG*)MCcutgFile.Get("mcHe6");
		GCutmcPPAngAng = (TCutG*)MCcutgFile.Get("mcPPAngAng");
		GCutmcPPLeneLang = (TCutG*)MCcutgFile.Get("mcPPLeneLang");
		GCutmcDDAngAng = (TCutG*)MCcutgFile.Get("mcDDAngAng");
		GCutmcDDLeneLang = (TCutG*)MCcutgFile.Get("mcDDLeneLang");
		return df.Define("mcHe6", [](Double_t sqrang, Double_t sqlang){return GCutmcHe6->IsInside(sqrang,sqlang);}, {"sqretot","sqrde"})
				 .Define("mcPPAngAng", [](Double_t sqrang, Double_t sqlang){return GCutmcPPAngAng->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
				 .Define("mcPPLeneLang", [](Double_t sqlde, Double_t sqletot, Double_t sqlang){return GCutmcPPLeneLang->IsInside(sqlang, sqlde+sqletot);}, {"sqlde","sqletot", "sqlang"})
				 .Define("mcPP", "mcPPAngAng && mcPPLeneLang && mcHe6")
				 .Define("mcDDAngAng", [](Double_t sqrang, Double_t sqlang){return GCutmcDDAngAng->IsInside(sqrang,sqlang);}, {"sqrang","sqlang"})
				 .Define("mcDDLeneLang", [](Double_t sqlde, Double_t sqletot, Double_t sqlang){return GCutmcDDLeneLang->IsInside(sqlang, sqlde+sqletot);}, {"sqlde","sqletot", "sqlang"})
				 .Define("mcDD", "mcDDAngAng && mcDDLeneLang && mcHe6")
				 .Define("sumCut", "mcPP || mcDD");
	}

	else if (cutsType==dt)
	{
		TFile dtCutFile("/home/zalewski/aku_dec/analysis/dtCutFile.root","READ");
		GCutHe4 = (TCutG*)dtCutFile.Get("dehe4");
		GCutdtLeneLang = (TCutG*)dtCutFile.Get("dtLeneLang");
		GCutdtLeneRang = (TCutG*)dtCutFile.Get("dtLeneRang");		
		GCnoProtDeut = (TCutG*)dtCutFile.Get("noProtDeut");
		GCutdtEneEne = (TCutG*)dtCutFile.Get("dtEneEne");
		GCutdtAngAng = (TCutG*)dtCutFile.Get("dtAngAng");
		GCbHe6 = (TCutG*)dtCutFile.Get("bHe6");
		GCdtTT = (TCutG*)dtCutFile.Get("dtTT");
		GCdtLTl = (TCutG*)dtCutFile.Get("dtLTl");
		GCdtLTr = (TCutG*)dtCutFile.Get("dtLTr");
		GCdtRTl = (TCutG*)dtCutFile.Get("dtRTl");
		GCdtRTr = (TCutG*)dtCutFile.Get("dtRTr");
		GCdts1 = (TCutG*)dtCutFile.Get("dts1");
		GCdts2 = (TCutG*)dtCutFile.Get("dts2");
		GCdtmm3H_N = (TCutG*)dtCutFile.Get("dt_mm3H_mmN");
		GCdtmm4He_N = (TCutG*)dtCutFile.Get("dt_mm4He_mmN");

		return df.Define("he4", [](Double_t sqretot, Double_t sqrde){return GCutHe4->IsInside(sqretot,sqrde);}, {"sqretot","sqrde"})
				 .Define("dtLeneLang", [](Double_t sqlang, Double_t sqlde, Double_t sqletot){return GCutdtLeneLang->IsInside(sqlang, sqlde+sqletot);}, {"sqlang", "sqlde", "sqletot"})
				 .Define("dtLeneRang", [](Double_t sqrang, Double_t sqlde, Double_t sqletot){return GCutdtLeneRang->IsInside(sqrang, sqlde+sqletot);}, {"sqrang", "sqlde", "sqletot"})
				 .Define("dtAngAng", [](Double_t sqrang, Double_t sqlang){return GCutdtAngAng->IsInside(sqrang, sqlang);}, {"sqrang", "sqlang"})
				 .Define("dtEneEne", [](Double_t sqlde, Double_t sqletot, Double_t sqrde, Double_t sqretot){return GCutdtEneEne->IsInside(sqlde+sqletot, sqrde+sqretot);}, 
				 																													{"sqlde", "sqletot", "sqrde", "sqretot"})
				 .Define("noProtDeut", [](Double_t sqlde, Double_t sqletot){return (!GCnoProtDeut->IsInside(sqletot, sqlde));}, {"sqlde", "sqletot"})
				 .Define("dtTT", [](int sqrtime, int sqltime){return GCdtTT->IsInside(sqrtime, sqltime);}, {"sqrtime", "sqltime"})
				 .Define("bHe6", [](Double_t tof, Double_t aF5){return (GCbHe6->IsInside(tof, aF5) ||  (aF5==0 && tof>168 && tof<181));}, {"tof", "aF5"})
				 .Define("dtLTl", [](ROOT::RVec<unsigned short> &tdcF5, int sqltime){return GCdtLTl->IsInside(tdcF5[0]-sqltime, sqltime);}, {"tdcF5", "sqltime"})
				 .Define("dtLTr", [](ROOT::RVec<unsigned short> &tdcF5, int sqltime){return GCdtLTr->IsInside(tdcF5[0]-sqltime, sqltime);}, {"tdcF5", "sqltime"})
				 .Define("dtRTl", [](ROOT::RVec<unsigned short> &tdcF5, int sqrtime){return GCdtRTl->IsInside(tdcF5[0]-sqrtime, sqrtime);}, {"tdcF5", "sqrtime"})
				 .Define("dtRTr", [](ROOT::RVec<unsigned short> &tdcF5, int sqrtime){return GCdtRTr->IsInside(tdcF5[0]-sqrtime, sqrtime);}, {"tdcF5", "sqrtime"})
				 .Define("timeCut","(dtLTl || dtLTr) && (dtRTl || dtRTr)")
				 .Define("dts2", [](Double_t sqlang, Double_t sqlde){return GCdts2->IsInside(sqlang, sqlde);}, {"sqlang", "sqlde"})
				 .Define("dtmm3H_N", [](Double_t mm3H, Double_t mmN){return GCdtmm3H_N->IsInside(mm3H, mmN);}, {"mm3H", "mmN"})
				 .Define("dtmm4He_N", [](Double_t mm4He, Double_t mmN){return GCdtmm4He_N->IsInside(mm4He, mmN);}, {"mm4He", "mmN"})
				 .Define("shortcut", "he4 && noProtDeut && dtTT && bHe6 && timeCut && sqlde<23")
				 .Define("sumCut", "he4 && noProtDeut && dtTT && timeCut && bHe6 && sqlde<23 && dtmm3H_N && dtmm4He_N && dtAngAng");
	}

	else if (cutsType==mcDT)
	{
		TFile dtCutFile("/home/zalewski/aku_dec/analysis/dtCutFile.root","READ");
		GCmcHe4 = (TCutG*)dtCutFile.Get("mcHe4");
		GCmcLeneLang = (TCutG*)dtCutFile.Get("mcLeneLang");
		GCutdtAngAng = (TCutG*)dtCutFile.Get("dtAngAng");

		return df.Define("mcHe4", [](Double_t sqretot, Double_t sqrde){return GCmcHe4->IsInside(sqretot, sqrde);}, {"sqretot", "sqrde"})
				 .Define("dtAngAng", [](Double_t sqrang, Double_t sqlang){return GCutdtAngAng->IsInside(sqrang, sqlang);}, {"sqrang", "sqlang"})
				 .Define("mcLeneLang", [](Double_t sqlang, Double_t sqlde, Double_t sqletot){return GCmcLeneLang->IsInside(sqlang, sqlde+sqletot);}, {"sqlang", "sqlde", "sqletot"});
	}
	
	else return df;
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

void loadCorrectionParameters(Int_t m_runNo)
{
	parameters.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/aligning_v18125/chosen3.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("Failed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	float tmpContainer;
	for (int iii = 0; iii < 12; iii++)
	{
		outStreamGenerated>>tmpContainer;
		if (iii!=0)
		{
			parameters.push_back(tmpContainer);
			printf("%s = %f\t", parNames[iii-1].c_str(), parameters[iii-1]);
		}
	}
	std::cout<<std::endl;
}

void loadCalibrationParameters()
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

		std::ifstream instream(fileName.Data());

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
						//.Filter(filterToF,{"tdcF3", "tdcF5"})
						//.Filter([sqlThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_L){return filterSQ(SQX_L, sqlThreshold);},{"SQX_L"})
						//.Filter([sqrThreshold](ROOT::VecOps::RVec<unsigned short> &SQX_R){return filterSQ(SQX_R, sqrThreshold);},{"SQX_R"})
						.Define("tof", calculateToF, {"tdcF3", "tdcF5"});

	outDF.Snapshot("cleaned", "/home/zalewski/Desktop/6He/thesis/cal_geo1wToF.root"/*outFilename.Data()*/);
}

struct Calibrator
{
	int columnID;
	//create object containing ID of the column which is to be calibrated
	Calibrator(int columnToCalibrate) 
	{
		columnID = columnToCalibrate;
	};

	ROOT::RVecD operator()(ROOT::VecOps::RVec<unsigned short> &inputArray)
	{		
		auto mySize = inputArray.size();
		ROOT::RVecD outputArray(mySize);
		for (size_t iii = 0; iii < mySize; iii++)
		{
			//std::cout<<vecAcoef.at(columnID).at(iii)<<" "<<"My ID is: "<<vecCalibratedColumns.at(columnID)<<std::endl;
			outputArray.at(iii) = 0.0;
			if (inputArray.at(iii)>0.0)
			{
				outputArray.at(iii) = vecAcoef.at(columnID).at(iii) + vecBcoef.at(columnID).at(iii) * inputArray.at(iii);
			}	
			
		}

		return outputArray;
	}
};

ROOT::RDF::RNode ApplyDefines(ROOT::RDF::RNode df, const std::vector<std::string> &colNames, int iii = 0)
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
				 .Define("geo", "cs::runNo");
	}

	std::string inputColumn = colNames[iii];
	inputColumn.insert(0, "cal_");
	//std::cout<<colNames[iii]<<"\t"<<inputColumn<<std::endl;
	//calling 'Calibrator' so I can pass other info
	return ApplyDefines(df.Define(inputColumn, Calibrator(iii), {colNames[iii]}), colNames, iii + 1);
}

void calibratorCaller(TString inFileName)
{	
	inFileName.ReplaceAll("simp","cln");
	std::cout<<"Calibrating "<<inFileName<<std::endl;
	ROOT::RDataFrame inDF("cleaned", inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cln","cal");
	// Print columns' names
	auto dfWithDefines = ApplyDefines(inDF, vecCalibratedColumns);
	dfWithDefines.Snapshot("calibrated", outFilename.Data(), dfWithDefines.GetColumnNames());
}

void analysis(TString inFileName, Int_t geoID=0)
{
	inFileName.ReplaceAll("simp","cal");
	TString treeName = "calibrated";
	ROOT::RDataFrame inDF(treeName.Data(), inFileName.Data());
	TString outFilename = inFileName.ReplaceAll("cal","dE").ReplaceAll("mc","mc_out");
	std::cout<<"Analysing "<<inFileName<<std::endl;

	h1_Si = new ELC(1, 1, si_Nel, 2.33, si_A, si_Z, si_W, 200.,3000);
	h2_Si = new ELC(2, 1, si_Nel, 2.33, si_A, si_Z, si_W, 200.,3000);
	h3_Si = new ELC(3, 1, si_Nel, 2.33, si_A, si_Z, si_W, 200.,3000);
	he4_Si = new ELC(4, 2, si_Nel, 2.33, si_A, si_Z, si_W, 200.,3000);
	he6_Si = new ELC(6, 2, si_Nel, 2.33, si_A, si_Z, si_W, 300.,3000);
	/*
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
	siEloss6He.SetDeltaEtab(1500);
*/
	Double_t qDT_45 = cs::Qdt;
	Double_t qPT_75 = cs::Qpt;

	parameters[sMWPC_1_X] = -1.0;
	parameters[sMWPC_1_Y] = -2.1375;
	parameters[sMWPC_2_X] = 0.2; 
	parameters[sMWPC_2_Y] = -1.125;
	parameters[sTarPos] = 10.0;
	parameters[sLang1] = 0.3;
	parameters[sLang2] = 0.0;
	parameters[sLang3] = 0.0;
	parameters[sRang] = 0.0;
	parameters[sDistL] = -20.0;
	parameters[sDistR] = -30.0;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos]-2.0;

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + parameters[sDistL];
	sqrDist = 250.0 + parameters[sDistR];

	inDF.Filter(filterLeftDetector, {"cal_SQX_L"})
		.Filter(filterLeftDetector, {"cal_SQY_L"})
		.Filter(filterRightDetector, {"cal_SQX_R"})
		.Filter(filterRightDetector, {"cal_SQY_R"});

	auto outDF = inDF//.Range(100)
	.Define("kinE", getKineticEnergy, {"tof"})
					 //.Define("elo", getEloss, {"tof"})
					 .Define("X1",[](Double_t MWPC_1_X){return (MWPC_1_X + parameters[sMWPC_1_X]);}, {"MWPC_1_X"})
					 .Define("Y1",[](Double_t MWPC_1_Y){return (MWPC_1_Y + parameters[sMWPC_1_Y]);}, {"MWPC_1_Y"})
					 .Define("Z1",[](){return -816.0;})
					 .Define("X2",[](Double_t MWPC_2_X){return (MWPC_2_X + parameters[sMWPC_2_X]);}, {"MWPC_2_X"})
					 .Define("Y2",[](Double_t MWPC_2_Y){return (MWPC_2_Y + parameters[sMWPC_2_Y]);}, {"MWPC_2_Y"})
					 .Define("Z2",[](){return -270.0;})
					 .Define("MWPC", getMWPC, {"X1", "Y1", "Z1", "X2", "Y2", "Z2"})
					 .Define("vBeam", getBeamVector, {"MWPC", "kinE"})
					 .Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
					 .Define("tarVertex", [](ROOT::RVecD rvecMWPC, Int_t geo){return getTarVertex(rvecMWPC, geo);}, {"MWPC", "geo"})
					 .Define("SQ300_strip", getStripNumber, {"cal_SQ300"})
					 .Define("SQX_L_strip", getStripNumber, {"cal_SQX_L"})
					 .Define("SQY_L_strip", getStripNumber, {"cal_SQY_L"})
					 .Define("SQX_R_strip", getStripNumber, {"cal_SQX_R"})
					 .Define("SQY_R_strip", getStripNumber, {"cal_SQY_R"})
					 .Define("CsI_L_strip", getStripNumber, {"cal_CsI_L"})
					 .Define("CsI_R_strip", getStripNumber, {"cal_CsI_R"})
					 .Define("sq300", [](ROOT::RVecD &SQ300, Int_t stripNumber){return SQ300[stripNumber];}, {"cal_SQ300", "SQ300_strip"})
					 .Define("sqlde", getSQlde, {"cal_SQX_L"})
					 .Define("fsqlde1", [](Double_t m_sqlde){return h1_Si->GetE(m_sqlde, -(ionDeadLayer + 3.0));}, {"sqlde"})
					 .Define("fsqlde2", [](Double_t m_sqlde){return h2_Si->GetE(m_sqlde, -(ionDeadLayer + 3.0));}, {"sqlde"})
					 .Define("fsqlde3", [](Double_t m_sqlde){return h3_Si->GetE(m_sqlde, -(ionDeadLayer + 3.0));}, {"sqlde"})
					 .Define("sqrde", [](ROOT::RVecD &SQX_R, Int_t stripNumber){return SQX_R[stripNumber];}, {"cal_SQX_R", "SQX_R_strip"})
					 .Define("sqletot", [](ROOT::RVecD &CsI_L, Int_t stripNumber){return CsI_L[stripNumber];}, {"cal_CsI_L", "CsI_L_strip"})
					 .Define("sqretot", [](ROOT::RVecD &CsI_R, Int_t stripNumber){return CsI_R[stripNumber];}, {"cal_CsI_R", "CsI_R_strip"})
					 //.Define("sqletot2", calibrateCsIL_2H, {"CsI_L", "CsI_L_strip"})

					 .Define("sqltime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_L, Int_t stripNumber){return (int)tSQX_L[stripNumber];}, {"tSQX_L", "SQX_L_strip"})
					 .Define("sqrtime", [](ROOT::VecOps::RVec<unsigned short> &tSQX_R, Int_t stripNumber){return (int)tSQX_R[stripNumber];}, {"tSQX_R", "SQX_R_strip"})
					 .Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
					 .Define("leftLabVertex", getLeftDetPosition, {"geo"})
					 .Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
					 .Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip", "geo"})
					 .Define("rightLabVertex", getRightDetPosition, {"geo"})
					 .Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})
					 .Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
					 .Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
					 .Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
					 .Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})					 
					 .Define("lv1H", getLV1H, {"sqlde", "sqletot", "v2H"})
					 .Define("lv2H", getLV2H, {"sqlde", "sqletot", "v2H"})
					 .Define("lv3H", getLV3H, {"sqlde", "sqletot", "v2H"})
					 .Define("lv4He", getLV4He, {"sqrde", "sqretot", "v6He"})
					 .Define("lv5He", getLV5He, {"sqrde", "sqretot", "v6He"})
					 .Define("mm1H", getMissingMass1H, {"lv1H", "lvBeam"})
					 .Define("mm2H", getMissingMass2H, {"lv2H", "lvBeam"})
					 .Define("mm3H", getMissingMass3H, {"lv3H", "lvBeam"})
					 .Define("mm4He", getMissingMass4He, {"lv4He", "lvBeam"})
					 .Define("mmN", getMissingMassN, {"lv4He", "lv3H", "lvBeam"})
					 .Define("mc1H", getCMAngle1H, {"sqlang", "lvBeam"})
					 .Define("mc2H", getCMAngle2H, {"sqlang", "lvBeam"})
					 .Define("dtCM", getCMAngleDT, {"sqrang", "lvBeam"})
				 	 .Define("dtCMene", getCMAngleDTene, {"sqlde", "sqletot", "lvBeam"})
					 //kinematics generation
					 .Define("thetaCM", "rnd->Uniform(0.0,TMath::Pi())")
					 .Define("inangl", inelSQLang2H_18, {"thetaCM","lvBeam"})
					 .Define("inangr", inelSQRang2H_18, {"thetaCM","lvBeam"})
					 .Define("sqlangpp", getPPLang, {"thetaCM","lvBeam"})
					 .Define("sqrangpp", getPPRang, {"thetaCM","lvBeam"})
					 .Define("sqlepp", getEnePP, {"thetaCM","lvBeam"})
					 .Define("sqlangdd", getDDLang, {"thetaCM","lvBeam"})
					 .Define("sqrangdd", getDDRang, {"thetaCM","lvBeam"})
					 .Define("sqledd", getEneDD, {"thetaCM","lvBeam"})
					 .Define("sqlangdp", getDPLang, {"thetaCM","lvBeam"})
					 .Define("sqrangdp", getDPRang, {"thetaCM","lvBeam"})

					 .Define("sqlangdt", [qDT_45](Double_t thetaCM, TLorentzVector lvBeam){return getDTLang(thetaCM, lvBeam, qDT_45);}, {"thetaCM","lvBeam"})
					 .Define("sqrangdt", [qDT_45](Double_t thetaCM, TLorentzVector lvBeam){return getDTRang(thetaCM, lvBeam, qDT_45);}, {"thetaCM","lvBeam"})
					 .Define("sqlangpt", [qPT_75](Double_t thetaCM, TLorentzVector lvBeam){return getPTLang(thetaCM, lvBeam, qPT_75);}, {"thetaCM","lvBeam"})
					 .Define("sqrangpt", [qPT_75](Double_t thetaCM, TLorentzVector lvBeam){return getPTRang(thetaCM, lvBeam, qPT_75);}, {"thetaCM","lvBeam"});

					auto filteredDF = ApplyGraphicalCuts(outDF, dt);

	ROOT::RDF::RSnapshotOptions myOpts;
/*
	//small output file pp && dd
	TString smallRealFname = TString::Format("/home/zalewski/dataTmp/small/geo%d/small%d.root", cs::runNo, cs::runNo);
	filteredDF.Snapshot("smallReal", smallRealFname.Data(), columnList, myOpts);

*/
/*
	//small output file dt
	TString smallRealFnameDT = "/home/zalewski/dataTmp/small/small_dt.root";
	filteredDF.Snapshot("smallRealDT", smallRealFnameDT.Data(), columnListDT);
*/
/*
	//small elastic MC output file
	TString smallMCFname = outFilename;
	filteredDF.Snapshot("smallMC", smallMCFname.ReplaceAll("MC","small").Data(), MCcolumnListLong);
*/

//	//small dt MC output file
//	TString smallMCFname = outFilename;
//	filteredDF.Snapshot("smallMC", smallMCFname.ReplaceAll("MC","small").Data()/*, MCcolumnListDT*/);


	//full output file
	filteredDF.Snapshot("analyzed", outFilename.Data());

	delete h1_Si;
	delete h2_Si;
	delete h3_Si;
	delete he4_Si;
	delete he6_Si;
}

void miscTranlator()
{
	TString str_name;
	TString sourceDir = "/home/zalewski/dataTmp/misc/raw/";

	TSystemDirectory *dir_data = new TSystemDirectory("data",sourceDir.Data());

	TIter bluster = dir_data->GetListOfFiles();
	while (TObject *obj = bluster())
	{
		str_name = obj->GetName();
		//std::cout<<str_name<<std::endl;

		if ((str_name.Contains("root") && !str_name.Contains("li9")) ||
		(str_name.Contains("mc") && !str_name.Contains("out")))
		{
			std::cout<<str_name<<std::endl;
			TString inputFilePath = sourceDir + str_name;
			TString versionNumber(str_name(7,1));
			//printf("fName:\t%s\tversionNo: %d\n",str_name.Data(), versionNumber.Atoi());
			//translator(inputFilePath);
			//cleaner(inputFilePath);
			//calibratorCaller(inputFilePath);
			//loadCorrectionParameters(1);
			//analysis(inputFilePath, 1);
		}
	}

}

void ui()
{
	ROOT::EnableImplicitMT();
	TStopwatch stopwatch;
	
	beamDeadLayer = 500.0;
	ionDeadLayer = 20.0;
	if (cs::runNo>1)
	{
		ionDeadLayer*=2.0;
	}
	

	rnd = new TRandom3();
	selector = TSelector::GetSelector("/home/zalewski/aku_dec/analysis/simplifier123.C");
	std::cout<<"..."<<std::endl<<".."<<std::endl<<"."<<std::endl;


	if (cs::runNo==5)
	{
		selector = TSelector::GetSelector("/home/zalewski/aku_dec/analysis/simplifier5.C");
	}

	loadCalibrationParameters();
	//loadCorrectionParameters();
	TString str_name, sourceDir;

	if (cs::runNo==1 || cs::runNo==2 || cs::runNo==3)
	{
		sourceDir = "/home/zalewski/dataTmp/simp/geo" + std::to_string(cs::runNo) + "/";
	}

	else
	{
		printf("Wrong geometry (geo = %d)\tExiting now.\n", cs::runNo);

	}

	TSystemDirectory *dir_data = new TSystemDirectory("data",sourceDir.Data());

	
	TIter bluster = dir_data->GetListOfFiles();
	while (TObject *obj = bluster())
	{
		str_name = obj->GetName();
		//std::cout<<str_name<<std::endl;
		if (str_name.Contains("run00"))
		{
			ToFconstant = cs::tof_const + 7.7;
		}
		

		if ((str_name.Contains("simp_geo") && !str_name.Contains("li9")) ||
		(str_name.Contains("mc") && !str_name.Contains("out")))
		{
			TString inputFilePath = sourceDir + str_name;
			TString versionNumber(str_name(7,1));
			//printf("fName:\t%s\tversionNo: %d\n",str_name.Data(), versionNumber.Atoi());
			//translator(inputFilePath);
			cleaner(inputFilePath);
			calibratorCaller(inputFilePath);
			analysis(inputFilePath, 0);
			//makeSmallFile();
		}
	}

	//loadCorrectionParameters(6);

	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo1_1H.root", 1);
	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo2_1H.root", 2);
	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo3_1H.root", 3);
	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo1_2H.root", 1);
	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo2_2H.root", 2);
	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo3_2H.root", 3);
	

	//analysis("/home/zalewski/dataTmp/MC/v3/mc_geo2_2H_dt.root", 2);

	//joinDE();
	//miscTranlator();
	stopwatch.Print();
}