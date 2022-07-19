#ifndef ui_hh
#define ui_hh 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "TLorentzVector.h"
#include "TSelector.h"
#include "TCutG.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "Math/Vector4D.h"
#include "TSystemDirectory.h"

#include "/home/zalewski/aku/ELC/AELC.h"
#include "/home/zalewski/aku/ELC/ELC.h"
#include "csi_L_2.hh"
//#include "/home/zalewski/aku/TELoss/TELoss.h"

//global variables used by couple of functions
	Int_t numberOfThreads;
	TSelector *selector;	
	TString raw_data_path = "/home/zalewski/dataTmp/raw/geo" + std::to_string(cs::runNo) + "/";

	std::vector<std::string> vecCalibratedColumns{"SQX_L",
												  "SQY_L",
												  "CsI_L",
												  "SQX_R",
												  "SQY_R",
												  "CsI_R",
												  "tdcF3",
												  "tdcF5",
												  "tdcF6",
												  "tdcMWPC",
												  "SQ300"};

	std::vector<std::string> columnList{"SQX_L_strip",
										"SQX_R_strip",
										"SQY_L_strip",
										"SQY_R_strip",
										"sqlde",
										"sqletot",
										"sqlang",
										"sqrde",
										"sqretot",
										"sqrang",
										"geo",
										"MWPC_1_X",
										"MWPC_2_X",
										"MWPC_1_Y",
										"MWPC_2_Y",
										"kinE"};//,
										//"pp",
										//"dd"};

	std::vector<std::string> MCcolumnList{"geo",
										"fZ2H",
										"mcPP",
										"mcDD",
										"fevx",
										"fevy",
										"fThetaCM"};

std::vector<std::string> smallerTree{"kinE",
									 "elo",
									 "vBeam",
									 "lvBeam",
									 "tarVertex",
									 "sq300",
									 "sqlde",
									 "resqlde",
									 "sqrde",
									 "sqletot",
									 "sqretot",
									 "sqltime",
									 "sqrtime",
									 "v2H",
									 "v6He",
									 "sqlang",
									 "sqrang",
									 "mc1H",
									 "mc2H",
									 "he4",
									 "he6",
									 "p",
									 "d",
									 "t",
									 "pp",
									 "dd",
									 "pt",
									 "partPT",
									 "thetaCM",
									 "sqlangpp",
									 "sqrangpp",
									 "sqlepp",
									 "sqlangdd",
									 "sqrangdd",
									 "sqledd",
									 "sqlangdt_45",
									 "sqrangdt_45",
									 "sqlangdt_0",
									 "sqrangdt_0",
									 "sqlangpt_75",
									 "sqrangpt_75",
									 "sqlangpt_0",
									 "sqrangpt_0",
									 "mlv1H",
									 "mlv2H",
									 "lv2H.",
									 "mm1H",
									 "mm2H",
									 "remm1H",
									 "remm2H",
									 //"lv6He",
									 "geo"};

	enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sLang1, sLang2, sLang3, sRang, sTarPos, sDistL, sDistR};
	std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sLang1", "sLang2", "sLang3", "sRang", "sTarPos", "sDistL", "sDistR"};
	std::vector<Double_t> parameters(11);

	std::vector< std::vector< float > > vecAcoef;
	std::vector< std::vector< float > > vecBcoef;
	std::vector< float > Acoef;
	std::vector< float > Bcoef;

//variables used by cleaner
	Double_t ToFconstant = cs::tof_const;
	Double_t tdcBinning = 0.125;
	Double_t detDist = 250.0;
	

//variables used by analysis
	TRandom3 *rnd;
	const TLorentzVector lvTar1H(0.0, 0.0, 0.0, cs::mass1H);
	const TLorentzVector lvTar2H(0.0, 0.0, 0.0, cs::mass2H);

	//general cuts
	TCutG *GCutHe4, *GCutHe6, *GCutP, *GCutD, *GCutT;

	//elastic cuts
	TCutG *GCutangAngPP, *GCutangAngDD, *GCutlangEPP, *GCutlangEDD;
	TCutG *GCtimeCutR, *GCtimeCutL;

	//MC cuts
	TCutG *GCutmcPPAngAng, *GCutmcPPLeneLang, *GCutmcPPReneRang, *GCutmcHe6;
	TCutG *GCutmcDDAngAng, *GCutmcDDLeneLang, *GCutmcDDReneRang;

	//dt cuts
	TCutG *GCdtLT1, *GCdtLT2, *GCdts1, *GCdts2, *GCdts3, *GCdts4;
	TCutG *GCutdtLEneLAng, *GCutdtLEneRAng;
	TCutG *GCnoProtDeut, *GCbHe6, *GCdtTT;


	//TELoss siEloss1H, siEloss2H, siEloss3H, siEloss4He, siEloss6He;
	AELC *h1_Si;
	AELC *h2_Si;
	AELC *h3_Si;
	AELC *he4_Si;
	AELC *he6_Si;
	
	Double_t tarPos[3], leftAngle[3], tarAngle[3];
	Double_t sqlDist, sqrDist;
	Double_t rightAngle;

	Int_t myVal;
	int si_Nel=1;
	double si_A[1] = {28};
	double si_Z[1] = {14};
	double si_W[1] = {1};

	Double_t beamDeadLayer;
	Double_t ionDeadLayer;

	enum reactionType{elastic, MC, dt, all};

#endif