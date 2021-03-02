#ifndef ui_hh
#define ui_hh 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "TLorentzVector.h"
#include "TSelector.h"
#include "TCutG.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "Math/Vector4D.h"
#include "TSystemDirectory.h"

#include "/home/zalewski/aku/ELC/AELC.h"
#include "/home/zalewski/aku/ELC/ELC.h"
#include "/home/zalewski/aku/TELoss/TELoss.h"

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

	std::vector< std::vector< float > > vecAcoef;
	std::vector< std::vector< float > > vecBcoef;
	std::vector< float > Acoef;
	std::vector< float > Bcoef;

//variables used by cleaner
	Double_t ToFconstant = cs::tof_const;
	Double_t tdcBinning = 0.125;

//variables used by analysis
	TRandom3 *rnd;
	const TLorentzVector lvTar1H(0.0, 0.0, 0.0, cs::mass1H);
	const TLorentzVector lvTar2H(0.0, 0.0, 0.0, cs::mass2H);
	TCutG *GCutHe4, *GCutHe6, *GCutdE2H, *GCutdE3H, *GCutangAngPT, *GCutlangEPT, *GCutrangEPT, *GCutangAngPP, *GCutangAngDD, *GCutlangEPP, *GCutlangEDD;
	TELoss siEloss1H, siEloss2H, siEloss3H, siEloss4He, siEloss6He;
	
	//TVector3

	Double_t sqlAng, sqrAng;
	Double_t sqlDist,sqrDist;
	Double_t tarPos, tarThickness, tarAngle;
	int si_Nel=1;
	double si_A[1] = {28};
	double si_Z[1] = {14};
	double si_W[1] = {1};

#endif