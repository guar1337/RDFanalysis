#ifndef helloworld_hxx
#define helloworld_hxx 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "TRandomGen.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"


std::string smallFile = "/home/zalewski/dataTmp/small/small.root";
enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos1, sTarPos2, sTarPos3, sLang1, sLang2, sLang3, sRang};
std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sTarPos1", "sTarPos2", "sTarPos3", "sLang1", "sLang2", "sLang3", "sRang"};


Double_t tarPos[3], leftAngle[3], tarAngle[3];
Double_t sqlDist, sqrDist;
Double_t rightAngle;
Bool_t calculate;

const Int_t numberOfParameters = 11;
Double_t lowerParamBound[numberOfParameters], upperParamBound[numberOfParameters];
Int_t firstRun, lastRun;

std::string outGeneratedParams = "/home/zalewski/Desktop/6He/analysis/experimental/generated.txt";
std::string outObtainedParams = "/home/zalewski/Desktop/6He/analysis/experimental/obtained.txt";

Double_t myMemo[4] = {1000.0, 1000.0, 1000.0, 1000.0};
#endif
