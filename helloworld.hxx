#ifndef helloworld_hxx
#define helloworld_hxx 1

#include "/home/zalewski/aku/analysis/constants.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ROOT/RDataFrame.hxx"
#include "TVector3.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLorentzVector.h"


std::string smallFile = "/home/zalewski/dataTmp/small/small.root";
enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos1, sTarPos2, sTarPos3, sLang1, sLang2, sLang3, sRang, sLeftDetShift, sLeftDetDist, sRightDetShift, sRightDetDist};
Double_t tarPos[3], leftAngle[3], tarAngle[3];
Double_t sqlDist, sqrDist;
Double_t rightAngle;
TRandom3 *rnd;
#endif
