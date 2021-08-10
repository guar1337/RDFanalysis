#ifndef modifiedAnalysis_hh
#define modifiedAnalysis_hh 1

#include "constants.h"
#include "ROOT/RDataFrame.hxx"

enum sPar {sMWPC_1_X, sMWPC_1_Y, sMWPC_2_X, sMWPC_2_Y, sTarPos, sLang1, sLang2, sLang3, sRang};
std::vector<std::string> parNames = {"sMWPC_1_X", "sMWPC_1_Y", "sMWPC_2_X", "sMWPC_2_Y", "sTarPos", "sLang1", "sLang2", "sLang3", "sRang"};
std::vector<Double_t> parameters;

enum sCuts {mX1, X1, mY1, Y1, mX2, X2, mY2, Y2, mX3, X3, mY3, Y3};
std::vector<std::string> corrNames = {"-X1", "X1", "-Y1", "Y1", "-X2", "X2", "-Y2", "Y2", "-X3", "X3", "-Y3", "Y3"};
std::vector<Int_t> vCut;
Double_t tarPos[3], leftAngle[3], tarAngle[3];
Double_t sqlDist, sqrDist;
Double_t rightAngle;
Bool_t verbosity;


#endif