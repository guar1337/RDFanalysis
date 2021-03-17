// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//	 this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//	 this list of conditions and the following disclaimer in the documentation
//	 and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//	 used to endorse or promote products derived from this software without
//	 specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


#include "helloworld.hxx"


TVector3 getBeamVector(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Double_t m_kinE)
{
	Double_t m_MWPC_1_Z = -816.0;
	Double_t m_MWPC_2_Z = -270.0;
	TVector3 m_beamVector(m_MWPC_2_X - m_MWPC_1_X, m_MWPC_2_Y - m_MWPC_1_Y, m_MWPC_2_Z - m_MWPC_1_Z);
	Double_t m_eneBeam = cs::mass6He + m_kinE;
	Double_t m_momBeam = sqrt(m_eneBeam*m_eneBeam - cs::mass6He*cs::mass6He);
	m_beamVector.SetMag(m_momBeam);
	return m_beamVector;
}

TVector3 getTarVertex(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y, Int_t m_geo, const Double_t *m_pars)
{
	Double_t m_tarPos(0.0);
	Double_t m_tarAngle;

	switch (m_geo)
	{
	case 1:
		m_tarAngle = 45.0*TMath::DegToRad();
		m_tarPos = tarPos[m_geo-1] + m_pars[sTarPos1];
		break;

	case 2:
		m_tarAngle = 6.0*TMath::DegToRad();
		m_tarPos = tarPos[m_geo-1] + m_pars[sTarPos2];
		break;

	case 3:
		m_tarAngle = 0.0*TMath::DegToRad();
		m_tarPos = tarPos[m_geo-1] + m_pars[sTarPos3];
		break;
	
	default:
		break;
	}

	Double_t m_MWPC_1_Z = -816.0;
	Double_t m_MWPC_2_Z = -270.0;
	Double_t m_dX = m_MWPC_2_X - m_MWPC_1_X;
	Double_t m_dY = m_MWPC_2_Y - m_MWPC_1_Y;
	Double_t m_dZ = m_MWPC_2_Z - m_MWPC_1_Z;
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, m_tarPos);
	TVector3 m_beamPoint(m_MWPC_2_X, m_MWPC_2_Y, m_MWPC_2_Z);
	TVector3 m_tarPerpendicular(sin(m_tarAngle), 0.0, cos(m_tarAngle));
	Double_t m_dCoeff = ((m_tarPoint-m_beamPoint).Dot(m_tarPerpendicular))/(m_vBeam.Dot(m_tarPerpendicular));
		
	Double_t m_evX = m_MWPC_2_X + m_dX * m_dCoeff;
	Double_t m_evY = m_MWPC_2_Y + m_dY * m_dCoeff;
	Double_t m_evZ = m_MWPC_2_Z + m_dZ * m_dCoeff;

	return TVector3(m_evX, m_evY, m_evZ);
}

TVector3 getLeftDetVertex(Int_t m_xStrip, Int_t m_yStrip, Int_t m_geo)
{
	m_xStrip+=rnd->Uniform(0.0,1.0)-0.5;
	m_yStrip+=rnd->Uniform(0.0,1.0)-0.5;

	// coordinates of hit in LAB system	
	Double_t X2hDet = -cs::widthStripX * m_xStrip * cos(leftAngle[m_geo-1]);
	Double_t Y2hDet = cs::widthStripY * m_yStrip;
	Double_t Z2hDet = cs::widthStripX * m_xStrip * sin(leftAngle[m_geo-1]);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip)
{
	m_xStrip+=rnd->Uniform(0.0,1.0)-0.5;
	m_yStrip+=rnd->Uniform(0.0,1.0)-0.5;

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * m_xStrip * cos(rightAngle);
	Double_t Y6HeDet = cs::widthStripY * m_yStrip;
	Double_t Z6HeDet = cs::widthStripX * m_xStrip * sin(rightAngle);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

TVector3 getLeftDetPosition(Int_t m_geo, const Double_t *m_pars)
{
	Double_t X2Hlab = sqlDist*sin(leftAngle[m_geo-1]) + (cs::sqlXzero + m_pars[sLeftDetShift]) * cos(leftAngle[m_geo-1]);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(leftAngle[m_geo-1]) - (cs::sqlXzero + m_pars[sLeftDetShift]) * sin(leftAngle[m_geo-1]);
	return TVector3(X2Hlab, Y2Hlab, Z2Hlab);
}

TVector3 getRightDetPosition(const Double_t *m_pars)
{
	Double_t X6Helab = sqrDist*sin(-rightAngle) - (cs::sqrXzero + m_pars[sRightDetShift]) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart + cs::widthStripY;
	Double_t Z6Helab = sqrDist*cos(rightAngle) - (cs::sqrXzero + m_pars[sRightDetShift]) * sin(rightAngle);
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

Double_t recoPT(Double_t m_sqlang, TLorentzVector m_lvBeam)
{
	//generate pt reaction (lv3H actually) with:
	//			theta in range 0:3.14
	//			beam vector
	//			q of the reaction
	Int_t multi = (rnd->Integer(2)==1) ? 1 : -1;
	multi = 1;
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

struct CostFunctor
{
	bool operator()(const double *par, double *residual) const
	{
/*		#3 different datasets:
		#pp, dd, pt
		#same pool of par: 3x target shift, 4x2 MWPC, 1+3x angle, 2x dist, 2x det shift
		create column (obtained-Expected)**2 value for each dataset - i.e. sqrang, esqrang(sqlang, vBeam)
		calcualte Chi2 - Sum("newCol")
*/
		Double_t tarMass1H = cs::mass1H;
		Double_t tarMass2H = cs::mass2H;

		tarPos[0] = 10.0 + par[sTarPos1];
		tarPos[1] = 10.0 + par[sTarPos2];
		tarPos[2] = 10.0 + par[sTarPos3];

		tarAngle[0] = 45.0 * TMath::DegToRad();
		tarAngle[1] = 6.0 * TMath::DegToRad();
		tarAngle[2] = 0.0 * TMath::DegToRad();

		leftAngle[0] = (65.0 + par[sLang1]) * TMath::DegToRad();
		leftAngle[1] = (50.0 + par[sLang2]) * TMath::DegToRad();
		leftAngle[2] = (35.0 + par[sLang3]) * TMath::DegToRad();
		rightAngle = (15.0 + par[sRang]) * TMath::DegToRad();

		sqlDist = 170.0 + par[sLeftDetDist];
		sqrDist = 250.0 + par[sRightDetDist];

		ROOT::RDataFrame smallDF("small", smallFile);

		auto newDF = smallDF.Define("X1",[&par](Double_t MWPC_1_X){return (MWPC_1_X  + par[sMWPC_1_X]);}, {"MWPC_1_X"})
				.Define("Y1",[&par](Double_t MWPC_1_Y){return (MWPC_1_Y  + par[sMWPC_1_Y]);}, {"MWPC_1_Y"})
				.Define("X2",[&par](Double_t MWPC_2_X){return (MWPC_2_X  + par[sMWPC_2_X]);}, {"MWPC_2_X"})
				.Define("Y2",[&par](Double_t MWPC_2_Y){return (MWPC_2_Y  + par[sMWPC_2_Y]);}, {"MWPC_2_Y"})
				.Define("vBeam", getBeamVector, {"X1", "Y1", "X2", "Y2", "kinE"})
				.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
				.Define("tarVertex", [par](Double_t X1, Double_t Y1, Double_t X2, Double_t Y2, Int_t geo){return getTarVertex(X1, Y1, X2, Y2, geo, par);}, {"X1", "Y1", "X2", "Y2", "geo"})

				.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
				.Define("leftLabVertex", [par](Int_t geo){return getLeftDetPosition(geo, par);}, {"geo"})
				.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
				.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
				.Define("rightLabVertex", [par](){return getRightDetPosition(par);})
				.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

				.Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
				.Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
				.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
				.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"});


auto ppDF = newDF.Filter("pp").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
auto ddDF = newDF.Filter("dd").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
auto ptDF = newDF.Filter("pt && sqrang>10").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});

// for pp and dd sqlang describes the ang-ang relation. For the pr reaction I might have to divide the data into two ranges
auto nEventsPP = ppDF.Count().GetValue()/(1000);
auto sumPP = ppDF.Define("resqrang",[tarMass1H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass1H);}, {"sqlang", "lvBeam"})
	 				  .Define("difSqrang","pow(resqrang-sqrang,2)")
					  .Sum<double>("difSqrang");

auto nEventsDD = ddDF.Count().GetValue()/(1000);
auto sumDD = ddDF.Define("resqrang",[tarMass2H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass2H);}, {"sqlang", "lvBeam"})
	 				  .Define("difSqrang","pow(resqrang-sqrang,2)")
					  .Sum<double>("difSqrang");

double nEventsPT = ptDF.Count().GetValue();
auto sumPT = ptDF.Define("resqrang", recoPT, {"sqlang", "lvBeam"})
					  .Define("difSqrang", "pow(resqrang-sqrang,2)");
					  
auto tmp = sumPT.Filter("difSqrang>0").Sum<double>("difSqrang");

printf("chi2PP: %f\tchi2DD: %f\tchi2PT: %f\n", sumPP.GetValue()/nEventsPP, sumDD.GetValue()/nEventsDD, tmp.GetValue()/100.0);


		residual[0] = sumPP.GetValue()/nEventsPP;
		residual[1] = sumDD.GetValue()/nEventsDD;
		residual[2] = tmp.GetValue()/100.0;
		return true;
	}
};

void drawResults(double *finalPars)
{
/*	
	#3 different datasets:
	#pp, dd, pt
	#same pool of par: 3x target shift, 4x2 MWPC, 1+3x angle, 2x dist, 2x det shift
	create column (obtained-Expected)**2 value for each dataset - i.e. sqrang, esqrang(sqlang, vBeam)
	calcualte Chi2 - Sum("newCol")
*/

	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;

	tarPos[0] = 10.0 + finalPars[sTarPos1];
	tarPos[1] = 10.0 + finalPars[sTarPos2];
	tarPos[2] = 10.0 + finalPars[sTarPos3];

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + finalPars[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + finalPars[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + finalPars[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + finalPars[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + finalPars[sLeftDetDist];
	sqrDist = 250.0 + finalPars[sRightDetDist];

	ROOT::RDataFrame smallDF("small", smallFile);

	auto newDF = smallDF.Define("X1",[&finalPars](Double_t MWPC_1_X){return (MWPC_1_X  + finalPars[sMWPC_1_X]);}, {"MWPC_1_X"})
							.Define("Y1",[&finalPars](Double_t MWPC_1_Y){return (MWPC_1_Y  + finalPars[sMWPC_1_Y]);}, {"MWPC_1_Y"})
							.Define("X2",[&finalPars](Double_t MWPC_2_X){return (MWPC_2_X  + finalPars[sMWPC_2_X]);}, {"MWPC_2_X"})
							.Define("Y2",[&finalPars](Double_t MWPC_2_Y){return (MWPC_2_Y  + finalPars[sMWPC_2_Y]);}, {"MWPC_2_Y"})
							.Define("vBeam", getBeamVector, {"X1", "Y1", "X2", "Y2", "kinE"})
							.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
							.Define("tarVertex", [finalPars](Double_t X1, Double_t Y1, Double_t X2, Double_t Y2, Int_t geo){return getTarVertex(X1, Y1, X2, Y2, geo, finalPars);}, {"X1", "Y1", "X2", "Y2", "geo"})

							.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
							.Define("leftLabVertex", [finalPars](Int_t geo){return getLeftDetPosition(geo, finalPars);}, {"geo"})
							.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
							.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip"})
							.Define("rightLabVertex", [finalPars](){return getRightDetPosition(finalPars);})
							.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

							.Define("v2H", [](TVector3 leftGlobVertex, TVector3 tarVertex){return TVector3(leftGlobVertex-tarVertex);}, {"leftGlobVertex", "tarVertex"})
							.Define("v6He", [](TVector3 rightGlobVertex, TVector3 tarVertex){return TVector3(rightGlobVertex-tarVertex);}, {"rightGlobVertex", "tarVertex"})
							.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
							.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"});


	auto ppDF = newDF.Filter("pp").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
	auto ddDF = newDF.Filter("dd").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});
	auto ptDF = newDF.Filter("pt").Cache<Double_t, Double_t, TLorentzVector>({"sqlang", "sqrang", "lvBeam"});

	auto nEventsPP = ppDF.Count().GetValue();
	auto sumPP = ppDF.Define("resqrang",[tarMass1H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass1H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");
						
	auto nEventsDD = ddDF.Count().GetValue();
	auto sumDD = ddDF.Define("resqrang",[tarMass2H](Double_t sqlang, TLorentzVector lvBeam){return myAngAngFit(sqlang, lvBeam, tarMass2H);}, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");

	auto nEventsPT = ptDF.Count().GetValue();
	auto sumPT = ptDF.Define("resqrang", recoPT, {"sqlang", "lvBeam"})
						  .Define("difSqrang","resqrang-sqrang")
						  .Define("difSqrang2", "pow(resqrang-sqrang,2)");

	TGraph ppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	TGraph ptAngAng = sumPT.Graph<Double_t, Double_t>("sqlang", "sqrang").GetValue();
	
	TGraph fppAngAng = sumPP.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	TGraph fddAngAng = sumDD.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();
	TGraph fptAngAng = sumPT.Graph<Double_t, Double_t>("sqlang", "resqrang").GetValue();

	TGraph difPPSqrang = sumPP.Graph<Double_t, Double_t>("sqlang", "difSqrang").GetValue();
	TGraph difDDSqrang = sumDD.Graph<Double_t, Double_t>("sqlang", "difSqrang").GetValue();
	TGraph difPTSqrang = sumPT.Graph<Double_t, Double_t>("sqlang", "difSqrang").GetValue();

	TCanvas myCanvas("myCanvas", "Minimize results", 1200, 800);
	myCanvas.Divide(2,2);

	myCanvas.cd(1);
	ppAngAng.Draw("AP");
	fppAngAng.SetMarkerColor(kRed);
	fppAngAng.Draw("P,same");

	myCanvas.cd(2);
	ddAngAng.Draw("AP");
	fddAngAng.SetMarkerColor(kRed);
	fddAngAng.Draw("P,same");

	myCanvas.cd(3);
	ptAngAng.Draw("AP");
	fptAngAng.SetMarkerColor(kRed);
	fptAngAng.Draw("P,same");

	myCanvas.cd(4);
	difDDSqrang.Draw("AP");

	myCanvas.Print("/home/zalewski/Desktop/myData.cvs");

/*

myCanvas->SetBatch();
ppGraph->Draw("ACL");

//gPad->Print("/home/zalewski/Desktop/myPP.jpeg");
*/
}


int main(int argc, char** argv)
{
	google::InitGoogleLogging(argv[0]);
	rnd = new TRandom3();
	// The variable to solve for with its initial value. It will be
	// mutated in place by the solver.
	double x[15];

	// Build the problem.
	ceres::Problem problem;

	// Set up the only cost function (also known as residual). This uses
	// auto-differentiation to obtain the derivative (jacobian).
	ceres::NumericDiffOptions diffOpts;
	diffOpts.relative_step_size = 0.01;

	ceres::CostFunction* cost_function =
		new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, 3, 15>(new CostFunctor, ceres::DO_NOT_TAKE_OWNERSHIP, 2, diffOpts);
//                                                           |        |   |
//                               Finite Differencing Scheme -+        |   |
//                               Dimension of residual ---------------+   |
//                               Dimension of x --------------------------+
	problem.AddResidualBlock(cost_function, nullptr, x);

	problem.SetParameterLowerBound(x, sMWPC_1_X, -50.0);
	problem.SetParameterUpperBound(x, sMWPC_1_X, 50.0);
	problem.SetParameterLowerBound(x, sMWPC_1_Y, -50.0);
	problem.SetParameterUpperBound(x, sMWPC_1_Y, 50.0);
	problem.SetParameterLowerBound(x, sMWPC_2_X, -50.0);
	problem.SetParameterUpperBound(x, sMWPC_2_X, 50.0);
	problem.SetParameterLowerBound(x, sMWPC_2_Y, -50.0);
	problem.SetParameterUpperBound(x, sMWPC_2_Y, 50.0);
	problem.SetParameterLowerBound(x, sTarPos1, -20.0);
	problem.SetParameterUpperBound(x, sTarPos1, 20.0);
	problem.SetParameterLowerBound(x, sTarPos2, -20.0);
	problem.SetParameterUpperBound(x, sTarPos2, 20.0);
	problem.SetParameterLowerBound(x, sTarPos3, -20.0);
	problem.SetParameterUpperBound(x, sTarPos3, 20.0);
	problem.SetParameterLowerBound(x, sLang1, -3.0);
	problem.SetParameterUpperBound(x, sLang1, 3.0);
	problem.SetParameterLowerBound(x, sLang2, -3.0);
	problem.SetParameterUpperBound(x, sLang2, 3.0);
	problem.SetParameterLowerBound(x, sLang3, -3.0);
	problem.SetParameterUpperBound(x, sLang3, 3.0);
	problem.SetParameterLowerBound(x, sRang, -3.0);
	problem.SetParameterUpperBound(x, sRang, 3.0);
	problem.SetParameterLowerBound(x, sLeftDetShift, -10.0);
	problem.SetParameterUpperBound(x, sLeftDetShift, 10.0);
	problem.SetParameterLowerBound(x, sLeftDetDist, -10.0);
	problem.SetParameterUpperBound(x, sLeftDetDist, 10.0);
	problem.SetParameterLowerBound(x, sRightDetShift, -10.0);
	problem.SetParameterUpperBound(x, sRightDetShift, 10.0);
	problem.SetParameterLowerBound(x, sRightDetDist, -10.0);
	problem.SetParameterUpperBound(x, sRightDetDist, 10.0);

	for (int iii = 0; iii < 15; iii++)
	{
		x[iii] = 10;
	}
	

	// Run the solver!
   ceres::Solver::Options options;
	options.trust_region_strategy_type = ceres::DOGLEG;
	options.dogleg_type = ceres::SUBSPACE_DOGLEG;
   options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = true;
	options.num_threads = 20;
	options.max_num_iterations = 200;
	ceres::Solver::Summary summary;
	Solve(options, &problem, &summary);
	drawResults(x);
	std::cout << summary.FullReport() << "\n";
	return 0;
}

