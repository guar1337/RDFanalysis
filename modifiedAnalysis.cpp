#include "modifiedAnalysis.hh"

void trimHistogramLeft(TH1F *tmpHisto, Int_t noBins)
{
	Int_t minimumBin = tmpHisto->FindFirstBinAbove(0.1);
	for (int iii = minimumBin; iii < minimumBin + noBins; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}

}

void trimHistogramRight(TH1F *tmpHisto, Int_t noBins)
{
	Int_t maximumBin = tmpHisto->FindLastBinAbove(0.1);
	for (int iii = maximumBin - noBins; iii <= maximumBin; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}
}

void makeGraph(Int_t m_runNo)
{
	TCanvas *canvasPP = new TCanvas("canvasPP","canvasPP",10,10,900,700);
	TString mcDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut/MCDataCMHisto.root";
	TFile mcDataFile(mcDataPath.Data(), "UPDATE");

	TString histName = "ppMC_geo1_CM_set" + std::to_string(m_runNo);
	TH1F *mcPPgeo1Hist = (TH1F*)mcDataFile.Get(histName.Data());
	TH1F *mcPPgeo2Hist = (TH1F*)mcDataFile.Get(histName.ReplaceAll("geo1","geo2").Data());
	TH1F *mcPPgeo3Hist = (TH1F*)mcDataFile.Get(histName.ReplaceAll("geo2","geo3").Data());
	
	histName.ReplaceAll("pp", "dd").ReplaceAll("geo3","geo1");
	TH1F *mcDDgeo1Hist = (TH1F*)mcDataFile.Get(histName.Data());
	TH1F *mcDDgeo2Hist = (TH1F*)mcDataFile.Get(histName.ReplaceAll("geo1","geo2").Data());
	TH1F *mcDDgeo3Hist = (TH1F*)mcDataFile.Get(histName.ReplaceAll("geo2","geo3").Data());

	TString realDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut/realDataCMHisto.root";
	TFile realDataFile(realDataPath.Data(), "READ");
	histName = "ppReal_geo1_CM_set" + std::to_string(m_runNo);
	TH1F *realPPgeo1Hist = (TH1F*)realDataFile.Get(histName.Data());
	TH1F *realPPgeo2Hist = (TH1F*)realDataFile.Get(histName.ReplaceAll("geo1","geo2").Data());
	TH1F *realPPgeo3Hist = (TH1F*)realDataFile.Get(histName.ReplaceAll("geo2","geo3").Data());
	histName.ReplaceAll("pp", "dd").ReplaceAll("geo3","geo1");
	TH1F *realDDgeo1Hist = (TH1F*)realDataFile.Get(histName.Data());
	TH1F *realDDgeo2Hist = (TH1F*)realDataFile.Get(histName.ReplaceAll("geo1","geo2").Data());
	TH1F *realDDgeo3Hist = (TH1F*)realDataFile.Get(histName.ReplaceAll("geo2","geo3").Data());

	TString ppGraphTitle = "Elastic scattering 6He+p ver. " + std::to_string(m_runNo);
	TString ddGraphTitle = "Elastic scattering 6He+d ver. " + std::to_string(m_runNo);
	
	THStack *ppxSectStack = new THStack("PP X-section stack", ppGraphTitle.Data());
	TH1F *wolski = new TH1F("wolski","wolski",180,0,180);
	#include "/home/zalewski/Desktop/6He/analysis/dataOut/source/wolski.hh"
	canvasPP->SetLogy();
	Int_t normFactor = 2000;
	realPPgeo1Hist->Divide(mcPPgeo1Hist);
	realPPgeo1Hist->Scale(normFactor);
	realPPgeo1Hist->SetMarkerStyle(20);
	realPPgeo1Hist->SetMarkerSize(1);
	realPPgeo1Hist->SetMarkerColor(kRed);

	realPPgeo2Hist->Divide(mcPPgeo2Hist);
	realPPgeo2Hist->Scale(normFactor/3.72);
	realPPgeo2Hist->SetMarkerStyle(20);
	realPPgeo2Hist->SetMarkerSize(1);
	realPPgeo2Hist->SetMarkerColor(kGreen);

	realPPgeo3Hist->Divide(mcPPgeo3Hist);
	realPPgeo3Hist->Scale(normFactor/4.0);
	realPPgeo3Hist->SetMarkerStyle(20);
	realPPgeo3Hist->SetMarkerSize(1);
	realPPgeo3Hist->SetMarkerColor(kBlue);

	wolski->SetLineStyle(1);


	trimHistogramLeft(realPPgeo1Hist, 3);
	trimHistogramRight(realPPgeo1Hist, 4);
	trimHistogramLeft(realPPgeo2Hist, 5);
	trimHistogramRight(realPPgeo2Hist, 4);
	trimHistogramLeft(realPPgeo3Hist, 4);
	trimHistogramRight(realPPgeo3Hist, 1);

	ppxSectStack->Add(realPPgeo1Hist);
	ppxSectStack->Add(realPPgeo2Hist);
	ppxSectStack->Add(realPPgeo3Hist);
	ppxSectStack->Add(wolski, "L");

	ppxSectStack->SetMinimum(0.01);
	ppxSectStack->SetMaximum(1000);
	ppxSectStack->Draw();
	ppxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ppxSectStack->GetXaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ppxSectStack->Draw("HIST, nostack,p");
	TString outPPFileName = "/home/zalewski/Desktop/6He/analysis/dataOut/pp_v" + std::to_string(m_runNo) + ".png";
	canvasPP->Print(outPPFileName.Data());
	delete ppxSectStack;
	delete canvasPP;

	TCanvas *canvasDD = new TCanvas("canvasDD","canvasDD",10,10,1200,700);
	THStack *ddxSectStack = new THStack("DD X-section stack", ddGraphTitle.Data());
	TH1F *pereyPerey = new TH1F("perey","pereyPerey",90,0,180);
	#include "/home/zalewski/Desktop/6He/analysis/dataOut/source/pereyPerey.hh"
	TH1F *daehnick = new TH1F("dae","daehnick",90,0,180);
	#include "/home/zalewski/Desktop/6He/analysis/dataOut/source/daehnick.hh"
	canvasDD->SetLogy();
	normFactor = 5500;
	realDDgeo1Hist->Divide(mcDDgeo1Hist);
	realDDgeo1Hist->Scale(normFactor);
	realDDgeo1Hist->SetMarkerStyle(20);
	realDDgeo1Hist->SetMarkerSize(1);
	realDDgeo1Hist->SetMarkerColor(kRed);

	realDDgeo2Hist->Divide(mcDDgeo2Hist);
	realDDgeo2Hist->Scale(normFactor/3.72);
	realDDgeo2Hist->SetMarkerStyle(20);
	realDDgeo2Hist->SetMarkerSize(1);
	realDDgeo2Hist->SetMarkerColor(kGreen);

	realDDgeo3Hist->Divide(mcDDgeo3Hist);
	realDDgeo3Hist->Scale(normFactor/4.0);
	realDDgeo3Hist->SetMarkerStyle(20);
	realDDgeo3Hist->SetMarkerSize(1);
	realDDgeo3Hist->SetMarkerColor(kBlue);

	pereyPerey->SetLineStyle(1);
	pereyPerey->SetLineColor(kBlack);

	daehnick->SetLineStyle(10);
	daehnick->SetLineColor(kBlack);

	trimHistogramLeft(realDDgeo1Hist, 3);
	trimHistogramRight(realDDgeo1Hist, 6);
	trimHistogramLeft(realDDgeo2Hist, 5);
	trimHistogramRight(realDDgeo2Hist, 5);
	trimHistogramLeft(realDDgeo3Hist, 6);
	trimHistogramRight(realDDgeo3Hist, 1);

	ddxSectStack->Add(realDDgeo1Hist);
	ddxSectStack->Add(realDDgeo2Hist);
	ddxSectStack->Add(realDDgeo3Hist);
	ddxSectStack->Add(daehnick, "L");
	ddxSectStack->Add(pereyPerey, "L");

	ddxSectStack->SetMinimum(0.01);
	ddxSectStack->SetMaximum(10000);
	ddxSectStack->Draw();
	ddxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ddxSectStack->GetXaxis()->CenterTitle();
	ddxSectStack->GetXaxis()->SetRangeUser(0,140);
	ddxSectStack->GetYaxis()->CenterTitle();
	ddxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ddxSectStack->Draw("HIST, nostack,p");

	TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
	legend->AddEntry(realDDgeo1Hist,"Geometry 1","p");
	legend->AddEntry(realDDgeo2Hist,"Geometry 2","p");
	legend->AddEntry(realDDgeo3Hist,"Geometry 3","p");
	legend->AddEntry(daehnick,"Daehnick, PRC 6, 2253 (1980)","L");
	legend->AddEntry(pereyPerey,"Perey-Perey","L");
	legend->Draw();

	TString outDDFileName = "/home/zalewski/Desktop/6He/analysis/dataOut/dd_v" + std::to_string(m_runNo) + ".png";
	canvasDD->Print(outDDFileName.Data());

	delete canvasDD;
	delete ddxSectStack;
}

void makeSmallFile(Int_t m_runNo)
{
	TString chainName = "smallMC";
	TChain smallChain(chainName.Data());
	TString nameModifier = std::to_string(m_runNo);
	smallChain.Add("/home/zalewski/dataTmp/small/geo1/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo1/mc_out_2H_v" + nameModifier + ".root");

	smallChain.Add("/home/zalewski/dataTmp/small/geo2/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo2/mc_out_2H_v" + nameModifier + ".root");

	smallChain.Add("/home/zalewski/dataTmp/small/geo3/mc_out_1H_v" + nameModifier + ".root");
	smallChain.Add("/home/zalewski/dataTmp/small/geo3/mc_out_2H_v" + nameModifier + ".root");

	ROOT::RDataFrame smallDF(smallChain);
	TString outFName = "/home/zalewski/dataTmp/small/smallMC_v" + nameModifier + ".root";
	smallDF.Snapshot("small", outFName.Data());
}

void loadGeometryCorrectionParameters(Int_t m_runNo)
{
	parameters.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/experimental2/chosen.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("loadGeometryCorrectionParameters:\tFailed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	float tmpContainer;
    if (verbosity==true) printf("%d.\t", m_runNo);
	for (int iii = 0; iii < 10; iii++)
	{
		outStreamGenerated>>tmpContainer;
		if (iii!=0)
		{
			parameters.push_back(tmpContainer);
			if (verbosity==true) printf("%s = %f\t", parNames[iii-1].c_str(), parameters[iii-1]);
		}
	}
	if (verbosity==true) std::cout<<std::endl;
}

void loadCutsCorrectionParameters(Int_t m_runNo)
{
	vCut.clear();
	std::string line;
	std::string fName = "/home/zalewski/Desktop/6He/analysis/dataOut/tarVertexCuts.txt";

	std::ifstream outStreamGenerated(fName, std::ios::in);
	if (!outStreamGenerated)
	{
		printf("loadCutsCorrectionParameters:\tFailed to open file: %s\n", fName.c_str());
	}

	int jumpTo = m_runNo;
	for (int iii = 0; iii<jumpTo; iii++)
	{
		std::getline(outStreamGenerated, line, ';');
	}
	
	//printf("%s\n", line.c_str());

	Int_t tmpContainer;
    if (verbosity==true) printf("%d.\t", m_runNo);
	for (int iii = 0; iii < 12; iii++)
	{
		outStreamGenerated>>tmpContainer;
		vCut.push_back(tmpContainer);
		if (verbosity==true) printf("%s = %d\t", corrNames[iii].c_str(), vCut[iii]);
	}
	if (verbosity==true) std::cout<<std::endl;
	
}

std::vector<Double_t> getMWPC(Double_t m_MWPC_1_X, Double_t m_MWPC_1_Y, Double_t m_MWPC_2_X, Double_t m_MWPC_2_Y)
{
	std::vector<Double_t> rvecMWPC(6);
	rvecMWPC.at(0) = m_MWPC_1_X + parameters[sMWPC_1_X];
	rvecMWPC.at(1) = m_MWPC_1_Y + parameters[sMWPC_1_Y];
	rvecMWPC.at(2) = -816.0;
	rvecMWPC.at(3) = m_MWPC_2_X + parameters[sMWPC_2_X];
	rvecMWPC.at(4) = m_MWPC_2_Y + parameters[sMWPC_2_Y];
	rvecMWPC.at(5) = -270.0;

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

TVector3 getTarVertex(std::vector<Double_t> rvecMWPC, Int_t m_geo)
{
	Double_t m_tarPos(0.0);
	Double_t m_tarAngle;

	switch (m_geo)
	{
	case 1:
		m_tarAngle = 45.0*TMath::DegToRad();
		m_tarPos = parameters[sTarPos];
		break;

	case 2:
		m_tarAngle = 6.0*TMath::DegToRad();
		m_tarPos = parameters[sTarPos];
		break;

	case 3:
		m_tarAngle = 0.0*TMath::DegToRad();
		m_tarPos = parameters[sTarPos];
		break;
	
	default:
		break;
	}

	Double_t m_dX = rvecMWPC[3] - rvecMWPC[0];
	Double_t m_dY = rvecMWPC[4] - rvecMWPC[1];
	Double_t m_dZ = rvecMWPC[5] - rvecMWPC[2];
	TVector3 m_vBeam(m_dX, m_dY, m_dZ);
		
	TVector3 m_tarPoint(0.0, 0.0, parameters[sTarPos]);
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
	Double_t X2hDet = -cs::widthStripX * m_xStrip * cos(leftAngle[m_geo-1]);
	Double_t Y2hDet = cs::widthStripY * m_yStrip;
	Double_t Z2hDet = cs::widthStripX * m_xStrip * sin(leftAngle[m_geo-1]);
	return TVector3(X2hDet, Y2hDet, Z2hDet);
}

TVector3 getRightDetVertex(Int_t m_xStrip, Int_t m_yStrip, Int_t m_geo)
{
	Double_t xStrip = m_xStrip + (gRandom->Uniform(0.0,1.0)-0.5);
	Double_t yStrip = m_yStrip + (gRandom->Uniform(0.0,1.0)-0.5);

	// coordinates of hit in LAB system
	Double_t X6HeDet = cs::widthStripX * m_xStrip * cos(rightAngle);
	Double_t Y6HeDet = cs::widthStripY * m_yStrip;
	Double_t Z6HeDet = cs::widthStripX * m_xStrip * sin(rightAngle);
	return TVector3(X6HeDet, Y6HeDet, Z6HeDet);
}

TVector3 getLeftDetPosition(Int_t m_geo)
{
	Double_t X2Hlab = sqlDist*sin(leftAngle[m_geo-1]) + (cs::sqlXzero) * cos(leftAngle[m_geo-1]);
	Double_t Y2Hlab = cs::sqlYstart + cs::widthStripY;
	Double_t Z2Hlab = sqlDist*cos(leftAngle[m_geo-1]) - (cs::sqlXzero) * sin(leftAngle[m_geo-1]);
	return TVector3(X2Hlab, Y2Hlab, Z2Hlab);
}

TVector3 getRightDetPosition(Int_t m_geo)
{
	Double_t X6Helab = sqrDist*sin(-rightAngle) - (cs::sqlXzero) * cos(rightAngle);
	Double_t Y6Helab = cs::sqrYstart + cs::widthStripY;
	Double_t Z6Helab = sqrDist*cos(rightAngle) - (cs::sqlXzero) * sin(rightAngle);
	return TVector3(X6Helab, Y6Helab, Z6Helab);
}

Double_t getMCAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi() - (acos(m_cosLeftAng)-m_ThetaCM2H))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getMCAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	TVector3 m_vBoost = m_lvCM.BoostVector();

	Double_t m_ThetaCM2H = m_lvCM.Theta();
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi() - (acos(m_cosLeftAng)-m_ThetaCM2H))*TMath::RadToDeg();
	return m_sqlangCM;
}

void realDataAnalyzer(Int_t setNo)
{
	loadGeometryCorrectionParameters(setNo);
	loadCutsCorrectionParameters(setNo);
	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos];

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0;
	sqrDist = 250.0;

	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small/small.root";
	ROOT::RDataFrame smallDF("small", smallFileName.Data());

	auto newDF = smallDF.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
						.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
						.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
						.Define("tarVertex", getTarVertex, {"MWPC", "geo"})

						.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
						.Define("leftLabVertex", [](Int_t geo){return getLeftDetPosition(geo);}, {"geo"})
						.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
						.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip", "geo"})
						.Define("rightLabVertex", [](Int_t geo){return getRightDetPosition(geo);}, {"geo"})
						.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

						.Define("v2H", "leftGlobVertex-tarVertex")
						.Define("v6He", "rightGlobVertex-tarVertex")
						.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
                        .Define("mc1H", getMCAngle1H, {"sqlang", "lvBeam"})
					    .Define("mc2H", getMCAngle2H, {"sqlang", "lvBeam"})
						.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut1", "tarVertex.X()>vCut[mX1] && tarVertex.X()<vCut[X1] && tarVertex.Y()>vCut[mY1] && tarVertex.Y()<vCut[Y1]")
						.Define("tarCut2", "tarVertex.X()>vCut[mX2] && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]")
						.Define("tarCut3", "tarVertex.X()>vCut[mX3] && tarVertex.X()<vCut[X3] && tarVertex.Y()>vCut[mY3] && tarVertex.Y()<vCut[Y3]");

	TString ppHistoName1 = "ppReal_geo1_CM_set" + std::to_string(setNo);
	TString ppHistoName2 = "ppReal_geo2_CM_set" + std::to_string(setNo);
	TString ppHistoName3 = "ppReal_geo3_CM_set" + std::to_string(setNo);
	auto ppDF = newDF.Filter("pp").Cache<Double_t, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc1H", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ppCMHist1 = ppDF.Filter("geo==1 && tarCut1").Histo1D({ppHistoName1.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist2 = ppDF.Filter("geo==2 && tarCut2").Histo1D({ppHistoName2.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist3 = ppDF.Filter("geo==3 && tarCut3").Histo1D({ppHistoName3.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName1 = "ddReal_geo1_CM_set" + std::to_string(setNo);
	TString ddHistoName2 = "ddReal_geo2_CM_set" + std::to_string(setNo);
	TString ddHistoName3 = "ddReal_geo3_CM_set" + std::to_string(setNo);
	auto ddDF = newDF.Filter("dd").Cache<Double_t, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc2H", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ddCMHist1 = ddDF.Filter("geo==1 && tarCut1").Histo1D({ddHistoName1.Data(), "Real elastic scattering on deuterons", 90,0,180}, {"mc2H"});
	auto ddCMHist2 = ddDF.Filter("geo==2 && tarCut2").Histo1D({ddHistoName2.Data(), "Real elastic scattering on deuterons", 90,0,180}, {"mc2H"});
	auto ddCMHist3 = ddDF.Filter("geo==3 && tarCut3").Histo1D({ddHistoName3.Data(), "Real elastic scattering on deuterons", 90,0,180}, {"mc2H"});

	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut/realDataCMHisto.root", "UPDATE");
	ppCMHist1->Write(ppCMHist1->GetName(), 1);
	ppCMHist2->Write(ppCMHist2->GetName(), 1);
	ppCMHist3->Write(ppCMHist3->GetName(), 1);
	ddCMHist1->Write(ddCMHist1->GetName(), 1);
	ddCMHist2->Write(ddCMHist2->GetName(), 1);
	ddCMHist3->Write(ddCMHist3->GetName(), 1);
	histOutputFile.Close();
}

void MCDataAnalyzer(Int_t setNo)
{
	loadGeometryCorrectionParameters(setNo);
	loadCutsCorrectionParameters(setNo);
	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos];

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0;
	sqrDist = 250.0;

	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small/smallMC_v" + std::to_string(setNo) + ".root";
	ROOT::RDataFrame smallDF("small", smallFileName.Data());

	auto newDF = smallDF.Define("MWPC", getMWPC, {"MWPC_1_X", "MWPC_1_Y", "MWPC_2_X", "MWPC_2_Y"})
						.Define("vBeam", getBeamVector, {"MWPC", "kinE"})
						.Define("lvBeam", [](TVector3 vBeam, Double_t kinE){return TLorentzVector(vBeam,cs::mass6He+kinE);}, {"vBeam", "kinE"})
						.Define("tarVertex", getTarVertex, {"MWPC", "geo"})

						.Define("leftDetVertex", getLeftDetVertex, {"SQX_L_strip", "SQY_L_strip", "geo"})
						.Define("leftLabVertex", [](Int_t geo){return getLeftDetPosition(geo);}, {"geo"})
						.Define("leftGlobVertex", [](TVector3 leftDetVertex, TVector3 leftLabVertex){return TVector3(leftDetVertex+leftLabVertex);}, {"leftDetVertex", "leftLabVertex"})
						.Define("rightDetVertex", getRightDetVertex, {"SQX_R_strip", "SQY_R_strip", "geo"})
						.Define("rightLabVertex", [](Int_t geo){return getRightDetPosition(geo);}, {"geo"})
						.Define("rightGlobVertex", [](TVector3 rightDetVertex, TVector3 rightLabVertex){return TVector3(rightDetVertex+rightLabVertex);}, {"rightDetVertex", "rightLabVertex"})

						.Define("v2H", "leftGlobVertex-tarVertex")
						.Define("v6He", "rightGlobVertex-tarVertex")
						.Define("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
                        .Define("mc1H", getMCAngle1H, {"sqlang", "lvBeam"})
					    .Define("mc2H", getMCAngle2H, {"sqlang", "lvBeam"})
						.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut1", "tarVertex.X()>vCut[mX1] && tarVertex.X()<vCut[X1] && tarVertex.Y()>vCut[mY1] && tarVertex.Y()<vCut[Y1]")
						.Define("tarCut2", "tarVertex.X()>vCut[mX2] && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]")
						.Define("tarCut3", "tarVertex.X()>vCut[mX3] && tarVertex.X()<vCut[X3] && tarVertex.Y()>vCut[mY3] && tarVertex.Y()<vCut[Y3]");

	TString ppHistoName1 = "ppMC_geo1_CM_set" + std::to_string(setNo);
	TString ppHistoName2 = "ppMC_geo2_CM_set" + std::to_string(setNo);
	TString ppHistoName3 = "ppMC_geo3_CM_set" + std::to_string(setNo);
	auto ppDF = newDF.Filter("mcPP").Cache<Double_t, TLorentzVector, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc1H", "lv2H_CM.", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ppCMHist1 = ppDF.Filter("geo==1 && tarCut1").Histo1D({ppHistoName1.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist2 = ppDF.Filter("geo==2 && tarCut2").Histo1D({ppHistoName2.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist3 = ppDF.Filter("geo==3 && tarCut3").Histo1D({ppHistoName3.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName1 = "ddMC_geo1_CM_set" + std::to_string(setNo);
	TString ddHistoName2 = "ddMC_geo2_CM_set" + std::to_string(setNo);
	TString ddHistoName3 = "ddMC_geo3_CM_set" + std::to_string(setNo);
	auto ddDF = newDF.Filter("mcDD").Cache<Double_t, TLorentzVector, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc2H", "lv2H_CM.", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ddCMHist1 = ddDF.Filter("geo==1 && tarCut1").Histo1D({ddHistoName1.Data(), "MC elastic scattering on deuterons", 90,0,180}, {"mc2H"});
	auto ddCMHist2 = ddDF.Filter("geo==2 && tarCut2").Histo1D({ddHistoName2.Data(), "MC elastic scattering on deuterons", 90,0,180}, {"mc2H"});
	auto ddCMHist3 = ddDF.Filter("geo==3 && tarCut3").Histo1D({ddHistoName3.Data(), "MC elastic scattering on deuterons", 90,0,180}, {"mc2H"});

	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut/MCDataCMHisto.root", "UPDATE");
	ppCMHist1->Write(ppCMHist1->GetName(), 1);
	ppCMHist2->Write(ppCMHist2->GetName(), 1);
	ppCMHist3->Write(ppCMHist3->GetName(), 1);
	ddCMHist1->Write(ddCMHist1->GetName(), 1);
	ddCMHist2->Write(ddCMHist2->GetName(), 1);
	ddCMHist3->Write(ddCMHist3->GetName(), 1);
	histOutputFile.Close();
}

void modifiedAnalysis()
{
	ROOT::EnableImplicitMT();
	TStopwatch *stopwatch = new TStopwatch();
	verbosity=false;

    for(int iii=1; iii<=6; iii++)
    {
        //makeSmallFile(iii);
		//MCDataAnalyzer(iii);
		//realDataAnalyzer(iii);
		makeGraph(iii);
    }
	//makeSmallFile(5);
    //realDataAnalyzer(1);

}