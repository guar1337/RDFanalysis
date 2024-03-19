#include "modifiedAnalysis.hh"

TH1F makeGlobHist(TH1F *xSectHistGeo1, TH1F *xSectHistGeo2, TH1F *xSectHistGeo3)
{
	TString histName =  xSectHistGeo1->GetName();
	histName.ReplaceAll("_geo1_", "").ReplaceAll("_stack_1", "");

	TString histTitle =  xSectHistGeo1->GetTitle();
	histTitle.ReplaceAll(" in geo1 [CM deg]", "");

	TH1F xSectFinalHist(histName.Data(),histTitle.Data(), 90,0,180);
	xSectFinalHist.GetXaxis()->SetTitle(xSectHistGeo1->GetXaxis()->GetTitle());
	xSectFinalHist.GetYaxis()->SetTitle(xSectHistGeo1->GetYaxis()->GetTitle());

	Bool_t overlappingBinsFlag=false;
	for (int iii = 0; iii < 180; iii++)
	{
		if (xSectHistGeo1->GetBinContent(iii)>0.0001)
		{
			xSectFinalHist.SetBinContent(iii, xSectHistGeo1->GetBinContent(iii));
			xSectFinalHist.SetBinError(iii, xSectHistGeo1->GetBinError(iii));
			overlappingBinsFlag=true;
		}

		if (xSectHistGeo2->GetBinContent(iii)>0.0001)
		{
			xSectFinalHist.SetBinContent(iii, xSectHistGeo2->GetBinContent(iii));
			xSectFinalHist.SetBinError(iii, xSectHistGeo2->GetBinError(iii));
			if (overlappingBinsFlag==true)  throw std::invalid_argument( "overlapping bins in joining geoHists into global hist" );
		}
		
		if (xSectHistGeo3->GetBinContent(iii)>0.0001)
		{
			xSectFinalHist.SetBinContent(iii, xSectHistGeo3->GetBinContent(iii));
			xSectFinalHist.SetBinError(iii, xSectHistGeo3->GetBinError(iii));
			if (overlappingBinsFlag==true)  throw std::invalid_argument( "overlapping bins in joining geoHists into global hist" );
		}
		overlappingBinsFlag=false;
	}
	return xSectFinalHist;
}

void makeTGraph(TH1F inHist)
{
	//extract data from histogram
	std::vector<Double_t> xAxis;
	std::vector<Double_t> yAxis;
	std::vector<Double_t> yAxisError;
	std::vector<Double_t> xAxisError;
	TString outFilePath = TString::Format("/home/zalewski/Desktop/6He/analysis/dataOut/%s.csv", inHist.GetName());
	std::ofstream outStream(outFilePath.Data(), std::ios::trunc);

	outStream<<"CM angle\t"<<"dSigma/dPhi\t"<<"error"<<"\n";
	for (int iii = 1; iii < 90; iii++)
	{
		outStream<<inHist.GetBinCenter(iii)<<"\t"<<inHist.GetBinContent(iii)<<"\t";
		if (inHist.GetBinContent(iii)>0.0001)
		{
			outStream<<inHist.GetBinError(iii)<<"\n";
		}
		
		else
		outStream<<0.0<<"\n";
		//printf("bin: %d\tsigma: %f\tbinCenter: %f\n", iii, inHist.GetBinContent(iii), inHist.GetBinCenter(iii));
		if (inHist.GetBinContent(iii)>0.0001)
		{
			xAxis.push_back(inHist.GetBinCenter(iii));
			yAxis.push_back(inHist.GetBinContent(iii));
			yAxisError.push_back(inHist.GetBinError(iii));
			xAxisError.push_back(0.5);
			
			//printf("angle(CM): %f\tdSigma/dTheta: %f\terror: %f\n", inHist.GetBinCenter(iii), inHist.GetBinContent(iii), inHist.GetBinError(iii));
		}		
	}

	TGraphErrors myTGraph(xAxis.size(), &xAxis[0], &yAxis[0], &xAxisError[0], &yAxisError[0]);
	myTGraph.SetTitle(inHist.GetTitle());
	myTGraph.GetXaxis()->SetTitle(inHist.GetXaxis()->GetTitle());
	myTGraph.GetYaxis()->SetTitle(inHist.GetYaxis()->GetTitle());
	myTGraph.GetXaxis()->CenterTitle();
	myTGraph.GetYaxis()->CenterTitle();
	myTGraph.SetMinimum(0.001);
	myTGraph.SetMaximum(1000);
	TString histName = inHist.GetName();
	//pp reaction has larger xSections
	printf("%s\n", histName.Data());
	if (!histName.CompareTo("xSect1H"))
	{
		myTGraph.SetMinimum(0.1);
	}

	if (!histName.CompareTo("xSectDT"))
	{
		myTGraph.SetMinimum(0.05);
		myTGraph.SetMaximum(50);
	}

	TCanvas *myCanvas = new TCanvas("myCanvas","myCanvas",10,10,1200,700);
	myCanvas->SetLogy();
	myTGraph.Draw("AP");
	TString canvName = TString::Format("/home/zalewski/Desktop/6He/analysis/dataOut/%s.png", inHist.GetName());
	myCanvas->SaveAs(canvName.Data());
	delete myCanvas;	
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

Double_t getCMAngleDTsimple(TLorentzVector m_lv6He_CM)
{
	return m_lv6He_CM.Theta()*TMath::RadToDeg();
}

Bool_t filterTarget(Double_t m_fevx, Double_t m_fevy, Int_t m_geo)
{
	Bool_t tarCut = false;
	if (m_geo==1)
	{
		tarCut = (m_fevx>vCut[mX1] && m_fevx<vCut[X1] && m_fevy>vCut[mY1] && m_fevy<vCut[Y1]);
		if (trigger)
		{
			printf("vCut[mX%d] = %d\tvCut[X%d] = %d\tvCut[mY%d] = %d\tvCut[Y%d] = %d\n", m_geo, vCut[mX1], m_geo, vCut[X1], m_geo, vCut[mY1], m_geo, vCut[Y1]);
			trigger=false;
		}
		
	}

	else if (m_geo==2)
	{
		tarCut = (m_fevx>vCut[mX2] && m_fevx<vCut[X2] && m_fevy>vCut[mY2] && m_fevy<vCut[Y2]);
		if (trigger)
		{
			printf("vCut[mX%d] = %d\tvCut[X%d] = %d\tvCut[mY%d] = %d\tvCut[Y%d] = %d\n", m_geo, vCut[mX2], m_geo, vCut[X2], m_geo, vCut[mY2], m_geo, vCut[Y2]);
			trigger=false;
		}
	}
	
	else
	{
		tarCut = (m_fevx>vCut[mX3] && m_fevx<vCut[X3] && m_fevy>vCut[mY3] && m_fevy<vCut[Y3]);
		if (trigger)
		{
			printf("vCut[mX%d] = %d\tvCut[X%d] = %d\tvCut[mY%d] = %d\tvCut[Y%d] = %d\n", m_geo, vCut[mX3], m_geo, vCut[X3], m_geo, vCut[mY3], m_geo, vCut[Y3]);
			trigger=false;
		}
	}
	
	return tarCut;
}

void loadCutsCorrectionParameters(Int_t m_runNo=1)
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

Double_t getCMAngle1H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass1H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi()-acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t getCMAngle2H(Double_t m_LAng, TLorentzVector m_lvBeam)
{
	TLorentzVector m_lvTar(0.0, 0.0, 0.0, cs::mass2H);
	TLorentzVector m_lvCM = m_lvTar + m_lvBeam;
	Double_t m_gammaSquare2H = m_lvCM.Gamma()*m_lvCM.Gamma();
	Double_t m_tanSquare = tan(m_LAng*TMath::DegToRad()) * tan(m_LAng*TMath::DegToRad());
	Double_t m_cosLeftAng = (1.0 - m_gammaSquare2H*m_tanSquare)/(1.0 + m_gammaSquare2H*m_tanSquare);
	Double_t m_sqlangCM = (TMath::Pi()-acos(m_cosLeftAng))*TMath::RadToDeg();
	return m_sqlangCM;
}

Double_t recoCMAngle(TLorentzVector m_lv2H_CM, TLorentzVector m_flvBeam)
{
	return m_lv2H_CM.Vect().Angle(m_flvBeam.Vect())*TMath::RadToDeg();
}

void getStatisticalError(Int_t geo, Int_t mass)
{	

	Long_t nBeam[3] = {3521461760, 7627521280, 20136056320};

	Double_t targetError[2][3] = {	{0.2555, 0.1646, 0.1755},
									{0.2238, 0.1732, 0.2268}};
	//Statistical error of each bin = Err(number of exp events) + Err(MC scattered) + Err(MC observed) + Err(tarThicnkess) + Err(nBeam)

	//error of exp events
	TFile realHitsFile("/home/zalewski/Desktop/6He/analysis/dataOut/realDataCMHisto.root", "READ");
	TString realHistName = TString::Format("real_geo%d_%dH", geo, mass);
	TH1F *realHist = (TH1F*)realHitsFile.Get(realHistName.Data());

	TString realHistErrorName = TString::Format("realError_geo%d_%dH", geo, mass);
	TH1F realHistError(realHistErrorName.Data(), realHistErrorName.Data(), 90,0,180);
	
	for (int iii = 0; iii < 90; iii++)
	{
		if (realHist->GetBinContent(iii)>0)
		{
			realHistError.SetBinContent(iii, sqrt(realHist->GetBinContent(iii))/realHist->GetBinContent(iii));
		}
	}

	//error of MC scattered
	TFile MCScatteredEventsFile("/home/zalewski/Desktop/6He/analysis/dataOut/efficiency.root", "READ");
	TString MCScatteredHistName = TString::Format("scattered_geo%d_%dH", geo, mass);
	TH1F *MCScatteredHist = (TH1F*)MCScatteredEventsFile.Get(MCScatteredHistName.Data());

	TString MCScatteredHistErrorName = TString::Format("scatteredError_geo%d_%dH", geo, mass);
	TH1F  MCScatteredHistError(MCScatteredHistErrorName.Data(), MCScatteredHistErrorName.Data(), 90,0,180);
	
	for (int iii = 0; iii < 90; iii++)
	{
		if (MCScatteredHist->GetBinContent(iii)>0)
		{
			MCScatteredHistError.SetBinContent(iii, sqrt(MCScatteredHist->GetBinContent(iii))/MCScatteredHist->GetBinContent(iii));
			//printf("Exp events: %f\tErr: %f\n", MCScatteredHist->GetBinContent(iii), sqrt(MCScatteredHist->GetBinContent(iii))/MCScatteredHist->GetBinContent(iii));
		}		
	}	


	//error of MC observed
	TFile MCObservedEventsFile("/home/zalewski/Desktop/6He/analysis/dataOut/efficiency.root", "READ");
	TString MCObservedHistName = TString::Format("observed_geo%d_%dH", geo, mass);
	TH1F *MCObservedHist = (TH1F*)MCObservedEventsFile.Get(MCObservedHistName.Data());

	TString MCObservedHistErrorName = TString::Format("observedError_geo%d_%dH", geo, mass);
	TH1F  MCObservedHistError(MCObservedHistErrorName.Data(), MCObservedHistErrorName.Data(), 90,0,180);
	
	for (int iii = 0; iii < 90; iii++)
	{
		if (MCObservedHist->GetBinContent(iii)>0)
		{
			MCObservedHistError.SetBinContent(iii, sqrt(MCObservedHist->GetBinContent(iii))/MCObservedHist->GetBinContent(iii));
			//printf("Exp events: %f\tErr: %f\n", MCObservedHist->GetBinContent(iii), sqrt(MCObservedHist->GetBinContent(iii))/MCObservedHist->GetBinContent(iii));
		}
	}	

	//error of Beam and Tar
	TString beamTargetHistErrorName = TString::Format("beamTargetError_geo%d_%dH", geo, mass);
	TH1F  beamTargetHistError(beamTargetHistErrorName.Data(), beamTargetHistErrorName.Data(), 90,0,180);
	for (int iii = 0; iii < 90; iii++)
	{
		beamTargetHistError.SetBinContent(iii, sqrt(nBeam[geo-1])/nBeam[geo-1] + targetError[mass-1][geo-1]);
		//printf("Exp events: %f\tErr: %f\n", sqrt(nBeam[geo-1])/nBeam[geo-1], sqrt(nBeam[geo-1])/nBeam[geo-1] + targetError[mass-1][geo-1]);
	}	

	TString totalHistErrorName = TString::Format("totalError_geo%d_%dH", geo, mass);
	TH1F totalErrorHist(totalHistErrorName.Data(), totalHistErrorName.Data(), 90,0,180);

	totalErrorHist.Add(&realHistError);
	//totalErrorHist.Add(&MCScatteredHistError);
	//totalErrorHist.Add(&MCObservedHistError);
	//totalErrorHist.Add(&beamTargetHistError);

	for (int iii = 0; iii < 90; iii++)
	{
		printf("Exp events: %d\tErr: %f\n", iii, totalErrorHist.GetBinContent(iii));
	}
	
	TFile errorOutFile("/home/zalewski/Desktop/6He/analysis/dataOut/error.root", "UPDATE");
	totalErrorHist.Write(totalErrorHist.GetName(), 1);
}

void getAngleError(Int_t geo, Int_t mass)
{
	//TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut/tmp.root", "UPDATE");

	TString histoName = TString::Format("angError_geo%d_%dH", geo, mass);	
	TString histoTitle = TString::Format("angError of theta angle scattered of {}^{%d}H in geo%d [CM deg]", mass, geo);
	TString inFName = TString::Format("/home/zalewski/dataTmp/small/v3/mc_out_geo%d_%dH_big.root", geo, mass);
	ROOT::RDataFrame inDF("smallMC", inFName.Data());
	if (mass==1)
	{
		auto error = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
		 				 .Filter([geo](Double_t fZ2H){return (fZ2H>lCut || geo>1);}, {"fZ2H"})
						 .Filter("sqlde>0 && sqrde>0 && mcHe6 && mcPP")
						 .Define("thetaCMAngle", getCMAngle1H, {"sqlang","lvBeam"})
						 .Define("reTheta", recoCMAngle, {"lv6He_CM.", "flvBeam."})
						 .Define("diff", "thetaCMAngle-reTheta");
						 //.Histo1D<Double_t>({"errorHist", "errorHist", 100,-5,5}, {"diff"});
		for (int iii = 0; iii < 90; iii++)
		{
			Double_t err = error.Filter([iii](Double_t thetaCM){return (thetaCM>iii*2&&thetaCM<(iii+1)*2);}, {"reTheta"}).StdDev("diff").GetValue();
			printf("geo: %d\tion: %dH\terror[%d] = %f\n", geo, mass, iii, err);
		}
	}

	else if (mass==2)
	{
		auto error = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
		 				 .Filter([geo](Double_t fZ2H){return (fZ2H>lCut || geo>1);}, {"fZ2H"})
						 .Filter("sqlde>0 && sqrde>0 && mcHe6 && mcDD")
						 .Define("thetaCMAngle", getCMAngle2H, {"sqlang","lvBeam"})
						 .Define("reTheta", recoCMAngle, {"lv6He_CM.", "flvBeam."})
						 .Define("diff", "thetaCMAngle-reTheta");
						 //.Histo1D<Double_t>({"errorHist", "errorHist", 100,-5,5}, {"diff"});
		for (int iii = 0; iii < 90; iii++)
		{
			Double_t err = error.Filter([iii](Double_t thetaCM){return (thetaCM>iii*2&&thetaCM<(iii+1)*2);}, {"reTheta"}).StdDev("diff").GetValue();
			printf("geo: %d\tion: %dH\terror[%d] = %f\n", geo, mass, iii, err);
		}
	}
}

std::vector<Double_t> getRelativeError(TH1F *my1DHist)
{
	//prepare relative errors for real data
	Int_t noBins = my1DHist->GetXaxis()->GetNbins();
	std::vector<Double_t> relativeError(noBins);
	for (int iii = 0; iii < noBins; iii++)
	{
		if (my1DHist->GetBinContent(iii)>0.0)
		{
			relativeError.at(iii) = my1DHist->GetBinError(iii)/my1DHist->GetBinContent(iii);
			//printf("Bin: %d\tbin value: %f\tbin error: %f\t%f\n", iii, my1DHist->GetBinContent(iii),my1DHist->GetBinError(iii), relativeError.at(iii)*100.0);
		}	

	}
	return relativeError;
}

void setRelativeError(TH1F *myHist, std::vector<Double_t> relativeError, Double_t additionalError = 0.0)
{
	for (int iii = 0; iii < myHist->GetXaxis()->GetNbins(); iii++)
	{
		//realPPgeo1HistError.at(iii) = realPPgeo1Hist->GetBinError(iii)/realPPgeo1Hist->GetBinContent(iii);
		//printf("Error at bin[%d] = %f\n", iii, realPPgeo1HistError.at(iii));
		myHist->SetBinError(iii, myHist->GetBinContent(iii)*(relativeError.at(iii) + additionalError));
	}
}

void calculateXsection(Int_t geo, Int_t mass, Double_t beamTargetIntegral, Double_t additionalError = 0.0)
{
	//efficiency histogram
	TString efficiencyDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut/efficiency.root";
	TFile efficiencyDataFile(efficiencyDataPath.Data(), "READ");	
	TString efficiencyHistName = TString::Format("efficiency_geo%d_%dH", geo, mass);
	TH1F *efficiencyHist = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());

	//theta integral histogram
	TH1F thetaIntegral("thetaIntegral", "thetaIntegral", 90,0,180);
	for (int iii = 0; iii < 90; iii++)
	{
		thetaIntegral.SetBinContent(iii, cos(2.0*iii*TMath::DegToRad())-cos(2.0*(iii+1)*TMath::DegToRad()));
	}

	//real data histogram
	TString realDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut/realDataCMHisto.root";
	TFile realDataFile(realDataPath.Data(), "READ");
	TString realHistName = TString::Format("real_geo%d_%dH", geo, mass);
	TH1F *realHist = (TH1F*)realDataFile.Get(realHistName.Data());

	std::vector<Double_t> relativeError = getRelativeError(realHist);

	TString xSectionHistName = TString::Format("xSect_geo%d_%dH", geo, mass);
	TString xSectionHistTitle = TString::Format("Elastic scattering cross section for {}^{6}He + {}^{%d}H reaction in geo%d [CM deg]", mass, geo);
	TH1F xSection(*realHist);
	xSection.SetNameTitle(xSectionHistName.Data(), xSectionHistTitle.Data());
	xSection.Divide(&thetaIntegral);
	xSection.Divide(efficiencyHist);
	xSection.Scale(beamTargetIntegral/(2.0*TMath::Pi()));
	setRelativeError(&xSection, relativeError, additionalError);


	xSection.SetMarkerStyle(20);
	xSection.SetMarkerSize(1);
	xSection.SetMarkerColor(kRed);
	if (geo==2)
	{
		xSection.SetMarkerColor(kGreen);
	}

	else if (geo==3)
	{
		xSection.SetMarkerColor(kBlue);
	}
	
	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/dataOut/xSection.root", "UPDATE");
	xSection.Write(xSection.GetName(), 1);

	getRelativeError(&xSection);
}

void calculateXsectionDT()
{
	Double_t beamTargetIntegral = 1e-4*1.514;
	Double_t additionalError = 0.0;
	//efficiency histogram
	TString efficiencyDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut_dt/efficiency.root";
	TFile efficiencyDataFile(efficiencyDataPath.Data(), "READ");	
	TString efficiencyHistName = TString::Format("efficiency_dt");
	TH1F *efficiencyHist = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());

	//theta integral histogram
	TH1F thetaIntegral("thetaIntegral", "thetaIntegral", 90,0,180);
	for (int iii = 0; iii < 90; iii++)
	{
		thetaIntegral.SetBinContent(iii, cos(2.0*iii*TMath::DegToRad())-cos(2.0*(iii+1)*TMath::DegToRad()));
	}

	thetaIntegral.SaveAs("/home/zalewski/Desktop/thetaIntegral.png");
	//real data histogram
	TString realDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut_dt/realDataCMHisto.root";
	TFile realDataFile(realDataPath.Data(), "READ");
	TString realHistName = TString::Format("real_dt");
	TH1F *realHist = (TH1F*)realDataFile.Get(realHistName.Data());

	std::vector<Double_t> relativeError = getRelativeError(realHist);

	TString xSectionHistName = TString::Format("xSect_dt");
	TString xSectionHistTitle = TString::Format("{}^{6}He({}^{2}H,{}^{3}H){}^{5}He reaction cross section[CM deg]");
	TH1F xSection(*realHist);
	xSection.SetNameTitle(xSectionHistName.Data(), xSectionHistTitle.Data());
	xSection.Divide(&thetaIntegral);
	xSection.Multiply(efficiencyHist);
	xSection.Scale(beamTargetIntegral/(2.0*TMath::Pi()));
	setRelativeError(&xSection, relativeError, additionalError);


	xSection.SetMarkerStyle(20);
	xSection.SetMarkerSize(1);
	xSection.SetMarkerColor(kBlack);
	
	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/dataOut_dt/xSection.root", "UPDATE");
	xSection.Write(xSection.GetName(), 1);

	getRelativeError(&xSection);
}

void calculateEfficiency(Int_t geo, Int_t mass)
{
	//trigger=true;
	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut/efficiency.root", "UPDATE");
	
	TString inFName = TString::Format("/home/zalewski/dataTmp/small/v3/mc_out_geo%d_%dH.root", geo, mass);
	TString histoName = TString::Format("scattered_geo%d_%dH", geo, mass);	
	TString histoTitle = TString::Format("Theta angle scattered of {}^{%d}H in geo%d [CM deg]", mass, geo);
	TString tarCutFilter = TString::Format("fevx>vCut[mX%d] && fevx<vCut[X%d] && fevy>vCut[mY%d] && fevy<vCut[Y%d]", geo, geo, geo, geo);

	ROOT::RDataFrame inDF("smallMC", inFName.Data());
	if (mass==1)
	{
		auto scattered = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
							 .Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"fThetaCM"});
		scattered->Write(scattered->GetName(), 1);

		//pp in geo1 observed
		histoName.ReplaceAll("scattered", "observed");
		histoTitle.ReplaceAll("scattered", "observed");
		auto observed = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
		 					.Filter([geo](Double_t fZ2H){return (fZ2H>lCut || geo>1);}, {"fZ2H"})
							.Filter("sqlde>0 && sqrde>0 && mcPP")
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"fThetaCM"});
		observed->Write(observed->GetName(), 1);

		TH1D efficiency(*observed);
		efficiency.Divide(&(scattered.GetValue()));
		histoName.ReplaceAll("observed", "efficiency");
		histoTitle = TString::Format("Efficiency of {}^{%d}H detection in geo%d [CM deg]", mass, geo);
		efficiency.SetNameTitle(histoName.Data(), histoTitle.Data());
		efficiency.SetMarkerStyle(2);
		efficiency.SetMarkerColor(kRed);

		if (geo==2)
		{
			efficiency.SetMarkerColor(kGreen);
		}

		else if (geo==3)
		{
			efficiency.SetMarkerColor(kBlue);
		}

		efficiency.Write(efficiency.GetName(), 1);
	}

	else if (mass==2)
	{
		auto scattered = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
							 .Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"fThetaCM"});
		scattered->Write(scattered->GetName(), 1);

		//dd in geo1 observed
		histoName.ReplaceAll("scattered", "observed");
		histoTitle.ReplaceAll("scattered", "observed");
		auto observed = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
							.Filter([geo](Double_t fZ2H){return (fZ2H>lCut || geo>1);}, {"fZ2H"})
							.Filter("sqlde>0 && sqrde>0 && mcDD")
							.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"fThetaCM"});
		observed->Write(observed->GetName(), 1);

		TH1D efficiency(*observed);
		efficiency.Divide(&(scattered.GetValue()));
		histoName.ReplaceAll("observed", "efficiency");
		histoTitle = TString::Format("Efficiency of {}^{%d}H detection in geo%d [CM deg]", mass, geo);
		efficiency.SetNameTitle(histoName.Data(), histoTitle.Data());
		efficiency.SetMarkerStyle(2);
		efficiency.SetMarkerColor(kRed);

		if (geo==2)
		{
			efficiency.SetMarkerColor(kGreen);
		}

		else if (geo==3)
		{
			efficiency.SetMarkerColor(kBlue);
		}

		efficiency.Write(efficiency.GetName(), 1);
	}
	
	histOutputFile.Close();
}

void calculateEfficiencyDT()
{
	Int_t geo=2;
	//trigger=true;
	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut_dt/efficiency.root", "UPDATE");
	
	TString inFName = TString::Format("/home/zalewski/dataTmp/small/v3/mc_out_geo2_2H_dt.root");
	TString histoName = TString::Format("scattered_dt");	
	TString histoTitle = TString::Format("Theta angle scattered of {}^{3}H in geo2 [CM deg]");
	TString tarCutFilter = TString::Format("fevx>vCut[mX2] && fevx<vCut[X2] && fevy>vCut[mY2] && fevy<vCut[Y2]");

	ROOT::RDataFrame inDF("smallMC", inFName.Data());

	auto scattered = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
						 //.Define("thetaCMAngle", getCMAngleDT,{"fsqrang","lvBeam"})
						 //.Define("thetaCMAngle", getCMAngleDT, {"sqrang", "lvBeam"})
						 .Define("thetaCMAngle", "180.0 - fThetaCM")
						 .Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"thetaCMAngle"});
	scattered->Write(scattered->GetName(), 1);

	histoName.ReplaceAll("scattered", "observed");
	histoTitle.ReplaceAll("scattered", "observed");
	auto observed = inDF.Filter([geo](Double_t fevx, Double_t fevy){return filterTarget(fevx, fevy, geo);}, {"fevx", "fevy"})
						.Filter("mcLeneLang && mcHe4")
						.Define("thetaCMAngle", getCMAngleDT, {"sqrang", "lvBeam"})
						.Histo1D<Double_t>({histoName.Data(), histoTitle.Data(), 90,0,180}, {"dtCM"});
	observed->Write(observed->GetName(), 1);

	TH1D efficiency(*scattered);
	efficiency.Divide(&(observed.GetValue()));
	histoName.ReplaceAll("observed", "efficiency");
	histoTitle = TString::Format("Efficiency of dt reaction detection in geo2 [CM deg]");
	efficiency.SetNameTitle(histoName.Data(), histoTitle.Data());
	efficiency.SetMarkerStyle(2);
	efficiency.SetMarkerColor(kBlack);
	efficiency.Write(efficiency.GetName(), 1);

	histOutputFile.Close();
}

void trimHistogramLeft(TH1F *tmpHisto, Int_t noBins)
{
	Int_t minimumBin = tmpHisto->FindFirstBinAbove(0.001);
	for (int iii = minimumBin; iii < minimumBin + noBins; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}
}

void trimHistogramRight(TH1F *tmpHisto, Int_t noBins)
{
	Int_t maximumBin = tmpHisto->FindLastBinAbove(0.001);
	for (int iii = maximumBin - noBins; iii <= maximumBin; iii++)
	{
		tmpHisto->SetBinContent(iii, 0.0);
	}
}

void makeGraph()
{
	//load histograms
	TString efficiencyDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut/efficiency.root";
	TFile efficiencyDataFile(efficiencyDataPath.Data(), "READ");	
	TString efficiencyHistName = "efficiency_geo1_1H";
	TH1F *eff_geo1_1H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());
	TH1F *eff_geo2_1H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.ReplaceAll("geo1","geo2").Data());
	TH1F *eff_geo3_1H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.ReplaceAll("geo2","geo3").Data());

	efficiencyHistName.ReplaceAll("1H", "2H").ReplaceAll("geo3","geo1");
	TH1F *eff_geo1_2H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());
	TH1F *eff_geo2_2H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.ReplaceAll("geo1","geo2").Data());
	TH1F *eff_geo3_2H = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.ReplaceAll("geo2","geo3").Data());

	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/dataOut/xSection.root", "UPDATE");
	TString histoName = "xSect_geo1_1H";
	
	TH1F *xSect_geo1_1H = (TH1F*)xSectionFile.Get(histoName.Data());
	TH1F *xSect_geo2_1H = (TH1F*)xSectionFile.Get(histoName.ReplaceAll("geo1","geo2").Data());
	TH1F *xSect_geo3_1H = (TH1F*)xSectionFile.Get(histoName.ReplaceAll("geo2","geo3").Data());

	histoName.ReplaceAll("1H", "2H").ReplaceAll("geo3","geo1");
	TH1F *xSect_geo1_2H = (TH1F*)xSectionFile.Get(histoName.Data());
	TH1F *xSect_geo2_2H = (TH1F*)xSectionFile.Get(histoName.ReplaceAll("geo1","geo2").Data());
	TH1F *xSect_geo3_2H = (TH1F*)xSectionFile.Get(histoName.ReplaceAll("geo2","geo3").Data());

	TFile theoHistFile("/home/zalewski/Desktop/6He/analysis/dataOut/source/theoHists.root", "READ");
	TH1F *wolskiExp = (TH1F*)theoHistFile.Get("wolskiExp");
	TH1F *wolskiOM = (TH1F*)theoHistFile.Get("wolskiOM");
	TH1F *expPP = (TH1F*)theoHistFile.Get("expPP");
	TH1F *expPP1 = (TH1F*)theoHistFile.Get("expPP1");
	TH1F *perey = (TH1F*)theoHistFile.Get("perey");
	TH1F *daehnick = (TH1F*)theoHistFile.Get("daehnick");
	TH1F *Li6d = (TH1F*)theoHistFile.Get("Li6d");
	TH1F *expDD = (TH1F*)theoHistFile.Get("expDD");
	TH1F *expDD1 = (TH1F*)theoHistFile.Get("expDD1");

/*
	//prepare relative errors for real data
	std::vector<Double_t> realPPgeo1HistError = getRelativeError(realPPgeo1Hist);
	std::vector<Double_t> realPPgeo2HistError = getRelativeError(realPPgeo2Hist);
	std::vector<Double_t> realPPgeo3HistError = getRelativeError(realPPgeo3Hist);

	std::vector<Double_t> realDDgeo1HistError = getRelativeError(realDDgeo1Hist);
	std::vector<Double_t> realDDgeo2HistError = getRelativeError(realDDgeo2Hist);
	std::vector<Double_t> realDDgeo3HistError = getRelativeError(realDDgeo3Hist);
*/
	TCanvas *canvasPP = new TCanvas("canvasPP","canvasPP",10,10,1200,700);
	TString ppGraphTitle = "Elastic scattering 6He+p";
	THStack *ppxSectStack = new THStack("PP X-section stack", ppGraphTitle.Data());
	canvasPP->SetLogy();

	eff_geo1_1H->Scale(1e-2);
	eff_geo2_1H->Scale(1e-2);
	eff_geo3_1H->Scale(1e-2);
	eff_geo1_2H->Scale(1e-2);
	eff_geo2_2H->Scale(1e-2);
	eff_geo3_2H->Scale(1e-2);

	trimHistogramLeft(xSect_geo1_1H, 2);
	trimHistogramRight(xSect_geo1_1H, 3);
	trimHistogramLeft(xSect_geo2_1H, 9);
	trimHistogramRight(xSect_geo2_1H, 6);
	trimHistogramLeft(xSect_geo3_1H, 8);
	trimHistogramRight(xSect_geo3_1H, 3);

	ppxSectStack->Add(xSect_geo1_1H, "E1");
	ppxSectStack->Add(xSect_geo2_1H, "E1");
	ppxSectStack->Add(xSect_geo3_1H, "E1");
	//ppxSectStack->Add(eff_geo1_1H);
	//ppxSectStack->Add(eff_geo2_1H);
	//ppxSectStack->Add(eff_geo3_1H);
	ppxSectStack->Add(wolskiExp, "E1");
	trimHistogramLeft(expPP1, 4);
	ppxSectStack->Add(expPP1, "E1");
	ppxSectStack->Add(wolskiOM, "L");

	ppxSectStack->SetMinimum(0.01);
	ppxSectStack->SetMaximum(1000);
	ppxSectStack->Draw();
	ppxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ppxSectStack->GetXaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->CenterTitle();
	ppxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ppxSectStack->Draw("HIST, nostack,p");

	TLegend *legendPP = new TLegend(0.7,0.7,0.9,0.9);
	legendPP->AddEntry(xSect_geo1_1H,"Geometry 1","p");
	legendPP->AddEntry(xSect_geo2_1H,"Geometry 2","p");
	legendPP->AddEntry(xSect_geo3_1H,"Geometry 3","p");
	legendPP->AddEntry(expPP1, "gas target data", "p");
	legendPP->AddEntry(wolskiExp,"Wolski data","p");
	legendPP->AddEntry(wolskiOM,"Wolski OM","L");
	legendPP->Draw();

	TString outPPFileName = TString::Format("/home/zalewski/Desktop/6He/analysis/dataOut/pp.png");
	canvasPP->Print(outPPFileName.Data());
	if (saveHistogram)
	{
		canvasPP->SaveAs("/home/zalewski/Desktop/6He/analysis/dataOut/pp_2.C");
	}
	
	delete ppxSectStack;
	delete canvasPP;

	TCanvas *canvasDD = new TCanvas("canvasDD","canvasDD",10,10,1200,700);
	TString ddGraphTitle = "Elastic scattering 6He+d";
	THStack *ddxSectStack = new THStack("DD X-section stack", ddGraphTitle.Data());
	canvasDD->SetLogy();

	trimHistogramLeft(xSect_geo1_2H, 4);
	trimHistogramRight(xSect_geo1_2H, 10);
	trimHistogramLeft(xSect_geo2_2H, 10);
	trimHistogramRight(xSect_geo2_2H, 6);
	trimHistogramLeft(xSect_geo3_2H, 10);
	trimHistogramRight(xSect_geo3_2H, 0);

	ddxSectStack->Add(xSect_geo1_2H, "E1");
	ddxSectStack->Add(xSect_geo2_2H, "E1");
	ddxSectStack->Add(xSect_geo3_2H, "E1");
	//ddxSectStack->Add(eff_geo1_2H);
	//ddxSectStack->Add(eff_geo2_2H);
	//ddxSectStack->Add(eff_geo3_2H);
	trimHistogramLeft(expDD1, 4);
	ddxSectStack->Add(expDD1, "E1");
	ddxSectStack->Add(daehnick, "L");
	//ddxSectStack->Add(perey, "L");
	//ddxSectStack->Add(Li6d);

	//ddxSectStack->SetMinimum(0.01);
	ddxSectStack->SetMaximum(10000);
	ddxSectStack->Draw();
	ddxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	ddxSectStack->GetXaxis()->CenterTitle();
	ddxSectStack->GetXaxis()->SetRangeUser(0,140);
	ddxSectStack->GetYaxis()->CenterTitle();
	ddxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	ddxSectStack->Draw("HIST, nostack,p");

	TLegend *legendDD = new TLegend(0.7,0.7,0.9,0.9);
	legendDD->AddEntry(xSect_geo1_2H,"Geometry 1","p");
	legendDD->AddEntry(xSect_geo2_2H,"Geometry 2","p");
	legendDD->AddEntry(xSect_geo3_2H,"Geometry 3","p");
	legendDD->AddEntry(daehnick,"Daehnick, PRC 6, 2253 (1980)","L");
	legendDD->AddEntry(perey,"Perey-Perey","L");
	legendDD->Draw();

	TString outDDFileName = TString::Format("/home/zalewski/Desktop/6He/analysis/dataOut/dd.png");
	canvasDD->Print(outDDFileName.Data());
	if (saveHistogram)
	{
		canvasDD->SaveAs("/home/zalewski/Desktop/6He/analysis/dataOut/dd_2.C");
	}

	delete canvasDD;
	delete ddxSectStack;

	TH1F ddFinalHist1H(makeGlobHist(xSect_geo1_1H, xSect_geo2_1H, xSect_geo3_1H));
	ddFinalHist1H.GetXaxis()->SetTitle("CM angle [deg]");
	ddFinalHist1H.GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	makeTGraph(ddFinalHist1H);

	TH1F ddFinalHist2H(makeGlobHist(xSect_geo1_2H, xSect_geo2_2H, xSect_geo3_2H));
	ddFinalHist2H.GetXaxis()->SetTitle("CM angle [deg]");
	ddFinalHist2H.GetYaxis()->SetTitle("Differential cross section [mb/sr]");	
	makeTGraph(ddFinalHist2H);
}



void makeGraphDT()
{
	//load histograms
	TString efficiencyDataPath = "/home/zalewski/Desktop/6He/analysis/dataOut_dt/efficiency.root";
	TFile efficiencyDataFile(efficiencyDataPath.Data(), "READ");	
	TString efficiencyHistName = "efficiency_dt";
	TH1F *eff_dt = (TH1F*)efficiencyDataFile.Get(efficiencyHistName.Data());

	TFile xSectionFile("/home/zalewski/Desktop/6He/analysis/dataOut_dt/xSection.root", "UPDATE");
	TString histoName = "xSect_dt";

	TFile theoHistFile("/home/zalewski/Desktop/6He/analysis/dataOut/source/theoHists.root", "READ");
	TH1F *dtTheo = (TH1F*)theoHistFile.Get("dtTheo");
	
	TH1F *xSect_dt = (TH1F*)xSectionFile.Get(histoName.Data());

	TCanvas *canvasDT = new TCanvas("canvasDT","canvasDT",10,10,1200,700);
	TString dtGraphTitle = "{}^{6}He({}^{2}H,{}^{3}H){}^{5}He reaction cross section";
	THStack *dtxSectStack = new THStack("DT X-section stack", dtGraphTitle.Data());
	canvasDT->SetLogy();

	eff_dt->Scale(1e-2);
	trimHistogramLeft(xSect_dt, 1);
	trimHistogramRight(xSect_dt, 12);

	dtxSectStack->Add(xSect_dt, "E1");
	//dtxSectStack->Add(eff_dt);
	dtxSectStack->Add(dtTheo, "L");

	dtxSectStack->SetMinimum(0.01);
	dtxSectStack->SetMaximum(1000);
	
	dtxSectStack->Draw();
	dtxSectStack->GetXaxis()->SetRangeUser(0, 140);
	dtxSectStack->GetXaxis()->SetTitle("CM angle [deg]");
	dtxSectStack->GetXaxis()->CenterTitle();
	dtxSectStack->GetYaxis()->CenterTitle();
	dtxSectStack->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	dtxSectStack->Draw("HIST, nostack,p");

/*
	TLegend *legendDT = new TLegend(0.7,0.7,0.9,0.9);
	legendDT->AddEntry(xSect_dt,"Geometry 1","p");
	legendDT->Draw();
*/

	TString outDTFileName = TString::Format("/home/zalewski/Desktop/6He/analysis/dataOut_dt/dt.png");
	canvasDT->Print(outDTFileName.Data());
	if (saveHistogram)
	{
		canvasDT->SaveAs("/home/zalewski/Desktop/6He/analysis/dataOut_dt/dt.C");
	}
	
	delete dtxSectStack;
	delete canvasDT;
	xSect_dt->SetName("xSectDT");
	xSect_dt->GetXaxis()->SetTitle("CM angle [deg]");
	xSect_dt->GetYaxis()->SetTitle("Differential cross section [mb/sr]");
	makeTGraph(*xSect_dt);
}

void makeSmallFile()
{
	TString chainName = "smallMC";
	TChain smallChain(chainName.Data());
	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo1_1H.root");
	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo1_2H.root");

	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo2_1H.root");
	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo2_2H.root");

	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo3_1H.root");
	smallChain.Add("/home/zalewski/dataTmp/small/v3/mc_out_geo3_2H.root");

	ROOT::RDataFrame smallDF(smallChain);
	TString outFName = "/home/zalewski/dataTmp/small/smallMC.root";
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

void realDataAnalyzer()
{
	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;

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
	tarPos[2] = parameters[sTarPos] - 2.0;

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + parameters[sDistL];
	sqrDist = 250.0 + parameters[sDistR];

	
	ROOT::EnableImplicitMT();
	//TString smallFileName = "/home/zalewski/dataTmp/small/small_2.root";
	TString smallFileName = "/home/zalewski/dataTmp/small/small.root";
	ROOT::RDataFrame smallDF("smallReal", smallFileName.Data());

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
						.Redefine("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
                        .Define("mc1H", getCMAngle1H, {"sqlang", "lvBeam"})
					    .Define("mc2H", getCMAngle2H, {"sqlang", "lvBeam"})
						.Define("lCut", "(v2H.X()>56.0 || geo>1)")
						.Redefine("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut1", "tarVertex.X()>vCut[mX1] && tarVertex.X()<vCut[X1] && tarVertex.Y()>vCut[mY1] && tarVertex.Y()<vCut[Y1]")
						.Define("tarCut2", "tarVertex.X()>vCut[mX2] && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]")
						.Define("tarCut3", "tarVertex.X()>vCut[mX3] && tarVertex.X()<vCut[X3] && tarVertex.Y()>vCut[mY3] && tarVertex.Y()<vCut[Y3]");

	TString ppHistoName1 = "real_geo1_1H";
	TString ppHistoName2 = "real_geo2_1H";
	TString ppHistoName3 = "real_geo3_1H";
	auto ppDF = newDF.Filter("pp").Cache<Double_t, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc1H", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ppCMHist1 = ppDF.Filter("geo==1 && tarCut1").Histo1D({ppHistoName1.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist2 = ppDF.Filter("geo==2 && tarCut2").Histo1D({ppHistoName2.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist3 = ppDF.Filter("geo==3 && tarCut3").Histo1D({ppHistoName3.Data(), "Real elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName1 = "real_geo1_2H";
	TString ddHistoName2 = "real_geo2_2H";
	TString ddHistoName3 = "real_geo3_2H";
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

void realDataAnalyzerDT()
{
	Double_t tarMass2H = cs::mass2H;

	parameters[sMWPC_1_X] = -1.0;
	parameters[sMWPC_1_Y] = -2.1375;
	parameters[sMWPC_2_X] = 0.2; 
	parameters[sMWPC_2_Y] = -1.125;
	parameters[sTarPos] = 10.0;
	parameters[sLang1] = 0.0;
	parameters[sLang2] = 0.0;
	parameters[sLang3] = 0.0;
	parameters[sRang] = 0.0;
	parameters[sDistL] = -20.0;
	parameters[sDistR] = -30.0;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos] - 2.0;

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + parameters[sDistL];
	sqrDist = 250.0 + parameters[sDistR];

	
	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small/small_dt.root";
	ROOT::RDataFrame smallDF("smallRealDT", smallFileName.Data());

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
						.Redefine("sqlang", [](TVector3 v2H, TVector3 vBeam){return v2H.Angle(vBeam)*TMath::RadToDeg();}, {"v2H", "vBeam"})
						.Redefine("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("dtCM", getCMAngleDT, {"sqrang", "lvBeam"})
						.Define("tarCut", "tarVertex.X()>vCut[mX2] && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]");

	TString dtHistoName = "real_dt";
	auto dtDF = newDF.Filter("tarCut").Filter("sumCut").Cache<Double_t, TVector3>({"dtCM", "tarVertex"});
	auto dtCMHist = dtDF.Filter("1").Histo1D({dtHistoName.Data(), "Real dt reaction", 90,0,180}, {"dtCM"});

	TFile histOutputFile("/home/zalewski/Desktop/6He/analysis/dataOut_dt/realDataCMHisto.root", "UPDATE");
	dtCMHist->Write(dtCMHist->GetName(), 1);
	histOutputFile.Close();
}

void MCDataAnalyzer()
{
	Double_t tarMass1H = cs::mass1H;
	Double_t tarMass2H = cs::mass2H;
	parameters[sMWPC_1_X] = -1.0;
	parameters[sMWPC_1_Y] = -2.1375;
	parameters[sMWPC_2_X] = 0.2; 
	parameters[sMWPC_2_Y] = -1.125;
	parameters[sTarPos] = 10.0;
	parameters[sLang1] = 0.0;
	parameters[sLang2] = 0.0;
	parameters[sLang3] = -1.0;
	parameters[sRang] = 0.0;
	parameters[sDistL] = -15.0;
	parameters[sDistR] = -25.0;

	tarPos[0] = parameters[sTarPos];
	tarPos[1] = parameters[sTarPos];
	tarPos[2] = parameters[sTarPos] + 2.0;

	tarAngle[0] = 45.0 * TMath::DegToRad();
	tarAngle[1] = 6.0 * TMath::DegToRad();
	tarAngle[2] = 0.0 * TMath::DegToRad();

	leftAngle[0] = (65.0 + parameters[sLang1]) * TMath::DegToRad();
	leftAngle[1] = (50.0 + parameters[sLang2]) * TMath::DegToRad();
	leftAngle[2] = (35.0 + parameters[sLang3]) * TMath::DegToRad();
	rightAngle = (15.0 + parameters[sRang]) * TMath::DegToRad();

	sqlDist = 170.0 + parameters[sDistL];
	sqrDist = 250.0 + parameters[sDistR];

	ROOT::EnableImplicitMT();
	TString smallFileName = "/home/zalewski/dataTmp/small/smallMC.root";
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
                        .Define("mc1H", getCMAngle1H, {"sqlang", "lvBeam"})
					    .Define("mc2H", getCMAngle2H, {"sqlang", "lvBeam"})
						.Define("sqrang", [](TVector3 v6He, TVector3 vBeam){return v6He.Angle(vBeam)*TMath::RadToDeg();}, {"v6He", "vBeam"})
						.Define("tarCut1", "tarVertex.X()>(vCut[mX1]) && tarVertex.X()<vCut[X1] && tarVertex.Y()>vCut[mY1] && tarVertex.Y()<vCut[Y1]")
						.Define("tarCut2", "tarVertex.X()>(vCut[mX2]) && tarVertex.X()<vCut[X2] && tarVertex.Y()>vCut[mY2] && tarVertex.Y()<vCut[Y2]")
						.Define("tarCut3", "tarVertex.X()>(vCut[mX3]) && tarVertex.X()<vCut[X3] && tarVertex.Y()>vCut[mY3] && tarVertex.Y()<vCut[Y3]");

	TString ppHistoName1 = "ppMC_geo1_CM";
	TString ppHistoName2 = "ppMC_geo2_CM";
	TString ppHistoName3 = "ppMC_geo3_CM";
	auto ppDF = newDF.Filter("mcPP").Cache<Double_t, TLorentzVector, Int_t, TVector3, Bool_t, Bool_t, Bool_t>({"mc1H", "lv2H_CM.", "geo", "tarVertex", "tarCut1", "tarCut2", "tarCut3"});
	auto ppCMHist1 = ppDF.Filter("geo==1 && tarCut1").Histo1D({ppHistoName1.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist2 = ppDF.Filter("geo==2 && tarCut2").Histo1D({ppHistoName2.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});
	auto ppCMHist3 = ppDF.Filter("geo==3 && tarCut3").Histo1D({ppHistoName3.Data(), "MC elastic scattering on protons", 90,0,180}, {"mc1H"});

	TString ddHistoName1 = "ddMC_geo1_CM";
	TString ddHistoName2 = "ddMC_geo2_CM";
	TString ddHistoName3 = "ddMC_geo3_CM";
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
	saveHistogram = true;
	loadCutsCorrectionParameters();

	//makeSmallFile();
	
	//MCDataAnalyzer();
	//realDataAnalyzerDT();
	calculateEfficiencyDT();
	calculateXsectionDT();
	makeGraphDT();



	//geo, mass

	//realDataAnalyzer();
/*
	calculateEfficiency(geo1, protium);
	calculateEfficiency(geo2, protium);
	calculateEfficiency(geo3, protium);

	calculateEfficiency(geo1, deuterium);
	calculateEfficiency(geo2, deuterium);
	calculateEfficiency(geo3, deuterium);


	calculateXsection(geo1, protium, 1e-3*1.484, 0.0);
	calculateXsection(geo2, protium, 1e-4*3.514, 0.0);
	calculateXsection(geo3, protium, 1e-4*1.392, 0.0);

	calculateXsection(geo1, deuterium, 1e-4*5.32, 0.0);
	calculateXsection(geo2, deuterium, 1e-4*1.514, 0.0);
	calculateXsection(geo3, deuterium, 1e-5*3.805, 0.0);

	
	makeGraph();
*/
/*
	getAngleError(geo1, protium);
	getAngleError(geo2, protium);
	getAngleError(geo3, protium);

	getAngleError(geo1, deuterium);
	getAngleError(geo2, deuterium);
	getAngleError(geo3, deuterium);

	getStatisticalError(geo1, protium);
	getStatisticalError(geo2, protium);
	getStatisticalError(geo3, protium);

	getStatisticalError(geo1, deuterium);
	getStatisticalError(geo2, deuterium);
	getStatisticalError(geo3, deuterium);
*/
}