int MCCutMaker()
{
   TFile mcCuts("/home/zalewski/aku/analysis/mcCuts.root", "RECREATE");

   TCutG *mcPP = new TCutG("mcPP",13);
   mcPP->SetVarX("sqrang");
   mcPP->SetVarY("sqlang");
   mcPP->SetTitle("Graph");
   mcPP->SetFillStyle(1000);
   mcPP->SetPoint(0,9.44293,22.4098);
   mcPP->SetPoint(1,11.2545,31.7516);
   mcPP->SetPoint(2,10.4242,51.1783);
   mcPP->SetPoint(3,7.72569,65.7219);
   mcPP->SetPoint(4,3.30993,81.2208);
   mcPP->SetPoint(5,0.422705,88.4395);
   mcPP->SetPoint(6,0.517059,84.724);
   mcPP->SetPoint(7,4.29121,69.6497);
   mcPP->SetPoint(8,7.53699,54.2569);
   mcPP->SetPoint(9,8.51827,39.1826);
   mcPP->SetPoint(10,7.76344,26.9745);
   mcPP->SetPoint(11,7.49925,22.4098);
   mcPP->SetPoint(12,9.44293,22.4098);

   TCutG *mcDD = new TCutG("mcDD",22);
   mcDD->SetVarX("sqrang");
   mcDD->SetVarY("sqlang");
   mcDD->SetTitle("Graph");
   mcDD->SetFillStyle(1000);
   mcDD->SetPoint(0,17.2177,17.6327);
   mcDD->SetPoint(1,19.9728,25.0637);
   mcDD->SetPoint(2,21.0862,34.8301);
   mcDD->SetPoint(3,19.1048,47.4628);
   mcDD->SetPoint(4,14.5758,62.8556);
   mcDD->SetPoint(5,4.42331,84.4055);
   mcDD->SetPoint(6,1.93237,88.2272);
   mcDD->SetPoint(7,0.403834,90.7749);
   mcDD->SetPoint(8,0.120773,88.9703);
   mcDD->SetPoint(9,0.215127,87.9087);
   mcDD->SetPoint(10,0.460447,88.6518);
   mcDD->SetPoint(11,1.51721,86.2102);
   mcDD->SetPoint(12,3.10235,82.4947);
   mcDD->SetPoint(13,5.93297,75.1699);
   mcDD->SetPoint(14,12.1037,60.9448);
   mcDD->SetPoint(15,16.3496,48.5244);
   mcDD->SetPoint(16,18.3877,38.8641);
   mcDD->SetPoint(17,18.482,33.8747);
   mcDD->SetPoint(18,17.5008,27.8238);
   mcDD->SetPoint(19,16.6138,24.5329);
   mcDD->SetPoint(20,15.6514,19.6497);
   mcDD->SetPoint(21,17.2177,17.6327);

   TCutG *mcHe6 = new TCutG("mcHe6",20);
   mcHe6->SetVarX("sqretot");
   mcHe6->SetVarY("sqrde");
   mcHe6->SetTitle("Graph");
   mcHe6->SetFillStyle(1000);
   mcHe6->SetPoint(0,173.45,14.8639);
   mcHe6->SetPoint(1,152.511,16.6194);
   mcHe6->SetPoint(2,94.4263,22.4713);
   mcHe6->SetPoint(3,78.3635,24.8957);
   mcHe6->SetPoint(4,67.7506,26.9857);
   mcHe6->SetPoint(5,49.3931,31.4164);
   mcHe6->SetPoint(6,30.1751,38.8567);
   mcHe6->SetPoint(7,17.1241,46.297);
   mcHe6->SetPoint(8,3.78623,55.9944);
   mcHe6->SetPoint(9,1.77838,54.824);
   mcHe6->SetPoint(10,9.37953,48.1361);
   mcHe6->SetPoint(11,24.8687,38.7731);
   mcHe6->SetPoint(12,37.4894,33.0048);
   mcHe6->SetPoint(13,52.9786,27.6545);
   mcHe6->SetPoint(14,75.3518,22.6385);
   mcHe6->SetPoint(15,97.4381,18.2914);
   mcHe6->SetPoint(16,120.098,15.6998);
   mcHe6->SetPoint(17,150.79,13.4427);
   mcHe6->SetPoint(18,172.302,12.5231);
   mcHe6->SetPoint(19,173.45,14.8639);


   mcPP->Write();
   mcDD->Write();
   mcHe6->Write();

   mcCuts.Write();
   mcCuts.Close();
   return 5;
}