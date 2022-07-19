#include <vector>

int counter = 1;
int ID=0;
TCanvas myCanv("boxB","boxB long name",10,10,4*297,4*210);

void drawIV(std::vector<int> current, std::string boxID, std::string detNo, std::string comment, int detClass)
{
	
	std::vector<int> diodeVoltage =	{10, 20, 50, 80, 100, 120};
	TGraph ivGraph(diodeVoltage.size(), &diodeVoltage[0], &current[0]);
	ivGraph.GetXaxis()->SetLimits(0.0, 140.0);
	ivGraph.GetYaxis()->SetRangeUser(0.0, 1.1*current.at(5));
	ivGraph.GetXaxis()->SetTitle("voltage [V]");
	ivGraph.GetYaxis()->SetTitle("current [nA]");
	TString graphTitle = TString::Format("Box %s det %s class: %d", boxID.c_str(), detNo.c_str(), detClass);

	ivGraph.SetTitle(graphTitle.Data());
	ivGraph.SetMarkerStyle(55);
	ivGraph.SetMarkerSize(2);
	ivGraph.SetMarkerColor(kRed);
	ivGraph.GetYaxis()->SetTitleOffset(-1.0);
	ivGraph.GetYaxis()->CenterTitle();

	myCanv.cd(counter);
	printf("%d\n",counter);
	/*(counter==1) ? ivGraph.Draw("APL") : */ivGraph.DrawClone("APL,same");
//	myCanv.Modified();
//	myCanv.Update();

	if (counter==6 || (counter==5 && ID==10))
	{
		TString picName = TString::Format("/home/zalewski/Desktop/PIN/w%d.png", ID++);
		myCanv.Print(picName.Data());
		counter=1;
	}

	else
	{
		counter++;
	}

}

int readDiodeCurrent()
{
	std::string parPath = "/home/zalewski/Desktop/pin.csv";
	std::string boxID, detNo, comment;
	int detClass;

	std::ifstream in(parPath);
	std::vector<int> current(6);

	while (true)
	{	
		in >> boxID >> detNo >> detClass >> current.at(0) >> current.at(1) >> current.at(2) >> current.at(3) >> current.at(4) >> current.at(5) >> comment;
		if( in.eof() ) break;
		std::replace(comment.begin(), comment.end(), '-',' '); // replace all '-' to ' '
		printf("Box: %s\tdet: %s\t%d\t%d\t%s\n", boxID.c_str(), detNo.c_str(), current.at(0), current.at(5), comment.c_str());
		drawIV(current, boxID, detNo, comment, detClass);
	}

	return 1;
}

void clipboard()
{
	myCanv.Divide(3,2);
	//drawIV();
	readDiodeCurrent();
}