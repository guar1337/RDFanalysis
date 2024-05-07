#define simplifier_cxx
// The class definition in simplifier.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//	 Begin():		  called every time a loop on the tree starts,
//						  a convenient place to create your histograms.
//	 SlaveBegin():	called after Begin(),	when on PROOF called only on the
//						  slave servers.
//	 Process():		called for each event,	in this function you decide what
//						  to read and fill your histograms.
//	 SlaveTerminate: called at the end of the loop on the tree,	when on PROOF
//						  called only on the slave servers.
//	 Terminate():	 called at the end of the loop on the tree,
//						  a convenient place to draw/fit your histograms.
//
// To use this file,	try the following session on your Tree T:
//
// root> T->Process("simplifier.C")
// root> T->Process("simplifier.C","some options")
// root> T->Process("simplifier.C+")
//


#include "simplifier.h"
#include <TH2.h>
#include <TStyle.h>

void simplifier::Begin(TTree *inTree)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString outFileName = GetOption();

	outFile = new TFile(outFileName.Data(),	"RECREATE");
	outTree = new TTree("simplified","simplified tree");
	outTree->Branch("SQX_L",	SQX_L,	"SQX_L[32]/s");
	outTree->Branch("SQY_L",	SQY_L,	"SQY_L[16]/s");
	outTree->Branch("CsI_L",	CsI_L,	"CsI_L[16]/s");
	outTree->Branch("tSQX_L",	tSQX_L,"tSQX_L[32]/s");
	outTree->Branch("tSQY_L",	tSQY_L,"tSQY_L[16]/s");
	outTree->Branch("tCsI_L",	tCsI_L,	"tCsI_L[16]/s");
	outTree->Branch("SQX_R",	SQX_R,	"SQX_R[32]/s");
	outTree->Branch("SQY_R",	SQY_R,	"SQY_R[16]/s");
	outTree->Branch("tSQX_R",	tSQX_R,"tSQX_R[32]/s");
	outTree->Branch("tSQY_R",	tSQY_R,"tSQY_R[16]/s");
	outTree->Branch("CsI_R",	CsI_R,	"CsI_R[16]/s");
	outTree->Branch("tCsI_R",	tCsI_R,	"tCsI_R[16]/s");
	outTree->Branch("SQ300",	SQ300,	"SQ300[16]/s");
	outTree->Branch("tSQ300",	tSQ300,	"tSQ300[16]/s");
	outTree->Branch("F3",		F3,		"F3[4]/s");
	outTree->Branch("tdcF3",	tF3,	"tdcF3[4]/s");
	outTree->Branch("F5",		F5,		"F5[4]/s");
	outTree->Branch("tdcF5",	tF5,	"tdcF5[4]/s");
	outTree->Branch("F6",		F6,		"F6[4]/s");
	outTree->Branch("tdcF6",	tF6,	"tdcF6[4]/s");
	outTree->Branch("tdcMWPC",	tMWPC, 	"tdcMWPC[4]/s");
	outTree->Branch("nx1",		&nx1,	"nx1/s");
	outTree->Branch("ny1",		&ny1,	"ny1/s");
	outTree->Branch("nx2",		&nx2,	"nx2/s");
	outTree->Branch("ny2",		&ny2,	"ny2/s");
	outTree->Branch("x1",		x1,		"x1[32]/s");
	outTree->Branch("y1",		y1,		"y1[32]/s");
	outTree->Branch("x2",		x2,		"x2[32]/s");
	outTree->Branch("y2",		y2,		"y2[32]/s");
	outTree->Branch("trigger",	&trigger,"trigger/I");
  
}

void simplifier::SlaveBegin(TTree * /*tree*/)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();

}

Bool_t simplifier::Process(Long64_t entry)
{
	// The Process() function is called for each entry in the tree (or possibly
	// keyed object in the case of PROOF) to be processed. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	// When processing keyed objects with PROOF,	the object is already loaded
	// and is available via the fObject pointer.
	//
	// This function should contain the \"body\" of the analysis. It can contain
	// simple or elaborate selection criteria,	run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	//
	// Use fStatus to set the return value of TTree::Process().
	//
	// The return value is currently not used.

	fReader.SetLocalEntry(entry);


   nx1 = *NeEvent_nx1;
   ny1 = *NeEvent_ny1;
   nx2 = *NeEvent_nx2;
   ny2 = *NeEvent_ny2;
   trigger = *NeEvent_trigger;

	for (int iii = 0; iii < 4; iii++)
	{
		F3[iii]		= NeEvent_F3[iii];
		F5[iii]		= NeEvent_F5[iii];
		F6[iii]		= NeEvent_F6[iii];
		tF3[iii]		= NeEvent_tF3[iii];
		tF5[iii]		= NeEvent_tF5[iii];
		tF6[iii]		= NeEvent_tF6[iii];
		tMWPC[iii]	= NeEvent_tMWPC[iii];
	}
	

	for (int iii = 0; iii < 16; iii++)
	{
		SQY_L[iii]	=	NeEvent_SQY_L[iii];
		CsI_L[iii]	=	NeEvent_CsI_L[iii];
		tSQY_L[iii]	=	NeEvent_tSQY_L[iii];
		tCsI_L[iii]	=	NeEvent_tCsI_L[iii];
		SQY_R[iii]	=	NeEvent_SQY_R[iii];
		tSQY_R[iii]	=	NeEvent_tSQY_R[iii];
		CsI_R[iii]	=	NeEvent_CsI_R[iii];
		tCsI_R[iii]	=	NeEvent_tCsI_R[iii];
		SQ300[iii]	=	0;//NeEvent_SQ300[iii];
		tSQ300[iii]	=	0;//NeEvent_tSQ300[iii];
	}
	

	for (int iii = 0; iii < 32; iii++)
	{
		SQX_L[iii] = NeEvent_SQX_L[iii];
		SQX_R[iii] = NeEvent_SQX_R[iii];
		tSQX_L[iii] = NeEvent_tSQX_L[iii];
		tSQX_R[iii] = NeEvent_tSQX_R[iii];
		x1[iii] = NeEvent_x1[iii];
		y1[iii] = NeEvent_y1[iii];
		x2[iii] = NeEvent_x2[iii];
		y2[iii] = NeEvent_y2[iii];
	}
	outTree->Fill();

	return kTRUE;
}

void simplifier::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
}

void simplifier::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client,	it can be used to presenttranslator
	delete outFile;

}