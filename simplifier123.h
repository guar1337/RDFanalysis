//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct  4 13:35:28 2020 by ROOT version 6.22/00
// from TTree AnalysisxTree/Go4FileStore
// found on file: he6_01.root
//////////////////////////////////////////////////////////

#ifndef simplifier123_h
#define simplifier123_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class simplifier123 : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   TTree *outTree;
   TFile *outFile;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<unsigned short> NeEvent_CsI_L = {fReader, "NeEvent.CsI_L[16]"};
   TTreeReaderArray<unsigned short> NeEvent_tCsI_L = {fReader, "NeEvent.tCsI_L[16]"};
   TTreeReaderArray<unsigned short> NeEvent_CsI_R = {fReader, "NeEvent.CsI_R[16]"};
   TTreeReaderArray<unsigned short> NeEvent_tCsI_R = {fReader, "NeEvent.tCsI_R[16]"};
   TTreeReaderArray<unsigned short> NeEvent_SQX_L = {fReader, "NeEvent.SQX_L[32]"};
   TTreeReaderArray<unsigned short> NeEvent_tSQX_L = {fReader, "NeEvent.tSQX_L[32]"};
   TTreeReaderArray<unsigned short> NeEvent_SQY_L = {fReader, "NeEvent.SQY_L[16]"};
   TTreeReaderArray<unsigned short> NeEvent_tSQY_L = {fReader, "NeEvent.tSQY_L[16]"};
   TTreeReaderArray<unsigned short> NeEvent_SQX_R = {fReader, "NeEvent.SQX_R[32]"};
   TTreeReaderArray<unsigned short> NeEvent_tSQX_R = {fReader, "NeEvent.tSQX_R[32]"};
   TTreeReaderArray<unsigned short> NeEvent_SQY_R = {fReader, "NeEvent.SQY_R[16]"};
   TTreeReaderArray<unsigned short> NeEvent_SQ300 = {fReader, "NeEvent.SQ300[16]"};
   TTreeReaderArray<unsigned short> NeEvent_tSQ300 = {fReader, "NeEvent.tSQ300[16]"};
   TTreeReaderArray<unsigned short> NeEvent_tSQY_R = {fReader, "NeEvent.tSQY_R[16]"};
   TTreeReaderArray<unsigned short> NeEvent_F3 = {fReader, "NeEvent.F3[4]"};
   TTreeReaderArray<unsigned short> NeEvent_tF3 = {fReader, "NeEvent.tF3[4]"};
   TTreeReaderArray<unsigned short> NeEvent_F5 = {fReader, "NeEvent.F5[4]"};
   TTreeReaderArray<unsigned short> NeEvent_tF5 = {fReader, "NeEvent.tF5[4]"};
   TTreeReaderArray<unsigned short> NeEvent_F6 = {fReader, "NeEvent.F6[4]"};
   TTreeReaderArray<unsigned short> NeEvent_tF6 = {fReader, "NeEvent.tF6[4]"};
   TTreeReaderArray<unsigned short> NeEvent_tMWPC = {fReader, "NeEvent.tMWPC[4]"};
   TTreeReaderValue<unsigned short> NeEvent_nx1 = {fReader, "NeEvent.nx1"};
   TTreeReaderValue<unsigned short> NeEvent_ny1 = {fReader, "NeEvent.ny1"};
   TTreeReaderValue<unsigned short> NeEvent_nx2 = {fReader, "NeEvent.nx2"};
   TTreeReaderValue<unsigned short> NeEvent_ny2 = {fReader, "NeEvent.ny2"};
   TTreeReaderArray<unsigned short> NeEvent_x1 = {fReader, "NeEvent.x1[32]"};
   TTreeReaderArray<unsigned short> NeEvent_y1 = {fReader, "NeEvent.y1[32]"};
   TTreeReaderArray<unsigned short> NeEvent_x2 = {fReader, "NeEvent.x2[32]"};
   TTreeReaderArray<unsigned short> NeEvent_y2 = {fReader, "NeEvent.y2[32]"};
   TTreeReaderValue<Int_t> NeEvent_trigger = {fReader, "NeEvent.trigger"};
   

   simplifier123(TTree * /*tree*/ =0) { }
   virtual ~simplifier123() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   unsigned short nx1, ny1, nx2, ny2;
   Int_t trigger;

   unsigned short F3[4], F5[4], F6[4];
   unsigned short tF3[4], tF5[4], tF6[4];
   unsigned short tMWPC[4];

   unsigned short SQY_L[16], CsI_L[16], tSQY_L[16];
   unsigned short tCsI_L[16], SQY_R[16], tSQY_R[16];
	unsigned short CsI_R[16], tCsI_R[16];
   unsigned short SQ300[16], tSQ300[16];

   unsigned short x1[32], y1[32], x2[32], y2[32];
   unsigned short SQX_L[32], SQX_R[32];
   unsigned short tSQX_L[32], tSQX_R[32];

   ClassDef(simplifier123,0);

};

#endif

#ifdef simplifier123_cxx
void simplifier123::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t simplifier123::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef simplifier123_cxx
