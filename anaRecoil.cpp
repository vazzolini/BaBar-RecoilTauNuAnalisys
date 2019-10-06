//
//  $Id: anaRecoil.cpp,v 1.4 2003/02/13 18:49:12 cartaro Exp $
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "recoilTaunu.hh"


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  TString progName  = argv[0]; progName.ReplaceAll("bin/", "");
  TString writeName, fileName, dirName, cutFile("");

  TString treeName("h1");
  Int_t dump(0), write(0);
  Int_t file(0), lun(0);
  Int_t isMC(0), isVerbose(0), switchOff(0);
  Int_t dirspec(0);
  Int_t nevents(0), skipevents(0);
  Int_t blind(0), runGammas(0), runKlongs(0), runCategories(0),runSplitOff(0), filterK0s(-1), newFormat(1); 
  Int_t pidKilling(0), pidKillingKaon(0), pidKillingEl(0), pidKillingMu(0);
  Int_t isTauNu(0), isVcb(0), isDstar(0);
  Int_t doUnBrem(0), vDotagstudy(0), generatemodelist(0), trackcutstudy(0);
  Int_t smearTracks(0);
  Int_t smearNeut(0);
  Int_t oneProng(0);
  Int_t randomSeed(processID);
  Bool_t selectB0( false );
  Bool_t smearNeutrals( false );
  Bool_t genericBBbar ( false );

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-b")  == 0) {blind      = 1;}                                             // skip events with Mx(fit) < 1.6 GeV
    if(strcmp(argv[i],"-C")  == 0) {cutFile    = TString(argv[++i]); }                           // file with cuts
    if(strcmp(argv[i],"-c")  == 0) {fileName   = TString(argv[++i]); file = 0; }                 // file with chain definition
    if(strcmp(argv[i],"-D")  == 0) {dirName     = TString(argv[++i]);  dirspec = 1; }            // where to put the output
    if(strcmp(argv[i],"-d")  == 0) dump        = atoi(argv[++i]);                                // printout
    if(strcmp(argv[i],"-Dstar")  == 0) isDstar = 1;
    if(strcmp(argv[i],"-f")  == 0) {fileName   = TString(argv[++i]); file = 1; }                 // single file instead of chain
    if(strcmp(argv[i],"-F")  == 0) {filterK0s  = atoi(argv[++i]);}                               // require 0 or x K0S->pi0pi0 decays
    if(strcmp(argv[i],"-g")  == 0) {runGammas  = 1;}                                             // run gamma study
    if(strcmp(argv[i],"-K")  == 0) {runKlongs  = 1;}                                             // run Klong study
    if(strcmp(argv[i],"-l")  == 0) {runCategories = 1;}                                          // run Categories
    if(strcmp(argv[i],"-MC") == 0)  isMC = 1;                                                    // flag to avoid warnings
    if(strcmp(argv[i],"-mc") == 0)  isMC = 1;                                                    // flag to avoid warnings
    if(strcmp(argv[i],"-n")  == 0) nevents     = atoi(argv[++i]);                                // there is still a small bug in here
    if(strcmp(argv[i],"-p")  == 0) pidKilling  = 1;                                              // redo PidKilling on tuple (overall)
    if(strcmp(argv[i],"-pe")  == 0) pidKillingEl  = 1;                                           // redo PidKilling on tuple for el
    if(strcmp(argv[i],"-pk")  == 0) pidKillingKaon  = 1;                                         // redo PidKilling on tuple for kaons
    if(strcmp(argv[i],"-pm")  == 0) pidKillingMu  = 1;                                           // redo PidKilling on tuple for muons
    if(strcmp(argv[i],"-r")  == 0) randomSeed  = atoi(argv[++i]);                                // set seed for random number generator
    if(strcmp(argv[i],"-TAUNU")  == 0) isTauNu = 1;
    if(strcmp(argv[i],"-B0")  == 0) selectB0 = true;

    if(strcmp(argv[i],"-vcb")  == 0) isVcb = 1;
    if(strcmp(argv[i],"-vcbtagstudy")  == 0) vDotagstudy = 1;
    if(strcmp(argv[i],"-unbrem")  == 0) doUnBrem = 1;
    if(strcmp(argv[i],"-modelist")  == 0) generatemodelist = 1;
    if(strcmp(argv[i],"-trackcutstudy") == 0) trackcutstudy = 1;

    if(strcmp(argv[i],"-t")  == 0) treeName    = TString(argv[++i]);                             // ...
    if(strcmp(argv[i],"-ts") == 0)  lun = atoi(argv[++i]);                                       // timestamp business
    if(strcmp(argv[i],"-s")  == 0) skipevents  = atoi(argv[++i]);                                // skip events at the beginning
    if(strcmp(argv[i],"-ST") == 0) smearTracks = 1;                                              // smear tracks
    if(strcmp(argv[i],"-SN") == 0) smearNeutrals = 1;                                            // smear neut
    if(strcmp(argv[i],"-intpur") == 0) genericBBbar = 1;                                        // int purity histos dump

    if(strcmp(argv[i],"-SOff")  == 0) {runSplitOff = 1;}                                         // study splitoffs
    if(strcmp(argv[i],"-v")  == 0)  isVerbose = 1;                                               // be verbose
    if(strcmp(argv[i],"-w")  == 0) {writeName = TString(argv[++i]); write = 1;}                  // write event list into red. file
    if(strcmp(argv[i],"-y")  == 0)  newFormat = atoi(argv[++i]);                                       // flag to indicate the tree version
                                                                                                 // 0: Feb02, 1: apr 02, 2: aug02
    if(strcmp(argv[i],"-z")  == 0)  switchOff = 1;                                               // flag to speed up reading of tree
    if(strcmp(argv[i],"-1P")  == 0) oneProng  = 1;                                               // flag for one prong tau decays
  }

  cout << "Setting gRandom seed: " << randomSeed << endl;
  gRandom->SetSeed(processID);


  // -- Rootfile with output histograms
  TString  barefile(fileName), chainFile, meta, histfile;
  if (file == 0) {
    if (barefile.Contains("chains/")) {
      meta = barefile;
      //      histfile = "outp/" + barefile.ReplaceAll("chains/", "") + ".root";
      histfile = barefile.ReplaceAll("chains/", "") + ".root";
      if(dirspec) histfile = dirName + "/" + histfile;
    } else {
      //      meta = "chains/" + barefile;
      meta = barefile;
      histfile =  barefile + ".root";
      if(dirspec) histfile = dirName + "/" + histfile;
    }

    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla + ".root";
    if(dirspec) histfile = dirName + "/" + histfile;
  }  else if (file == 1) {
    // -- The following lines strip everything from the string up to and including the last '/'
    int fl = barefile.Last('/');
    TString bla(barefile);
    bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
    histfile =  bla;
    if(dirspec) histfile = dirName + "/" + histfile;
  }

  // use the following line to dump the output into some other place: "scratch/anaRecoil/bla.root"
  TString histFile = histfile;
  cout << "Opening " << histFile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.Data() << " for input" << endl;


  // -- Set up chain
  TChain *chain = new TChain(treeName);
  cout << "Chaining ... " << treeName << endl;
  if (file == 0) {
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      cout << meta << endl; 
      chain->Add(meta); 
    }
    is.close();
  }
  else if (file == 1) {
    cout << fileName << endl;
    chain->Add(fileName);
  }

  recoilTaunu a(chain,isMC, newFormat, selectB0,smearNeutrals,genericBBbar);
  if (pidKilling) {
    // do nothing 
  }
    
    if (strcmp(cutFile.Data(), "")) a.readCuts(cutFile, 0);
    a.dumpCuts();
    if (!(write)) a.openHistFile(histFile);
    
    // if (write) a.makeEventList(1); 
    
    if (!(write)) a.bookHist(dump);
    cout << "====================" << endl;
    cout << "isVerbose      = " << isVerbose << endl;
    cout << "runGammas      = " << runGammas << endl;
    cout << "runKlongs      = " << runKlongs << endl;
    cout << "PidKilling     = " << pidKilling << endl;
    cout << "PidKillingKaon = " << pidKillingKaon << endl;
    cout << "TrkSmearing    = " << smearTracks << endl;
    cout << "NeuSmearing    = " << smearNeut << endl;
    cout << "Filter K0S     = " << filterK0s << endl;
    cout << "categories     = " << runCategories << endl;
    cout << "splitoff       = " << runSplitOff   << endl;
    cout << "blind          = " << blind << endl;
    cout << "write          = " << write << endl;
    cout << "====================" << endl;

    
    if (!(write)) a.Loop(nevents, skipevents, isVerbose, lun);
    if (write) {
      //  a.Skim(.8,nevents, skipevents,isVerbose,"");
      //a.dumpEventList(writeName.Data());
    }
    if (!(write)) a.closeHistFile();
  
  return 0;

}

