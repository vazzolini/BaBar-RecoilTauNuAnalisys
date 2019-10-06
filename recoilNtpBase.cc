#include "recoilNtpBase.hh"

#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include "TH2.h"
#include "TH1.h"
#include <TMap.h>
#include <TObjString.h>
#include "TLorentzVector.h"
#include "TRandom.h"

// ----------------------------------------------------------------------
recoilNtpBase::recoilNtpBase(TTree *tree, int isMC, int newFormat) {
  fNewFormat = newFormat; 
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    if (!f) {
      f = new TFile("/u/ec/ursl/d/output/breco-110201/mc-b0cocktail-20.root");
    }
    tree = (TTree*)gDirectory->Get("h9");
    
  }
  Init(tree,isMC);
  fToBeCopied = new TEventList("toBeCopied", "Events to be copied", 1000);
  initRest();
}

// ----------------------------------------------------------------------
recoilNtpBase::~recoilNtpBase() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


// ----------------------------------------------------------------------
recoilNtpBase::recoilNtpBase(TString filename, TString treename,int isMC, int newFormat) {
  fNewFormat = newFormat;
  TFile *f = new TFile(filename);
  TTree *tree = (TTree*)f->Get(treename);
  if (!tree) { 
    cout << "Did not find " << treename << " in file " << filename << endl;
    f->ls();
  } else {
    Init(tree,isMC);
  }
  initRest();
}


// ----------------------------------------------------------------------
Bool_t recoilNtpBase::isAncestor(int ancestor, int cand) {
  int mom(cand); 
  while (mom > 0) {
    //    cout<<TMath::Abs(idMc[mom])<<endl;
    mom = mothMc[mom]-1; 
    //    cout << "   now looking at " << mom << " which is a " << idMc[mom] << endl;
    if (mom == ancestor) return true;
  }
  return false;
}

// ----------------------------------------------------------------------
Int_t recoilNtpBase::isRecoed(int imc) {
  int result(-1), i(0); 
  for (i = 0; i < nTrk; ++i) {
    if (imc == (IndexTrk[i]-1)) {
      result = i;
      break;
    }
  }
  if (result > -1) return result;
  for (i = 0; i < nGam; ++i) {
    if (imc == (IndexGam[i]-1)) {
      result = i;
      break;
    }
  }
  return result;
}

// ----------------------------------------------------------------------
void  recoilNtpBase::printLorentz(const TLorentzVector &p) {
  char line[200]; 
  sprintf(line, "p = (%7.5f, %7.5f, %7.5f, %7.5f), m = %7.5f",  p.X(), p.Y(), p.Z(), p.E(), p.Mag()); 
  cout << line;
}

// ----------------------------------------------------------------------

// #include "util.icc"

void recoilNtpBase::timestamp(const char * lun) {

  //  This function prints out eventIDs and run numbers in the official
  //  BaBar format (e.g. eventID in hex as aa:bb:xxxxxx/yyyyyyyy:Z, 
  //  run number as an 8-digit decimal number following "run=" token
  //  Output is written to LUN <lun>

  int k;
  int temp;
  int iplat;
  Int_t ipart;
  char cplat[2];
  char cpart[8],tspart[8];
  int nplat,npart;
  
  const char * array[17]={"G","H","J","K","L","M","N","P","Q","R","S",
			  "T","U","V","W","X","Y"};

  if(lun == "") return;

  iplat=platform;
  ipart = partition;

  //      iplat=127
  //      partition=8388607             ! 0x7fffff

  // implement wildcard behavior of platform
  if ( iplat == 0 ) {
    sprintf(cplat,"%s","\0");
    nplat=1;
  }  else {
    sprintf(cplat,"%X%s",iplat,"\0");
    nplat=2;
  }

  // implement wildcard behavior of partitionMask

  if ( partition == 0 ) {
    sprintf(cpart,"%s","*\0");
    sprintf(tspart,"%s","0\0");
    npart=1;
    //   Need to redo twos complement interpretation (MC has negative partition id)
    //  i.e. map -2 ->FFFFFFFE , -1 -> FFFFFFFF, etc.
  } else if (partition < 0) {
    temp = (Int_t) (pow(2,32) + partition);
    sprintf(cpart,"%X%s", temp,"\0");
    sprintf(tspart,"%X%s", temp,"\0");
    npart=8;
  }  else {
    sprintf(cpart,"%X", partition);
    sprintf(tspart,"%X", partition);
    npart=6;
  }

  // checksum
  k= (TMath::Abs(iplat+ipart+upperID+lowerID)%17) +1;

  // now print
  //  if(fIntPurity < 30) cout<<" Ts file must be set up "<<fIntPurity<<endl;

  char name[100];
  char name2[100];

  sprintf(name,"%s",lun);
  sprintf(name2,"%s%s",lun,"MI");

  //  cout<<name<<" "<<name2<<endl;

  ofstream outFile(name,ios::app);
  ofstream outFile2(name2,ios::app);
  outFile.setf(ios::uppercase);
  outFile2.setf(ios::uppercase);
  int pur(0);// tmp fix
  outFile<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<cpart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<endl;
  outFile2<<hex<<setw(nplat)<<setfill('0')<<cplat<<":"<<setw(npart)
	 <<setfill('0')<<tspart<<":"<<setw(8)<<setfill('0')<<upperID
	 <<"/"<<setw(8)<<setfill('0')<<lowerID<<":"<<array[k-1]<<"  run="<<dec<<setw(8)
	 <<setfill('0')<<runNumber<<"  "<<dec<<pur<<" "
	 <<fBmode<<" "<<endl;
}


// ----------------------------------------------------------------------
bool recoilNtpBase::checkPurity(int chbcand) {

  if (fVerbose) cout << "Checking purity" << endl;

  Int_t Brectrktmp, Brecgamtmp;
  Int_t nB(0);
  double tmppurB(0.), tmpIpurB(0.);
  bool check = true;  
  //Assign the correct Lund to the B reco

  Int_t tmpblund;
  tmpblund = B0LUND;
  nB = nB0;
  int modeB = modeB0[indexbestB];
  if(bestB0==0) {
    modeB = modeChB[indexbestB]; 
    tmpblund = CHBLUND;
    nB = nChB;
  }  
  tmpIpurB = brecointpur[modeB-10000];
  tmppurB = brecosig[modeB-10000]/(brecosig[modeB-10000]+brecobkg[modeB-10000]);
  
  fPurity = tmppurB;
  fIntPurity = tmpIpurB;

  if (tmppurB < 0.1) {
    if (fReturnLog[1] == 0) { 
      fReturnString[1] = TString(Form("purity too low"));
      return false;
    }
  }
  
  if (fVerbose) cout << "Survived cut on breco purity" << endl;

  if ((fSeedMode == 0) && (tmppurB < IPURDC)) {
    if (fReturnLog[2] == 0) fReturnString[2] = TString(Form("Dc, ipur too low"));
    fReturnLog[2]++;
    return false;
  }

  if ((fSeedMode == 1) && (tmppurB < IPURDSTAR)) {
    if (fReturnLog[3] == 0) fReturnString[3] = TString(Form("D*, ipur too low"));
    fReturnLog[3]++;
    return false;
  }
  if ((fSeedMode == 2) && (tmppurB < IPURD0)) {
    if (fReturnLog[4] == 0) fReturnString[4] = TString(Form("D0, ipur too low"));
    fReturnLog[4]++;
    return false;
  }
  if ((fSeedMode == 3) && (tmppurB < IPURDSTAR0)) {
    if (fReturnLog[5] == 0) fReturnString[5] = TString(Form("D*0, ipur too low"));
    fReturnLog[5]++;
    return false;
  }

  if (fVerbose) cout << "Survived cuts on breco integrated purities" << endl;
  return true;
} 
 
// ----------------------------------------------------------------------
void recoilNtpBase::dumpOneB(int b1) {
  char line[200];
  for (int i = 0; i < nMc; ++i) {
    if ((i == b1) || isAncestor(b1, i)) {
      sprintf(line, "%3d %+6d mom(%3d) ndau(%3d) p=%5.3f, t=%5.3f f=%+5.3f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, nDauMc[i],
	      pMc[i], thetaMc[i], phiMc[i],
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  }
}

// ----------------------------------------------------------------------
void recoilNtpBase::dumpGeneratorBlock(int b1, int b2) {
  char line[200];
  if (b1 < 0) {
    for (int i = 0; i < nMc; ++i) {
      sprintf(line, "%3d %+6d mom(%3d) ndau(%3d) p=%5.3f, t=%5.3f f=%+5.3f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, nDauMc[i],
	      pMc[i], thetaMc[i], phiMc[i],
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  } else {
    int overlap(-1);
    char recoed[2];
    for (int i = 0; i < nMc; ++i) {
      if (isAncestor(b1, i)) {
	overlap = 1;
      } else if (isAncestor(b2, i)) {
	overlap = 2;
      } else {
	overlap = 0;
      }
      if (isRecoed(i) > -1) {
	sprintf(recoed, "r");
      } else {
	sprintf(recoed, " ");
      }
      sprintf(line, "%3d %+6d mom(%3d) flag(%1d%s) p=%5.3f, t=%5.3f f=%+5.3f v=(%+7.3f,%+7.3f,%+7.3f)", 
	      i, idMc[i], mothMc[i]-1, overlap, recoed,
	      pMc[i], thetaMc[i], phiMc[i],
	      xMc[i], yMc[i], zMc[i]);
      cout << line << endl;
    }
  }
}

// ----------------------------------------------------------------------
void recoilNtpBase::Loop(Int_t maxEvent, Int_t startEvent, Int_t isVerbose, Int_t lun) {
  int step(1000);
  int ischB(0);
  
  fVerbose = isVerbose; 

  double tmpMassPB, tmpMassThetaB, tmpMassPhiB ;
  double tmpPB, tmpThetaB, tmpPhiB ;
  double tmpPgen, tmpThetagen, tmpPhigen, tmpMassgen ;
  double tmpMB, tmpBevM;
  //


  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  if (startEvent > 0) {
    cout << "Will start at event " << startEvent << endl;
    if (startEvent+maxEvent >  nentries) {
      cout << "Requested " << maxEvent << " events, but will run only to end of chain"  << endl;
      maxEvent = nentries - startEvent; 
    }
  }

  Int_t nbytes = 0, nb = 0;
  int nvxbevt(0); 
  int nk0sEvents(0); 
  const char *pChar; 

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
    if (fReturnLog[0] == 0) fReturnString[0] = TString("Loop event counter");
    fReturnLog[0]++;
    if (fVerbose) cout << "->  new event " << endl;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) {
      cout << "File " << fChain->GetCurrentFile()->GetName(); 
      fFileChanged = 1; // file has changed
      if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2000")) {
	fRunRange = TString("Run 1"); 
	cout << " Run 1: 2000" << endl;
      } else if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2001")) {
	fRunRange = TString("Run 2a"); 
	cout << " Run 2a: 2001" << endl;
      } else if (pChar = strstr(fChain->GetCurrentFile()->GetName(), "2002")) {
	fRunRange = TString("Run 2b"); 
	cout << " Run 2b: 2002" << endl;
      }
      else {
	fRunRange = TString("undefined"); 
	cout << " Runrange ?" << endl;
      }
    } else {
      fFileChanged = 0;  // staying in same file
    }

    

    // -- Initialize event
    initVariables();
    findbestB();

    
    brecoOverlap = 1; // this is flag for overlap with BRECO candidate
    if(indexbestB == 1) brecoOverlap = 2;

    if (skipBadBreco()) {
      if (fReturnLog[6] == 0) fReturnString[6] = TString("Bad B-D flavor correlation");
      fReturnLog[6]++;
      continue;
    }
    processEvent();
  }
}

void recoilNtpBase::processEvent()
{
  // selectTracks();
  // doSplitOffStudy();
  // selectPhotons();

  int ischB = 0;
  if (bestB0 == 0) {
    ischB = 1;
  }

  // -- MonteCarlo Truth
  if (fIsMC) {
    mcTruth(ischB);
  }
  Double_t tmpMassPB, tmpMassThetaB, tmpMassPhiB;
  Double_t tmpPB, tmpBevM, tmpThetaB, tmpPhiB, tmpMB;
  // -- Reco quantities
  if(ischB == 0) {
    tmpMassPB =MassPB0[indexbestB];
    tmpMassThetaB =MassThetaB0[indexbestB];
    tmpMassPhiB =MassPhiB0[indexbestB];
    tmpPB = pB0[indexbestB];
    tmpBevM = massB0[indexbestB];
    tmpThetaB =thB0[indexbestB];
    tmpPhiB =phiB0[indexbestB];
    tmpMB = BZMASS;
  } else {
    tmpMassPB =MassPChB[indexbestB];
    tmpMassThetaB =MassThetaChB[indexbestB];
    tmpMassPhiB =MassPhiChB[indexbestB];
    tmpPB = pChB[indexbestB];
    tmpBevM = massChB[indexbestB];
    tmpThetaB =thChB[indexbestB];
    tmpPhiB =phiChB[indexbestB];
    tmpMB = BPMASS;
  }
  mk4Vector(p4Breco, tmpMassPB , tmpMassThetaB, tmpMassPhiB, tmpMB); 
  mk4Vector(p4BrecoNC, tmpPB , tmpThetaB, tmpPhiB, tmpBevM); 
  
  p4Upsilon = TLorentzVector(pxUps, pyUps, pzUps, eUps); 
  upsBoost = TVector3(pxUps, pyUps, pzUps);
  upsBoost.SetMag(upsBoost.Mag()/eUps);

//  Changed sign of p4Brecoil with respect to original
  p4Brecoil = p4Upsilon - p4Breco; 
  cmsBoost = p4Brecoil.Vect();
  cmsBoost.SetMag(cmsBoost.Mag()/p4Brecoil.E());
  
  TLorentzVector p4t = p4Breco; 
  p4t.Boost(-upsBoost);
  fPcmsBreco = p4t.Vect().Mag(); 

  
  if (fVerbose) cout << "CALLING breco()" << endl;
  breco(ischB);
  
  bool pureEnough = checkPurity( ischB );
  if ( pureEnough ) {
    if (fVerbose) cout << "CALLING recoil()" << endl;
    recoil(ischB);
  }
  
}

/*
// ----------------------------------------------------------------------
void recoilNtpBase::LoopKill(const char *Ifile, Int_t maxEvent, Int_t startEvent, Int_t isVerbose, const char *ITSfile) {

  int step(1000);
  fDupcount = 0;
  
  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntries());
  if (maxEvent == 0) maxEvent = nentries;
  if (nentries < 1) {
    cout << "Found no entries in " << fChain->GetName() << endl;
  } else {
    cout << "Found " << nentries << " entries in tree " << fChain->GetName() << endl;
  }

  if (startEvent > 0) {
    cout << "Will start at event " << startEvent << endl;
    if (startEvent+maxEvent >  nentries) {
      cout << "Requested " << maxEvent << " events, but will run only to end of chain"  << endl;
      maxEvent = nentries - startEvent; 
    }
  }

  Int_t nbytes = 0, nb = 0;

  if(strlen(Ifile) > 0)  read_killTab(Ifile);

  for (Int_t jentry = startEvent; jentry < startEvent+maxEvent; jentry++) {
    if (isVerbose) cout << "->  new event " << endl;
    fEvent = jentry;
    // in case of a TChain, ientry is the entry number in the current file
    tsdump = kFALSE;
    Int_t ientry = LoadTree(jentry); 
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%step == 0) cout << "Event " << jentry << endl;
    if (ientry == 0) cout << "File " << fChain->GetCurrentFile()->GetName() << endl;

    fisDuplicate = kFALSE;
    cout<<upperID<< " Up ; Low "<<lowerID<<" In loopkill"<<endl;
    if(ITSfile != "") timestamp(ITSfile);

    //    if(strlen(Ifile) > 0)  fisDuplicate = kill_dupli(upper,lower);

    if (fOptMakeEventList) {
      if (!fisDuplicate) {
	//	cout << "--> Copying " << jentry << endl;
	fToBeCopied->Enter(jentry);
      }
    }

    fIVal = 0;

    const char *Arstr[64] = { 
      "1", "2" ,"3", "4",       
      "11", "12" ,"13", "14",   
      "21", "22" ,"23", "24",   
      "31", "32" ,"33", "34",
      "41", "42" ,"43", "44",
      "51", "52" ,"53", "54",
      "61", "62" ,"63", "64",
      "71", "72" ,"73", "74",
      "81", "82" ,"83", "84",
      "91", "92" ,"93", "94",
      "101", "102" ,"103", "104",
      "111", "112" ,"113", "114",
      "121", "122" ,"123", "124",
      "131", "132" ,"133", "134",
      "141", "142" ,"143", "144",
      "151", "152" ,"153", "154" };
    
    if(fisDuplicate) {
      cout<<" Duplicated event!!!!   Recoil part will not be computed"<<endl;
      for (int iVal= 0; iVal < 64; iVal++) {
	//	cout<<Arstr[iVal]<<" "<<fValMap<<endl;
	if (strcmp(fValMap,Arstr[iVal]) == 0) {
	  fIVal = iVal + 10;
	}
      }
      fDupcount++;
    }


  }
  
  if(Ifile != " ") cout<<fDupcount<<" Num. of duplicates"<<endl;

}
*/ 
// ----------------------------------------------------------------------
TFile* recoilNtpBase::openHistFile(TString name) {
  cout << name << endl;
  fHistFile = new TFile(name.Data(), "RECREATE");
  fHistFile->cd();
  cout << "Opened " << fHistFile->GetName() << endl;
  return fHistFile;
}


// ----------------------------------------------------------------------
void recoilNtpBase::closeHistFile() {
  cout << "Return Log: " << endl;
  for (int i = 0; i < 100; ++i) {
    if (fReturnLog[i] > 0) cout << Form("fReturnLog[%3d] = %6d   %s", i, fReturnLog[i], (fReturnString[i]).Data()) << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;

  cout << "Writing " << fHistFile->GetName() << endl;
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;

}

Int_t recoilNtpBase::isMcKs( int iks ) 
{
  Int_t i1 = IndexTrk[d1KsIndex[iks]-1];
  Int_t i2 = IndexTrk[d2KsIndex[iks]-1];
  if (i1>0 && i2>0 && i1<201 && i2<201) {
    Int_t ipar = mothMc[i1-1];
    if ( ipar == mothMc[i2-1] )
      if ( idMc[ipar-1] == 310 )
        return (ipar-1);
  }
  return -1;
}

void recoilNtpBase::findbestB()
{
    double tmpintpur=-999;
    indexbestB=-999;
    bestB0=1;
    
    for (int iB0=0; iB0<nB0; iB0++){
      int mode = modeB0[iB0]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iB0;    
        tmpintpur = brecointpur[mode];
//	tmpintpur = intpurB0[iB0];   // old 
      }
    }
    
    for (int iChB=0; iChB<nChB; iChB++){
      int mode = modeChB[iChB]-10000;
      if(brecointpur[mode]>tmpintpur) {
	indexbestB = iChB;    
        tmpintpur = brecointpur[mode];
//	tmpintpur = intpurChB[iChB];  // old
	bestB0=0;
      }
    }
}

void recoilNtpBase::breco(int chbcand) 
{
  // #trx #pi0 #ks
  nnpi0=0; nntrk=0; nnks=0; nnpar=0;
  //B0
  if(d1B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d1B0Lund[indexbestB]==310) {nnks++;}
  else if(d1B0Lund[indexbestB]!=0) {nntrk++;}
  if(d2B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d2B0Lund[indexbestB]==310) {nnks++;}
  else if(d2B0Lund[indexbestB]!=0) {nntrk++;}
  if(d3B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d3B0Lund[indexbestB]==310) {nnks++;}
  else if(d3B0Lund[indexbestB]!=0) {nntrk++;}
  if(d4B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d4B0Lund[indexbestB]==310) {nnks++;}
  else if(d4B0Lund[indexbestB]!=0) {nntrk++;}
  if(d5B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d5B0Lund[indexbestB]==310) {nnks++;}
  else if(d5B0Lund[indexbestB]!=0) {nntrk++;}
  if(d6B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d6B0Lund[indexbestB]==310) {nnks++;}
  else if(d6B0Lund[indexbestB]!=0) {nntrk++;}
  if(d7B0Lund[indexbestB]==111) {nnpi0++;}
  else if(d7B0Lund[indexbestB]==310) {nnks++;}
  else if(d7B0Lund[indexbestB]!=0) {nntrk++;}
  //ChB
  if(d1ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d1ChBLund[indexbestB]==310) {nnks++;}
  else if(d1ChBLund[indexbestB]!=0) {nntrk++;}
 if(d2ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d2ChBLund[indexbestB]==310) {nnks++;}
  else if(d2ChBLund[indexbestB]!=0) {nntrk++;}
  if(d3ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d3ChBLund[indexbestB]==310) {nnks++;}
  else if(d3ChBLund[indexbestB]!=0) {nntrk++;}
  if(d4ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d4ChBLund[indexbestB]==310) {nnks++;}
  else if(d4ChBLund[indexbestB]!=0) {nntrk++;}
  if(d5ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d5ChBLund[indexbestB]==310) {nnks++;}
  else if(d5ChBLund[indexbestB]!=0) {nntrk++;}
  if(d6ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d6ChBLund[indexbestB]==310) {nnks++;}
  else if(d6ChBLund[indexbestB]!=0) {nntrk++;}
  if(d7ChBLund[indexbestB]==111) {nnpi0++;}
  else if(d7ChBLund[indexbestB]==310) {nnks++;}
  else if(d7ChBLund[indexbestB]!=0) {nntrk++;}
  nnpar=nnpi0+nnks+nntrk;
  
  fBrecoMc = MCB0[indexbestB];
  fMes = mseB0[indexbestB];
  fDeltaE = deltaeB0[indexbestB];

  // Thrust of the B0
  fThru = ThruB0[indexbestB];
  fthThru = thThruB0[indexbestB];
  fphThru = phiThruB0[indexbestB];

  // cos(T_BReco^T_RestOfEvent)
  fcosTBR = cosTBB0[indexbestB];

  int tmpdauB = d1B0Lund[indexbestB];
  int tmpmodeB = modeB0[indexbestB];
  if(chbcand) {
    fBrecoMc = MCChB[indexbestB];
    fMes = mseChB[indexbestB];
    fDeltaE = deltaeChB[indexbestB];
    
    // Thrust of the B0
    fThru = ThruChB[indexbestB];
    fthThru = thThruChB[indexbestB];
    fphThru = phiThruChB[indexbestB];
  
    // cos(T_BReco^T_RestOfEvent)
    fcosTBR = cosTBChB[indexbestB];
    
    tmpdauB = d1ChBLund[indexbestB];
    tmpmodeB = modeChB[indexbestB];
  }

  fBmode = tmpmodeB;
   // BetaCoreTools/BtaExclusiveDecayList.hh
  // the following is an arbitrary definition
  // fSeedMode = 0 dc
  // fSeedMode = 1 dstar
  // fSeedMode = 2 d0
  // fSeedMode = 3 dstar0
  if ((11000 <= tmpmodeB) &&  (tmpmodeB < 12000)) {
    fSeedMode = 2;
  } else if  ((12000 <= tmpmodeB) &&  (tmpmodeB < 13000)) {
    fSeedMode = 0;
  } else if  ((13000 <= tmpmodeB) &&  (tmpmodeB < 14000)) {
    fSeedMode = 1;
  } else if  ((14000 <= tmpmodeB) &&  (tmpmodeB < 16000)) {
    fSeedMode = 3;
  } else {
    fSeedMode = -1;
  }
 fBrecoCharge= 0;
  fBrecoFlavor = 1;  

  if (chbcand) {
    if(d2ChBLund[indexbestB]!=0&&d2ChBLund[indexbestB]!=111&&d2ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d2ChBLund[indexbestB])/d2ChBLund[indexbestB]);}	
     if(d3ChBLund[indexbestB]!=0&&d3ChBLund[indexbestB]!=111&&d3ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d3ChBLund[indexbestB])/d3ChBLund[indexbestB]);}
     if(d4ChBLund[indexbestB]!=0&&d4ChBLund[indexbestB]!=111&&d4ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d4ChBLund[indexbestB])/d4ChBLund[indexbestB]);}
     if(d5ChBLund[indexbestB]!=0&&d5ChBLund[indexbestB]!=111&&d5ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d5ChBLund[indexbestB])/d5ChBLund[indexbestB]);}
     if(d6ChBLund[indexbestB]!=0&&d6ChBLund[indexbestB]!=111&&d6ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d6ChBLund[indexbestB])/d6ChBLund[indexbestB]);}
     if(d7ChBLund[indexbestB]!=0&&d7ChBLund[indexbestB]!=111&&d7ChBLund[indexbestB]!=310){ fBrecoCharge=fBrecoCharge+(TMath::Abs(d7ChBLund[indexbestB])/d7ChBLund[indexbestB]);}	
     fBrecoFlavor=fBrecoCharge;	
  }else{
    fBrecoCharge=0;
    //     cout << indexbestB <<  " " << bestB0 << " " << nB0 << " " << d1B0Lund[indexbestB] << " " << fMes << " " << fIntPurity << endl;
    fBrecoFlavor=-1*(TMath::Abs(d1B0Lund[indexbestB])/d1B0Lund[indexbestB]);
  }  
}

Bool_t recoilNtpBase::kill_dupli(int Up, int Low)
{
  // char Bts[100];
  char iUpLo[100];
  char low[100];
  bool tmpDupli = kFALSE;

  fValMap = " ";
  sprintf(low,"%X",Low);
  int hmany =0;

  if(strlen(low)<8) {
    hmany = 8 - strlen(low);	
  }

  if(hmany == 0) {
    sprintf(iUpLo,"%s%X%s%X","00",Up,"/",Low); 
  } else if(hmany == 1) {
    sprintf(iUpLo,"%s%X%s%X","00",Up,"/0",Low); 
  } else if(hmany == 2) {
    sprintf(iUpLo,"%s%X%s%X","00",Up,"/00",Low); 
  } else if(hmany == 3) {
    sprintf(iUpLo,"%s%X%s%X","00",Up,"/000",Low); 
  } else {
    cout<<"Something Fishy"<<endl;
  }

  if(strlen(iUpLo)<17) {
    cout<<"Small string"<<endl;
  }

  TMapIter *IterMap = new TMapIter(map);  TObject *NextItem;
  IterMap->Reset();
  NextItem = IterMap->Next();

//  cout<<NextItem<<" Poiter to the item"<<endl;

  while (NextItem != 0) {

//  cout<<iUpLo<<" <-- Still debugging "<<endl;

    if( NextItem != 0)  {
//      cout<<iUpLo<<" "<<((TObjString*)(map->GetValue(NextItem)))->GetString()<<endl;
      if(  strcmp(iUpLo,
		  ((TObjString*)(map->GetValue(NextItem)))->GetString())
	   == 0 ) {
	fValMap = ((TObjString*)(map2->GetValue(NextItem)))->GetString();
	tmpDupli = kTRUE;
	break;
      }
      NextItem = IterMap->Next();
    }
  }
  return tmpDupli;
}

void recoilNtpBase::initRest()
{
   
  PURITY     = 0.;
  INTPURITY  = 0.;
  IPURDSTAR  = 0.;
  IPURDC     = 0.;
  IPURDSTAR0 = 0.;
  IPURD0     = 0.;

  KAMOMLO = 0.1;

  // -- CUTS
  DESIGNALLO = -0.1; 
  DESIGNALHI =  0.1; 

  MESSIGNALLO   =  5.27; 
  MESSIGNALHI   =  5.29; 

  MESSIDEBANDLO = 5.20; 
  MESSIDEBANDHI = 5.26; 

  MESSIGNALBANDLO = 5.27; 
  MESSIGNALBANDHI = 5.29; 

  fHistFile = 0;
  fDump = 0;
  fVerbose = 0;

  readintpur();
}

void recoilNtpBase::initVariables()
{
  TLorentzVector p4Null(0., 0., 0., 0.);
  fRunnumber = fLower = fUpper = -99;
  p4Brecoil = p4Breco = p4BrecoNC = p4BrecoGen = p4BrecoGenSmear = p4Null;
  lundB = fSeedMode = fBmode = -999999;
  
  fBrecoCharge = 100; 
  fBrecoFlavor = 100;
  fBrecoMc = -99;
  fPurity = fIntPurity = -99.;
  fMes = fDeltaE = fBmass = -99.;
 
  // Thrust of the BReco
  fThru = fthThru = fphThru = -99.;

  // cos(T_BReco^T_RestOfEvent)
  fcosTBR = -99.;

 // -- Some event quantities
  fRunnumber = runNumber;
  fLower = lowerID; 
  fUpper = upperID;
  fR2 = FoxWol2; 
  fThrust = thrust; 
 
}
// ----------------------------------------------------------------------
void recoilNtpBase::readintpur(){    

   char tableName[1000], buffer[200], fname[100];
   float bmode, dmode, sig, bkg, pur, sb;
   ifstream is("tablepurity.dat");
   int mode;
   while (is.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     sscanf(buffer, "%f %f %f %f %f %f", &bmode, &dmode, &sig, &bkg, &pur, &sb);
     mode = (int) ( (dmode+100) * 100 + bmode-10000 );
     brecosig[mode] = sig;	
     brecobkg[mode] = bkg;
     brecointpur[mode] = pur; 
   }
}

// ----------------------------------------------------------------------
void recoilNtpBase::Init(TTree *tree, int isMC) {
//   Set branch addresses
  
  if (isMC) { 
    fIsMC = kTRUE; 
  } else {
    fIsMC = kFALSE; 
  }
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event",&event);
   fChain->SetBranchAddress("runNumber",&runNumber);
   fChain->SetBranchAddress("platform",&platform);
   fChain->SetBranchAddress("partition",&partition);
   fChain->SetBranchAddress("upperID",&upperID);
   fChain->SetBranchAddress("lowerID",&lowerID);
   //just for Killing
//   fChain->SetBranchAddress("lower",&lower);
//   fChain->SetBranchAddress("upper",&upper);
//   fChain->SetBranchAddress("run",&run); 
//   fChain->SetBranchAddress("pur",&pur);
//   fChain->SetBranchAddress("mode",&mode);

   fChain->SetBranchAddress("beamSX",&beamSX);
   fChain->SetBranchAddress("beamSY",&beamSY);
   fChain->SetBranchAddress("beamSZ",&beamSZ);
   fChain->SetBranchAddress("beamSCovXX",&beamSCovXX);
   fChain->SetBranchAddress("beamSCovYY",&beamSCovYY);
   fChain->SetBranchAddress("beamSCovZZ",&beamSCovZZ);
   fChain->SetBranchAddress("beamSCovXZ",&beamSCovXZ);
   fChain->SetBranchAddress("pxUps",&pxUps);
   fChain->SetBranchAddress("pyUps",&pyUps);
   fChain->SetBranchAddress("pzUps",&pzUps);
   fChain->SetBranchAddress("eUps",&eUps);
   fChain->SetBranchAddress("nTrkTot",&nTrkTot);
   fChain->SetBranchAddress("W2",&W2);
   fChain->SetBranchAddress("FoxWol2",&FoxWol2);
   fChain->SetBranchAddress("FoxWol2Neu",&FoxWol2Neu);
   fChain->SetBranchAddress("thrust",&thrust);
   fChain->SetBranchAddress("thrustNeu",&thrustNeu);
   //MC block
   if(isMC) { 
     fChain->SetBranchAddress("nMc",&nMc);
     fChain->SetBranchAddress("pMc",pMc);
     fChain->SetBranchAddress("massMc",massMc);
     fChain->SetBranchAddress("thetaMc",thetaMc);
     fChain->SetBranchAddress("phiMc",phiMc);
     fChain->SetBranchAddress("idMc",idMc);
     fChain->SetBranchAddress("mothMc",mothMc);
     fChain->SetBranchAddress("nDauMc",nDauMc);
     fChain->SetBranchAddress("xMc",xMc);
     fChain->SetBranchAddress("yMc",yMc);
     fChain->SetBranchAddress("zMc",zMc);
   }
   fChain->SetBranchAddress("nB0",&nB0);
   fChain->SetBranchAddress("massB0",massB0);
   fChain->SetBranchAddress("pB0",pB0);
   fChain->SetBranchAddress("thB0",thB0);
   fChain->SetBranchAddress("phiB0",phiB0);
   fChain->SetBranchAddress("errMassB0",errMassB0);
   fChain->SetBranchAddress("m0B0",m0B0);
   fChain->SetBranchAddress("xB0",xB0);
   fChain->SetBranchAddress("yB0",yB0);
   fChain->SetBranchAddress("zB0",zB0);
   fChain->SetBranchAddress("s2xB0",s2xB0);
   fChain->SetBranchAddress("s2yB0",s2yB0);
   fChain->SetBranchAddress("s2zB0",s2zB0);
   fChain->SetBranchAddress("chi2B0",chi2B0);
   fChain->SetBranchAddress("dofB0",dofB0);
   fChain->SetBranchAddress("stB0",stB0);
   fChain->SetBranchAddress("ndauB0",ndauB0);
   if(isMC)   fChain->SetBranchAddress("MCB0",MCB0);
   fChain->SetBranchAddress("mseB0",mseB0);
   fChain->SetBranchAddress("mHatB0",mHatB0);
   fChain->SetBranchAddress("deltaeB0",deltaeB0);
   fChain->SetBranchAddress("ThruB0",ThruB0);
   fChain->SetBranchAddress("thThruB0",thThruB0);
   fChain->SetBranchAddress("phiThruB0",phiThruB0);
   fChain->SetBranchAddress("cosTBB0",cosTBB0);
   fChain->SetBranchAddress("d1B0Index",d1B0Index);
   fChain->SetBranchAddress("d1B0Lund",d1B0Lund);
   fChain->SetBranchAddress("d2B0Index",d2B0Index);
   fChain->SetBranchAddress("d2B0Lund",d2B0Lund);
   fChain->SetBranchAddress("d3B0Index",d3B0Index);
   fChain->SetBranchAddress("d3B0Lund",d3B0Lund);
   fChain->SetBranchAddress("d4B0Index",d4B0Index);
   fChain->SetBranchAddress("d4B0Lund",d4B0Lund);
   fChain->SetBranchAddress("d5B0Index",d5B0Index);
   fChain->SetBranchAddress("d5B0Lund",d5B0Lund);
   fChain->SetBranchAddress("d6B0Index",d6B0Index);
   fChain->SetBranchAddress("d6B0Lund",d6B0Lund);
   fChain->SetBranchAddress("d7B0Index",d7B0Index);
   fChain->SetBranchAddress("d7B0Lund",d7B0Lund);
   fChain->SetBranchAddress("modeB0",modeB0);
   fChain->SetBranchAddress("purB0",purB0);
   fChain->SetBranchAddress("intpurB0",intpurB0);
   fChain->SetBranchAddress("VtxXLepB0",VtxXLepB0);
   fChain->SetBranchAddress("VtxYLepB0",VtxYLepB0);
   fChain->SetBranchAddress("VtxZLepB0",VtxZLepB0);
   fChain->SetBranchAddress("VtxCovXXLepB0",VtxCovXXLepB0);
   fChain->SetBranchAddress("VtxCovYYLepB0",VtxCovYYLepB0);
   fChain->SetBranchAddress("VtxCovXYLepB0",VtxCovXYLepB0);
   fChain->SetBranchAddress("VtxCovZZLepB0",VtxCovZZLepB0);
   fChain->SetBranchAddress("VtxCovXZLepB0",VtxCovXZLepB0);
   fChain->SetBranchAddress("VtxCovYZLepB0",VtxCovYZLepB0);
   fChain->SetBranchAddress("VtxChiSqLepB0",VtxChiSqLepB0);
   fChain->SetBranchAddress("VtxNDofLepB0",VtxNDofLepB0);
   fChain->SetBranchAddress("VtxStatLepB0",VtxStatLepB0);
   fChain->SetBranchAddress("VtxNUsedLepB0",VtxNUsedLepB0);
   fChain->SetBranchAddress("DocaLepB0",DocaLepB0);
   fChain->SetBranchAddress("DocaErrLepB0",DocaErrLepB0);
   fChain->SetBranchAddress("VtxXXB0",VtxXXB0);
   fChain->SetBranchAddress("VtxYXB0",VtxYXB0);
   fChain->SetBranchAddress("VtxZXB0",VtxZXB0);
   fChain->SetBranchAddress("VtxCovXXXB0",VtxCovXXXB0);
   fChain->SetBranchAddress("VtxCovYYXB0",VtxCovYYXB0);
   fChain->SetBranchAddress("VtxCovXYXB0",VtxCovXYXB0);
   fChain->SetBranchAddress("VtxCovZZXB0",VtxCovZZXB0);
   fChain->SetBranchAddress("VtxCovXZXB0",VtxCovXZXB0);
   fChain->SetBranchAddress("VtxCovYZXB0",VtxCovYZXB0);
   fChain->SetBranchAddress("VtxChiSqXB0",VtxChiSqXB0);
   fChain->SetBranchAddress("VtxNDofXB0",VtxNDofXB0);
   fChain->SetBranchAddress("VtxStatXB0",VtxStatXB0);
   fChain->SetBranchAddress("VtxNUsedXB0",VtxNUsedXB0);
   fChain->SetBranchAddress("VtxPXB0",VtxPXB0);
   fChain->SetBranchAddress("VtxPhiXB0",VtxPhiXB0);
   fChain->SetBranchAddress("VtxThetaXB0",VtxThetaXB0);
   fChain->SetBranchAddress("ThrustXB0",ThrustXB0);
   fChain->SetBranchAddress("ThrustXPhiB0",ThrustXPhiB0);
   fChain->SetBranchAddress("ThrustXThetaB0",ThrustXThetaB0);
   fChain->SetBranchAddress("MassPB0",MassPB0);
   fChain->SetBranchAddress("MassPhiB0",MassPhiB0);
   fChain->SetBranchAddress("MassThetaB0",MassThetaB0);
   fChain->SetBranchAddress("Cov00B0",Cov00B0);
   fChain->SetBranchAddress("Cov10B0",Cov10B0);
   fChain->SetBranchAddress("Cov11B0",Cov11B0);
   fChain->SetBranchAddress("Cov20B0",Cov20B0);
   fChain->SetBranchAddress("Cov21B0",Cov21B0);
   fChain->SetBranchAddress("Cov22B0",Cov22B0);
   fChain->SetBranchAddress("Cov30B0",Cov30B0);
   fChain->SetBranchAddress("Cov31B0",Cov31B0);
   fChain->SetBranchAddress("Cov32B0",Cov32B0);
   fChain->SetBranchAddress("Cov33B0",Cov33B0);
   fChain->SetBranchAddress("nChB",&nChB);
   fChain->SetBranchAddress("massChB",massChB);
   fChain->SetBranchAddress("pChB",pChB);
   fChain->SetBranchAddress("thChB",thChB);
   fChain->SetBranchAddress("phiChB",phiChB);
   fChain->SetBranchAddress("errMassChB",errMassChB);
   fChain->SetBranchAddress("m0ChB",m0ChB);
   fChain->SetBranchAddress("xChB",xChB);
   fChain->SetBranchAddress("yChB",yChB);
   fChain->SetBranchAddress("zChB",zChB);
   fChain->SetBranchAddress("s2xChB",s2xChB);
   fChain->SetBranchAddress("s2yChB",s2yChB);
   fChain->SetBranchAddress("s2zChB",s2zChB);
   fChain->SetBranchAddress("chi2ChB",chi2ChB);
   fChain->SetBranchAddress("dofChB",dofChB);
   fChain->SetBranchAddress("stChB",stChB);
   fChain->SetBranchAddress("ndauChB",ndauChB);
   if(isMC)   fChain->SetBranchAddress("MCChB",MCChB);
   fChain->SetBranchAddress("mseChB",mseChB);
   fChain->SetBranchAddress("mHatChB",mHatChB);
   fChain->SetBranchAddress("deltaeChB",deltaeChB);
   fChain->SetBranchAddress("ThruChB",ThruChB);
   fChain->SetBranchAddress("thThruChB",thThruChB);
   fChain->SetBranchAddress("phiThruChB",phiThruChB);
   fChain->SetBranchAddress("cosTBChB",cosTBChB);
   fChain->SetBranchAddress("d1ChBIndex",d1ChBIndex);
   fChain->SetBranchAddress("d1ChBLund",d1ChBLund);
   fChain->SetBranchAddress("d2ChBIndex",d2ChBIndex);
   fChain->SetBranchAddress("d2ChBLund",d2ChBLund);
   fChain->SetBranchAddress("d3ChBIndex",d3ChBIndex);
   fChain->SetBranchAddress("d3ChBLund",d3ChBLund);
   fChain->SetBranchAddress("d4ChBIndex",d4ChBIndex);
   fChain->SetBranchAddress("d4ChBLund",d4ChBLund);
   fChain->SetBranchAddress("d5ChBIndex",d5ChBIndex);
   fChain->SetBranchAddress("d5ChBLund",d5ChBLund);
   fChain->SetBranchAddress("d6ChBIndex",d6ChBIndex);
   fChain->SetBranchAddress("d6ChBLund",d6ChBLund);
   fChain->SetBranchAddress("d7ChBIndex",d7ChBIndex);
   fChain->SetBranchAddress("d7ChBLund",d7ChBLund);
   fChain->SetBranchAddress("modeChB",modeChB);
   fChain->SetBranchAddress("purChB",purChB);
   fChain->SetBranchAddress("intpurChB",intpurChB);
   fChain->SetBranchAddress("VtxXLepChB",VtxXLepChB);
   fChain->SetBranchAddress("VtxYLepChB",VtxYLepChB);
   fChain->SetBranchAddress("VtxZLepChB",VtxZLepChB);
   fChain->SetBranchAddress("VtxCovXXLepChB",VtxCovXXLepChB);
   fChain->SetBranchAddress("VtxCovYYLepChB",VtxCovYYLepChB);
   fChain->SetBranchAddress("VtxCovXYLepChB",VtxCovXYLepChB);
   fChain->SetBranchAddress("VtxCovZZLepChB",VtxCovZZLepChB);
   fChain->SetBranchAddress("VtxCovXZLepChB",VtxCovXZLepChB);
   fChain->SetBranchAddress("VtxCovYZLepChB",VtxCovYZLepChB);
   fChain->SetBranchAddress("VtxChiSqLepChB",VtxChiSqLepChB);
   fChain->SetBranchAddress("VtxNDofLepChB",VtxNDofLepChB);
   fChain->SetBranchAddress("VtxStatLepChB",VtxStatLepChB);
   fChain->SetBranchAddress("VtxNUsedLepChB",VtxNUsedLepChB);
   fChain->SetBranchAddress("DocaLepChB",DocaLepChB);
   fChain->SetBranchAddress("DocaErrLepChB",DocaErrLepChB);
   fChain->SetBranchAddress("VtxXXChB",VtxXXChB);
   fChain->SetBranchAddress("VtxYXChB",VtxYXChB);
   fChain->SetBranchAddress("VtxZXChB",VtxZXChB);
   fChain->SetBranchAddress("VtxCovXXXChB",VtxCovXXXChB);
   fChain->SetBranchAddress("VtxCovYYXChB",VtxCovYYXChB);
   fChain->SetBranchAddress("VtxCovXYXChB",VtxCovXYXChB);
   fChain->SetBranchAddress("VtxCovZZXChB",VtxCovZZXChB);
   fChain->SetBranchAddress("VtxCovXZXChB",VtxCovXZXChB);
   fChain->SetBranchAddress("VtxCovYZXChB",VtxCovYZXChB);
   fChain->SetBranchAddress("VtxChiSqXChB",VtxChiSqXChB);
   fChain->SetBranchAddress("VtxNDofXChB",VtxNDofXChB);
   fChain->SetBranchAddress("VtxStatXChB",VtxStatXChB);
   fChain->SetBranchAddress("VtxNUsedXChB",VtxNUsedXChB);
   fChain->SetBranchAddress("VtxPXChB",VtxPXChB);
   fChain->SetBranchAddress("VtxPhiXChB",VtxPhiXChB);
   fChain->SetBranchAddress("VtxThetaXChB",VtxThetaXChB);
   fChain->SetBranchAddress("ThrustXChB",ThrustXChB);
   fChain->SetBranchAddress("ThrustXPhiChB",ThrustXPhiChB);
   fChain->SetBranchAddress("ThrustXThetaChB",ThrustXThetaChB);
   fChain->SetBranchAddress("MassPChB",MassPChB);
   fChain->SetBranchAddress("MassPhiChB",MassPhiChB);
   fChain->SetBranchAddress("MassThetaChB",MassThetaChB);
   fChain->SetBranchAddress("Cov00ChB",Cov00ChB);
   fChain->SetBranchAddress("Cov10ChB",Cov10ChB);
   fChain->SetBranchAddress("Cov11ChB",Cov11ChB);
   fChain->SetBranchAddress("Cov20ChB",Cov20ChB);
   fChain->SetBranchAddress("Cov21ChB",Cov21ChB);
   fChain->SetBranchAddress("Cov22ChB",Cov22ChB);
   fChain->SetBranchAddress("Cov30ChB",Cov30ChB);
   fChain->SetBranchAddress("Cov31ChB",Cov31ChB);
   fChain->SetBranchAddress("Cov32ChB",Cov32ChB);
   fChain->SetBranchAddress("Cov33ChB",Cov33ChB);
   fChain->SetBranchAddress("nDstar",&nDstar);
   fChain->SetBranchAddress("massDstar",massDstar);
   fChain->SetBranchAddress("pDstar",pDstar);
   fChain->SetBranchAddress("thDstar",thDstar);
   fChain->SetBranchAddress("phiDstar",phiDstar);
   fChain->SetBranchAddress("errMassDstar",errMassDstar);
   fChain->SetBranchAddress("m0Dstar",m0Dstar);
   fChain->SetBranchAddress("xDstar",xDstar);
   fChain->SetBranchAddress("yDstar",yDstar);
   fChain->SetBranchAddress("zDstar",zDstar);
   fChain->SetBranchAddress("s2xDstar",s2xDstar);
   fChain->SetBranchAddress("s2yDstar",s2yDstar);
   fChain->SetBranchAddress("s2zDstar",s2zDstar);
   fChain->SetBranchAddress("chi2Dstar",chi2Dstar);
   fChain->SetBranchAddress("dofDstar",dofDstar);
   fChain->SetBranchAddress("stDstar",stDstar);
   fChain->SetBranchAddress("ndauDstar",ndauDstar);
   if(isMC)   fChain->SetBranchAddress("MCDstar",MCDstar);
   fChain->SetBranchAddress("d1DstarIndex",d1DstarIndex);
   fChain->SetBranchAddress("d1DstarLund",d1DstarLund);
   fChain->SetBranchAddress("d2DstarIndex",d2DstarIndex);
   fChain->SetBranchAddress("d2DstarLund",d2DstarLund);
   fChain->SetBranchAddress("nDstarBS",&nDstarBS);
   fChain->SetBranchAddress("massDstarBS",massDstarBS);
   fChain->SetBranchAddress("chi2DstarBS",chi2DstarBS);
   fChain->SetBranchAddress("dofDstarBS",dofDstarBS);
   //     fChain->SetBranchAddress("spixDstarBS",spixDstarBS);
   //     fChain->SetBranchAddress("spiyDstarBS",spiyDstarBS);
   //     fChain->SetBranchAddress("spizDstarBS",spizDstarBS);
   fChain->SetBranchAddress("nDstar0",&nDstar0);
   fChain->SetBranchAddress("massDstar0",massDstar0);
   fChain->SetBranchAddress("pDstar0",pDstar0);
   fChain->SetBranchAddress("thDstar0",thDstar0);
   fChain->SetBranchAddress("phiDstar0",phiDstar0);
   fChain->SetBranchAddress("errMassDstar0",errMassDstar0);
   fChain->SetBranchAddress("m0Dstar0",m0Dstar0);
   fChain->SetBranchAddress("xDstar0",xDstar0);
   fChain->SetBranchAddress("yDstar0",yDstar0);
   fChain->SetBranchAddress("zDstar0",zDstar0);
   fChain->SetBranchAddress("s2xDstar0",s2xDstar0);
   fChain->SetBranchAddress("s2yDstar0",s2yDstar0);
   fChain->SetBranchAddress("s2zDstar0",s2zDstar0);
   fChain->SetBranchAddress("chi2Dstar0",chi2Dstar0);
   fChain->SetBranchAddress("dofDstar0",dofDstar0);
   fChain->SetBranchAddress("stDstar0",stDstar0);
   fChain->SetBranchAddress("ndauDstar0",ndauDstar0);
   if(isMC)   fChain->SetBranchAddress("MCDstar0",MCDstar0);
   fChain->SetBranchAddress("d1Dstar0Index",d1Dstar0Index);
   fChain->SetBranchAddress("d1Dstar0Lund",d1Dstar0Lund);
   fChain->SetBranchAddress("d2Dstar0Index",d2Dstar0Index);
   fChain->SetBranchAddress("d2Dstar0Lund",d2Dstar0Lund);
   fChain->SetBranchAddress("nD0",&nD0);
   fChain->SetBranchAddress("massD0",massD0);
   fChain->SetBranchAddress("pD0",pD0);
   fChain->SetBranchAddress("thD0",thD0);
   fChain->SetBranchAddress("phiD0",phiD0);
   fChain->SetBranchAddress("errMassD0",errMassD0);
   fChain->SetBranchAddress("m0D0",m0D0);
   fChain->SetBranchAddress("xD0",xD0);
   fChain->SetBranchAddress("yD0",yD0);
   fChain->SetBranchAddress("zD0",zD0);
   fChain->SetBranchAddress("s2xD0",s2xD0);
   fChain->SetBranchAddress("s2yD0",s2yD0);
   fChain->SetBranchAddress("s2zD0",s2zD0);
   fChain->SetBranchAddress("chi2D0",chi2D0);
   fChain->SetBranchAddress("dofD0",dofD0);
   fChain->SetBranchAddress("stD0",stD0);
   fChain->SetBranchAddress("ndauD0",ndauD0);
   if(isMC)   fChain->SetBranchAddress("MCD0",MCD0);
   fChain->SetBranchAddress("d1D0Index",d1D0Index);
   fChain->SetBranchAddress("d1D0Lund",d1D0Lund);
   fChain->SetBranchAddress("d2D0Index",d2D0Index);
   fChain->SetBranchAddress("d2D0Lund",d2D0Lund);
   fChain->SetBranchAddress("d3D0Index",d3D0Index);
   fChain->SetBranchAddress("d3D0Lund",d3D0Lund);
   fChain->SetBranchAddress("d4D0Index",d4D0Index);
   fChain->SetBranchAddress("d4D0Lund",d4D0Lund);
   fChain->SetBranchAddress("nChD",&nChD);
   fChain->SetBranchAddress("massChD",massChD);
   fChain->SetBranchAddress("pChD",pChD);
   fChain->SetBranchAddress("thChD",thChD);
   fChain->SetBranchAddress("phiChD",phiChD);
   fChain->SetBranchAddress("errMassChD",errMassChD);
   fChain->SetBranchAddress("m0ChD",m0ChD);
   fChain->SetBranchAddress("xChD",xChD);
   fChain->SetBranchAddress("yChD",yChD);
   fChain->SetBranchAddress("zChD",zChD);
   fChain->SetBranchAddress("s2xChD",s2xChD);
   fChain->SetBranchAddress("s2yChD",s2yChD);
   fChain->SetBranchAddress("s2zChD",s2zChD);
   fChain->SetBranchAddress("chi2ChD",chi2ChD);
   fChain->SetBranchAddress("dofChD",dofChD);
   fChain->SetBranchAddress("stChD",stChD);
   fChain->SetBranchAddress("ndauChD",ndauChD);
   if(isMC)   fChain->SetBranchAddress("MCChD",MCChD);
   fChain->SetBranchAddress("d1ChDIndex",d1ChDIndex);
   fChain->SetBranchAddress("d1ChDLund",d1ChDLund);
   fChain->SetBranchAddress("d2ChDIndex",d2ChDIndex);
   fChain->SetBranchAddress("d2ChDLund",d2ChDLund);
   fChain->SetBranchAddress("d3ChDIndex",d3ChDIndex);
   fChain->SetBranchAddress("d3ChDLund",d3ChDLund);
   fChain->SetBranchAddress("d4ChDIndex",d4ChDIndex);
   fChain->SetBranchAddress("d4ChDLund",d4ChDLund);
   fChain->SetBranchAddress("nKs",&nKs);
   fChain->SetBranchAddress("massKs",massKs);
   fChain->SetBranchAddress("pKs",pKs);
   fChain->SetBranchAddress("thKs",thKs);
   fChain->SetBranchAddress("phiKs",phiKs);
   fChain->SetBranchAddress("errMassKs",errMassKs);
   fChain->SetBranchAddress("m0Ks",m0Ks);
   fChain->SetBranchAddress("xKs",xKs);
   fChain->SetBranchAddress("yKs",yKs);
   fChain->SetBranchAddress("zKs",zKs);
   fChain->SetBranchAddress("s2xKs",s2xKs);
   fChain->SetBranchAddress("s2yKs",s2yKs);
   fChain->SetBranchAddress("s2zKs",s2zKs);
   fChain->SetBranchAddress("chi2Ks",chi2Ks);
   fChain->SetBranchAddress("dofKs",dofKs);
   fChain->SetBranchAddress("stKs",stKs);
   fChain->SetBranchAddress("ndauKs",ndauKs);
   if(isMC)   fChain->SetBranchAddress("MCKs",MCKs);
   fChain->SetBranchAddress("d1KsIndex",d1KsIndex);
   fChain->SetBranchAddress("d1KsLund",d1KsLund);
   fChain->SetBranchAddress("d2KsIndex",d2KsIndex);
   fChain->SetBranchAddress("d2KsLund",d2KsLund);
   fChain->SetBranchAddress("nPi0",&nPi0);
   fChain->SetBranchAddress("massPi0",massPi0);
   fChain->SetBranchAddress("pPi0",pPi0);
   fChain->SetBranchAddress("thPi0",thPi0);
   fChain->SetBranchAddress("phiPi0",phiPi0);
   fChain->SetBranchAddress("errMassPi0",errMassPi0);
   fChain->SetBranchAddress("m0Pi0",m0Pi0);
   fChain->SetBranchAddress("xPi0",xPi0);
   fChain->SetBranchAddress("yPi0",yPi0);
   fChain->SetBranchAddress("zPi0",zPi0);
   fChain->SetBranchAddress("s2xPi0",s2xPi0);
   fChain->SetBranchAddress("s2yPi0",s2yPi0);
   fChain->SetBranchAddress("s2zPi0",s2zPi0);
   fChain->SetBranchAddress("chi2Pi0",chi2Pi0);
   fChain->SetBranchAddress("dofPi0",dofPi0);
   fChain->SetBranchAddress("stPi0",stPi0);
   fChain->SetBranchAddress("ndauPi0",ndauPi0);
   if(isMC)   fChain->SetBranchAddress("MCPi0",MCPi0);
   fChain->SetBranchAddress("d1Pi0Index",d1Pi0Index);
   fChain->SetBranchAddress("d1Pi0Lund",d1Pi0Lund);
   fChain->SetBranchAddress("d2Pi0Index",d2Pi0Index);
   fChain->SetBranchAddress("d2Pi0Lund",d2Pi0Lund);
   fChain->SetBranchAddress("nGConv",&nGConv);
   fChain->SetBranchAddress("massGConv",massGConv);
   fChain->SetBranchAddress("pGConv",pGConv);
   fChain->SetBranchAddress("thGConv",thGConv);
   fChain->SetBranchAddress("phiGConv",phiGConv);
   fChain->SetBranchAddress("errMassGConv",errMassGConv);
   fChain->SetBranchAddress("m0GConv",m0GConv);
   fChain->SetBranchAddress("xGConv",xGConv);
   fChain->SetBranchAddress("yGConv",yGConv);
   fChain->SetBranchAddress("zGConv",zGConv);
   fChain->SetBranchAddress("s2xGConv",s2xGConv);
   fChain->SetBranchAddress("s2yGConv",s2yGConv);
   fChain->SetBranchAddress("s2zGConv",s2zGConv);
   fChain->SetBranchAddress("chi2GConv",chi2GConv);
   fChain->SetBranchAddress("dofGConv",dofGConv);
   fChain->SetBranchAddress("stGConv",stGConv);
   fChain->SetBranchAddress("ndauGConv",ndauGConv);
   if(isMC)   fChain->SetBranchAddress("MCGConv",MCGConv);
   fChain->SetBranchAddress("d1GConvIndex",d1GConvIndex);
   fChain->SetBranchAddress("d1GConvLund",d1GConvLund);
   fChain->SetBranchAddress("d2GConvIndex",d2GConvIndex);
   fChain->SetBranchAddress("d2GConvLund",d2GConvLund);
   if (fNewFormat) {
     fChain->SetBranchAddress("nDalitz",&nDalitz);
     fChain->SetBranchAddress("massDalitz",massDalitz);
     fChain->SetBranchAddress("pDalitz",pDalitz);
     fChain->SetBranchAddress("thDalitz",thDalitz);
     fChain->SetBranchAddress("phiDalitz",phiDalitz);
     fChain->SetBranchAddress("errMassDalitz",errMassDalitz);
     fChain->SetBranchAddress("m0Dalitz",m0Dalitz);
     fChain->SetBranchAddress("xDalitz",xDalitz);
     fChain->SetBranchAddress("yDalitz",yDalitz);
     fChain->SetBranchAddress("zDalitz",zDalitz);
     fChain->SetBranchAddress("s2xDalitz",s2xDalitz);
     fChain->SetBranchAddress("s2yDalitz",s2yDalitz);
     fChain->SetBranchAddress("s2zDalitz",s2zDalitz);
     fChain->SetBranchAddress("chi2Dalitz",chi2Dalitz);
     fChain->SetBranchAddress("dofDalitz",dofDalitz);
     fChain->SetBranchAddress("stDalitz",stDalitz);
     fChain->SetBranchAddress("ndauDalitz",ndauDalitz);
     if(isMC)   fChain->SetBranchAddress("MCDalitz",MCDalitz);
     fChain->SetBranchAddress("d1DalitzIndex",d1DalitzIndex);
     fChain->SetBranchAddress("d1DalitzLund",d1DalitzLund);
     fChain->SetBranchAddress("d2DalitzIndex",d2DalitzIndex);
     fChain->SetBranchAddress("d2DalitzLund",d2DalitzLund);
     fChain->SetBranchAddress("nJpsi",&nJpsi);
     fChain->SetBranchAddress("massJpsi",massJpsi);
     fChain->SetBranchAddress("pJpsi",pJpsi);
     fChain->SetBranchAddress("thJpsi",thJpsi);
     fChain->SetBranchAddress("phiJpsi",phiJpsi);
     fChain->SetBranchAddress("errMassJpsi",errMassJpsi);
     fChain->SetBranchAddress("m0Jpsi",m0Jpsi);
     fChain->SetBranchAddress("xJpsi",xJpsi);
     fChain->SetBranchAddress("yJpsi",yJpsi);
     fChain->SetBranchAddress("zJpsi",zJpsi);
     fChain->SetBranchAddress("s2xJpsi",s2xJpsi);
     fChain->SetBranchAddress("s2yJpsi",s2yJpsi);
     fChain->SetBranchAddress("s2zJpsi",s2zJpsi);
     fChain->SetBranchAddress("chi2Jpsi",chi2Jpsi);
     fChain->SetBranchAddress("dofJpsi",dofJpsi);
     fChain->SetBranchAddress("stJpsi",stJpsi);
     fChain->SetBranchAddress("ndauJpsi",ndauJpsi);
     if(isMC)   fChain->SetBranchAddress("MCJpsi",MCJpsi);
     fChain->SetBranchAddress("d1JpsiIndex",d1JpsiIndex);
     fChain->SetBranchAddress("d1JpsiLund",d1JpsiLund);
     fChain->SetBranchAddress("d2JpsiIndex",d2JpsiIndex);
     fChain->SetBranchAddress("d2JpsiLund",d2JpsiLund);
   }   
   fChain->SetBranchAddress("nTrk",&nTrk);
   fChain->SetBranchAddress("IfrLayTrk",IfrLayTrk);
   fChain->SetBranchAddress("IfrNsTrk",IfrNsTrk);
   fChain->SetBranchAddress("IfrInnerTrk",IfrInnerTrk);
   fChain->SetBranchAddress("IfrBarrelTrk",IfrBarrelTrk);
   fChain->SetBranchAddress("IfrFWDTrk",IfrFWDTrk);
   fChain->SetBranchAddress("IfrBWDTrk",IfrBWDTrk);
   fChain->SetBranchAddress("IfrMeasIntLenTrk",IfrMeasIntLenTrk);
   fChain->SetBranchAddress("IfrFirstHitTrk",IfrFirstHitTrk);
   fChain->SetBranchAddress("IfrLastHitTrk",IfrLastHitTrk);
   fChain->SetBranchAddress("lMomTrk",lMomTrk);
   fChain->SetBranchAddress("ZMom42Trk",ZMom42Trk);
   fChain->SetBranchAddress("ecalTrk",ecalTrk);
   fChain->SetBranchAddress("ecalXTrk",ecalXTrk);
   fChain->SetBranchAddress("ecalYTrk",ecalYTrk);
   fChain->SetBranchAddress("ecalZTrk",ecalZTrk);
   fChain->SetBranchAddress("nCryTrk",nCryTrk);
   fChain->SetBranchAddress("nBumpTrk",nBumpTrk);
   fChain->SetBranchAddress("ZMom20Trk",ZMom20Trk);
   fChain->SetBranchAddress("secMomTrk",secMomTrk);
   fChain->SetBranchAddress("s1s9Trk",s1s9Trk);
   fChain->SetBranchAddress("s9s25Trk",s9s25Trk);
   fChain->SetBranchAddress("erawTrk",erawTrk);
   fChain->SetBranchAddress("phiClusterTrk",phiClusterTrk);
   fChain->SetBranchAddress("thetaClusterTrk",thetaClusterTrk);
   fChain->SetBranchAddress("covEETrk",covEETrk);
   fChain->SetBranchAddress("covTTTrk",covTTTrk);
   fChain->SetBranchAddress("covPPTrk",covPPTrk);
   fChain->SetBranchAddress("covRRTrk",covRRTrk);
   fChain->SetBranchAddress("phicMatTrk",phicMatTrk);
   fChain->SetBranchAddress("trkcMatTrk",trkcMatTrk);
   fChain->SetBranchAddress("nPidTrk",nPidTrk);
   fChain->SetBranchAddress("emcStatusTrk",emcStatusTrk);
   fChain->SetBranchAddress("phiAtEMCTrk",phiAtEMCTrk);
   fChain->SetBranchAddress("thetaAtEMCTrk",thetaAtEMCTrk);
   fChain->SetBranchAddress("isvtTrk",isvtTrk);
   fChain->SetBranchAddress("nsvtTrk",nsvtTrk);
   if (fNewFormat) {
     fChain->SetBranchAddress("fhitTrk",fhitTrk);
     fChain->SetBranchAddress("ndchTrk",ndchTrk);
     fChain->SetBranchAddress("lhitTrk",lhitTrk);
     fChain->SetBranchAddress("tLenTrk",tLenTrk);
     fChain->SetBranchAddress("ntdofTrk",ntdofTrk);
     fChain->SetBranchAddress("tproTrk",tproTrk);
     fChain->SetBranchAddress("tChi2Trk",tChi2Trk);
     fChain->SetBranchAddress("cPidTrk",cPidTrk);
     fChain->SetBranchAddress("sfRangeTrk",sfRangeTrk);
     fChain->SetBranchAddress("trkFitStatusTrk",trkFitStatusTrk);
   }
   fChain->SetBranchAddress("chargeTrk",chargeTrk);
   fChain->SetBranchAddress("momentumTrk",momentumTrk);
   fChain->SetBranchAddress("ppcov00",ppcov00);
   fChain->SetBranchAddress("ppcov10",ppcov10);
   fChain->SetBranchAddress("ppcov11",ppcov11);
   fChain->SetBranchAddress("ppcov20",ppcov20);
   fChain->SetBranchAddress("ppcov21",ppcov21);
   fChain->SetBranchAddress("ppcov22",ppcov22);
   fChain->SetBranchAddress("xPocaTrk",xPocaTrk);
   fChain->SetBranchAddress("yPocaTrk",yPocaTrk);
   fChain->SetBranchAddress("zPocaTrk",zPocaTrk);
   fChain->SetBranchAddress("thetaTrk",thetaTrk);
   fChain->SetBranchAddress("phiTrk",phiTrk);
   fChain->SetBranchAddress("muonIdTrk",muonIdTrk);
   fChain->SetBranchAddress("elecIdTrk",elecIdTrk);
   fChain->SetBranchAddress("kaonIdTrk",kaonIdTrk);
   fChain->SetBranchAddress("pionIdTrk",pionIdTrk);
   if (isMC) {
     fChain->SetBranchAddress("idTrk",idTrk);
     fChain->SetBranchAddress("IndexTrk",IndexTrk);
     fChain->SetBranchAddress("IndexNtTrk",IndexNtTrk);
   }
   fChain->SetBranchAddress("B0RecTrk",B0RecTrk);
   fChain->SetBranchAddress("chBRecTrk",chBRecTrk);
   fChain->SetBranchAddress("nGam",&nGam);
   fChain->SetBranchAddress("IfrLayGam",IfrLayGam);
   fChain->SetBranchAddress("IfrNsGam",IfrNsGam);
   fChain->SetBranchAddress("IfrInnerGam",IfrInnerGam);
   fChain->SetBranchAddress("IfrBarrelGam",IfrBarrelGam);
   fChain->SetBranchAddress("IfrFWDGam",IfrFWDGam);
   fChain->SetBranchAddress("IfrBWDGam",IfrBWDGam);
   fChain->SetBranchAddress("IfrMeasIntLenGam",IfrMeasIntLenGam);
   fChain->SetBranchAddress("IfrFirstHitGam",IfrFirstHitGam);
   fChain->SetBranchAddress("IfrLastHitGam",IfrLastHitGam);
   fChain->SetBranchAddress("IfrExpIntLenGam",IfrExpIntLenGam);
   fChain->SetBranchAddress("IfrIntLenBeforeIronGam",IfrIntLenBeforeIronGam);
   fChain->SetBranchAddress("IfrTrkMatchGam",IfrTrkMatchGam);
   fChain->SetBranchAddress("IfrEmcMatchGam",IfrEmcMatchGam);
   fChain->SetBranchAddress("IfrLastBarrelGam",IfrLastBarrelGam);
   fChain->SetBranchAddress("IfrCLFitChi2Gam",IfrCLFitChi2Gam);
   fChain->SetBranchAddress("IfrStrips0",IfrStrips0);
   fChain->SetBranchAddress("IfrStrips1",IfrStrips1);
   fChain->SetBranchAddress("IfrStrips2",IfrStrips2);
   fChain->SetBranchAddress("IfrStrips3",IfrStrips3);
   fChain->SetBranchAddress("IfrStrips4",IfrStrips4);
   fChain->SetBranchAddress("IfrStrips5",IfrStrips5);
   fChain->SetBranchAddress("IfrStrips6",IfrStrips6);
   fChain->SetBranchAddress("IfrStrips7",IfrStrips7);
   fChain->SetBranchAddress("IfrStrips8",IfrStrips8);
   fChain->SetBranchAddress("IfrStrips9",IfrStrips9);
   fChain->SetBranchAddress("IfrStrips10",IfrStrips10);
   fChain->SetBranchAddress("IfrStrips11",IfrStrips11);
   fChain->SetBranchAddress("IfrStrips12",IfrStrips12);
   fChain->SetBranchAddress("IfrStrips13",IfrStrips13);
   fChain->SetBranchAddress("IfrStrips14",IfrStrips14);
   fChain->SetBranchAddress("IfrStrips15",IfrStrips15);
   fChain->SetBranchAddress("IfrStrips16",IfrStrips16);
   fChain->SetBranchAddress("IfrStrips17",IfrStrips17);
   fChain->SetBranchAddress("IfrStrips18",IfrStrips18);
   fChain->SetBranchAddress("IfrStrips19",IfrStrips19);
   fChain->SetBranchAddress("lMomGam",lMomGam);
   fChain->SetBranchAddress("ZMom42Gam",ZMom42Gam);
   fChain->SetBranchAddress("ecalGam",ecalGam);
   fChain->SetBranchAddress("ecalXGam",ecalXGam);
   fChain->SetBranchAddress("ecalYGam",ecalYGam);
   fChain->SetBranchAddress("ecalZGam",ecalZGam);
   fChain->SetBranchAddress("nCryGam",nCryGam);
   fChain->SetBranchAddress("nBumpGam",nBumpGam);
   fChain->SetBranchAddress("ZMom20Gam",ZMom20Gam);
   fChain->SetBranchAddress("secMomGam",secMomGam);
   fChain->SetBranchAddress("s1s9Gam",s1s9Gam);
   fChain->SetBranchAddress("s9s25Gam",s9s25Gam);
   fChain->SetBranchAddress("erawGam",erawGam);
   fChain->SetBranchAddress("phiClusterGam",phiClusterGam);
   fChain->SetBranchAddress("thetaClusterGam",thetaClusterGam);
   fChain->SetBranchAddress("covEEGam",covEEGam);
   fChain->SetBranchAddress("covTTGam",covTTGam);
   fChain->SetBranchAddress("covPPGam",covPPGam);
   fChain->SetBranchAddress("covRRGam",covRRGam);
   fChain->SetBranchAddress("emcStatusGam",emcStatusGam);
   fChain->SetBranchAddress("thetaGam",thetaGam);
   fChain->SetBranchAddress("phiGam",phiGam);
   fChain->SetBranchAddress("energyGam",energyGam);
   if(isMC){
     fChain->SetBranchAddress("idGam",idGam);
     fChain->SetBranchAddress("IndexGam",IndexGam);
     fChain->SetBranchAddress("IndexNtGam",IndexNtGam);
   }
   fChain->SetBranchAddress("B0RecGam",B0RecGam);
   fChain->SetBranchAddress("chBRecGam",chBRecGam);

   Notify();
}
// ----------------------------------------------------------------------
Bool_t recoilNtpBase::Notify() {
//   called when loading a new file
//   get branch pointers
   b_event = fChain->GetBranch("event");
   b_runNumber = fChain->GetBranch("runNumber");
   b_platform = fChain->GetBranch("platform");
   b_partition = fChain->GetBranch("partition");
   b_upperID = fChain->GetBranch("upperID");
   b_lowerID = fChain->GetBranch("lowerID");
//Just for Killing
   // b_upper = fChain->GetBranch("upper");
   // b_lower = fChain->GetBranch("lower");
   //b_run   = fChain->GetBranch("run");
   //b_pur   = fChain->GetBranch("pur");
   //b_mode  = fChain->GetBranch("mode");
//Just for Killing
   b_beamSX = fChain->GetBranch("beamSX");
   b_beamSY = fChain->GetBranch("beamSY");
   b_beamSZ = fChain->GetBranch("beamSZ");
   b_beamSCovXX = fChain->GetBranch("beamSCovXX");
   b_beamSCovYY = fChain->GetBranch("beamSCovYY");
   b_beamSCovZZ = fChain->GetBranch("beamSCovZZ");
   b_beamSCovXZ = fChain->GetBranch("beamSCovXZ");
   b_pxUps = fChain->GetBranch("pxUps");
   b_pyUps = fChain->GetBranch("pyUps");
   b_pzUps = fChain->GetBranch("pzUps");
   b_eUps = fChain->GetBranch("eUps");
   b_nTrkTot = fChain->GetBranch("nTrkTot");
   b_W2 = fChain->GetBranch("W2");
   b_FoxWol2 = fChain->GetBranch("FoxWol2");
   b_FoxWol2Neu = fChain->GetBranch("FoxWol2Neu");
   b_thrust = fChain->GetBranch("thrust");
   b_thrustNeu = fChain->GetBranch("thrustNeu");
   b_nMc = fChain->GetBranch("nMc");
   b_pMc = fChain->GetBranch("pMc");
   b_massMc = fChain->GetBranch("massMc");
   b_thetaMc = fChain->GetBranch("thetaMc");
   b_phiMc = fChain->GetBranch("phiMc");
   b_idMc = fChain->GetBranch("idMc");
   b_mothMc = fChain->GetBranch("mothMc");
   b_nDauMc = fChain->GetBranch("nDauMc");
   b_xMc = fChain->GetBranch("xMc");
   b_yMc = fChain->GetBranch("yMc");
   b_zMc = fChain->GetBranch("zMc");
   b_nB0 = fChain->GetBranch("nB0");
   b_massB0 = fChain->GetBranch("massB0");
   b_pB0 = fChain->GetBranch("pB0");
   b_thB0 = fChain->GetBranch("thB0");
   b_phiB0 = fChain->GetBranch("phiB0");
   b_errMassB0 = fChain->GetBranch("errMassB0");
   b_m0B0 = fChain->GetBranch("m0B0");
   b_xB0 = fChain->GetBranch("xB0");
   b_yB0 = fChain->GetBranch("yB0");
   b_zB0 = fChain->GetBranch("zB0");
   b_s2xB0 = fChain->GetBranch("s2xB0");
   b_s2yB0 = fChain->GetBranch("s2yB0");
   b_s2zB0 = fChain->GetBranch("s2zB0");
   b_chi2B0 = fChain->GetBranch("chi2B0");
   b_dofB0 = fChain->GetBranch("dofB0");
   b_stB0 = fChain->GetBranch("stB0");
   b_ndauB0 = fChain->GetBranch("ndauB0");
   b_MCB0 = fChain->GetBranch("MCB0");
   b_mseB0 = fChain->GetBranch("mseB0");
   b_mHatB0 = fChain->GetBranch("mHatB0");
   b_deltaeB0 = fChain->GetBranch("deltaeB0");
   b_ThruB0 = fChain->GetBranch("ThruB0");
   b_thThruB0 = fChain->GetBranch("thThruB0");
   b_phiThruB0 = fChain->GetBranch("phiThruB0");
   b_cosTBB0 = fChain->GetBranch("cosTBB0");
   b_d1B0Index = fChain->GetBranch("d1B0Index");
   b_d1B0Lund = fChain->GetBranch("d1B0Lund");
   b_d2B0Index = fChain->GetBranch("d2B0Index");
   b_d2B0Lund = fChain->GetBranch("d2B0Lund");
   b_d3B0Index = fChain->GetBranch("d3B0Index");
   b_d3B0Lund = fChain->GetBranch("d3B0Lund");
   b_d4B0Index = fChain->GetBranch("d4B0Index");
   b_d4B0Lund = fChain->GetBranch("d4B0Lund");
   b_d5B0Index = fChain->GetBranch("d5B0Index");
   b_d5B0Lund = fChain->GetBranch("d5B0Lund");
   b_d6B0Index = fChain->GetBranch("d6B0Index");
   b_d6B0Lund = fChain->GetBranch("d6B0Lund");
   b_d7B0Index = fChain->GetBranch("d7B0Index");
   b_d7B0Lund = fChain->GetBranch("d7B0Lund");
   b_modeB0 = fChain->GetBranch("modeB0");
   b_purB0 = fChain->GetBranch("purB0");
   b_intpurB0 = fChain->GetBranch("intpurB0");
   b_VtxXLepB0 = fChain->GetBranch("VtxXLepB0");
   b_VtxYLepB0 = fChain->GetBranch("VtxYLepB0");
   b_VtxZLepB0 = fChain->GetBranch("VtxZLepB0");
   b_VtxCovXXLepB0 = fChain->GetBranch("VtxCovXXLepB0");
   b_VtxCovYYLepB0 = fChain->GetBranch("VtxCovYYLepB0");
   b_VtxCovXYLepB0 = fChain->GetBranch("VtxCovXYLepB0");
   b_VtxCovZZLepB0 = fChain->GetBranch("VtxCovZZLepB0");
   b_VtxCovXZLepB0 = fChain->GetBranch("VtxCovXZLepB0");
   b_VtxCovYZLepB0 = fChain->GetBranch("VtxCovYZLepB0");
   b_VtxChiSqLepB0 = fChain->GetBranch("VtxChiSqLepB0");
   b_VtxNDofLepB0 = fChain->GetBranch("VtxNDofLepB0");
   b_VtxStatLepB0 = fChain->GetBranch("VtxStatLepB0");
   b_VtxNUsedLepB0 = fChain->GetBranch("VtxNUsedLepB0");
   b_DocaLepB0 = fChain->GetBranch("DocaLepB0");
   b_DocaErrLepB0 = fChain->GetBranch("DocaErrLepB0");
   b_VtxXXB0 = fChain->GetBranch("VtxXXB0");
   b_VtxYXB0 = fChain->GetBranch("VtxYXB0");
   b_VtxZXB0 = fChain->GetBranch("VtxZXB0");
   b_VtxCovXXXB0 = fChain->GetBranch("VtxCovXXXB0");
   b_VtxCovYYXB0 = fChain->GetBranch("VtxCovYYXB0");
   b_VtxCovXYXB0 = fChain->GetBranch("VtxCovXYXB0");
   b_VtxCovZZXB0 = fChain->GetBranch("VtxCovZZXB0");
   b_VtxCovXZXB0 = fChain->GetBranch("VtxCovXZXB0");
   b_VtxCovYZXB0 = fChain->GetBranch("VtxCovYZXB0");
   b_VtxChiSqXB0 = fChain->GetBranch("VtxChiSqXB0");
   b_VtxNDofXB0 = fChain->GetBranch("VtxNDofXB0");
   b_VtxStatXB0 = fChain->GetBranch("VtxStatXB0");
   b_VtxNUsedXB0 = fChain->GetBranch("VtxNUsedXB0");
   b_VtxPXB0 = fChain->GetBranch("VtxPXB0");
   b_VtxPhiXB0 = fChain->GetBranch("VtxPhiXB0");
   b_VtxThetaXB0 = fChain->GetBranch("VtxThetaXB0");
   b_ThrustXB0 = fChain->GetBranch("ThrustXB0");
   b_ThrustXPhiB0 = fChain->GetBranch("ThrustXPhiB0");
   b_ThrustXThetaB0 = fChain->GetBranch("ThrustXThetaB0");
   b_MassPB0 = fChain->GetBranch("MassPB0");
   b_MassPhiB0 = fChain->GetBranch("MassPhiB0");
   b_MassThetaB0 = fChain->GetBranch("MassThetaB0");
   b_Cov00B0 = fChain->GetBranch("Cov00B0");
   b_Cov10B0 = fChain->GetBranch("Cov10B0");
   b_Cov11B0 = fChain->GetBranch("Cov11B0");
   b_Cov20B0 = fChain->GetBranch("Cov20B0");
   b_Cov21B0 = fChain->GetBranch("Cov21B0");
   b_Cov22B0 = fChain->GetBranch("Cov22B0");
   b_Cov30B0 = fChain->GetBranch("Cov30B0");
   b_Cov31B0 = fChain->GetBranch("Cov31B0");
   b_Cov32B0 = fChain->GetBranch("Cov32B0");
   b_Cov33B0 = fChain->GetBranch("Cov33B0");
   b_nChB = fChain->GetBranch("nChB");
   b_massChB = fChain->GetBranch("massChB");
   b_pChB = fChain->GetBranch("pChB");
   b_thChB = fChain->GetBranch("thChB");
   b_phiChB = fChain->GetBranch("phiChB");
   b_errMassChB = fChain->GetBranch("errMassChB");
   b_m0ChB = fChain->GetBranch("m0ChB");
   b_xChB = fChain->GetBranch("xChB");
   b_yChB = fChain->GetBranch("yChB");
   b_zChB = fChain->GetBranch("zChB");
   b_s2xChB = fChain->GetBranch("s2xChB");
   b_s2yChB = fChain->GetBranch("s2yChB");
   b_s2zChB = fChain->GetBranch("s2zChB");
   b_chi2ChB = fChain->GetBranch("chi2ChB");
   b_dofChB = fChain->GetBranch("dofChB");
   b_stChB = fChain->GetBranch("stChB");
   b_ndauChB = fChain->GetBranch("ndauChB");
   b_MCChB = fChain->GetBranch("MCChB");
   b_mseChB = fChain->GetBranch("mseChB");
   b_mHatChB = fChain->GetBranch("mHatChB");
   b_deltaeChB = fChain->GetBranch("deltaeChB");
   b_ThruChB = fChain->GetBranch("ThruChB");
   b_thThruChB = fChain->GetBranch("thThruChB");
   b_phiThruChB = fChain->GetBranch("phiThruChB");
   b_cosTBChB = fChain->GetBranch("cosTBChB");
   b_d1ChBIndex = fChain->GetBranch("d1ChBIndex");
   b_d1ChBLund = fChain->GetBranch("d1ChBLund");
   b_d2ChBIndex = fChain->GetBranch("d2ChBIndex");
   b_d2ChBLund = fChain->GetBranch("d2ChBLund");
   b_d3ChBIndex = fChain->GetBranch("d3ChBIndex");
   b_d3ChBLund = fChain->GetBranch("d3ChBLund");
   b_d4ChBIndex = fChain->GetBranch("d4ChBIndex");
   b_d4ChBLund = fChain->GetBranch("d4ChBLund");
   b_d5ChBIndex = fChain->GetBranch("d5ChBIndex");
   b_d5ChBLund = fChain->GetBranch("d5ChBLund");
   b_d6ChBIndex = fChain->GetBranch("d6ChBIndex");
   b_d6ChBLund = fChain->GetBranch("d6ChBLund");
   b_d7ChBIndex = fChain->GetBranch("d7ChBIndex");
   b_d7ChBLund = fChain->GetBranch("d7ChBLund");
   b_modeChB = fChain->GetBranch("modeChB");
   b_purChB = fChain->GetBranch("purChB");
   b_intpurChB = fChain->GetBranch("intpurChB");
   b_VtxXLepChB = fChain->GetBranch("VtxXLepChB");
   b_VtxYLepChB = fChain->GetBranch("VtxYLepChB");
   b_VtxZLepChB = fChain->GetBranch("VtxZLepChB");
   b_VtxCovXXLepChB = fChain->GetBranch("VtxCovXXLepChB");
   b_VtxCovYYLepChB = fChain->GetBranch("VtxCovYYLepChB");
   b_VtxCovXYLepChB = fChain->GetBranch("VtxCovXYLepChB");
   b_VtxCovZZLepChB = fChain->GetBranch("VtxCovZZLepChB");
   b_VtxCovXZLepChB = fChain->GetBranch("VtxCovXZLepChB");
   b_VtxCovYZLepChB = fChain->GetBranch("VtxCovYZLepChB");
   b_VtxChiSqLepChB = fChain->GetBranch("VtxChiSqLepChB");
   b_VtxNDofLepChB = fChain->GetBranch("VtxNDofLepChB");
   b_VtxStatLepChB = fChain->GetBranch("VtxStatLepChB");
   b_VtxNUsedLepChB = fChain->GetBranch("VtxNUsedLepChB");
   b_DocaLepChB = fChain->GetBranch("DocaLepChB");
   b_DocaErrLepChB = fChain->GetBranch("DocaErrLepChB");
   b_VtxXXChB = fChain->GetBranch("VtxXXChB");
   b_VtxYXChB = fChain->GetBranch("VtxYXChB");
   b_VtxZXChB = fChain->GetBranch("VtxZXChB");
   b_VtxCovXXXChB = fChain->GetBranch("VtxCovXXXChB");
   b_VtxCovYYXChB = fChain->GetBranch("VtxCovYYXChB");
   b_VtxCovXYXChB = fChain->GetBranch("VtxCovXYXChB");
   b_VtxCovZZXChB = fChain->GetBranch("VtxCovZZXChB");
   b_VtxCovXZXChB = fChain->GetBranch("VtxCovXZXChB");
   b_VtxCovYZXChB = fChain->GetBranch("VtxCovYZXChB");
   b_VtxChiSqXChB = fChain->GetBranch("VtxChiSqXChB");
   b_VtxNDofXChB = fChain->GetBranch("VtxNDofXChB");
   b_VtxStatXChB = fChain->GetBranch("VtxStatXChB");
   b_VtxNUsedXChB = fChain->GetBranch("VtxNUsedXChB");
   b_VtxPXChB = fChain->GetBranch("VtxPXChB");
   b_VtxPhiXChB = fChain->GetBranch("VtxPhiXChB");
   b_VtxThetaXChB = fChain->GetBranch("VtxThetaXChB");
   b_ThrustXChB = fChain->GetBranch("ThrustXChB");
   b_ThrustXPhiChB = fChain->GetBranch("ThrustXPhiChB");
   b_ThrustXThetaChB = fChain->GetBranch("ThrustXThetaChB");
   b_MassPChB = fChain->GetBranch("MassPChB");
   b_MassPhiChB = fChain->GetBranch("MassPhiChB");
   b_MassThetaChB = fChain->GetBranch("MassThetaChB");
   b_Cov00ChB = fChain->GetBranch("Cov00ChB");
   b_Cov10ChB = fChain->GetBranch("Cov10ChB");
   b_Cov11ChB = fChain->GetBranch("Cov11ChB");
   b_Cov20ChB = fChain->GetBranch("Cov20ChB");
   b_Cov21ChB = fChain->GetBranch("Cov21ChB");
   b_Cov22ChB = fChain->GetBranch("Cov22ChB");
   b_Cov30ChB = fChain->GetBranch("Cov30ChB");
   b_Cov31ChB = fChain->GetBranch("Cov31ChB");
   b_Cov32ChB = fChain->GetBranch("Cov32ChB");
   b_Cov33ChB = fChain->GetBranch("Cov33ChB");
   b_nDstar = fChain->GetBranch("nDstar");
   b_massDstar = fChain->GetBranch("massDstar");
   b_pDstar = fChain->GetBranch("pDstar");
   b_thDstar = fChain->GetBranch("thDstar");
   b_phiDstar = fChain->GetBranch("phiDstar");
   b_errMassDstar = fChain->GetBranch("errMassDstar");
   b_m0Dstar = fChain->GetBranch("m0Dstar");
   b_xDstar = fChain->GetBranch("xDstar");
   b_yDstar = fChain->GetBranch("yDstar");
   b_zDstar = fChain->GetBranch("zDstar");
   b_s2xDstar = fChain->GetBranch("s2xDstar");
   b_s2yDstar = fChain->GetBranch("s2yDstar");
   b_s2zDstar = fChain->GetBranch("s2zDstar");
   b_chi2Dstar = fChain->GetBranch("chi2Dstar");
   b_dofDstar = fChain->GetBranch("dofDstar");
   b_stDstar = fChain->GetBranch("stDstar");
   b_ndauDstar = fChain->GetBranch("ndauDstar");
   b_MCDstar = fChain->GetBranch("MCDstar");
   b_d1DstarIndex = fChain->GetBranch("d1DstarIndex");
   b_d1DstarLund = fChain->GetBranch("d1DstarLund");
   b_d2DstarIndex = fChain->GetBranch("d2DstarIndex");
   b_d2DstarLund = fChain->GetBranch("d2DstarLund");
   b_nDstarBS = fChain->GetBranch("nDstarBS");
   b_massDstarBS = fChain->GetBranch("massDstarBS");
   b_chi2DstarBS = fChain->GetBranch("chi2DstarBS");
   b_dofDstarBS = fChain->GetBranch("dofDstarBS");
   //     b_spixDstarBS = fChain->GetBranch("spixDstarBS");
   //     b_spiyDstarBS = fChain->GetBranch("spiyDstarBS");
   //     b_spizDstarBS = fChain->GetBranch("spizDstarBS");
   b_nDstar0 = fChain->GetBranch("nDstar0");
   b_massDstar0 = fChain->GetBranch("massDstar0");
   b_pDstar0 = fChain->GetBranch("pDstar0");
   b_thDstar0 = fChain->GetBranch("thDstar0");
   b_phiDstar0 = fChain->GetBranch("phiDstar0");
   b_errMassDstar0 = fChain->GetBranch("errMassDstar0");
   b_m0Dstar0 = fChain->GetBranch("m0Dstar0");
   b_xDstar0 = fChain->GetBranch("xDstar0");
   b_yDstar0 = fChain->GetBranch("yDstar0");
   b_zDstar0 = fChain->GetBranch("zDstar0");
   b_s2xDstar0 = fChain->GetBranch("s2xDstar0");
   b_s2yDstar0 = fChain->GetBranch("s2yDstar0");
   b_s2zDstar0 = fChain->GetBranch("s2zDstar0");
   b_chi2Dstar0 = fChain->GetBranch("chi2Dstar0");
   b_dofDstar0 = fChain->GetBranch("dofDstar0");
   b_stDstar0 = fChain->GetBranch("stDstar0");
   b_ndauDstar0 = fChain->GetBranch("ndauDstar0");
   b_MCDstar0 = fChain->GetBranch("MCDstar0");
   b_d1Dstar0Index = fChain->GetBranch("d1Dstar0Index");
   b_d1Dstar0Lund = fChain->GetBranch("d1Dstar0Lund");
   b_d2Dstar0Index = fChain->GetBranch("d2Dstar0Index");
   b_d2Dstar0Lund = fChain->GetBranch("d2Dstar0Lund");
   b_nD0 = fChain->GetBranch("nD0");
   b_massD0 = fChain->GetBranch("massD0");
   b_pD0 = fChain->GetBranch("pD0");
   b_thD0 = fChain->GetBranch("thD0");
   b_phiD0 = fChain->GetBranch("phiD0");
   b_errMassD0 = fChain->GetBranch("errMassD0");
   b_m0D0 = fChain->GetBranch("m0D0");
   b_xD0 = fChain->GetBranch("xD0");
   b_yD0 = fChain->GetBranch("yD0");
   b_zD0 = fChain->GetBranch("zD0");
   b_s2xD0 = fChain->GetBranch("s2xD0");
   b_s2yD0 = fChain->GetBranch("s2yD0");
   b_s2zD0 = fChain->GetBranch("s2zD0");
   b_chi2D0 = fChain->GetBranch("chi2D0");
   b_dofD0 = fChain->GetBranch("dofD0");
   b_stD0 = fChain->GetBranch("stD0");
   b_ndauD0 = fChain->GetBranch("ndauD0");
   b_MCD0 = fChain->GetBranch("MCD0");
   b_d1D0Index = fChain->GetBranch("d1D0Index");
   b_d1D0Lund = fChain->GetBranch("d1D0Lund");
   b_d2D0Index = fChain->GetBranch("d2D0Index");
   b_d2D0Lund = fChain->GetBranch("d2D0Lund");
   b_d3D0Index = fChain->GetBranch("d3D0Index");
   b_d3D0Lund = fChain->GetBranch("d3D0Lund");
   b_d4D0Index = fChain->GetBranch("d4D0Index");
   b_d4D0Lund = fChain->GetBranch("d4D0Lund");
   b_nChD = fChain->GetBranch("nChD");
   b_massChD = fChain->GetBranch("massChD");
   b_pChD = fChain->GetBranch("pChD");
   b_thChD = fChain->GetBranch("thChD");
   b_phiChD = fChain->GetBranch("phiChD");
   b_errMassChD = fChain->GetBranch("errMassChD");
   b_m0ChD = fChain->GetBranch("m0ChD");
   b_xChD = fChain->GetBranch("xChD");
   b_yChD = fChain->GetBranch("yChD");
   b_zChD = fChain->GetBranch("zChD");
   b_s2xChD = fChain->GetBranch("s2xChD");
   b_s2yChD = fChain->GetBranch("s2yChD");
   b_s2zChD = fChain->GetBranch("s2zChD");
   b_chi2ChD = fChain->GetBranch("chi2ChD");
   b_dofChD = fChain->GetBranch("dofChD");
   b_stChD = fChain->GetBranch("stChD");
   b_ndauChD = fChain->GetBranch("ndauChD");
   b_MCChD = fChain->GetBranch("MCChD");
   b_d1ChDIndex = fChain->GetBranch("d1ChDIndex");
   b_d1ChDLund = fChain->GetBranch("d1ChDLund");
   b_d2ChDIndex = fChain->GetBranch("d2ChDIndex");
   b_d2ChDLund = fChain->GetBranch("d2ChDLund");
   b_d3ChDIndex = fChain->GetBranch("d3ChDIndex");
   b_d3ChDLund = fChain->GetBranch("d3ChDLund");
   b_d4ChDIndex = fChain->GetBranch("d4ChDIndex");
   b_d4ChDLund = fChain->GetBranch("d4ChDLund");
   b_nKs = fChain->GetBranch("nKs");
   b_massKs = fChain->GetBranch("massKs");
   b_pKs = fChain->GetBranch("pKs");
   b_thKs = fChain->GetBranch("thKs");
   b_phiKs = fChain->GetBranch("phiKs");
   b_errMassKs = fChain->GetBranch("errMassKs");
   b_m0Ks = fChain->GetBranch("m0Ks");
   b_xKs = fChain->GetBranch("xKs");
   b_yKs = fChain->GetBranch("yKs");
   b_zKs = fChain->GetBranch("zKs");
   b_s2xKs = fChain->GetBranch("s2xKs");
   b_s2yKs = fChain->GetBranch("s2yKs");
   b_s2zKs = fChain->GetBranch("s2zKs");
   b_chi2Ks = fChain->GetBranch("chi2Ks");
   b_dofKs = fChain->GetBranch("dofKs");
   b_stKs = fChain->GetBranch("stKs");
   b_ndauKs = fChain->GetBranch("ndauKs");
   b_MCKs = fChain->GetBranch("MCKs");
   b_d1KsIndex = fChain->GetBranch("d1KsIndex");
   b_d1KsLund = fChain->GetBranch("d1KsLund");
   b_d2KsIndex = fChain->GetBranch("d2KsIndex");
   b_d2KsLund = fChain->GetBranch("d2KsLund");
   b_nPi0 = fChain->GetBranch("nPi0");
   b_massPi0 = fChain->GetBranch("massPi0");
   b_pPi0 = fChain->GetBranch("pPi0");
   b_thPi0 = fChain->GetBranch("thPi0");
   b_phiPi0 = fChain->GetBranch("phiPi0");
   b_errMassPi0 = fChain->GetBranch("errMassPi0");
   b_m0Pi0 = fChain->GetBranch("m0Pi0");
   b_xPi0 = fChain->GetBranch("xPi0");
   b_yPi0 = fChain->GetBranch("yPi0");
   b_zPi0 = fChain->GetBranch("zPi0");
   b_s2xPi0 = fChain->GetBranch("s2xPi0");
   b_s2yPi0 = fChain->GetBranch("s2yPi0");
   b_s2zPi0 = fChain->GetBranch("s2zPi0");
   b_chi2Pi0 = fChain->GetBranch("chi2Pi0");
   b_dofPi0 = fChain->GetBranch("dofPi0");
   b_stPi0 = fChain->GetBranch("stPi0");
   b_ndauPi0 = fChain->GetBranch("ndauPi0");
   b_MCPi0 = fChain->GetBranch("MCPi0");
   b_d1Pi0Index = fChain->GetBranch("d1Pi0Index");
   b_d1Pi0Lund = fChain->GetBranch("d1Pi0Lund");
   b_d2Pi0Index = fChain->GetBranch("d2Pi0Index");
   b_d2Pi0Lund = fChain->GetBranch("d2Pi0Lund");
   b_nGConv = fChain->GetBranch("nGConv");
   b_massGConv = fChain->GetBranch("massGConv");
   b_pGConv = fChain->GetBranch("pGConv");
   b_thGConv = fChain->GetBranch("thGConv");
   b_phiGConv = fChain->GetBranch("phiGConv");
   b_errMassGConv = fChain->GetBranch("errMassGConv");
   b_m0GConv = fChain->GetBranch("m0GConv");
   b_xGConv = fChain->GetBranch("xGConv");
   b_yGConv = fChain->GetBranch("yGConv");
   b_zGConv = fChain->GetBranch("zGConv");
   b_s2xGConv = fChain->GetBranch("s2xGConv");
   b_s2yGConv = fChain->GetBranch("s2yGConv");
   b_s2zGConv = fChain->GetBranch("s2zGConv");
   b_chi2GConv = fChain->GetBranch("chi2GConv");
   b_dofGConv = fChain->GetBranch("dofGConv");
   b_stGConv = fChain->GetBranch("stGConv");
   b_ndauGConv = fChain->GetBranch("ndauGConv");
   b_MCGConv = fChain->GetBranch("MCGConv");
   b_d1GConvIndex = fChain->GetBranch("d1GConvIndex");
   b_d1GConvLund = fChain->GetBranch("d1GConvLund");
   b_d2GConvIndex = fChain->GetBranch("d2GConvIndex");
   b_d2GConvLund = fChain->GetBranch("d2GConvLund");
   if (fNewFormat) {
     b_nDalitz = fChain->GetBranch("nDalitz");
     b_massDalitz = fChain->GetBranch("massDalitz");
     b_pDalitz = fChain->GetBranch("pDalitz");
     b_thDalitz = fChain->GetBranch("thDalitz");
     b_phiDalitz = fChain->GetBranch("phiDalitz");
     b_errMassDalitz = fChain->GetBranch("errMassDalitz");
     b_m0Dalitz = fChain->GetBranch("m0Dalitz");
     b_xDalitz = fChain->GetBranch("xDalitz");
     b_yDalitz = fChain->GetBranch("yDalitz");
     b_zDalitz = fChain->GetBranch("zDalitz");
     b_s2xDalitz = fChain->GetBranch("s2xDalitz");
     b_s2yDalitz = fChain->GetBranch("s2yDalitz");
     b_s2zDalitz = fChain->GetBranch("s2zDalitz");
     b_chi2Dalitz = fChain->GetBranch("chi2Dalitz");
     b_dofDalitz = fChain->GetBranch("dofDalitz");
     b_stDalitz = fChain->GetBranch("stDalitz");
     b_ndauDalitz = fChain->GetBranch("ndauDalitz");
     b_MCDalitz = fChain->GetBranch("MCDalitz");
     b_d1DalitzIndex = fChain->GetBranch("d1DalitzIndex");
     b_d1DalitzLund = fChain->GetBranch("d1DalitzLund");
     b_d2DalitzIndex = fChain->GetBranch("d2DalitzIndex");
     b_d2DalitzLund = fChain->GetBranch("d2DalitzLund");
     b_nJpsi = fChain->GetBranch("nJpsi");
     b_massJpsi = fChain->GetBranch("massJpsi");
     b_pJpsi = fChain->GetBranch("pJpsi");
     b_thJpsi = fChain->GetBranch("thJpsi");
     b_phiJpsi = fChain->GetBranch("phiJpsi");
     b_errMassJpsi = fChain->GetBranch("errMassJpsi");
     b_m0Jpsi = fChain->GetBranch("m0Jpsi");
     b_xJpsi = fChain->GetBranch("xJpsi");
     b_yJpsi = fChain->GetBranch("yJpsi");
     b_zJpsi = fChain->GetBranch("zJpsi");
     b_s2xJpsi = fChain->GetBranch("s2xJpsi");
     b_s2yJpsi = fChain->GetBranch("s2yJpsi");
     b_s2zJpsi = fChain->GetBranch("s2zJpsi");
     b_chi2Jpsi = fChain->GetBranch("chi2Jpsi");
     b_dofJpsi = fChain->GetBranch("dofJpsi");
     b_stJpsi = fChain->GetBranch("stJpsi");
     b_ndauJpsi = fChain->GetBranch("ndauJpsi");
     b_MCJpsi = fChain->GetBranch("MCJpsi");
     b_d1JpsiIndex = fChain->GetBranch("d1JpsiIndex");
     b_d1JpsiLund = fChain->GetBranch("d1JpsiLund");
     b_d2JpsiIndex = fChain->GetBranch("d2JpsiIndex");
     b_d2JpsiLund = fChain->GetBranch("d2JpsiLund");
   }
   b_nTrk = fChain->GetBranch("nTrk");
   b_IfrLayTrk = fChain->GetBranch("IfrLayTrk");
   b_IfrNsTrk = fChain->GetBranch("IfrNsTrk");
   b_IfrInnerTrk = fChain->GetBranch("IfrInnerTrk");
   b_IfrBarrelTrk = fChain->GetBranch("IfrBarrelTrk");
   b_IfrFWDTrk = fChain->GetBranch("IfrFWDTrk");
   b_IfrBWDTrk = fChain->GetBranch("IfrBWDTrk");
   b_IfrMeasIntLenTrk = fChain->GetBranch("IfrMeasIntLenTrk");
   b_IfrFirstHitTrk = fChain->GetBranch("IfrFirstHitTrk");
   b_IfrLastHitTrk = fChain->GetBranch("IfrLastHitTrk");
   b_lMomTrk = fChain->GetBranch("lMomTrk");
   b_ZMom42Trk = fChain->GetBranch("ZMom42Trk");
   b_ecalTrk = fChain->GetBranch("ecalTrk");
   b_ecalXTrk = fChain->GetBranch("ecalXTrk");
   b_ecalYTrk = fChain->GetBranch("ecalYTrk");
   b_ecalZTrk = fChain->GetBranch("ecalZTrk");
   b_nCryTrk = fChain->GetBranch("nCryTrk");
   b_nBumpTrk = fChain->GetBranch("nBumpTrk");
   b_ZMom20Trk = fChain->GetBranch("ZMom20Trk");
   b_secMomTrk = fChain->GetBranch("secMomTrk");
   b_s1s9Trk = fChain->GetBranch("s1s9Trk");
   b_s9s25Trk = fChain->GetBranch("s9s25Trk");
   b_erawTrk = fChain->GetBranch("erawTrk");
   b_phiClusterTrk = fChain->GetBranch("phiClusterTrk");
   b_thetaClusterTrk = fChain->GetBranch("thetaClusterTrk");
   b_covEETrk = fChain->GetBranch("covEETrk");
   b_covTTTrk = fChain->GetBranch("covTTTrk");
   b_covPPTrk = fChain->GetBranch("covPPTrk");
   b_covRRTrk = fChain->GetBranch("covRRTrk");
   b_phicMatTrk = fChain->GetBranch("phicMatTrk");
   b_trkcMatTrk = fChain->GetBranch("trkcMatTrk");
   b_nPidTrk = fChain->GetBranch("nPidTrk");
   b_emcStatusTrk = fChain->GetBranch("emcStatusTrk");
   b_phiAtEMCTrk = fChain->GetBranch("phiAtEMCTrk");
   b_thetaAtEMCTrk = fChain->GetBranch("thetaAtEMCTrk");
   b_isvtTrk = fChain->GetBranch("isvtTrk");
   b_nsvtTrk = fChain->GetBranch("nsvtTrk");
   b_fhitTrk = fChain->GetBranch("fhitTrk");
   b_ndchTrk = fChain->GetBranch("ndchTrk");
   b_lhitTrk = fChain->GetBranch("lhitTrk");
   b_tLenTrk = fChain->GetBranch("tLenTrk");
   //     b_ntdofTrk = fChain->GetBranch("ntdofTrk");
   //     b_tproTrk = fChain->GetBranch("tproTrk");
   //     b_tChi2Trk = fChain->GetBranch("tChi2Trk");
   //     b_cPidTrk = fChain->GetBranch("cPidTrk");
   //     b_sfRangeTrk = fChain->GetBranch("sfRangeTrk");
   //     b_trkFitStatusTrk = fChain->GetBranch("trkFitStatusTrk");
   b_chargeTrk = fChain->GetBranch("chargeTrk");
   b_momentumTrk = fChain->GetBranch("momentumTrk");
   b_ppcov00 = fChain->GetBranch("ppcov00");
   b_ppcov10 = fChain->GetBranch("ppcov10");
   b_ppcov11 = fChain->GetBranch("ppcov11");
   b_ppcov20 = fChain->GetBranch("ppcov20");
   b_ppcov21 = fChain->GetBranch("ppcov21");
   b_ppcov22 = fChain->GetBranch("ppcov22");
   b_xPocaTrk = fChain->GetBranch("xPocaTrk");
   b_yPocaTrk = fChain->GetBranch("yPocaTrk");
   b_zPocaTrk = fChain->GetBranch("zPocaTrk");
   b_thetaTrk = fChain->GetBranch("thetaTrk");
   b_phiTrk = fChain->GetBranch("phiTrk");
   b_muonIdTrk = fChain->GetBranch("muonIdTrk");
   b_elecIdTrk = fChain->GetBranch("elecIdTrk");
   b_kaonIdTrk = fChain->GetBranch("kaonIdTrk");
   b_pionIdTrk = fChain->GetBranch("pionIdTrk");
   b_idTrk = fChain->GetBranch("idTrk");
   b_IndexTrk = fChain->GetBranch("IndexTrk");
   b_IndexNtTrk = fChain->GetBranch("IndexNtTrk");
   b_B0RecTrk = fChain->GetBranch("B0RecTrk");
   b_chBRecTrk = fChain->GetBranch("chBRecTrk");
   b_nGam = fChain->GetBranch("nGam");
   b_IfrLayGam = fChain->GetBranch("IfrLayGam");
   b_IfrNsGam = fChain->GetBranch("IfrNsGam");
   b_IfrInnerGam = fChain->GetBranch("IfrInnerGam");
   b_IfrBarrelGam = fChain->GetBranch("IfrBarrelGam");
   b_IfrFWDGam = fChain->GetBranch("IfrFWDGam");
   b_IfrBWDGam = fChain->GetBranch("IfrBWDGam");
   b_IfrMeasIntLenGam = fChain->GetBranch("IfrMeasIntLenGam");
   b_IfrFirstHitGam = fChain->GetBranch("IfrFirstHitGam");
   b_IfrLastHitGam = fChain->GetBranch("IfrLastHitGam");
   b_IfrExpIntLenGam = fChain->GetBranch("IfrExpIntLenGam");
   b_IfrIntLenBeforeIronGam = fChain->GetBranch("IfrIntLenBeforeIronGam");
   b_IfrTrkMatchGam = fChain->GetBranch("IfrTrkMatchGam");
   b_IfrEmcMatchGam = fChain->GetBranch("IfrEmcMatchGam");
   b_IfrLastBarrelGam = fChain->GetBranch("IfrLastBarrelGam");
   b_IfrCLFitChi2Gam = fChain->GetBranch("IfrCLFitChi2Gam");
   b_IfrStrips0 = fChain->GetBranch("IfrStrips0");
   b_IfrStrips1 = fChain->GetBranch("IfrStrips1");
   b_IfrStrips2 = fChain->GetBranch("IfrStrips2");
   b_IfrStrips3 = fChain->GetBranch("IfrStrips3");
   b_IfrStrips4 = fChain->GetBranch("IfrStrips4");
   b_IfrStrips5 = fChain->GetBranch("IfrStrips5");
   b_IfrStrips6 = fChain->GetBranch("IfrStrips6");
   b_IfrStrips7 = fChain->GetBranch("IfrStrips7");
   b_IfrStrips8 = fChain->GetBranch("IfrStrips8");
   b_IfrStrips9 = fChain->GetBranch("IfrStrips9");
   b_IfrStrips10 = fChain->GetBranch("IfrStrips10");
   b_IfrStrips11 = fChain->GetBranch("IfrStrips11");
   b_IfrStrips12 = fChain->GetBranch("IfrStrips12");
   b_IfrStrips13 = fChain->GetBranch("IfrStrips13");
   b_IfrStrips14 = fChain->GetBranch("IfrStrips14");
   b_IfrStrips15 = fChain->GetBranch("IfrStrips15");
   b_IfrStrips16 = fChain->GetBranch("IfrStrips16");
   b_IfrStrips17 = fChain->GetBranch("IfrStrips17");
   b_IfrStrips18 = fChain->GetBranch("IfrStrips18");
   b_IfrStrips19 = fChain->GetBranch("IfrStrips19");
   b_lMomGam = fChain->GetBranch("lMomGam");
   b_ZMom42Gam = fChain->GetBranch("ZMom42Gam");
   b_ecalGam = fChain->GetBranch("ecalGam");
   b_ecalXGam = fChain->GetBranch("ecalXGam");
   b_ecalYGam = fChain->GetBranch("ecalYGam");
   b_ecalZGam = fChain->GetBranch("ecalZGam");
   b_nCryGam = fChain->GetBranch("nCryGam");
   b_nBumpGam = fChain->GetBranch("nBumpGam");
   b_ZMom20Gam = fChain->GetBranch("ZMom20Gam");
   b_secMomGam = fChain->GetBranch("secMomGam");
   b_s1s9Gam = fChain->GetBranch("s1s9Gam");
   b_s9s25Gam = fChain->GetBranch("s9s25Gam");
   b_erawGam = fChain->GetBranch("erawGam");
   b_phiClusterGam = fChain->GetBranch("phiClusterGam");
   b_thetaClusterGam = fChain->GetBranch("thetaClusterGam");
   b_covEEGam = fChain->GetBranch("covEEGam");
   b_covTTGam = fChain->GetBranch("covTTGam");
   b_covPPGam = fChain->GetBranch("covPPGam");
   b_covRRGam = fChain->GetBranch("covRRGam");
   b_emcStatusGam = fChain->GetBranch("emcStatusGam");
   b_thetaGam = fChain->GetBranch("thetaGam");
   b_phiGam = fChain->GetBranch("phiGam");
   b_energyGam = fChain->GetBranch("energyGam");
   b_idGam = fChain->GetBranch("idGam");
   b_IndexGam = fChain->GetBranch("IndexGam");
   b_IndexNtGam = fChain->GetBranch("IndexNtGam");
   b_B0RecGam = fChain->GetBranch("B0RecGam");
   b_chBRecGam = fChain->GetBranch("chBRecGam");
   return kTRUE;
}

// ----------------------------------------------------------------------
void recoilNtpBase::Show(Int_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ----------------------------------------------------------------------
Int_t recoilNtpBase::Cut(Int_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

// ----------------------------------------------------------------------
Int_t recoilNtpBase::LoadTree(Int_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}
