#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include "recoilTaunu.hh"
#include "TH2.h"
#include "TH1.h"
#include <TMap.h>
#include <TObjString.h>
#include "TLorentzVector.h"
#include <assert.h>
#include "TRandom.h"
#include <vector>
#include <algorithm>


class MomentumComparator
{
public:
  MomentumComparator( Float_t * momenta ) : 
    _momenta ( momenta ) { }
  bool operator()( int i1, int i2 )
  { return _momenta[ i1 ] > _momenta[ i2 ]; }
private:
  Float_t * _momenta;
};



// ----------------------------------------------------------------------
recoilTaunu::recoilTaunu(TTree *tree, int isMC, int newFormat, bool selectB0, bool smearNeutrals, bool genericBBbar, int oneProng) : 
  recoilNtpBase(tree, isMC, newFormat) ,
  _selectB0( selectB0 ),
  _smearNeutrals( smearNeutrals ),
  _genericBB( genericBBbar )
{  
}

// ----------------------------------------------------------------------
recoilTaunu::recoilTaunu(TString filename, TString treename,int isMC, int newFormat, bool selectB0,  bool smearNeutrals, bool genericBBbar, int oneProng) :
  recoilNtpBase(filename, treename, isMC, newFormat),
  _selectB0( selectB0 ),
  _smearNeutrals( smearNeutrals )

{
  cout <<"recoilTaunu constructor filename treename isMC"<< endl;
}
// ---------------------------------------------------------------------- 
  

// ----------------------------------------------------------------------
void
recoilTaunu::recoil(int chbcand )
{
  maskKs();
  fPresel = false;
  threeProngs( chbcand );
  if (fPresel) fTree->Fill();
}


// ----------------------------------------------------------------------
void recoilTaunu::fillRecoilHist(const char *dir, int chbcand) 
{ } 


// ----------------------------------------------------------------------
void recoilTaunu::mcTruth( int chbcand ) {

  if (!fIsMC) return;

  TLorentzVector p4Null(0., 0., 0., 0.);
  Bool_t isBMc(kFALSE), isSemilep(kFALSE);
  Int_t ipB[] = {-1, -1};
  Int_t typB[] = {7, 7};
  Int_t cnt(-1),idpi(-999), idmu(-999), idel(-999);

  int tmpBLUND = B0LUND;
  if(chbcand) {
    tmpBLUND = CHBLUND;
  }

  fisEl = fisPi = fisMu = kFALSE;
  fB1Index = fB2Index = fBVxb = -99;

  // -- Find indices of B's
  for (Int_t imc = 0; imc < nMc; ++imc) {    
    if (TMath::Abs(idMc[imc]) == tmpBLUND) isBMc = kTRUE;
    if ((TMath::Abs(idMc[imc]) == 511) || (TMath::Abs(idMc[imc]) == 521)) {
      cnt++;
      ipB[cnt] = imc;
      if (fB1Index < 0) {
	fB1Index = imc; 
      } else {
	fB2Index = imc; 
      }
    }
    if (cnt == 1) break;
  }

  // -- Find number of K_L's
  fKlMc = 0; fKlMcEMC=0;
  for (Int_t imc = 0; imc < nMc; ++imc) {    
  // K_L counting
    if (idMc[imc] == 130) { 
      fKlMc++ ;
      // K_L's in the EMC acceptance
      if ((thetaMc[imc] > 0.410) &&(thetaMc[imc] < 2.409)) fKlMcEMC++;
    }
  }

  Int_t ib;
  fBVxbTyp = -1;

  // Meaning of fBVxbTyp:
  // -- Determine special modes of sl (or leptonic taunu) B decay: 
  //    1 UNUSED
  //    2 UNUSED
  //    3 UNUSED
  //    4 UNUSED
  //    5 TAU NU, TAU->E
  //    6 TAU NU, TAU->MU
  //    7 TAU NU, TAU->PI
  //    8 TAU NU, TAU->A1 (RHO0)
  //    9 TAU NU, TAU->A1 (RHO+-)
  //   10 TAU NU, TAU->RHO+-
  //   11 TAU NU, TAU->3 PI (+1 more particle)
  //  100 TAU NU, TAU->ALTRO

  // Move our code here because what we found did not work!
  for (ib = 0; ib < 2; ib++) {
    if (nDauMc[ipB[ib]]==2) {
      for (Int_t imc = 0; imc < nMc; ++imc) {
	if (isTruTau(imc)) {
	  if ((mothMc[imc]-1) == ipB[ib]) {
	    // tau -> 2-body decays 
	    if (nDauMc[imc]==2) {
	      for (Int_t p = 0; p < nMc; ++p) {
		// tau -> a1 nu.
		if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 20213)) {
		  for(Int_t r = 0; r < nMc; ++r) {
		    if( ((mothMc[r]-1) == p) && (TMath::Abs(idMc[r]) == 113 ) ) {
		      fBVxbTyp=8;
		    }
		    else if ( ((mothMc[r]-1) == p) &&  (TMath::Abs(idMc[r]) == 213)) {
		      fBVxbTyp = 9; // 1 prong ...
		    }
		  }
		}
		// tau -> pi nu.
		else if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 211)) {
		  fBVxbTyp=7;
		}
		// tau -> rho+- nu.
		else if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 213))
		  fBVxbTyp=10;
	      }
	    }
	    // tau -> 3-body decays 
	    else if (nDauMc[imc]==3) {
	      for (Int_t lep = 0; lep < nMc; ++lep) {
		// tau -> e nu nu
		if ( ((mothMc[lep]-1) == imc) &&  (TMath::Abs(idMc[lep]) == 11)) {
		  fBVxbTyp=5;
		}
		// tau -> mu nu nu
		if ( ((mothMc[lep]-1) == imc) &&  (TMath::Abs(idMc[lep]) == 13)) {
		  fBVxbTyp=6;
		}
	      } 
	    }
	    // tau -> 5-body decays 
	    else if( (nDauMc[imc]==5) ){
	      // tau -> pi pi pi  ... il quinto figlio qualunque (pi0 e' il piu' freq.)
	      Int_t primarypion = 0;
	      for (Int_t p = 0; p < nMc; ++p) {
		if ( ((mothMc[p]-1) == imc) && (TMath::Abs(idMc[p]) == 211) )
		  primarypion++;
	      }
	      if ( primarypion == 3) {
		fBVxbTyp = 11;
	      }
	    }
	  }
	}
      }
    }
  }
}



// ----------------------------------------------------------------------
void recoilTaunu::bookHist(int dump) {

  cout << " entering recoilTaunu bookHist " << endl;

  char name[100], title[100];

  if (!fHistFile) {
    cout << "Call recoilTaunu::openHistFile(...) before booking histograms" << endl;
  } else {
    fHistFile->cd();
  }    

  if (dump > 0) {
    cout << "Booking events tree" << endl;
    fDump = dump;
    fTree = new TTree("events", "events"); 
    fTree->Branch("run",    &fRunnumber, "run/I");
    fTree->Branch("lower",  &fLower, "lower/I");
    fTree->Branch("upper",  &fUpper, "upper/I");
    
    //CB Thrust and Fox-Wolfram

    fTree->Branch("thrust",     &fThrust, "thrust/D");
    fTree->Branch("R2",         &fR2, "R2/D");

    fTree->Branch("mycosth",    &fcosmyTh, "mycosth/D");


    // -- breco
    fTree->Branch("bmass",      &fBmass, "bmass/D");
    //    fTree->Branch("bmassfit",   &fBmassfit, "bmassfit/D");
    fTree->Branch("sbox",       &signalBox, "sbox/b");
    fTree->Branch("mes",        &fMes, "mes/D");
    fTree->Branch("de",         &fDeltaE, "de/D");
    fTree->Branch("pur",        &fPurity, "pur/D");
    fTree->Branch("intpur",     &fIntPurity, "intpur/D");
    fTree->Branch("brecoflav",  &fBrecoFlavor, "brecoflav/I");
    fTree->Branch("brecocharge",&fBrecoCharge , "brecocharge/I");
    fTree->Branch("brecomc",    &fBrecoMc,  "brecomc/I");
    fTree->Branch("MCB1",       &fB1Index,  "MCB1/I");
    fTree->Branch("MCB2",       &fB2Index,  "MCB2/I");
    fTree->Branch("ThruB",      &fThru, "ThruB/D");
    fTree->Branch("thThruB",    &fthThru, "thThruB/D");
    fTree->Branch("phThruB",    &fphThru, "phThruB/D");
    fTree->Branch("cosTBR",     &fcosTBR, "cosTBR/D");

    // MC info
    fTree->Branch("MBVxb",       &fBVxb, "MBVxb/I");
    fTree->Branch("MBVxbTyp",    &fBVxbTyp, "MBVxbTyp/I");
    fTree->Branch("nklmc",       &fKlMc, "nklmc/I");
    fTree->Branch("nklmcEMC",    &fKlMcEMC, "nklmcEMC/I");
    fTree->Branch("combibkd",    &fCombiBkd, "combibkd/I");
    fTree->Branch("bdaus",       &fBDaus, "bdaus/I");
    fTree->Branch("maxbdaus",    &fMaxBDaus, "maxbdaus/I");
    fTree->Branch("lundbdau",     fLundBDau,   "lundbdau[maxbdaus]/I");
    fTree->Branch("chamulti",    &fNFinalTrks, "chamulti/I");
    fTree->Branch("chamultiDCH", &fNFinalTrksDCH, "chamultiDCH/I");


    // Tina ntuple
    // -- TAU ANALYSIS BLOCK!
    // -- event prongs (max = 6)

    fTree->Branch("prongs",        &fProngs, "prongs/I");
    fTree->Branch("maxprongs",     &fMaxProngs, "maxprongs/I");

    fTree->Branch("nks",           &fNKshort, "nks/I");
    fTree->Branch("ThrTau",        &fThrTau,   "ThrTau/D");
    fTree->Branch("thThrTau",      &fthThrTau, "thThrTau/D");
    fTree->Branch("phThrTau",      &fphThrTau, "phThrTau/D");

    fTree->Branch("mprong",         ftrkM,   "mprong[maxprongs]/D");
    fTree->Branch("pprong",         ftrkP,   "pprong[maxprongs]/D");
    fTree->Branch("thetaprong",     ftrkTh,  "thetaprong[maxprongs]/D");
    fTree->Branch("phiprong",       ftrkPhi, "phiprong[maxprongs]/D");
    fTree->Branch("qprong",         ftrkCh,  "qprong[maxprongs]/I");

    fTree->Branch("mmiss",        &fMissM,  "mmiss/D");
    fTree->Branch("pmiss",        &fMissP,  "pmiss/D");
    fTree->Branch("thmiss",       &fMissTh, "thmiss/D");
    fTree->Branch("phmiss",       &fMissPh, "phmiss/D");

    //missing quantities with new neutral computation
    fTree->Branch("mmissLat",    &fMissMLat,  "mmissLat/D");
    fTree->Branch("pmissLat",    &fMissPLat,  "pmissLat/D");
    fTree->Branch("thmissLat",   &fMissThLat, "thmissLat/D");
    fTree->Branch("phmissLat",   &fMissPhLat, "phmissLat/D");
			       
    fTree->Branch("mmissAcc",    &fMissMAcc,  "mmissAcc/D");
    fTree->Branch("pmissAcc",    &fMissPAcc,  "pmissAcc/D");
    fTree->Branch("thmissAcc",   &fMissThAcc, "thmissAcc/D");
    fTree->Branch("phmissAcc",   &fMissPhAcc, "phmissAcc/D");
 
    fTree->Branch("mmissTight",  &fMissMTight,  "mmissTight/D");
    fTree->Branch("pmissTight",  &fMissPTight,  "pmissTight/D");
    fTree->Branch("thmissTight", &fMissThTight, "thmissTight/D");
    fTree->Branch("phmissTight", &fMissPhTight, "phmissTight/D");
    //----------------------------------------------------------//

    fTree->Branch("mprongs",      &fprongM,  "mprongs/D");
    fTree->Branch("pprongs",      &fprongP,  "pprongs/D");
    fTree->Branch("thprongs",     &fprongTh, "thprongs/D");
    fTree->Branch("phprongs",     &fprongPh, "phprongs/D");

    //--in the cms
 
    fTree->Branch("momprongcms",       ftrkcmsP,  "momprongcms[maxprongs]/D");
    fTree->Branch("thetaprongcms",     ftrkcmsTh, "thetaprongcms[maxprongs]/D");
    fTree->Branch("phiprongcms",       ftrkcmsPhi,"phiprongcms[maxprongs]/D");

    fTree->Branch("pprongscms",      &fcmsprongP,  "pprongscms/D");
    fTree->Branch("thprongscms",     &fcmsprongTh, "thprongscms/D");
    fTree->Branch("phprongscms",     &fcmsprongPh, "phprongscms/D");

    fTree->Branch("pmisscms",        &fcmsMissP,  "pmisscms/D");
    fTree->Branch("thmisscms",       &fcmsMissTh, "thmisscms/D");
    fTree->Branch("phmisscms",       &fcmsMissPh, "phmisscms/D");

    //missing quantities with new neutral computation
    fTree->Branch("pmisscmsLat",    &fcmsMissPLat,  "pmisscmsLat/D");
    fTree->Branch("thmisscmsLat",   &fcmsMissThLat, "thmisscmsLat/D");
    fTree->Branch("phmisscmsLat",   &fcmsMissPhLat, "phmisscmsLat/D");
				    
    fTree->Branch("pmisscmsAcc",    &fcmsMissPAcc,  "pmisscmsAcc/D");
    fTree->Branch("thmisscmsAcc",   &fcmsMissThAcc, "thmisscmsAcc/D");
    fTree->Branch("phmisscmsAcc",   &fcmsMissPhAcc, "phmisscmsAcc/D");

    fTree->Branch("pmisscmsTight",  &fcmsMissPTight,  "pmisscmsTight/D");
    fTree->Branch("thmisscmsTight", &fcmsMissThTight, "thmisscmsTight/D");
    fTree->Branch("phmisscmsTight", &fcmsMissPhTight, "phmisscmsTight/D");
    //-----------------------------------------------------------//

    //--PID bits
 
    fTree->Branch("elecId",      felecId, "elecId[maxprongs]/D");
    fTree->Branch("muonId",      fmuonId, "muonId[maxprongs]/D");
    fTree->Branch("kaonId",      fkaonId, "kaonId[maxprongs]/D");

    //-- gammas ---- with new cuts too -----
    fTree->Branch("maxgam",      &fMaxGam,      "maxgam/I");
    fTree->Branch("maskGamma",    fmaskGamma,   "maskGamma[maxgam]/I");
    fTree->Branch("ecalCluster",  fecalCluster, "ecalCluster[maxgam]/D");
    fTree->Branch("s9s25",        fs9s25,       "s9s25[maxgam]/D");
    fTree->Branch("lat",          flMom,        "lat[maxgam]/D");

    fTree->Branch("NeutMul",      &fNeutMult,      "NeutMul/I");
    fTree->Branch("ecalgam",      &fGamEcal,       "ecalgam/D");
    fTree->Branch("ecalgamnopi0", &fGamEcalNoPi0s, "ecalgamnopi0/D");

    fTree->Branch("NeutMulLat",      &fNeutMultLat,     "NeutMulLat/I");
    fTree->Branch("ecalgamLat",      &fGamEcalLat,      "ecalgamLat/D");
    fTree->Branch("ecalgamnopi0Lat", &fGamEcalNoPi0sLat,"ecalgamnopi0Lat/D");

    fTree->Branch("NeutMulAcc",      &fNeutMultAcc,     "NeutMulAcc/I");
    fTree->Branch("ecalgamAcc",      &fGamEcalAcc,      "ecalgamAcc/D");
    fTree->Branch("ecalgamnopi0Acc", &fGamEcalNoPi0sAcc,"ecalgamnopi0Acc/D");

    fTree->Branch("NeutMulTight",      &fNeutMultTight,     "NeutMulTight/I");
    fTree->Branch("ecalgamTight",      &fGamEcalTight,      "ecalgamTight/D");
    fTree->Branch("ecalgamnopi0Tight", &fGamEcalNoPi0sTight,"ecalgamnopi0Tight/D");

    fTree->Branch("mneu",         &fNeuM,  "mneu/D");
    fTree->Branch("pneu",         &fNeuP,  "pneu/D");
    fTree->Branch("thneu",        &fNeuTh, "thneu/D");
    fTree->Branch("phneu",        &fNeuPh, "phneu/D");

    fTree->Branch("mneuLat",      &fNeuMLat,  "mneuLat/D");
    fTree->Branch("pneuLat",      &fNeuPLat,  "pneuLat/D");
    fTree->Branch("thneuLat",     &fNeuThLat, "thneuLat/D");
    fTree->Branch("phneuLat",     &fNeuPhLat, "phneuLat/D");

    fTree->Branch("mneuAcc",      &fNeuMAcc,  "mneuAcc/D");
    fTree->Branch("pneuAcc",      &fNeuPAcc,  "pneuAcc/D");
    fTree->Branch("thneuAcc",     &fNeuThAcc, "thneuAcc/D");
    fTree->Branch("phneuAcc",     &fNeuPhAcc, "phneuAcc/D");

    fTree->Branch("mneuTight",    &fNeuMTight,  "mneuTight/D");
    fTree->Branch("pneuTight",    &fNeuPTight,  "pneuTight/D");
    fTree->Branch("thneuTight",   &fNeuThTight, "thneuTight/D");
    fTree->Branch("phneuTight",   &fNeuPhTight, "phneuTight/D");

    fTree->Branch("pneucms",      &fcmsNeuP,  "pneucms/D");
    fTree->Branch("thneucms",     &fcmsNeuTh, "thneucms/D");
    fTree->Branch("phneucms",     &fcmsNeuPh, "phneucms/D");

    fTree->Branch("pneucmsLat",   &fcmsNeuPLat,  "pneucmsLat/D");
    fTree->Branch("thneucmsLat",  &fcmsNeuThLat, "thneucmsLat/D");
    fTree->Branch("phneucmsLat",  &fcmsNeuPhLat, "phneucmsLat/D");

    fTree->Branch("pneucmsAcc",   &fcmsNeuPAcc,  "pneucmsAcc/D");
    fTree->Branch("thneucmsAcc",  &fcmsNeuThAcc, "thneucmsAcc/D");
    fTree->Branch("phneucmsAcc",  &fcmsNeuPhAcc, "phneucmsAcc/D");

    fTree->Branch("pneucmsTight", &fcmsNeuPTight,  "pneucmsTight/D");
    fTree->Branch("thneucmsTight",&fcmsNeuThTight, "thneucmsTight/D");
    fTree->Branch("phneucmsTight",&fcmsNeuPhTight, "phneucmsTight/D");

    //-- pi0's
    fTree->Branch("pi0s",        &fPi0s,   "pi0s/I");
    fTree->Branch("maxpi0s",     &fMaxPi0s,"maxpi0s/I");
    fTree->Branch("mpi0",        fpi0M,    "mpi0[maxpi0s]/D");
    fTree->Branch("ppi0",        fpi0P,    "ppi0[maxpi0s]/D");
    fTree->Branch("thetapi0",    fpi0Th,   "thetapi0[maxpi0s]/D");
    fTree->Branch("phipi0",      fpi0Ph,   "phipi0[maxpi0s]/D");
    
    fTree->Branch("ppi0cms",     fpi0cmsP,    "ppi0cms[maxpi0s]/D");
    fTree->Branch("thetapi0cms", fpi0cmsTh,   "thetapi0cms[maxpi0s]/D");
    fTree->Branch("phipi0cms",   fpi0cmsPh,   "phipi0cms[maxpi0s]/D");
    
    fTree->Branch("mpi0s",       &fpi0sM,    "mpi0s/D");
    fTree->Branch("ppi0s",       &fpi0sP,    "ppi0s/D");
    fTree->Branch("thetapi0s",   &fpi0sTh,   "thetapi0s/D");
    fTree->Branch("phipi0s",     &fpi0sPh,   "phipi0s/D");
    
    fTree->Branch("ppi0scms",     &fcmspi0sP,    "ppi0scms/D");
    fTree->Branch("thetapi0scms", &fcmspi0sTh,  "thetapi0scms/D");
    fTree->Branch("phipi0scms",   &fcmspi0sPh,   "phipi0scms/D");

  }
  
  fHistFile->cd();
}


// ----------------------------------------------------------------------
void recoilTaunu::threeProngs( int chbcand ) {

  // conditio sine qua non
  if ( chbcand == 0 && ! _selectB0 ) return;
  if ( chbcand != 0 && _selectB0 ) return;
  static bool first = true;
  if (first) {
    fHistFile->cd();
    fHistFile->mkdir("taunu", "taunu");
    fHistFile->cd("taunu");
    first = false;

    h_mES    =        new TH1F("mES","mes",50, 5.2, 5.3);
    h_mES030 =        new TH1F("mES030","mes030",50, 5.2, 5.3);
    h_mES040 =        new TH1F("mES040","mes040",50, 5.2, 5.3);
    h_mES045 =        new TH1F("mES045","mes045",50, 5.2, 5.3);
    h_mES050 =        new TH1F("mES050","mes050",50, 5.2, 5.3);
    h_mES055 =        new TH1F("mES055","mes055",50, 5.2, 5.3);
    h_mES060 =        new TH1F("mES060","mes060",50, 5.2, 5.3);
    h_mES065 =        new TH1F("mES065","mes065",50, 5.2, 5.3);
    h_mES070 =        new TH1F("mES070","mes070",50, 5.2, 5.3);
    h_mES080 =        new TH1F("mES080","mes080",50, 5.2, 5.3);
    h_mES090 =        new TH1F("mES090","mes090",50, 5.2, 5.3);

    h_deltaE =     new TH1F("deltaE","de",50, -0.12, 0.12);
    h_nProngs =    new TH1F("nProngs", "N. prongs",11, -0.5, 10.5);
    h_nPi0 =       new TH1F("nPi0", "N. pi0",11, -0.5, 10.5);
    h_nKs =        new TH1F("nKs", "N. Ks",11, -0.5, 10.5);
    h_nProngsSB =  new TH1F("nProngsSB", "N. prongs sideband",11, -0.5, 10.5);
    h_nPi0SB =     new TH1F("nPi0SB", "N. pi0 sideband",11, -0.5, 10.5);
    h_nKsSB =      new TH1F("nKsSB", "N. Ks sideband",11, -0.5, 10.5);
    h_mctruth =    new TH1F("mctruth", "mc truth code", 16, -0.5, 15.5);  
}
  fHistFile->cd();

  // inizializzazioni:
  // 4vettori  
  TLorentzVector p4Miss(0, 0, 0, 0);
  TLorentzVector p4MissLat(0, 0, 0, 0);
  TLorentzVector p4MissAcc(0, 0, 0, 0);
  TLorentzVector p4MissTight(0, 0, 0, 0);
  TLorentzVector p4Ups(0, 0, 0, 0);
  TLorentzVector p4BReco(0, 0, 0, 0);
  TLorentzVector p4Neu(0, 0, 0, 0);
  TLorentzVector p4NeuLat(0, 0, 0, 0);
  TLorentzVector p4NeuAcc(0, 0, 0, 0);
  TLorentzVector p4NeuTight(0, 0, 0, 0);
  TLorentzVector p4prongTot(0, 0, 0, 0);
  TLorentzVector p4Pi0Tot(0, 0, 0, 0);

  // mctruth
  fCombiBkd = 0 ;
  fBDaus = -1 ;
  fMaxBDaus = 0 ;
  fNFinalTrks = 0 ;
  fNFinalTrksDCH = 0 ;
  for (UInt_t kkh = 0; kkh < maxNBDau; ++kkh) {
    fLundBDau[kkh] = -77777 ;
  }

  // tracce cariche
  Int_t Brectrktmp;
  for (UInt_t jjh = 0; jjh < maxTracks; ++jjh)
  {
    ftrkM[jjh]      = -100.;
    ftrkP[jjh]      = -100.;
    ftrkTh[jjh]     = -100.;
    ftrkPhi[jjh]    = -100.0;
    ftrkCh[jjh]     = 0;
    ftrkcmsP[jjh]   = -100.;
    ftrkcmsTh[jjh]  = -100.;
    ftrkcmsPhi[jjh] = -100.0;
    
    felecId[jjh] = -100;
    fmuonId[jjh] = -100;
    fkaonId[jjh] = -100;
  }
  
  // pioni neutri  
  for (UInt_t iih = 0; iih < maxNPi0; ++iih)
  {
    fpi0M[iih]     = -100;
    fpi0P[iih]     = -100;
    fpi0Th[iih]    = -100;
    fpi0Ph[iih]    = -100;
    
    fpi0cmsP[iih]     = -100;
    fpi0cmsTh[iih]    = -100;
    fpi0cmsPh[iih]    = -100;
  }
  
  // cluster neutri e mask per i gamma
  int gam1 = 0;
  int gam2 = 0;

  for (Int_t gammaIndex = 0;  gammaIndex < maxGam; ++gammaIndex)
  {
    fmaskGamma[gammaIndex]=0;
    fecalCluster[gammaIndex]=0;
    fs9s25[gammaIndex]=0;
    flMom[gammaIndex]=0;
    gammaLock[gammaIndex]=0;
  }
  // globali
  p4BReco = p4Breco;
  p4Ups = p4Upsilon;
  fcosmyTh = -100.;

  // preselezione qualitativa sulle tracce del rinculo
  vector<Int_t> trkIndexes;
  for ( Int_t trkIndex = 0; trkIndex < nTrk; trkIndex++ )
  {
    Brectrktmp = chBRecTrk[trkIndex];
    if (Brectrktmp&brecoOverlap) continue;

    double theta = thetaTrk[ trkIndex ];
    if ( theta < 0.410 || theta > 2.54 ) continue;
    
    double p = momentumTrk[ trkIndex ];
    double pt = p * sin( theta );
    if ( pt < 0.1 || p > 10 ) continue;
    
    int ndch  = ndchTrk[trkIndex];
    if( ndch < 12 ) continue;
    
    double dz = fabs( zPocaTrk[ trkIndex ] );
    double dx = xPocaTrk[ trkIndex ];
    double dy = yPocaTrk[ trkIndex ];
    double dr = sqrt( dx * dx + dy * dy );
    if ( dz > 3 || dr > 1.5 ) continue;

    trkIndexes.push_back( trkIndex );
  }
  MomentumComparator c( momentumTrk );
  sort( trkIndexes.begin(), trkIndexes.end(), c );
  fProngs = trkIndexes.size();
  
  //  lista dei pi0 nel rinculo
  vector<Int_t> pi0Indexes;
  for( int pi0Index = 0; pi0Index < nPi0; pi0Index++ )
  {
    gam1 = d1Pi0Index[pi0Index] -1;
    gam2 = d2Pi0Index[pi0Index] -1;
    
    if ( (chBRecGam[gam1]&brecoOverlap) || (chBRecGam[gam2]&brecoOverlap) ) continue;
    bool overlap=false;
    for (vector<Int_t>::iterator it = pi0Indexes.begin(); it!=pi0Indexes.end(); it++)
    {
      if( fabs(pPi0[*it] - pPi0[pi0Index])<0.0001  && 
	  fabs(thPi0[*it] - thPi0[pi0Index])<0.0001  && 
	  fabs(phiPi0[*it] - phiPi0[pi0Index])<0.0001  )
      {
	overlap  = true;
	break;
      }
    }
    if (!overlap)
    pi0Indexes.push_back(pi0Index);
    // mask per i gamma
    gammaLock[gam1]=1;
    gammaLock[gam2]=1;
  }
  MomentumComparator p(pPi0);
  sort(pi0Indexes.begin(), pi0Indexes.end(), p );
  
  // lista dei gamma che non sono stati usati per fare i pi0
  /*
    vector<Int_t> gammaIndexes; 
    for (Int_t gammaIndex = 0; gammaIndex < nGam; ++gammaIndex) 
    if (  (gammaLock[gammaIndex] == 0)       && 
    (chBRecGam[gammaIndex] != brecoOv) &&
    (ecalGam[gammaIndex] > 0)          &&
    (energyGam[gammaIndex] >= 0.05)       )
    gammaIndexes.push_back(gammaIndex);
  */
  fPi0s = pi0Indexes.size();

  //************Charged tracks******************//
  //l'ordinamento per momento crescente potrebbe essere utile...
  fMaxProngs = 0;
  for ( vector<Int_t>::iterator it =  trkIndexes.begin(); it != trkIndexes.end() && fMaxProngs < maxTracks; it ++ )
  {
    TLorentzVector tmp;
    mk4Vector( tmp, momentumTrk[*it], thetaTrk[*it], phiTrk[*it], PIPMASS );
    
    double mass  = tmp.M2();
    mass         = ( mass > 0 ? sqrt( mass ) : -sqrt( - mass ) );
    ftrkM[fMaxProngs]   = mass;
    ftrkP[fMaxProngs]   = momentumTrk[*it];
    ftrkTh[fMaxProngs]  = thetaTrk[*it];
    ftrkPhi[fMaxProngs] = phiTrk[*it];
    ftrkCh[fMaxProngs]  = chargeTrk[*it];
    
    p4prongTot += tmp;
    tmp.Boost(-cmsBoost);
    
    ftrkcmsP[fMaxProngs]   = tmp.P();
    ftrkcmsTh[fMaxProngs]  = tmp.Theta();
    ftrkcmsPhi[fMaxProngs] = tmp.Phi();
    
    // store PID infos
    Int_t jbit;
    for (jbit = 0; jbit < 4; ++jbit)
      if (elecIdTrk[*it] & (0x1<<jbit)) felecId[fMaxProngs]=jbit;

    //    for (jbit = 0; jbit < 5; ++jbit)
    //      if (muonIdTrk[*it] & (0x1<<jbit)) fmuonId[fMaxProngs]=jbit;
    // patch for muonId in MC; verify only tight criterion    
    // comparing muonId bitmap with 01000 (=8) 
    // IF TIGHT ==> muonId = 1, IF LOOSE ==> muonId = 0,
    // IF VERYLOOSE ==> muonId = -2
 
    fmuonId[fMaxProngs] = 1;
    if ( ( (muonIdTrk[*it]) & ( 8 ) ) == 0 )  fmuonId[fMaxProngs]=0 ;
    if ( ( (muonIdTrk[*it]) & ( 4 ) ) == 0 )  fmuonId[fMaxProngs]=-1;
    if ( ( (muonIdTrk[*it]) & ( 2 ) ) == 0 )  fmuonId[fMaxProngs]=-2;
    for (jbit = 0; jbit < 4; ++jbit)
      if (kaonIdTrk[*it] & (0x1<<jbit)) fkaonId[fMaxProngs]=jbit;
    fMaxProngs++;
  }
  assert( fMaxProngs <= maxTracks );

  fThrTau   = ThruChB[0];
  fthThrTau = thThruChB[0];
  fphThrTau = phiThruChB[0] ;
  
  //**************Pi0s***********************//  
  fMaxPi0s = 0;
  for ( vector<Int_t>::iterator it = pi0Indexes.begin(); it != pi0Indexes.end() && fMaxPi0s < maxNPi0; it ++ )
  {
    TLorentzVector tmp;
    //    mk4Vector( tmp, pPi0[*it], thPi0[*it], phiPi0[*it], PIZMASS );
    // DO NOT WANT CONSTRAINT PI0 MASS
    mk4Vector( tmp, pPi0[*it], thPi0[*it], phiPi0[*it], m0Pi0[*it] );
    
    double mass  = tmp.M2();
    mass         = ( mass > 0 ? sqrt( mass ) : -sqrt( - mass ) );
    fpi0M[fMaxPi0s]  = mass;
    fpi0P[fMaxPi0s]  = pPi0[*it];
    fpi0Th[fMaxPi0s] = thPi0[*it];
    fpi0Ph[fMaxPi0s] = phiPi0[*it];
    
    p4Pi0Tot += tmp;
    tmp.Boost(-cmsBoost);

    fpi0cmsP[fMaxPi0s]  = tmp.P();
    fpi0cmsTh[fMaxPi0s] = tmp.Theta();
    fpi0cmsPh[fMaxPi0s] = tmp.Phi();
    fMaxPi0s++;
  }
  assert( fMaxPi0s <= maxNPi0 );

  //**************Gammas*****************//
  //MomentumComparator e( energyGam );
  //sort( gammaIndexes.begin(), gammaIndexes.end(), e );
  // for ( vector<Int_t>::iterator it = gammaIndexes.begin(); it < gammaIndexes.end(); it ++ ){
  //int i = *it;    
  // }
  fGamEcal = 0;
  fGamEcalNoPi0s = 0;
  fNeutMult = 0;
  TLorentzVector p4t(0,0,0,0);

  // neutral quantities with new lat and acceptance cuts
  fGamEcalLat = 0;
  fGamEcalNoPi0sLat = 0;
  fNeutMultLat = 0;
  TLorentzVector p4tLat(0,0,0,0);

  fGamEcalAcc = 0;
  fGamEcalNoPi0sAcc = 0;
  fNeutMultAcc = 0;
  TLorentzVector p4tAcc(0,0,0,0);

  fGamEcalTight = 0;
  fGamEcalNoPi0sTight = 0;
  fNeutMultTight = 0;
  TLorentzVector p4tTight(0,0,0,0);

  //---Neutral Energy Smearing-------//
  if(_smearNeutrals)
    smearNeutralEnergy();
  //---------------------------------//
  fMaxGam=0;
  //  cout << "nGam : " << nGam << endl;
  for ( int i = 0; (i<nGam) && (fMaxGam<maxGam); i++ )
  {    
    if ( (chBRecGam[i]&brecoOverlap) ) continue;
    if (ecalGam[i] < 0)              continue;

    //---all clusters info-----------//

    fmaskGamma[fMaxGam]=gammaLock[i];
    if(ecalGam[i]>=0)
      fecalCluster[fMaxGam]=ecalGam[i];
    else{
      fecalCluster[fMaxGam]=10000;
      //      cout<< "ecalGam[ " << i << " ] : " << ecalGam[i] <<endl;
      //      cout<< "s9s25Gam[ " << i << " ] : " << s9s25Gam[i] <<endl;
      //      cout<< "erawGam[ " << i << " ] : " << erawGam[i] <<endl;
    }
    fs9s25[fMaxGam]=s9s25Gam[i];
    flMom[fMaxGam]=lMomGam[i];
    fMaxGam++;

    //-------------------------------//
    if (energyGam[i] <= 0.05)        continue;
    fNeutMult++;
    fGamEcal += ecalGam[i];
    mk4Vector(p4t, energyGam[i], thetaGam[i], phiGam[i], 0.);
    p4Neu += p4t;
    if (  (gammaLock[i] == 0) ){
      fGamEcalNoPi0s += ecalGam[i];
    }

    //computing splitoff
    double phiMinNeu(1000.),thetaMinNeu(1000.);
    double thetaAtPhiMinNeu(1000.),phiAtThetaMinNeu(1000.);
    int splitOffGam = 0;
    for ( Int_t i1=0 ; i1<nTrk ; i1++ )
    {
      double dph=phiAtEMCTrk[i1]-phiGam[i];
      double dth=thetaAtEMCTrk[i1]-thetaGam[i];      
      if(TMath::Abs(dth) < thetaMinNeu)
      {
        thetaMinNeu = TMath::Abs(dth);
        phiAtThetaMinNeu = dph;
      }
      if(TMath::Abs(dph) < phiMinNeu)
      {
        phiMinNeu = TMath::Abs(dph);
        thetaAtPhiMinNeu = dth;
      }
      // select splitOff photons
      bool phicut = (chargeTrk[i1]>0 && dph>-0.03 && dph<0.07) || (chargeTrk[i1]<0 && dph>-0.07 && dph<0.03); 
      bool acc = thetaGam[i]> 0.41 && ecalTrk[i1] < 0; //track must not be matched
      if( phicut && acc && (TMath::Abs(dth)<0.03) )
	splitOffGam = 1;
    }

    // computing neutral quantities with semiexcl. recipe
    bool emcAccepted = false;   // emc acceptance + splitoff
    bool latCut = false;        // lat cut + s9s25 cut
    if ( lMomGam[i]>0.05 && lMomGam[i]<0.5 && s9s25Gam[i]> 0.9 ) 
      latCut = true;
    if ( (thetaGam[i] > 0.410) && (thetaGam[i] < 2.54) && splitOffGam==0 )
      emcAccepted = true;
    if ( latCut )
    {
      fNeutMultLat++;
      fGamEcalLat += ecalGam[i];
      mk4Vector(p4tLat, energyGam[i], thetaGam[i], phiGam[i], 0.);
      p4NeuLat += p4tLat;
      if (  (gammaLock[i] == 0) )
	fGamEcalNoPi0sLat += ecalGam[i];
    }      
    if ( emcAccepted )
    {
      fNeutMultAcc++;
      fGamEcalAcc += ecalGam[i];
      mk4Vector(p4tAcc, energyGam[i], thetaGam[i], phiGam[i], 0.);
      p4NeuAcc += p4tAcc;
      if (  (gammaLock[i] == 0) )
	fGamEcalNoPi0sAcc += ecalGam[i];
    }      
    if ( latCut && emcAccepted )
    {
      fNeutMultTight++;
      fGamEcalTight += ecalGam[i];
      mk4Vector(p4tTight, energyGam[i], thetaGam[i], phiGam[i], 0.);
      p4NeuTight += p4tTight;
      if (  (gammaLock[i] == 0) )
	fGamEcalNoPi0sTight += ecalGam[i];
    }
  }
  
  //  p4Miss = p4Ups - p4BReco - p4prongTot - p4Pi0Tot - p4Neu;
  p4Miss      = p4Ups - p4BReco - p4prongTot  - p4Neu;
  p4MissLat   = p4Ups - p4BReco - p4prongTot  - p4NeuLat;
  p4MissAcc   = p4Ups - p4BReco - p4prongTot  - p4NeuAcc;
  p4MissTight = p4Ups - p4BReco - p4prongTot  - p4NeuTight;
  
  fMissM  = p4Miss.M2();
  fMissM  = ( fMissM > 0 ? sqrt( fMissM) : -sqrt( - fMissM ) );
  fMissP  = p4Miss.P();
  fMissTh = p4Miss.Theta();
  fMissPh = p4Miss.Phi();

  //PMissing with new neutral computation
  fMissMLat  = p4MissLat.M2();
  fMissMLat  = ( fMissMLat > 0 ? sqrt( fMissMLat) : -sqrt( - fMissMLat ) );
  fMissPLat  = p4MissLat.P();
  fMissThLat = p4MissLat.Theta();
  fMissPhLat = p4MissLat.Phi();

  fMissMAcc  = p4MissAcc.M2();
  fMissMAcc  = ( fMissMAcc > 0 ? sqrt( fMissMAcc) : -sqrt( - fMissMAcc ) );
  fMissPAcc  = p4MissAcc.P();
  fMissThAcc = p4MissAcc.Theta();
  fMissPhAcc = p4MissAcc.Phi();

  fMissMTight  = p4MissTight.M2();
  fMissMTight  = ( fMissMTight > 0 ? sqrt( fMissMTight) : -sqrt( - fMissMTight ) );
  fMissPTight  = p4MissTight.P();
  fMissThTight = p4MissTight.Theta();
  fMissPhTight = p4MissTight.Phi();
  //------------------------------//

  fprongM  = p4prongTot.M2();
  fprongM  = ( fprongM > 0 ? sqrt( fprongM ) : -sqrt( - fprongM ) );
  fprongP  = p4prongTot.P();
  fprongTh = p4prongTot.Theta();
  fprongPh = p4prongTot.Phi();
  
  fNeuM  = p4Neu.M2();
  fNeuM  = ( fNeuM > 0 ? sqrt( fNeuM ) : -sqrt( - fNeuM ) );
  fNeuP  = p4Neu.P();
  fNeuTh = p4Neu.Theta();
  fNeuPh = p4Neu.Phi();

  // NEW neutral quantities
  fNeuMLat  = p4NeuLat.M2();
  fNeuMLat  = ( fNeuMLat > 0 ? sqrt( fNeuMLat ) : -sqrt( - fNeuMLat ) );
  fNeuPLat  = p4NeuLat.P();
  fNeuThLat = p4NeuLat.Theta();
  fNeuPhLat = p4NeuLat.Phi();

  fNeuMAcc  = p4NeuAcc.M2();
  fNeuMAcc  = ( fNeuMAcc > 0 ? sqrt( fNeuMAcc ) : -sqrt( - fNeuMAcc ) );
  fNeuPAcc  = p4NeuAcc.P();
  fNeuThAcc = p4NeuAcc.Theta();
  fNeuPhAcc = p4NeuAcc.Phi();

  fNeuMTight  = p4NeuTight.M2();
  fNeuMTight  = ( fNeuMTight > 0 ? sqrt( fNeuMTight ) : -sqrt( - fNeuMTight ) );
  fNeuPTight  = p4NeuTight.P();
  fNeuThTight = p4NeuTight.Theta();
  fNeuPhTight = p4NeuTight.Phi();
  //------------------------//

  frecoM  = p4BReco.M2();
  frecoM  = ( frecoM > 0 ? sqrt( frecoM ) : -sqrt( - frecoM ) );
  frecoP  = p4BReco.P();
  frecoTh = p4BReco.Theta();

  fpi0sM  = p4Pi0Tot.M2();
  fpi0sM  = ( fpi0M > 0 ? sqrt( fpi0sM ) : -sqrt( - fpi0sM ) );
  fpi0sP  = p4Pi0Tot.P();
  fpi0sTh = p4Pi0Tot.Theta();
  fpi0sPh = p4Pi0Tot.Phi();

  // compute cos(T_Breco^(p_prongs+p_pi0s))

  TVector3 thruVec(0,0,0);
  mk3Vector(thruVec, fThru, fthThru, fphThru);
  TLorentzVector p4PizPro(0,0,0,0);
  p4PizPro = p4prongTot + p4Pi0Tot;
  fcosmyTh = ( thruVec * (p4PizPro.Vect()) ) / (fThru * ( (p4PizPro.Vect()).Mag() ) );

  //boost the four vectors...
  p4prongTot.Boost(-cmsBoost);
  p4Neu.Boost(-cmsBoost);
  p4Miss.Boost(-cmsBoost);
  p4Pi0Tot.Boost(-cmsBoost);

  // boost also the new neutral quantity
  p4MissLat.Boost(-cmsBoost);
  p4MissAcc.Boost(-cmsBoost);
  p4MissTight.Boost(-cmsBoost);

  p4NeuLat.Boost(-cmsBoost);
  p4NeuAcc.Boost(-cmsBoost);
  p4NeuTight.Boost(-cmsBoost);

  // ...and store them.
  fcmsMissP  = p4Miss.P();
  fcmsMissTh = p4Miss.Theta();
  fcmsMissPh = p4Miss.Phi();

  //Pmissing with new neutral energy boost 
  fcmsMissPLat  = p4MissLat.P();
  fcmsMissThLat = p4MissLat.Theta();
  fcmsMissPhLat = p4MissLat.Phi();

  fcmsMissPAcc  = p4MissAcc.P();
  fcmsMissThAcc = p4MissAcc.Theta();
  fcmsMissPhAcc = p4MissAcc.Phi();

  fcmsMissPTight  = p4MissTight.P();
  fcmsMissThTight = p4MissTight.Theta();
  fcmsMissPhTight = p4MissTight.Phi();
  //----------------------------------//
      
  fcmsprongP  = p4prongTot.P();
  fcmsprongTh = p4prongTot.Theta();
  fcmsprongPh = p4prongTot.Phi();
      
  fcmsNeuP  = p4Neu.P();
  fcmsNeuTh = p4Neu.Theta();
  fcmsNeuPh = p4Neu.Phi();

  // store also the NEW neutral quantities
  fcmsNeuPLat  = p4NeuLat.P();
  fcmsNeuThLat = p4NeuLat.Theta();
  fcmsNeuPhLat = p4NeuLat.Phi();

  fcmsNeuPAcc  = p4NeuAcc.P();
  fcmsNeuThAcc = p4NeuAcc.Theta();
  fcmsNeuPhAcc = p4NeuAcc.Phi();

  fcmsNeuPTight  = p4NeuTight.P();
  fcmsNeuThTight = p4NeuTight.Theta();
  fcmsNeuPhTight = p4NeuTight.Phi();
  //-------------------------------//

  fcmspi0sP  = p4Pi0Tot.P();
  fcmspi0sTh = p4Pi0Tot.Theta();
  fcmspi0sPh = p4Pi0Tot.Phi();

  fPresel = true;
  if ( ( fProngs != 1 ) && ( fProngs != 3 ) ) 
    fPresel = false; 
  if ( fPi0s > 2 ) fPresel = false; 


  if ( fDeltaE > 0.08 || fDeltaE < -0.1 ) 
  {  
    fPresel = false;
    return;
  }

  // sbox : -0.1 < DE < 0.08 && mES > 5.27
  if ( fMes < 5.26 && fMes > 5.21 )
  {
    h_nProngsSB->Fill( fProngs );
    h_nPi0SB->Fill( fPi0s );
    h_nKsSB->Fill( fNKshort );
  }

  h_mES-> Fill( fMes );
  if(_genericBB)
  {
    if(fIntPurity>0.30)
      h_mES030-> Fill( fMes );
    if(fIntPurity>0.40)
      h_mES040-> Fill( fMes );
    if(fIntPurity>0.45)
      h_mES045-> Fill( fMes );
    if(fIntPurity>0.50)
      h_mES050-> Fill( fMes );
    if(fIntPurity>0.55)
      h_mES055-> Fill( fMes );
    if(fIntPurity>0.60)
      h_mES060-> Fill( fMes );
    if(fIntPurity>0.65)
      h_mES065-> Fill( fMes );
    if(fIntPurity>0.70)
      h_mES070-> Fill( fMes );
    if(fIntPurity>0.80)
      h_mES080-> Fill( fMes );
    if(fIntPurity>0.90)
      h_mES090-> Fill( fMes );
  }
  if ( fMes > 5.27 )
  {
    h_mctruth->Fill( fBVxbTyp );
    h_deltaE -> Fill( fDeltaE );
    h_nProngs->Fill( fProngs );
    h_nPi0->Fill( fPi0s );
    h_nKs->Fill( fNKshort );
  }

  //Let's study MC truth only if MC flag has been issued
  if(fIsMC) {

    // Study MCtruth B parent of the recoil
    // fCombiBkd = flag for combinatorial bkd 
    bool firstTrk = true ;
    Int_t indBProng = -1 ;
    // Study MC truth only for pre-selected events
    if (fPresel == true) {
      // loop on the tracks in the recoil
      for ( vector<Int_t>::iterator it =  trkIndexes.begin(); it != trkIndexes.end() ; it ++ ) {
	if ( fCombiBkd == -1 ) continue;
	
	// MC index of the track's mother
	Int_t mumInd = mothMc[IndexTrk[*it]-1] ;
	Int_t mumID = -7;
	//      cout << "mumID_INIT = " << mumID << endl;
	//      cout << " indBProng = " << indBProng << endl ;
	
	// iterate until we find the parent B
	bool foundB = false;
	Int_t iternum = 0 ;
	while ( foundB == false ) 
	  {
	    iternum ++;
	    assert(iternum < 100);
	    // lundID of the mother 
	    mumID = idMc[mumInd - 1] ;
	    //cout << "iteraz. a n. " << iternum << endl;
	    //cout << "mumID = " << mumID << endl ;
	    //cout << " mumInd = " << mumInd << endl ;
	    
	    // Sometimes the mother's index = 0... 
	    if (mumInd == 0) {
	      cout << " MC index of the track's mother = 0 !!!" << endl ;
	      fCombiBkd = -1 ;
	      foundB = true ;
	    } else {
	      // Is the mother a B+ or a B- ?
	      if(TMath::Abs(mumID) != 521) 
		{
		  // If not, go up another level
		  mumInd = mothMc[mumInd-1] ;
		} else {
		  // If so, stop going up the tree
		  foundB = true ;
		}
	    }
	  }
	
	//      cout << "indBProng = " << indBProng << endl ;
	if ((mumInd != indBProng) && (firstTrk == false) 
	    && (fCombiBkd != -1)) {
	  fCombiBkd = 1;
	}

	indBProng = mumInd ;
	firstTrk = false;
      }

      // Now, let's study the MC decay of the "correctly" reconstructed B recoil 
      // (i.e. not combinatorial bkd) 
      if (fCombiBkd == 0) {
	
	// number of daughters for correctly recontructed B recoil
	fBDaus = nDauMc[indBProng-1] ;      
	
	// lund ID of the first "maxNBDau" B daughters       
	for (Int_t imc = 0; imc < nMc && fMaxBDaus < maxNBDau ; ++imc) {    
	  if (mothMc[imc] == indBProng) { 
	    fLundBDau[fMaxBDaus] = idMc[imc] ;
	    fMaxBDaus++ ;
	  }
	  assert( fMaxBDaus <= maxNBDau );
	}
	
	// charged multiplicity of the decay
	for (Int_t imc = 0; imc < nMc ; ++imc) { 
	  if (fNFinalTrks==-1) continue;

	  // select only final state particles
	  if (nDauMc[imc]==0) {
	    
	    // select only charged pi, K or leptons 
	    Int_t idMCPart = TMath::Abs(idMc[imc]) ;
	    if ( (idMCPart == ELECLUND) || (idMCPart == MUONLUND) || 
		 (idMCPart == PIONLUND) || (idMCPart == KAONLUND) || 
		 (idMCPart == PROTLUND) ) {
	      
	      // Let's find the parent B of the final state trk
	      // We can replace this piece of code with a function (in the future...) 
	      bool foundB = false;
	      // MC index of the track's mother
	      Int_t mumInd = mothMc[imc] ;
	      Int_t mumID = -7;

	      // iterate until we find the parent B
	      Int_t iternum = 0 ;
	      while ( foundB == false ) 
		{
		  iternum ++;
		  assert(iternum < 100);
		  //cout << "iteraz. b n. " << iternum << endl;
		  //cout << "mumID = " << mumID << endl ;
		  //cout << " mumInd = " << mumInd << endl ;
		  
		  if (mumInd == 0) {
		    cout << " MC index of the track's mother = 0 !!!" << endl ;
		    fNFinalTrks = -1 ;
		    fNFinalTrksDCH = -1 ;
		    foundB = true ;
		  } else {
		    // lundID of the mother 
		    mumID = idMc[mumInd - 1] ;
		    // Is the mother a B+ or a B- ?
		    if(TMath::Abs(mumID) != 521) 
		      {
			// If not, go up another level
			mumInd = mothMc[mumInd-1] ;
		      } else {
			// If so, stop going up the tree
			foundB = true ;
		      }
		  }
		}
	      
	      //check that the parent B is the recoil B
	      if (mumInd == indBProng) {
		// charged multiplicity counting
		fNFinalTrks++ ;
		// charged multiplicity in the DCH acceptance
		if ( (thetaMc[imc] > 0.410 )  && ( thetaMc[imc] < 2.54 ) ) fNFinalTrksDCH++ ;
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }

  }

}

void
recoilTaunu::fillPidMaps() 
{
  for (Int_t i = 0; i < nTrk; ++i) {
    recKa[i] = 0; 
    recEl[i] = 0;
    recMu[i] = 0;
    // -- Kaons
    if (TMath::Abs(kaonIdTrk[i]) & IDKA) {
      recKa[i] = 1; 
      if (momentumTrk[i] < KAMOMLO) recKa[i] = 0; 
    } else {
      recKa[i] = 0; 
    }
    if (TMath::Abs(elecIdTrk[i]) & IDEL) {recEl[i] = 1; } 
    if ((TMath::Abs(muonIdTrk[i]) & IDMU) && ((TMath::Abs(kaonIdTrk[i]) & IDKA) == 0))  
      recMu[i] = 1;
  }
}
  
void
recoilTaunu::readCuts( TString filename, int dump )
{
  char  buffer[200];
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100], tablefile[1000];
  float CutValue;
  int ok(0);
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
     // -- breco 
    if (!strcmp(CutName, "mesSignalLo")) {MESSIGNALLO   = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalHi")) {MESSIGNALHI   = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSidebandLo")) {MESSIDEBANDLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSidebandHi")) {MESSIDEBANDHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalbandLo")) {MESSIGNALBANDLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "mesSignalbandHi")) {MESSIGNALBANDHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "deSignalLo")) {DESIGNALLO = CutValue; ok = 1;}
    if (!strcmp(CutName, "deSignalHi")) {DESIGNALHI = CutValue; ok = 1;}
    if (!strcmp(CutName, "purity")) {PURITY = CutValue; ok = 1;}
    if (!strcmp(CutName, "intPurity")) {INTPURITY = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDstar"))  {IPURDSTAR = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDc"))     {IPURDC = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurDstar0")) {IPURDSTAR0 = CutValue; ok = 1;}
    if (!strcmp(CutName, "ipurD0"))     {IPURD0 = CutValue; ok = 1;}
    if (dump == 1) dumpCuts();
  } 
  readintpur();
}


// ----------------------------------------------------------------------
void
recoilTaunu::dumpCuts() {
  cout << "====================================" << endl;
}

void recoilTaunu::breco( int chbcand ) 
{
  recoilNtpBase::breco( chbcand );

  if ((MESSIGNALLO < fMes) && (fMes < MESSIGNALHI)
      && (DESIGNALLO < fDeltaE) && (fDeltaE < DESIGNALHI)
      ) {
    signalBox = kTRUE;
  } else {
    signalBox = kFALSE;
  }

  
}

void
recoilTaunu::maskKs()
{
  fNKshort = 0;
  for (Int_t i = 0; i < nKs; ++i) {
    if (TMath::Abs(d1KsLund[i]) != 211) continue;
    if (TMath::Abs(d2KsLund[i]) != 211) continue;
    if ((d1KsIndex[i]-1 > nTrk) || (d2KsIndex[i]-1 > nTrk)) { 
      cout << "Daughter indices messed up: d1KsIndex = " << d1KsIndex[i]-1 << " d2KsIndex = " << d2KsIndex[i]-1 
	   << "  for nTrk = " << nTrk << endl;
      continue;
    }
    Int_t pi1 = d1KsIndex[i]-1;
    Int_t pi2 = d2KsIndex[i]-1;
    if (isRecEl(pi1)) continue;
    if (isRecEl(pi2)) continue;
    if (isRecMu(pi1)) continue;
    if (isRecMu(pi2)) continue;
    if (isRecKaon(pi1)) continue;
    if (isRecKaon(pi2)) continue;
    int overlap(0);
    Int_t Brectrktmp = B0RecTrk[pi1];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi1];
    if ((Brectrktmp&brecoOverlap)) {
      if(fVerbose) cout << "Rejecting Ks[" << i << "] due to pi1 at " << pi1 << " overlaps with BRECO " << endl;
      overlap++;
    }
    Brectrktmp = B0RecTrk[pi2];
    if(fBrecoCharge != 0) Brectrktmp = chBRecTrk[pi2];
    if ((Brectrktmp&brecoOverlap)) {
      if(fVerbose) cout << "Rejecting Ks[" << i << "] due to pi2 at " << pi2 << " overlaps with BRECO " << endl;
      overlap++;
    }
    if (overlap > 0) continue;
    fNKshort++;
  }
}

void
recoilTaunu::smearNeutralEnergy()
{
  int fOptSmearNeut = 0;
  if (fRunnumber < 100000) return; // no smearing for data
  if (fRunRange.Contains("Run 1")) {
    fOptSmearNeut = 1; 
    //    cout << "Smearing/ Shifting Neutrals with RUN1 recipe "  << endl;      
  } else if (fRunRange.Contains("Run 2a")) {
    fOptSmearNeut = 2; 
    //    cout << "Smearing/ Shifting Neutrals with RUN2 recipe "  << endl;      
  } else if (fRunRange.Contains("Run 2b")) {
    fOptSmearNeut = 2; 
    //    cout << "Smearing/ Shifting Neutrals with RUN2 recipe "  << endl;      
  } else {
    //      cout << "Warning: No run range determined. Taking RUN2 as default!" << endl;
      fOptSmearNeut = 2; 
  }
  double SIGMANEUT;

  for (int i = 0; i < nGam; ++i) 
  {
    SIGMANEUT = 0;
    if(energyGam[i]< 0.1 ) {
      SIGMANEUT =  0.03;
    } else if(energyGam[i]< 0.3 ) {
      SIGMANEUT =  0.026;
    } else if(energyGam[i]< 0.6 ) {
      SIGMANEUT =  (fOptSmearNeut==1) ? 0.024: 0.016;
    } else if(energyGam[i]< 1. ) {     
      SIGMANEUT =  (fOptSmearNeut==1) ? 0.020: 0.016;
    }
    if ( energyGam[i] > 0 )
    {
      double smear = gRandom->Gaus( 0, SIGMANEUT ); 
      energyGam[i] = clusterReCorrection(energyGam[i], thetaGam[i]) + smear * energyGam[i];
      ecalGam[i] = clusterReCorrection(ecalGam[i], thetaGam[i]) + smear * ecalGam[i];
    }  
  }

  return;
}

//------------------------------------------------------------------------------//
//
// Description: Function for the Correction of the Cluster Energy
//              => returns the corrected cluster energy
//
//   clusterCorrection(double rawEnergy, double clusterPostionTheta, bool newCorr=true)
//     => returns corrected energy with
//        - constants from March    2002 ("new") if newCorr=true (default)
//        - constants from November 1999 ("old") if newCorr=false
//
//   clusterReCorrection(oldCorrectedEnergy, clusterPostionTheta)
//     => convert old corrected energy
//             to new corrected energy
//
//              rawEnergy           - raw cluster energy (in GeV)
//              clusterPostionTheta - theta of the cluster positon (in rad)
//              oldCorrectedEnergy  - old corrected cluster energy (in GeV)
//
// Author : Enrico Maly  25 Apr 2002
//
//------------------------------------------------------------------------------//

double
recoilTaunu::clusterCorrection(const double rawEnergy,
                  const double clusterPositionTheta, 
                  const bool   newCorr) const
{

// constants from April 2002
 double constants[19] = {
   +9.046e-01,
   +1.243e-02,
   +7.785e-03,
   -2.178e-03,
   +1.620e-02,
   -5.686e-03,
   +2.063e-04,
   +1.408e-01,
   -1.133e-01,
   +2.214e-02,
   +8.409e-03,
   -1.650e-02,
   +5.301e-03,
   +4.998e-02,
   -2.964e-02,
   +9.960e-01,
   -9.097e-02,
   +4.351e-02,
   -5.891e-03
 };
 
 // constants from November 1999
 double oldConstants[19] = {
   +1.024e-00,
   -1.093e-01,
   +4.528e-02,
   -5.959e-03,
   +3.955e-04,
   +4.375e-04,
   -2.855e-04,
   +1.643e-02,
   -1.881e-02,
   +4.838e-03,
   +1.583e-02,
   -1.680e-02,
   +5.074e-03,
   +1.989e-02,
   -1.907e-02,
   +1.133e-00,
   -2.420e-01,
   +9.396e-02,
   -1.142e-02,
 };

  const double lgE = log10(1000*rawEnergy);
  const double cosTh = cos(clusterPositionTheta);

  double *c;
  if (newCorr) c=constants;
  else c=oldConstants;

  double eCorr=rawEnergy;

  // barrel section
  if (cosTh<=0.892)
    {

      const double p0 = c[0]+c[1]*lgE+c[2]*lgE*lgE+c[3]*pow(lgE,3);
      const double p1 = c[4]+c[5]*lgE+c[6]*lgE*lgE;
      const double p2 = c[7]+c[8]*lgE+c[9]*lgE*lgE;
      const double p3 = c[10]+c[11]*lgE+c[12]*lgE*lgE;
      const double p4 = c[13]+c[14]*lgE;

      const double correctionForBarrel = p0+p1*cosTh+p2*cosTh*cosTh+p3*pow(cosTh,3)+p4*pow(cosTh,4);
      eCorr = rawEnergy/correctionForBarrel;

    }
  // endcap section
  else if (cosTh>0.892)
    {
      const double correctionForEndcap = c[15]+c[16]*lgE+c[17]*lgE*lgE+c[18]*pow(lgE,3);
      eCorr = rawEnergy/correctionForEndcap;
    }

  return eCorr;
}

double
recoilTaunu::clusterReCorrection(const double oldCorrectedEnergy, 
                    const double clusterPositionTheta) const
{

  double eCorr=oldCorrectedEnergy;
  const double eRaw=eCorr*eCorr/clusterCorrection(eCorr,clusterPositionTheta,false);
  eCorr=clusterCorrection(eRaw,clusterPositionTheta,true);

  return eCorr;

}




