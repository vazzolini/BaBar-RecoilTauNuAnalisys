#ifndef recoilTaunu_h
#define recoilTaunu_h

#include "recoilNtpBase.hh"

class TH1F ; 
class recoilTaunu : public recoilNtpBase {


  Int_t fNKshort;
  Int_t fKlMc, fKlMcEMC;
  Int_t fNeutMult;
  Int_t fNeutMultLat;
  Int_t fNeutMultAcc;
  Int_t fNeutMultTight;

  Bool_t fPresel;
  Bool_t fPGoodChargeCorr, fPGoodEvent, fnoNeu, fnoChg; 
  Double_t fElab; 
  Bool_t fPion, fKaon;
  Int_t fNPions, fNKaons; 
  Int_t fTruTrk, fMomTrk;
  Double_t fcosmyTh;
  Double_t fPchpa, fTchpa, fFchpa, fEchpa, fMchpa;
  Double_t fThruLep;
  Int_t fPionCharge;
  Int_t fRecTrk;
  Bool_t fOnePion;
  Double_t fpnuRes, ftnuRes, fq2Res, fxmassRes;
  Int_t fBtaunu;
  Bool_t fisEl, fisPi, fisMu, fisOther; 
  Int_t fDstar0, fD0, fDc ,fDstar;
  TLorentzVector p4DauTau; 
  Int_t fNTauMC, fnElectron, fnMuon, fnPion;
  Double_t fPTauGen, fTTauGen, fFTauGen, fETauGen, fMTauGen, fPDauTau, fTDauTau, fFDauTau, fEDauTau, fMDauTau; 

  Double_t fPPlab, fPTlab, fPFlab, fPElab, fPPcms, fPTcms, fPFcms, fPEcms;
  Bool_t fPIsPrompt, fPIsCascade;
  // MC TRUTH 
  Int_t fB1Index, fB2Index;
  Int_t fBVxb, fBVxbTyp ;
  Int_t fCombiBkd ;
  Int_t fBDaus, fMaxBDaus;
  Int_t fNFinalTrks, fNFinalTrksDCH;
  static const UInt_t maxNBDau = 4 ;
  Int_t fLundBDau[maxNBDau];

  // Tina's Tree variables.
  static const UInt_t maxTracks = 3;

  UInt_t fProngs ; 
  UInt_t fMaxProngs;
 
  Double_t fThrTau, fthThrTau, fphThrTau;
  Double_t ftrkM[maxTracks];
  Double_t ftrkP[maxTracks];
  Double_t ftrkTh[maxTracks];
  Double_t ftrkPhi[maxTracks];
  Double_t ftrkcmsP[maxTracks];
  Double_t ftrkcmsTh[maxTracks];
  Double_t ftrkcmsPhi[maxTracks];
  Int_t ftrkCh[maxTracks];
  //  Double_t ftrkMass;
  Double_t felecId[maxTracks];
  Double_t fmuonId[maxTracks];
  Double_t fkaonId[maxTracks];

  static const Int_t maxGam = 1000;
  UInt_t fMaxGam;

  Double_t fecalCluster[maxGam];
  Double_t fs9s25[maxGam];
  Double_t flMom[maxGam];
  Int_t fmaskGamma[maxGam];
  Int_t gammaLock[maxGam];

  Double_t fGamEcal;
  Double_t fGamEcalLat;
  Double_t fGamEcalAcc;
  Double_t fGamEcalTight;
  Double_t fGamEcalNoPi0s;
  Double_t fGamEcalNoPi0sLat;
  Double_t fGamEcalNoPi0sAcc;
  Double_t fGamEcalNoPi0sTight;

  Double_t fMissM ; 
  Double_t fMissP ;
  Double_t fMissTh;
  Double_t fMissPh;

  Double_t fMissMLat ; 
  Double_t fMissPLat ;
  Double_t fMissThLat;
  Double_t fMissPhLat;
  Double_t fMissMAcc ; 
  Double_t fMissPAcc ;
  Double_t fMissThAcc;
  Double_t fMissPhAcc;
  Double_t fMissMTight ; 
  Double_t fMissPTight ;
  Double_t fMissThTight;
  Double_t fMissPhTight;

  Double_t fprongM ; 
  Double_t fprongP ;
  Double_t fprongTh;
  Double_t fprongPh;

  Double_t fNeuM ; 
  Double_t fNeuP ;
  Double_t fNeuTh;
  Double_t fNeuPh;

  Double_t fNeuMLat ; 
  Double_t fNeuPLat ;
  Double_t fNeuThLat;
  Double_t fNeuPhLat;
  Double_t fNeuMAcc ; 
  Double_t fNeuPAcc ;
  Double_t fNeuThAcc;
  Double_t fNeuPhAcc;
  Double_t fNeuMTight ; 
  Double_t fNeuPTight ;
  Double_t fNeuThTight;
  Double_t fNeuPhTight;

  Double_t frecoM ; 
  Double_t frecoP ;
  Double_t frecoTh;
  Double_t frecoPh;

  Double_t fcmsMissP ;
  Double_t fcmsMissTh;
  Double_t fcmsMissPh;

  Double_t fcmsMissPLat ;
  Double_t fcmsMissThLat;
  Double_t fcmsMissPhLat;
  Double_t fcmsMissPAcc ;
  Double_t fcmsMissThAcc;
  Double_t fcmsMissPhAcc;
  Double_t fcmsMissPTight ;
  Double_t fcmsMissThTight;
  Double_t fcmsMissPhTight;

  Double_t fcmsprongP ;
  Double_t fcmsprongTh;
  Double_t fcmsprongPh;
	       
  Double_t fcmsNeuP ;
  Double_t fcmsNeuTh;
  Double_t fcmsNeuPh;

  Double_t fcmsNeuPLat ;
  Double_t fcmsNeuThLat;
  Double_t fcmsNeuPhLat;
  Double_t fcmsNeuPAcc ;
  Double_t fcmsNeuThAcc;
  Double_t fcmsNeuPhAcc;
  Double_t fcmsNeuPTight ;
  Double_t fcmsNeuThTight;
  Double_t fcmsNeuPhTight;

  static const UInt_t maxNPi0 = 2; // troppi pi0 fake....
  UInt_t fPi0s;
  UInt_t fMaxPi0s;

  Double_t fpi0M[maxNPi0];
  Double_t fpi0P[maxNPi0];
  Double_t fpi0Th[maxNPi0];
  Double_t fpi0Ph[maxNPi0];
  
  Double_t fpi0cmsP[maxNPi0];
  Double_t fpi0cmsTh[maxNPi0];
  Double_t fpi0cmsPh[maxNPi0];
  
  Double_t fpi0sM;
  Double_t fpi0sP;
  Double_t fpi0sTh;
  Double_t fpi0sPh;
  
  Double_t fcmspi0sP;
  Double_t fcmspi0sTh;
  Double_t fcmspi0sPh;

  // histograms
  TH1F *h_mES,
    *h_mES030,
    *h_mES040,
    *h_mES045,
    *h_mES050,
    *h_mES055,
    *h_mES060,
    *h_mES065,
    *h_mES070,
    *h_mES080,
    *h_mES090,
    *h_deltaE,
    *h_nProngs,
    *h_nPi0,
    *h_nKs, 
    *h_nProngsSB,
    *h_nPi0SB,
    *h_nKsSB, 
    *h_mctruth;
  
  bool _selectB0;
  bool _smearNeutrals;
  bool _genericBB;
public :

  recoilTaunu(TTree *tree=0, int isMC =0, int newFormat=2, bool selectB0 = false,
	      bool smearNeutrals = false, bool genericBBbar = false, int oneProng = 0);
  recoilTaunu(TString, TString, int isMC=0, int newFormat=2, bool selectB0 = false,
	      bool smearNeutrals = false,  bool genericBBbar = false, int oneProng = 0);
  //  ~recoilTaunu();

  void      readCuts( TString filename, int dump = 1);
  void      dumpCuts();
  void      mcTruth( int chbcand );
  void      recoil( int chbcand );
  void      breco( int chbcand );
  void      fillRecoilHist(const char *, int chbcand);
  void      bookHist(int dump = 0);
  void      threeProngs( int chbcand );
  void      fillPidMaps();
  void      smearNeutralEnergy();
  double    clusterCorrection(const double rawEnergy,
			      const double clusterPositionTheta, 
			      const bool   newCorr) const;
  double    clusterReCorrection(const double oldCorrectedEnergy, 
				const double clusterPositionTheta) const;
  void maskKs();
  
};

#endif


