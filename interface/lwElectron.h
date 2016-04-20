#ifndef lwElectron_h
#define lwElectron_h

//
// light-weight electron candidate
//

#include "UserCode/diall/interface/particleBase.h"

class lwElectron : public particleBase {
 public:
  lwElectron();
  lwElectron(Double_t pt, Double_t eta, Double_t phi, Double_t m, Int_t id = -1);
  virtual ~lwElectron() {;}
  lwElectron(const lwElectron& obj); // copy constructor
  lwElectron& operator=(const lwElectron& other); // assignment

  Int_t       GetVetoId() {return fVetoId ; }
  Int_t       GetLooseId() {return fLooseId ; }
  Int_t       GetMediumId() {return fMediumId ; }
  Int_t       GetTightId() {return fTightId ; }
  Float_t     GetdEta() {return fdEtaAtVtx ; }
  Float_t     GetdPhi() {return fdPhiAtVtx ; }
  Float_t     GetHoverE() {return fHoverE ; }
  Float_t     GetD0() {return fDxy ; }
  Float_t     GetDZ() {return fDz ; }
  Float_t     GetSigmaIEtaIEta() {return fSigmaIEtaIEta ; }
  Float_t     GetEoverPInv() {return fEoverPInv ; }
  Int_t       GetMissHits() {return fMissHits ; }
  Int_t       GetConversionVeto() {return fConversionVeto ; }
  Float_t     GetPFChIso() {return fPFChIso ; }
  Float_t     GetPFPhoIso() {return fPFPhoIso ; }
  Float_t     GetPFNeuIso() {return fPFNeuIso ; }
  Float_t     GetPFPUIso() {return fPFPUIso ; }
  Float_t     GetPFChIso03() {return fPFChIso03 ; }
  Float_t     GetPFPhoIso03() {return fPFPhoIso03 ; }
  Float_t     GetPFNeuIso03() {return fPFNeuIso03 ; }
  Float_t     GetPFChIso04() {return fPFChIso04 ; }
  Float_t     GetPFPhoIso04() {return fPFPhoIso04 ; }
  Float_t     GetPFNeuIso04() {return fPFNeuIso04 ; }
  Float_t     GetEffAreaTimesRho() {return fEffAreaTimesRho ; }
  
  void        SetVetoId(int val) {fVetoId = val; }
  void        SetLooseId(int val) {fLooseId = val; }
  void        SetMediumId(int val) {fMediumId = val; }
  void        SetTightId(int val) {fTightId = val; }
  void        SetdEta(float val) {fdEtaAtVtx = val; }
  void        SetdPhi(float val) {fdPhiAtVtx = val; }
  void        SetHoverE(float val) {fHoverE = val; }
  void        SetD0(float val) {fDxy = val; }
  void        SetDZ(float val) {fDz = val; }
  void        SetSigmaIEtaIEta(float val) {fSigmaIEtaIEta = val; }
  void        SetEoverPInv(float val) {fEoverPInv = val; }
  void        SetMissHits(int val) {fMissHits = val; }
  void        SetConversionVeto(bool val) {fConversionVeto = val; }
  void        SetPFChIso(float val) {fPFChIso = val; }
  void        SetPFPhoIso(float val) {fPFPhoIso = val; }
  void        SetPFNeuIso(float val) {fPFNeuIso = val; }
  void        SetPFPUIso(float val) {fPFPUIso = val; }
  void        SetPFChIso03(float val) {fPFChIso03 = val; }
  void        SetPFPhoIso03(float val) {fPFPhoIso03 = val; }
  void        SetPFNeuIso03(float val) {fPFNeuIso03 = val; }
  void        SetPFChIso04(float val) {fPFChIso04 = val; }
  void        SetPFPhoIso04(float val) {fPFPhoIso04 = val; }
  void        SetPFNeuIso04(float val) {fPFNeuIso04 = val; }


  void        SetEffAreaTimesRho(float val) {fEffAreaTimesRho = val; }

  void        SetClosestGen(Int_t i)      { fMatchId1 = i; }
  void        SetClosestPF(Int_t i)       { fMatchId2  = i; }
  Int_t       GetClosestGen()   const     { return fMatchId1; }
  Int_t       GetClosestPF()    const     { return fMatchId2; }
  
 protected:
  
  Int_t                        fVetoId;
  Int_t                        fLooseId;
  Int_t                        fMediumId;
  Int_t                        fTightId;
  Float_t                      fdEtaAtVtx;   // Please add description!
  Float_t                      fdPhiAtVtx;   // 
  Float_t                      fSigmaIEtaIEta;
  Float_t                      fHoverE;     // 
  Float_t                      fDxy;        // 
  Float_t                      fDz;        //
  Float_t                      fEoverPInv;  // 
  Int_t                        fMissHits;   // 
  Int_t                        fConversionVeto;  // 
  Float_t                      fPFChIso;
  Float_t                      fPFPhoIso;
  Float_t                      fPFNeuIso;
  Float_t                      fPFPUIso;
  Float_t                      fPFChIso03;
  Float_t                      fPFPhoIso03;
  Float_t                      fPFNeuIso03;
  Float_t                      fPFChIso04;
  Float_t                      fPFPhoIso04;
  Float_t                      fPFNeuIso04;

  Float_t                      fEffAreaTimesRho;
  
  ClassDef(lwElectron,1)
    };
#endif
