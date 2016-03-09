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

  void     SetClosestGen(Int_t i)      { fMatchId1 = i; }
  void     SetClosestPF(Int_t i)       { fMatchId2  = i; }

  Float_t     GetdEta() {return fdEtaAtVtx ; }
  Float_t     GetdPhi() {return fdPhiAtVtx ; }
  Float_t     GetHoverE() {return fHoverE ; }
  Float_t     GetD0() {return fDxy ; }
  Float_t     GetDZ() {return fDz ; }
  Float_t     GetSigmaIEtaIEta() {return fSigmaIEtaIEta ; }
  Float_t     GetEoverPInv() {return fEoverPInv ; }
  Int_t       GetMissHits() {return fMissHits ; }
  Bool_t      GetConversionVeto() {return fConversionVeto ; }


  void     SetdEta(float val) {fdEtaAtVtx = val; }
  void     SetdPhi(float val) {fdPhiAtVtx = val; }
  void     SetHoverE(float val) {fHoverE = val; }
  void     SetD0(float val) {fDxy = val; }
  void     SetDZ(float val) {fDz = val; }
  void     SetSigmaIEtaIEta(float val) {fSigmaIEtaIEta = val; }
  void     SetEoverPInv(float val) {fEoverPInv = val; }
  void     SetMissHits(int val) {fMissHits = val; }
  void     SetConversionVeto(bool val) {fConversionVeto = val; }

  Int_t    GetClosestGen()   const     { return fMatchId1; }
  Int_t    GetClosestPF()    const     { return fMatchId2; }
  
 protected:
  
  Float_t                      fdEtaAtVtx;   // Please add description!
  Float_t                      fdPhiAtVtx;   // 
  Float_t                      fSigmaIEtaIEta;
  Float_t                      fHoverE;     // 
  Float_t                      fDxy;        // 
  Float_t                      fDz;        //
  Float_t                      fEoverPInv;  // 
  Int_t                        fMissHits;   // 
  Bool_t                       fConversionVeto;  // 
  
  ClassDef(lwElectron,1)
    };
#endif
