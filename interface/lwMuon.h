#ifndef lwMuon_h
#define lwMuon_h

//
// light-weight muon candidate
//

#include "UserCode/diall/interface/particleBase.h"

class lwMuon : public particleBase {
 public:
  lwMuon();
  lwMuon(Double_t pt, Double_t eta, Double_t phi, Double_t m, Int_t id = -1);
  virtual ~lwMuon() {;}
  lwMuon(const lwMuon& obj); // copy constructor
  lwMuon& operator=(const lwMuon& other); // assignment

  Int_t       GetClosestGen()   const     { return fMatchId1; }
  Int_t       GetClosestPF()    const     { return fMatchId2; }
  Float_t     GetPFChIso() {return fPFChIso ; }
  Float_t     GetPFPhoIso() {return fPFPhoIso ; }
  Float_t     GetPFNeuIso() {return fPFNeuIso ; }
  Float_t     GetPFPUIso() {return fPFPUIso ; }

  void        SetClosestGen(Int_t i)      { fMatchId1 = i; }
  void        SetClosestPF(Int_t i)       { fMatchId2  = i; }
  void        SetPFChIso(float val) {fPFChIso = val; }
  void        SetPFPhoIso(float val) {fPFPhoIso = val; }
  void        SetPFNeuIso(float val) {fPFNeuIso = val; }
  void        SetPFPUIso(float val) {fPFPUIso = val; }
  
 protected:
  
  Float_t                      fPFChIso;
  Float_t                      fPFPhoIso;
  Float_t                      fPFNeuIso;
  Float_t                      fPFPUIso;

  ClassDef(lwMuon,1)
};
#endif
