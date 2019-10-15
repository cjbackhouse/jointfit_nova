#pragma once

#include "CAFAna/Extrap/IExtrap.h"

namespace ana
{
  class Loaders;

  /// "Extrapolation" that simply returns the FD MC prediction
  class TrivialExtrap: public IExtrap
  {
  public:
    virtual OscillatableSpectrum NueSurvComponent()       override {return fNueSurv;}
    virtual OscillatableSpectrum AntiNueSurvComponent()   override {return fNueSurvAnti;}

    virtual OscillatableSpectrum NumuSurvComponent()      override {return fNumuSurv;}
    virtual OscillatableSpectrum AntiNumuSurvComponent()  override {return fNumuSurvAnti;}

    virtual OscillatableSpectrum NueAppComponent()        override {return fNueApp;}
    virtual OscillatableSpectrum AntiNueAppComponent()    override {return fNueAppAnti;}

    virtual OscillatableSpectrum NumuAppComponent()       override {return fNumuApp;}
    virtual OscillatableSpectrum AntiNumuAppComponent()   override {return fNumuAppAnti;}

    virtual OscillatableSpectrum TauFromEComponent()      override {return fTauFromE;}
    virtual OscillatableSpectrum AntiTauFromEComponent()  override {return fTauFromEAnti;}

    virtual OscillatableSpectrum TauFromMuComponent()     override {return fTauFromMu;}
    virtual OscillatableSpectrum AntiTauFromMuComponent() override {return fTauFromMuAnti;}

    virtual Spectrum NCTotalComponent() override {return fNCTot;}
    virtual Spectrum NCComponent()      override {return fNC;}
    virtual Spectrum NCAntiComponent()  override {return fNCAnti;}

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<TrivialExtrap> LoadFrom(TDirectory* dir);

  protected:
    TrivialExtrap()
      : fNueApp(0, {}, {}, 0, 0),    fNueAppAnti(0, {}, {}, 0, 0),
        fNumuSurv(0, {}, {}, 0, 0),  fNumuSurvAnti(0, {}, {}, 0, 0),
        fNumuApp(0, {}, {}, 0, 0),   fNumuAppAnti(0, {}, {}, 0, 0),
        fNueSurv(0, {}, {}, 0, 0),   fNueSurvAnti(0, {}, {}, 0, 0),
        fTauFromE(0, {}, {}, 0, 0),  fTauFromEAnti(0, {}, {}, 0, 0),
        fTauFromMu(0, {}, {}, 0, 0), fTauFromMuAnti(0, {}, {}, 0, 0),
        fNCTot(0, {}, {}, 0, 0),
        fNC(0, {}, {}, 0, 0),        fNCAnti(0, {}, {}, 0, 0)
    {}

    OscillatableSpectrum fNueApp,    fNueAppAnti;
    OscillatableSpectrum fNumuSurv,  fNumuSurvAnti;
    OscillatableSpectrum fNumuApp,   fNumuAppAnti;
    OscillatableSpectrum fNueSurv,   fNueSurvAnti;
    OscillatableSpectrum fTauFromE,  fTauFromEAnti;
    OscillatableSpectrum fTauFromMu, fTauFromMuAnti;
    Spectrum fNCTot;
    Spectrum fNC, fNCAnti;
  };
}
