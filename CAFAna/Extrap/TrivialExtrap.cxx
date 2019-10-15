#include "CAFAna/Extrap/TrivialExtrap.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  void TrivialExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("TrivialExtrap").Write("type");

    fNueApp.SaveTo(dir->mkdir("nue_app"));
    fNueAppAnti.SaveTo(dir->mkdir("nue_app_anti"));
    fNCTot.SaveTo(dir->mkdir("nc_tot"));
    fNC.SaveTo(dir->mkdir("nc"));
    fNCAnti.SaveTo(dir->mkdir("nc_anti"));
    fNumuSurv.SaveTo(dir->mkdir("numu_surv"));
    fNumuSurvAnti.SaveTo(dir->mkdir("numu_surv_anti"));
    fNumuApp.SaveTo(dir->mkdir("numu_app"));
    fNumuAppAnti.SaveTo(dir->mkdir("numu_app_anti"));
    fNueSurv.SaveTo(dir->mkdir("nue_surv"));
    fNueSurvAnti.SaveTo(dir->mkdir("nue_surv_anti"));
    fTauFromE.SaveTo(dir->mkdir("nutau_from_nue"));
    fTauFromEAnti.SaveTo(dir->mkdir("nutau_from_nue_anti"));
    fTauFromMu.SaveTo(dir->mkdir("nutau_from_numu"));
    fTauFromMuAnti.SaveTo(dir->mkdir("nutau_from_numu_anti"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TrivialExtrap> TrivialExtrap::LoadFrom(TDirectory* dir)
  {
    std::unique_ptr<TrivialExtrap> ret(new TrivialExtrap);

    // This is a lot of repetitive typing. Define some macros
#define LOAD_OSC(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *OscillatableSpectrum::LoadFrom(dir->GetDirectory(LABEL));
#define LOAD_SPECT(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *Spectrum::LoadFrom(dir->GetDirectory(LABEL));

    LOAD_OSC(fNueApp,        "nue_app");
    LOAD_OSC(fNueAppAnti,    "nue_app_anti");
    LOAD_OSC(fNumuSurv,      "numu_surv");
    LOAD_OSC(fNumuSurvAnti,  "numu_surv_anti");
    LOAD_OSC(fNumuApp,       "numu_app");
    LOAD_OSC(fNumuAppAnti,   "numu_app_anti");
    LOAD_OSC(fNueSurv,       "nue_surv");
    LOAD_OSC(fNueSurvAnti,   "nue_surv_anti");
    LOAD_OSC(fTauFromE,      "nutau_from_nue");
    LOAD_OSC(fTauFromEAnti,  "nutau_from_nue_anti");
    LOAD_OSC(fTauFromMu,     "nutau_from_numu");
    LOAD_OSC(fTauFromMuAnti, "nutau_from_numu_anti");

    LOAD_SPECT(fNCTot,  "nc_tot");
    LOAD_SPECT(fNC,     "nc");
    LOAD_SPECT(fNCAnti, "nc_anti");

    return ret;
  }
}
