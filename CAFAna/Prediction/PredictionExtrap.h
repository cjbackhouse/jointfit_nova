#pragma once

#include "TH1D.h"

#include "CAFAna/Prediction/IPrediction.h"

namespace ana
{
  class IExtrap;

  /// Take the output of an extrapolation and oscillate it as required
  class PredictionExtrap: public IPrediction
  {
  public:
    /// Takes ownership of \a extrap
    PredictionExtrap(std::unique_ptr<IExtrap>&& extrap);
    PredictionExtrap() = delete;

    // un-hide inherited method stubs so we don't get warnings from the compiler
    using IPrediction::Predict;
    using IPrediction::PredictComponent;
    using IPrediction::PredictSyst;

    Spectrum     Predict(osc::IOscCalculator* calc) const override;

    Spectrum PredictComponent(osc::IOscCalculator* calc,
                              Flavors::Flavors_t flav,
                              Current::Current_t curr,
                              Sign::Sign_t sign) const override;

    OscillatableSpectrum ComponentCC(int from, int to) const override;
    //nc
    Spectrum ComponentNCTotal()  const override;
    Spectrum ComponentNC()       const override;
    Spectrum ComponentNCAnti()   const override;
    //end nc
    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionExtrap> LoadFrom(TDirectory* dir);

    const IExtrap* GetExtrap() const {return fExtrap.get();}

  protected:
    std::unique_ptr<IExtrap> fExtrap;
  };

}
