#pragma once

#include "CAFAna/Experiment/IExperiment.h"

#include "CAFAna/Core/OscillatableSpectrum.h"

namespace ana
{
  class ToyExperiment: public IExperiment
  {
  public:
    ToyExperiment();
    double ChiSq(osc::IOscCalculatorAdjustable* osc,
                 const SystShifts& syst = SystShifts::Nominal()) const override;
  protected:
    OscillatableSpectrum fMC;
    Spectrum fData;
  };
}
