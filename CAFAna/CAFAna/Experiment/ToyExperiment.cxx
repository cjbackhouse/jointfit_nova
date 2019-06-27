#include "CAFAna/Experiment/ToyExperiment.h"

#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"

#include "TH1.h"
#include "TRandom3.h"

namespace ana
{
  osc::IOscCalculatorAdjustable* DefaultOscCalc()
  {
    osc::IOscCalculatorAdjustable* calc = new osc::OscCalculatorPMNSOpt;

    calc->SetL(810);
    calc->SetRho(2.84);
    calc->SetDmsq21(7.53e-5);
    calc->SetTh12(asin(sqrt(.307)));
    calc->SetTh23(M_PI/4);
    calc->SetDmsq32(2.45e-3);
    calc->SetTh13(asin(sqrt(2.10e-2)));
    calc->SetdCP(0);

    return calc;
  }

  ToyExperiment::ToyExperiment()
    : fMC("", 1, 1, Binning::Simple(20, 0, 5)),
      fData("", 1, 1, Binning::Simple(20, 0, 5))
  {
    for(int i = 0; i < 200; ++i){
      const double Etrue = gRandom->Gaus(2, 1);
      const double Ereco = Etrue*gRandom->Gaus(1, .1);
      fMC.Fill(Etrue, Ereco);
    }
    fData = fMC.Oscillated(DefaultOscCalc(), 14, 14).FakeData(1);
  }

  double ToyExperiment::ChiSq(osc::IOscCalculatorAdjustable* osc,
                              const SystShifts& syst) const
  {
    const Spectrum m = fMC.Oscillated(osc, 14, 14);
    TH1D* hm = m.ToTH1(1);
    TH1D* hd = fData.ToTH1(1);

    const double chisq = LogLikelihood(hm, hd);

    HistCache::Delete(hm);
    HistCache::Delete(hd);

    return chisq;
  }
}
