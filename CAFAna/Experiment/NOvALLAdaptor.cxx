#include "CAFAna/Experiment/NOvALLAdaptor.h"

#include "CAFAna/Experiment/ToyExperiment.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"

void NOvALLAdaptor::init()
{
  fExp = new ana::ToyExperiment;
  fOsc = new osc::OscCalculatorPMNSOpt;
  fSyst = new ana::SystShifts;
}

void NOvALLAdaptor::SetParameters(CovTypes iCov, std::vector<double> vals)
{
  if(vals.size() > 0) fSyst->SetShift(&ana::dummySyst, vals[0]);
}

void NOvALLAdaptor::SetOscParameters(OscPars p)
{
  fOsc->SetL(810);
  fOsc->SetRho(2.84);
  fOsc->SetDmsq21(p.dm21);
  fOsc->SetTh12(asin(p.sth12));
  fOsc->SetTh23(asin(p.sth23));
  fOsc->SetDmsq32(p.dm32);
  fOsc->SetTh13(asin(p.sth13));
  fOsc->SetdCP(p.dcp);
}

double NOvALLAdaptor::GetLikelihood()
{
  return fExp->ChiSq(fOsc, *fSyst);
}
