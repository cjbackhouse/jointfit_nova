// Copy-paste from DummyLLH

#include <vector>

enum CovTypes {//We should expand this to have all the T2K and NOvA classes
  kOsc=0,
  kSKdet=1,
  kBANFF=2,
  kXsec=3,
  kND280det=4,
  kNOvAdet=5,
};

class OscPars{
public:
  double dm32,dm21,sth13,sth12,sth23,dcp;
};

class TemplateLLHGetter {
public:
  TemplateLLHGetter(){};
  ~TemplateLLHGetter(){};

  virtual void init()=0;//Some sort of call once init fuction
  virtual void SetParameters(CovTypes iCov, std::vector<double> vals)=0;//Set systematic parameters
  virtual void SetOscParameters(OscPars oscpars)=0;//Set oscillation parameters
  virtual double GetLikelihood()=0;//Get likelihood

};


#include "CAFAna/Experiment/ToyExperiment.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"

class NOvAGetter: public TemplateLLHGetter
{
public:
  void init()
  {
    fExp = new ana::ToyExperiment;
    fOsc = new osc::OscCalculatorPMNSOpt;
  }
  void SetParameters(CovTypes iCov, std::vector<double> vals){}
  void SetOscParameters(OscPars p)
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
  double GetLikelihood()
  {
    return fExp->ChiSq(fOsc);
  }

protected:
  ana::IExperiment* fExp;
  osc::IOscCalculatorAdjustable* fOsc;
};


#include <iostream>

void test()
{
  class OscPars p;
  p.dm32 = 2.5e-3;
  p.dm21 = 7.5e-5;
  p.sth13 = 0.1;
  p.sth12 = 0.3;
  p.sth23 = 0.5;
  p.dcp = 0;

  NOvAGetter g;
  g.init();
  g.SetOscParameters(p);
  std::cout << "chisq = " << g.GetLikelihood() << std::endl;
}
