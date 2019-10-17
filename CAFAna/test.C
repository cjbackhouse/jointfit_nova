#include "CAFAna/Experiment/NOvALLAdaptor.h"

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

  TemplateLLHGetter* g = new NOvALLAdaptor;

  g->init();
  g->SetOscParameters(p);
  std::cout << "chisq = " << g->GetLikelihood() << std::endl;
  p.sth23 = 0.6;
  std::cout << "chisq = " << g->GetLikelihood() << std::endl;
  g->SetParameters(kNOvAdet, {1});
  std::cout << "chisq = " << g->GetLikelihood() << std::endl;
}
