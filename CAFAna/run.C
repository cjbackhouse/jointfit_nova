// Installed as part of the container
#include "/nova/bifrost/Bifrost.h"

// Installed as part of the container, but also inlined into NOvALLAdaptor
//#include "/nova/DummyLLH/src/TemplateLLHGetter.h"

#include "CAFAna/Experiment/NOvALLAdaptor.h"

Bifrost& operator>>(Bifrost& bf, OscPars& p)
{
  bf >> p.dm32 >> p.dm21 >> p.sth13 >> p.sth12 >> p.sth23 >> p.dcp;
  return bf;
}

void eval_loop(Bifrost& bf, TemplateLLHGetter* llh)
{
  OscPars p;
  std::vector<double> systs;

  while(true){
    bf >> p >> systs;
    llh->SetOscParameters(p);
    llh->SetParameters(kNOvAdet, systs);
    bf << llh->GetLikelihood();
  }
}

void run()
{
  TemplateLLHGetter* llh = new NOvALLAdaptor;
  llh->init();

  Bifrost& bf = *Bifrost::Inside();

  eval_loop(bf, llh);
}
