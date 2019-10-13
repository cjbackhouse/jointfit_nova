#include "CAFAna/Core/LoadFromFile.h"

#include "OscLib/func/OscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"

#include "TObjString.h"
#include "TH1.h"
#include "TVectorD.h"
#include <vector>

namespace ana
{
  //----------------------------------------------------------------------
  template<> std::unique_ptr<osc::IOscCalculator>
  LoadFrom<osc::IOscCalculator>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    const TString tag = ptag->GetString();

    if(tag == "NoOscillations") return std::unique_ptr<osc::IOscCalculator>(new osc::NoOscillations);

    osc::IOscCalculatorAdjustable* ret = 0;

    if(tag == "OscCalculator") ret = new osc::OscCalculator;
    if(tag == "OscCalculatorPMNSOpt") ret = new osc::OscCalculatorPMNSOpt;

    if(!ret){
      std::cout << "LoadFrom not implemented for " << tag << std::endl;
      abort();
    }

    TVectorD* params = (TVectorD*)dir->Get("params");
    assert(params);
    assert(params->GetNrows() == 8);

    ret->SetL     ((*params)[0]);
    ret->SetRho   ((*params)[1]);
    ret->SetDmsq21((*params)[2]);
    ret->SetDmsq32((*params)[3]);
    ret->SetTh12  ((*params)[4]);
    ret->SetTh13  ((*params)[5]);
    ret->SetTh23  ((*params)[6]);
    ret->SetdCP   ((*params)[7]);

    return std::unique_ptr<osc::IOscCalculatorAdjustable>(ret);
  }

  //----------------------------------------------------------------------
  template<> void SaveTo(const osc::IOscCalculator& x, TDirectory* dir)
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    if(dynamic_cast<const osc::NoOscillations*>(&x)){
      TObjString("NoOscillations").Write("type");
      tmp->cd();
      return;
    }

    const osc::IOscCalculatorAdjustable* y = dynamic_cast<const osc::IOscCalculatorAdjustable*>(&x);
    if(!y){
      std::cout << "Unknown calculator in SaveTo " << typeid(x).name() << std::endl;
      abort();
    }
    
    /* */if(dynamic_cast<const osc::OscCalculator*>(&x)) TObjString("OscCalculatorPMNS").Write("type");
    else if(dynamic_cast<const osc::OscCalculatorPMNSOpt*>(&x)) TObjString("OscCalculatorPMNSOpt").Write("type");
    else{
      std::cout << "Unimplemented calculator in SaveTo " << typeid(x).name() << std::endl;
      abort();
    }
    
    TVectorD params(8);

    params[0] = y->GetL();
    params[1] = y->GetRho();
    params[2] = y->GetDmsq21();
    params[3] = y->GetDmsq32();
    params[4] = y->GetTh12();
    params[5] = y->GetTh13();
    params[6] = y->GetTh23();
    params[7] = y->GetdCP();

    params.Write("params");

    tmp->cd();
  }
}
