#include "OscLib/func/IOscCalculator.h"

#include "Utilities/func/StanUtils.h"

namespace osc
{
  //---------------------------------------------------------------------------
  template<class T> _IOscCalculatorAdjustable<T>::~_IOscCalculatorAdjustable()
  {
  }

  //---------------------------------------------------------------------------
  template <typename T>
  TMD5* _IOscCalculatorAdjustable<T>::GetParamsHashDefault(const std::string& txt) const
  {
    TMD5* ret = new TMD5;
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 8;
    T buf[kNumParams] = {fRho, fL, fDmsq21, fDmsq32,
                         fTh12, fTh13, fTh23, fdCP};
    ret->Update((unsigned char*)buf, sizeof(T)*kNumParams);
    ret->Final();
    return ret;
  }

  template <typename T, typename U>
  void CopyParams(const osc::_IOscCalculatorAdjustable<T> * inCalc,
                  osc::_IOscCalculatorAdjustable<U> * outCalc)
  {
    assert (inCalc);
    assert (outCalc);

    outCalc->SetL(inCalc->GetL());
    outCalc->SetRho(inCalc->GetRho());

    outCalc->SetdCP(util::GetValAs<U>(inCalc->GetdCP()));
    outCalc->SetDmsq21(util::GetValAs<U>(inCalc->GetDmsq21()));
    outCalc->SetDmsq32(util::GetValAs<U>(inCalc->GetDmsq32()));
    outCalc->SetTh12(util::GetValAs<U>(inCalc->GetTh12()));
    outCalc->SetTh13(util::GetValAs<U>(inCalc->GetTh13()));
    outCalc->SetTh23(util::GetValAs<U>(inCalc->GetTh23()));
  }

  //---------------------------------------------------------------------------
  template class _IOscCalculatorAdjustable<double>;

#ifndef DARWINBUILD
  template class _IOscCalculatorAdjustable<stan::math::var>;
  template void CopyParams(const osc::_IOscCalculatorAdjustable<double> * inCalc,
                           osc::_IOscCalculatorAdjustable<stan::math::var> * outCalc);
  template void CopyParams(const osc::_IOscCalculatorAdjustable<stan::math::var> * inCalc,
                           osc::_IOscCalculatorAdjustable<double> * outCalc);
  template void CopyParams(const osc::_IOscCalculatorAdjustable<stan::math::var> * inCalc,
                           osc::_IOscCalculatorAdjustable<stan::math::var> * outCalc);
#endif

  template void CopyParams(const osc::_IOscCalculatorAdjustable<double> * inCalc,
                           osc::_IOscCalculatorAdjustable<double> * outCalc);
}
