#ifndef OSC_OSCCALCULATORPMNS_H
#define OSC_OSCCALCULATORPMNS_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file   OscCalculatorPMNS.h                                          //
//                                                                      //
// \brief  Adapt the PMNS calculator to standard interface              //
// \author <bckhouse@caltech.edu>					//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cassert>

#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/PMNS.h"

namespace osc
{
  /// Adapt the \ref PMNS calculator to standard interface
  template <typename T>
  class _OscCalculatorPMNS: public _IOscCalculatorAdjustable<T>
  {
    public:
      _OscCalculatorPMNS();
      virtual ~_OscCalculatorPMNS();

      virtual _IOscCalculatorAdjustable<T>* Copy() const override;

      virtual T P(int flavBefore, int flavAfter, double E) override;

      void SetL     (double   L     ) override {fPropDirty = true; this->fL      = L;}
      void SetRho   (double   rho   ) override {fPropDirty = true; this->fRho    = rho;}
      void SetDmsq21(const T& dmsq21) override {fDmDirty   = true; this->fDmsq21 = dmsq21;}
      void SetDmsq32(const T& dmsq32) override {fDmDirty   = true; this->fDmsq32 = dmsq32;}
      void SetTh12  (const T& th12  ) override {fMixDirty  = true; this->fTh12   = th12;}
      void SetTh13  (const T& th13  ) override {fMixDirty  = true; this->fTh13   = th13;}
      void SetTh23  (const T& th23  ) override {fMixDirty  = true; this->fTh23   = th23;}
      void SetdCP   (const T& dCP   ) override {fMixDirty  = true; this->fdCP    = dCP;}

      TMD5* GetParamsHash() const override
      {
        return _IOscCalculatorAdjustable<T>::GetParamsHashDefault("PMNS");
      }

    protected:
      _PMNS<T> fPMNS;

      bool fMixDirty;
      bool fDmDirty;
      bool fPropDirty;
      double fPrevE;
      int fPrevAnti;
  };

  typedef _OscCalculatorPMNS<double> OscCalculatorPMNS;

} // namespace

#endif
