#include "CAFAna/Core/ISyst.h"

#include "CAFAna/Core/SystRegistry.h"

#include "Utilities/func/MathUtil.h"

namespace ana
{
  //----------------------------------------------------------------------
  ISyst::ISyst(const std::string& shortName,
               const std::string& latexName)
    : fShortName(shortName), fLatexName(latexName)
  {
    SystRegistry::Register(this);
  }

  //----------------------------------------------------------------------
  ISyst::~ISyst()
  {
    // Normally ISysts should last for the life of the process, but in case one
    // is deleted it's best not to leave a dangling pointer in SystRegistry.
    SystRegistry::UnRegister(this);
  }

  //----------------------------------------------------------------------
  double ISyst::Penalty(double x) const
  {
    // Regular quadratic penalty term
    return x*x;
  }

  //----------------------------------------------------------------------
  double ISyst::PenaltyDerivative(double x) const
  {
    // Regular quadratic penalty term
    return 2*x;
  }
}
