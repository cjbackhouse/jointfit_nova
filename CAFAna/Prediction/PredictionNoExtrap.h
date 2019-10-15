#pragma once

#include "TDirectory.h"
#include "TObjString.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Extrap/TrivialExtrap.h"
#include "CAFAna/Prediction/PredictionExtrap.h"

#include "CAFAna/Extrap/TrivialExtrap.h"

namespace ana
{
  class Loaders;

  /// Prediction that just uses FD MC, with no extrapolation
  class PredictionNoExtrap: public PredictionExtrap
  {
    public:
      PredictionNoExtrap(std::unique_ptr<TrivialExtrap>&& extrap);

      virtual ~PredictionNoExtrap();

      static std::unique_ptr<PredictionNoExtrap> LoadFrom(TDirectory* dir);
      virtual void SaveTo(TDirectory* dir) const override;

  };

}
