#pragma once

#include "CAFAna/Experiment/SingleSampleExperiment.h"

namespace ana
{
  class ToyExperiment: public SingleSampleExperiment
  {
  public:
    ToyExperiment();
  };

  struct DummyNormSyst: public ISyst
  {
    DummyNormSyst();
  };

  extern const DummyNormSyst dummySyst;
}
