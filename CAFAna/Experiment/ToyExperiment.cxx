#include "CAFAna/Experiment/ToyExperiment.h"

namespace ana
{
  ToyExperiment::ToyExperiment()
    : SingleSampleExperiment(LoadFromFile<IPrediction>("state_v1.root",
                                                       "pred").release(),
                             *LoadFromFile<Spectrum>("state_v1.root",
                                                     "mock_data"))
  {
  }

  DummyNormSyst::DummyNormSyst() : ISyst("dummynorm", "Dummy Norm") {}

  const DummyNormSyst dummySyst;
}
