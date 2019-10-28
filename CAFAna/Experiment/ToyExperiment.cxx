#include "CAFAna/Experiment/ToyExperiment.h"

#include "CAFAna/Core/Utilities.h"

namespace ana
{
  ToyExperiment::ToyExperiment()
    : SingleSampleExperiment(LoadFromFile<IPrediction>(FindCAFAnaDir()+"/state_v1.root",
                                                       "pred").release(),
                             *LoadFromFile<Spectrum>(FindCAFAnaDir()+"/state_v1.root",
                                                     "mock_data"))
  {
  }

  DummyNormSyst::DummyNormSyst() : ISyst("dummynorm", "Dummy Norm") {}

  const DummyNormSyst dummySyst;
}
