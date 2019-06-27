#include "CAFAna/Experiment/IExperiment.h"

#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TObjString.h"

#include <cassert>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IExperiment> LoadFrom<IExperiment>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    std::cerr << "Unknown Experiment type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  void IExperiment::SaveTo(TDirectory* dir) const
  {
    assert(0 && "Not implemented");
  }
}
