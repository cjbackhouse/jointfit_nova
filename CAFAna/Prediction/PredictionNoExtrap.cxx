#include "CAFAna/Prediction/PredictionNoExtrap.h"

namespace ana
{
  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(std::unique_ptr<TrivialExtrap>&& extrap)
    : PredictionExtrap(std::move(extrap))
  {
  }

  //----------------------------------------------------------------------
  void PredictionNoExtrap::SaveTo(TDirectory *dir) const
  {
    TDirectory *tmp = gDirectory;

    dir->cd();

    TObjString("PredictionNoExtrap").Write("type");

    this->fExtrap->SaveTo(dir->mkdir("extrap"));

    tmp->cd();
  }


  //----------------------------------------------------------------------
  std::unique_ptr<PredictionNoExtrap> PredictionNoExtrap::LoadFrom(TDirectory *dir)
  {
    assert(dir->GetDirectory("extrap"));
    return std::make_unique<PredictionNoExtrap>(ana::LoadFrom<TrivialExtrap>(dir->GetDirectory("extrap")));
  }


  //----------------------------------------------------------------------
  PredictionNoExtrap::~PredictionNoExtrap()
  {
    // We created this in the constructor so it's our responsibility
    // delete fExtrap;
  }
}
