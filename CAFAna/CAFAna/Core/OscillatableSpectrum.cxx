#include "CAFAna/Core/OscillatableSpectrum.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/OscCurve.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/func/IOscCalculator.h"

#include "TDirectory.h"
#include "TH2.h"
#include "TMD5.h"
#include "TObjString.h"

#include <cassert>
#include <memory>

namespace ana
{
  //----------------------------------------------------------------------
  OscillatableSpectrum::OscillatableSpectrum(const std::string& label,
                                             const Binning& bins)
    : ReweightableSpectrum(label, bins),
      fCachedOsc(0, {}, {}, 0, 0),
      fCachedHash(0)
  {
    fTrueLabel = "True Energy (GeV)";

    DontAddDirectory guard;

    fPOT = 0;
    fLivetime = 0;

    fHist = HistCache::NewTH2D("", bins);
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum::OscillatableSpectrum(const std::string& label, double pot, double livetime,
                                             const Binning& bins)
    : ReweightableSpectrum(label, bins),
      fCachedOsc(0, {}, {}, 0, 0),
      fCachedHash(0)
  {
    fTrueLabel = "True Energy (GeV)";

    DontAddDirectory guard;

    fPOT = pot;
    fLivetime = livetime;

    fHist = HistCache::NewTH2D("", bins);
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum::OscillatableSpectrum(std::unique_ptr<TH2D> h,
                                             const std::vector<std::string>& labels,
                                             const std::vector<Binning>& bins,
                                             double pot, double livetime)
    : ReweightableSpectrum(std::move(h), labels, bins, pot, livetime),
      fCachedOsc(0, {}, {}, 0, 0),
      fCachedHash(0)
  {
    fTrueLabel = "True Energy (GeV)";
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum::~OscillatableSpectrum()
  {
    // Nulls fHist out, so it's safe that ~ReweightableSpectrum tries too
    HistCache::Delete(fHist, Bins1DX().ID(), kTrueEnergyBins.ID());

    delete fCachedHash;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum::OscillatableSpectrum(const OscillatableSpectrum& rhs)
    : ReweightableSpectrum(rhs.fLabels, rhs.fBins),
      fCachedOsc(0, {}, {}, 0, 0),
      fCachedHash(0)
  {
    DontAddDirectory guard;

    fHist = HistCache::Copy(rhs.fHist, rhs.Bins1DX(), kTrueEnergyBins);

    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;

    if(rhs.fCachedHash){
      fCachedOsc = rhs.fCachedOsc;
      fCachedHash = new TMD5(*rhs.fCachedHash);
    }
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum::OscillatableSpectrum(OscillatableSpectrum&& rhs)
    : ReweightableSpectrum(rhs.fLabels, rhs.fBins),
      fCachedOsc(0, {}, {}, 0, 0),
      fCachedHash(0)
  {
    DontAddDirectory guard;

    fHist = rhs.fHist;
    rhs.fHist = 0;

    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;

    if(rhs.fCachedHash){
      fCachedOsc = std::move(rhs.fCachedOsc);
      fCachedHash = rhs.fCachedHash;
      rhs.fCachedHash = 0;
    }
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum& OscillatableSpectrum::operator=(const OscillatableSpectrum& rhs)
  {
    if(this == &rhs) return *this;

    DontAddDirectory guard;

    if(fHist) HistCache::Delete(fHist, Bins1DX().ID());
    fHist = HistCache::Copy(rhs.fHist, rhs.Bins1DX(), kTrueEnergyBins);
    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;
    fLabels = rhs.fLabels;
    fBins = rhs.fBins;

    if(rhs.fCachedHash){
      fCachedOsc = rhs.fCachedOsc;
      delete fCachedHash;
      fCachedHash = new TMD5(*rhs.fCachedHash);
    }

    return *this;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum& OscillatableSpectrum::operator=(OscillatableSpectrum&& rhs)
  {
    if(this == &rhs) return *this;

    DontAddDirectory guard;

    if(fHist) HistCache::Delete(fHist, Bins1DX().ID());
    fHist = rhs.fHist;
    rhs.fHist = 0;
    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;
    fLabels = rhs.fLabels;
    fBins = rhs.fBins;

    if(rhs.fCachedHash){
      fCachedOsc = rhs.fCachedOsc;
      delete fCachedHash;
      fCachedHash = rhs.fCachedHash;
      rhs.fCachedHash = 0;
    }

    return *this;
  }

  //----------------------------------------------------------------------
  Spectrum OscillatableSpectrum::Oscillated(osc::IOscCalculator* calc,
                                            int from, int to) const
  {
    TMD5* hash = calc->GetParamsHash();
    if(hash && fCachedHash && *hash == *fCachedHash){
      delete hash;
      return fCachedOsc;
    }

    const OscCurve curve(calc, from, to);
    TH1D* Ps = curve.ToTH1();

    const Spectrum ret = WeightedBy(Ps);
    if(hash){
      fCachedOsc = ret;
      delete fCachedHash;
      fCachedHash = hash;
    }
    HistCache::Delete(Ps, kTrueEnergyBins.ID());
    return ret;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum& OscillatableSpectrum::operator+=(const OscillatableSpectrum& rhs)
  {
    if(rhs.fPOT){
      fHist->Add(rhs.fHist, fPOT/rhs.fPOT);
    }
    else{
      // How can it have events but no POT?
      assert(rhs.fHist->Integral() == 0);
    }

    delete fCachedHash;
    fCachedHash = 0; // Invalidate

    return *this;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum OscillatableSpectrum::operator+(const OscillatableSpectrum& rhs) const
  {
    OscillatableSpectrum ret = *this;
    ret += rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum& OscillatableSpectrum::operator-=(const OscillatableSpectrum& rhs)
  {
    if(rhs.fPOT){
      fHist->Add(rhs.fHist, -fPOT/rhs.fPOT);
    }
    else{
      // How can it have events but no POT?
      assert(rhs.fHist->Integral() == 0);
    }

    delete fCachedHash;
    fCachedHash = 0; // Invalidate

    return *this;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum OscillatableSpectrum::operator-(const OscillatableSpectrum& rhs) const
  {
    OscillatableSpectrum ret = *this;
    ret -= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  void OscillatableSpectrum::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TObjString("OscillatableSpectrum").Write("type");

    fHist->Write("hist");
    TH1D hPot("", "", 1, 0, 1);
    hPot.Fill(.5, fPOT);
    hPot.Write("pot");
    TH1D hLivetime("", "", 1, 0, 1);
    hLivetime.Fill(.5, fLivetime);
    hLivetime.Write("livetime");

    for(unsigned int i = 0; i < fBins.size(); ++i){
      TObjString(fLabels[i].c_str()).Write(TString::Format("label%d", i).Data());
      fBins[i].SaveTo(dir->mkdir(TString::Format("bins%d", i)));
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<OscillatableSpectrum> OscillatableSpectrum::LoadFrom(TDirectory* dir)
  {
    DontAddDirectory guard;

    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "OscillatableSpectrum");
    delete tag;

    TH2D* spect = (TH2D*)dir->Get("hist");
    assert(spect);
    TH1* hPot = (TH1*)dir->Get("pot");
    assert(hPot);
    TH1* hLivetime = (TH1*)dir->Get("livetime");
    assert(hLivetime);

    std::vector<std::string> labels;
    std::vector<Binning> bins;

    for(int i = 0; ; ++i){
      TDirectory* subdir = dir->GetDirectory(TString::Format("bins%d", i));
      if(!subdir) break;
      bins.push_back(*Binning::LoadFrom(subdir));
      TObjString* label = (TObjString*)dir->Get(TString::Format("label%d", i));
      labels.push_back(label ? label->GetString().Data() : "");
      delete subdir;
      delete label;
    }

    if(bins.empty() && labels.empty()){
      // Must be an old file. Make an attempt at backwards compatibility.
      bins.push_back(Binning::FromTAxis(spect->GetXaxis()));
      labels.push_back(spect->GetXaxis()->GetTitle());
    }

    auto ret = std::make_unique<OscillatableSpectrum>(std::unique_ptr<TH2D>(spect),
                                                      labels, bins,
                                                      hPot->GetBinContent(1),
                                                      hLivetime->GetBinContent(1));

    delete hPot;
    delete hLivetime;
    return ret;
  }
}
