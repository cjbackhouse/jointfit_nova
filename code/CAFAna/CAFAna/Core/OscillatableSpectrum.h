#pragma once

#include "CAFAna/Core/ReweightableSpectrum.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"

#include <string>

#include "TMD5.h"

class TH2;
class TH2D;

namespace osc{class IOscCalculator;}

namespace ana
{
  /// %Spectrum with true energy information, allowing it to be oscillated
  class OscillatableSpectrum: public ReweightableSpectrum
  {
  public:
    friend class SpectrumLoaderBase;
    friend class SpectrumLoader;
    friend class NullLoader;
    friend class MRCCLoader;
    friend class DUNERunPOTSpectrumLoader;

    OscillatableSpectrum(std::string label,
                         const Binning& bins,
                         SpectrumLoaderBase& loader,
                         const Var& var,
                         const Cut& cut,
                         const SystShifts& shift = kNoShift,
                         const Var& wei = kUnweighted,
                         int potRun = -1);

    OscillatableSpectrum(SpectrumLoaderBase& loader,
                         const HistAxis& axis,
                         const Cut& cut,
                         const SystShifts& shift = kNoShift,
                         const Var& wei = kUnweighted,
                         int potRun = -1);

    OscillatableSpectrum(std::string label, const Binning& bins);
    OscillatableSpectrum(std::string label, double pot, double livetime,
                         const Binning& bins);
    OscillatableSpectrum(TH2* h,
                         const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         double pot, double livetime);

    ~OscillatableSpectrum();

    /// Copy constructor
    OscillatableSpectrum(const OscillatableSpectrum& rhs);
    OscillatableSpectrum(OscillatableSpectrum&& rhs);
    /// Assignment operator
    OscillatableSpectrum& operator=(const OscillatableSpectrum& rhs);
    OscillatableSpectrum& operator=(OscillatableSpectrum&& rhs);

    // Expose these ones directly
    using ReweightableSpectrum::Fill;
    using ReweightableSpectrum::ToTH2;
    using ReweightableSpectrum::Clear;

    /// Rescale bins so that \ref TrueEnergy will return \a target
    using ReweightableSpectrum::ReweightToTrueSpectrum;
    /// Rescale bins so that \ref Unoscillated will return \a target
    using ReweightableSpectrum::ReweightToRecoSpectrum;

    // These under a different name
    Spectrum Unoscillated() const {return UnWeighted();}
    Spectrum TrueEnergy() const {return WeightingVariable();}

    Spectrum Oscillated(osc::IOscCalculator* calc, int from, int to) const;
    Spectrum Oscillated(osc::IOscCalculator* calc, int from, int to, int offLocation) const;

    OscillatableSpectrum& operator+=(const OscillatableSpectrum& rhs);
    OscillatableSpectrum operator+(const OscillatableSpectrum& rhs) const;

    OscillatableSpectrum& operator-=(const OscillatableSpectrum& rhs);
    OscillatableSpectrum operator-(const OscillatableSpectrum& rhs) const;

    void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<OscillatableSpectrum> LoadFrom(TDirectory* dir);

  protected:
    // Derived classes can be trusted take care of their own construction
    OscillatableSpectrum(const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         const Var& rwVar)
      : ReweightableSpectrum(labels, bins, rwVar),
        fOscCache(0, {}, {}, 0, 0), fOscHash(kUninitHash)
    {
    }

    OscillatableSpectrum(const std::string& label,
                         const Binning& bins,
                         const Var& rwVar)
      : ReweightableSpectrum(label, bins, rwVar),
        fOscCache(0, {}, {}, 0, 0), fOscHash(kUninitHash)
    {
    }

    const unsigned char kUninitHashData[16] = {0, 0, 0, 0,
                                               0, 0, 0, 0,
                                               0, 0, 0, 0,
                                               0, 0, 0, 0};
    const TMD5 kUninitHash = TMD5(kUninitHashData);
          
    mutable Spectrum fOscCache;
    mutable TMD5 fOscHash;
    mutable int fOscCacheFrom, fOscCacheTo;
  };
}
