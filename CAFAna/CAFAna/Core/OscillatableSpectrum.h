#pragma once

#include "CAFAna/Core/ReweightableSpectrum.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Spectrum.h"

#include <string>

class TH2;
class TH2D;
class TMD5;

namespace osc{class IOscCalculator;}

namespace ana
{
  /// %Spectrum with true energy information, allowing it to be oscillated
  class OscillatableSpectrum: public ReweightableSpectrum
  {
  public:
    OscillatableSpectrum(const std::string& label, const Binning& bins);
    OscillatableSpectrum(const std::string& label, double pot, double livetime,
                         const Binning& bins);

    OscillatableSpectrum(std::unique_ptr<TH2D> h,
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

    OscillatableSpectrum& operator+=(const OscillatableSpectrum& rhs);
    OscillatableSpectrum operator+(const OscillatableSpectrum& rhs) const;

    OscillatableSpectrum& operator-=(const OscillatableSpectrum& rhs);
    OscillatableSpectrum operator-(const OscillatableSpectrum& rhs) const;

    void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<OscillatableSpectrum> LoadFrom(TDirectory* dir);

  protected:
    // Derived classes can be trusted take care of their own construction
    OscillatableSpectrum(const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins)
      : ReweightableSpectrum(labels, bins),
        fCachedOsc(0, {}, {}, 0, 0),
        fCachedHash(0)
    {
    }

    mutable Spectrum fCachedOsc;
    mutable TMD5* fCachedHash;
  };
}
