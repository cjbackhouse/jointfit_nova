#pragma once

#include <list>
#include <string>

namespace ana
{
  class ISyst
  {
  public:
    ISyst(const std::string& shortName,
          const std::string& latexName);
    ISyst(const ISyst &) = delete;   // no copying.
    ISyst(ISyst && rhs) = delete;    // no moving either.
    virtual ~ISyst();

    void operator=(const ISyst &) = delete;  // still no copying.
    void operator=(ISyst &&)      = delete;  // etc.

    /// The name printed out to the screen
    virtual std::string ShortName() const final {return fShortName;}

    /// The name used on plots (ROOT's TLatex syntax)
    virtual std::string LatexName() const final {return fLatexName;}

    virtual double Penalty(double x) const;
    virtual double PenaltyDerivative(double x) const;

    /// PredictionInterp normally interpolates between spectra made at
    /// +/-1,2,3sigma. For some systematics that's overkill. Override this
    /// function to specify different behaviour for this systematic.
    virtual int PredInterpMaxNSigma() const
    {
      return 3;
    }

  private:
    std::string fShortName;
    std::string fLatexName;
    double fMin;
    double fMax;
  };

} // namespace
