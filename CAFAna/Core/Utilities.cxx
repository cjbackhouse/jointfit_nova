#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Ratio.h"

#include "Utilities/func/MathUtil.h"

#include "TArrayD.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TObjString.h"
#include "TString.h"
#include "TVector3.h"
#include "TVectorD.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include "sys/stat.h"
#include "wordexp.h"

namespace ana
{
  double LLPerBinFracSystErr::fgErr = -1;

  //----------------------------------------------------------------------
  std::string UniqueName()
  {
    static int N = 0;
    return TString::Format("cafanauniq%d", N++).Data();
  }

  //----------------------------------------------------------------------
  DontAddDirectory::DontAddDirectory()
  {
    fBackup = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
  }

  //----------------------------------------------------------------------
  DontAddDirectory::~DontAddDirectory()
  {
    TH1::AddDirectory(fBackup);
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::FloatingExceptionOnNaN(bool enable)
  {
    // Don't want any pending FPEs to trigger when we flip exceptions
    // on. Whoever had them off previously had a reason.
    feclearexcept(FE_INVALID);

    fegetexceptflag(&fBackup, FE_INVALID);

#ifndef DARWINBUILD
    if(enable)
      feenableexcept(FE_INVALID);
    else
      fedisableexcept(FE_INVALID);
#else
    std::cerr << "WARNING: CAFAna/Core/Utilities.cxx built on OS X, no feenableexcept available" << std::endl;
#endif
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::~FloatingExceptionOnNaN()
  {
    fesetexceptflag(&fBackup, FE_INVALID);
  }

  //----------------------------------------------------------------------
  double LogLikelihood(double e, double o)
  {
    // http://www.wolframalpha.com/input/?i=d%2Fds+m*(1%2Bs)+-d+%2B+d*ln(d%2F(m*(1%2Bs)))%2Bs%5E2%2FS%5E2%3D0
    // http://www.wolframalpha.com/input/?i=solve+-d%2F(s%2B1)%2Bm%2B2*s%2FS%5E2%3D0+for+s
    const double S = LLPerBinFracSystErr::GetError();
    if(S > 0){
      const double S2 = util::sqr(S);
      const double s = .25*(sqrt(8*o*S2+util::sqr(e*S2-2))-e*S2-2);
      e *= 1+s;
    }

    // With this value, negative expected events and one observed
    // event gives a chisq from this one bin of 182.
    const double minexp = 1e-40; // Don't let expectation go lower than this

    assert(o >= 0);
    if(e < minexp){
      if(!o) return 0;
      e = minexp;
    }

    if(o){
      /*
      const double x = (o-e)/e;
      if(fabs(x) < 1e-3){
        // For o/e very close to 1, the power expansion is much more stable
        // than the logarithm. With this many orders we're good to 1 part in
        // 10^21
        const double x2 = x*x;
        const double x3 = x2*x;
        const double x4 = x3*x;
        const double x5 = x4*x;
        const double x6 = x5*x;
        const double x7 = x6*x;
        return 2*(e-o) + 2*o*(x - x2/2 + x3/3 - x4/4 + x5/5 - x6/6 + x7/7);
      }
      */

      // This strange form is for numerical stability when e~o
      return 2*o*((e-o)/o + log1p((o-e)/e));
    }
    else{
      return 2*(e-o);
    }
  }

  //----------------------------------------------------------------------
  double LogLikelihood(const TH1* eh, const TH1* oh, bool useOverflow)
  {
    assert(eh->GetNbinsX() == oh->GetNbinsX());

    double chi = 0;

    int bufferBins = useOverflow? 2 : 1;

    for(int i = 0; i < eh->GetNbinsX()+bufferBins; ++i){
      const double e = eh->GetBinContent(i);
      const double o = oh->GetBinContent(i);

      chi += LogLikelihood(e, o);
    }

    return chi;
  }

  //----------------------------------------------------------------------
  // dLL/de
  double LogLikelihoodDerivative(double e, double o)
  {
    if(e == 0) return 2;
    return 2-2*o/e;
  }

  //----------------------------------------------------------------------
  // dLL/dx
  double LogLikelihoodDerivative(const TH1D* eh, const TH1D* oh,
                                 const std::vector<double>& dedx)
  {
    const double* ea = eh->GetArray();
    const double* oa = oh->GetArray();

    double ret = 0;
    for(unsigned int i = 0; i < dedx.size(); ++i){
      ret += LogLikelihoodDerivative(ea[i], oa[i]) * dedx[i];
    }
    return ret;
  }

  //----------------------------------------------------------------------
  TH2F* ExpandedHistogram(const std::string& title,
                          int nbinsx, double xmin, double xmax,
                          int nbinsy, double ymin, double ymax)
  {
    DontAddDirectory guard;

    // How wide the bins will be once we're done
    const double xwidth = (xmax-xmin)/(nbinsx-1);
    const double ywidth = (ymax-ymin)/(nbinsy-1);

    // Move the bin edges so that the limits occur at the centres
    xmin -= xwidth/2; ymin -= ywidth/2;
    xmax += xwidth/2; ymax += ywidth/2;

    return new TH2F(UniqueName().c_str(), title.c_str(),
                    nbinsx, xmin, xmax,
                    nbinsy, ymin, ymax);
  }

  // Helper functions for MakeTHND().
  namespace{
    // Eventually the bin parameters will all be unpacked and we just pass them
    // on to the regular constructor.
    template<class T, class... A> T* MakeHist(A... a)
    {
      DontAddDirectory guard;
      return new T(a...);
    }

    // This function consumes bins from the start of the argument list and
    // pushes their translations onto the list of arguments at the end.
    template<class T, class... A> T* MakeHist(const Binning& firstBin,
                                              A... args)
    {
      if(firstBin.IsSimple())
        return MakeHist<T>(args...,
                           firstBin.NBins(), firstBin.Min(), firstBin.Max());
      else
        return MakeHist<T>(args...,
                           firstBin.NBins(), &firstBin.Edges().front());
    }
  }

  // Concrete instantiations. MakeHist() requires us to put the bin arguments
  // first...
  //----------------------------------------------------------------------
  TH1D* MakeTH1D(const char* name, const char* title, const Binning& bins)
  {
    return MakeHist<TH1D>(bins, name, title);
  }

  //----------------------------------------------------------------------
  TH2D* MakeTH2D(const char* name, const char* title,
                 const Binning& binsx,
                 const Binning& binsy)
  {
    return MakeHist<TH2D>(binsx, binsy, name, title);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));
    return ToTH2Helper(std::move(h1), binsx, binsy, bintype);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Ratio& r,
             const Binning& binsx, const Binning& binsy)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());
    return ToTH2Helper(std::move(h1), binsx, binsy);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx, const Binning& binsy,
		   ana::EBinType bintype)
  {
    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == binsx.NBins()*binsy.NBins());

    TH2* ret = MakeTH2D("", UniqueName().c_str(), binsx, binsy);

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int ix = i / binsy.NBins();
      const int iy = i % binsy.NBins();

      ret->SetBinContent(ix+1, iy+1, val);
      ret->SetBinError  (ix+1, iy+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, const Binning& binsz,
	     ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz, bintype);
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Ratio& r,
             const Binning& binsx, const Binning& binsy, const Binning& binsz)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz);
  }

  //----------------------------------------------------------------------
  TH3* ToTH3Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx,
		   const Binning& binsy,
		   const Binning& binsz,
		   ana::EBinType bintype)
  {

    const int nx = binsx.NBins();
    const int ny = binsy.NBins();
    const int nz = binsz.NBins();

    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == nx*ny*nz);

    TH3* ret;

    // If all three axes are simple, we can call a simpler constructor
    if(binsx.IsSimple() && binsy.IsSimple() && binsz.IsSimple()){
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, binsx.Min(), binsx.Max(),
                     ny, binsy.Min(), binsy.Max(),
                     nz, binsz.Min(), binsz.Max());

      if(!binsx.IsSimple() || !binsy.IsSimple() || !binsz.IsSimple()){
        // TH3 doesn't have the constructors for mixed simple and custom
        std::cerr << "ToTH3: one or more axes is custom, but not all three. Applying Simple binning to all three axes" << std::endl;
      }
    }
    else{
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, &binsx.Edges().front(),
                     ny, &binsy.Edges().front(),
                     nz, &binsz.Edges().front());
    }

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int nynz = ny*nz;
      const int nmodnynz = i%nynz;
      const int ix = i/nynz;
      const int iy = nmodnynz/nz;
      const int iz = i%nz;

      ret->SetBinContent(ix+1, iy+1, iz+1, val);
      ret->SetBinError  (ix+1, iy+1, iz+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;

  }

  //----------------------------------------------------------------------
  std::vector<std::string> Wildcard(const std::string& wildcardString)
  {
    // Expand environment variables and wildcards like the shell would
    wordexp_t p;
    const int status = wordexp(wildcardString.c_str(), &p, WRDE_SHOWERR);

    if(status != 0){
      std::cerr << "Wildcard string '" << wildcardString
                << "' returned error " << status << " from wordexp()."
                << std::endl;
      return {};
    }

    std::vector<std::string> fileList;

    for(unsigned int i = 0; i < p.we_wordc; ++i){
      // Check the file exists before adding it
      struct stat sb;
      if(stat(p.we_wordv[i], &sb) == 0)
        fileList.push_back(p.we_wordv[i]);
    }

    wordfree(&p);

    return fileList;
  }

  //----------------------------------------------------------------------
  std::string FindCAFAnaDir()
  {
    const char* cafana = getenv("JOINTFIT_DIR");
    struct stat junk;
    if(!cafana || stat(cafana, &junk) != 0){
      std::cout << "Couldn't find CAFAna dir (using $JOINTFIT_DIR)" << std::endl;
      abort();
    }

    return cafana+std::string("/CAFAna/");
  }

  //----------------------------------------------------------------------
  std::vector<std::string> LoadFileList(const std::string& listfile)
  {
    std::vector<std::string> ret;

    std::ifstream is(listfile);
    if(!is.good()){
      std::cerr << "Can't open file list '" << listfile << "'. Aborting." << std::endl;
      abort();
    }

    while(!is.eof()){
      std::string fname;
      is >> fname;
      if(!fname.empty()) ret.push_back(fname);
    }
    return ret;
  }

  //----------------------------------------------------------------------
  bool AlmostEqual(double a, double b)
  {
    if(a == 0 && b == 0) return true;

    return fabs(a-b)/std::max(a, b) < .0001; // allow 0.01% error
  }

  //----------------------------------------------------------------------
  // Note that this does not work for 3D!
  TH1* GetMaskHist(const Spectrum& s, double xmin, double xmax, double ymin, double ymax)
  {
    if (s.GetBinnings().size() > 2){
      std::cout << "Error: unable to apply a mask in " << s.GetBinnings().size() << " dimensions" << std::endl;
      abort();
    }

    // The exposure isn't important here
    TH1* fMaskND  = s.ToTHX(s.POT());
    TH1D* fMask1D = s.ToTH1(s.POT());
    fMask1D->Reset();

    int ybins = fMaskND->GetNbinsY();

    for(int i = 0; i < fMask1D->GetNbinsX()+2; ++i){

      int ix = i / ybins;
      int iy = i % ybins;

      bool isMask = false;

      if (xmin < xmax){
	if (fMaskND->GetXaxis()->GetBinLowEdge(ix+1) < xmin) isMask=true;
	if (fMaskND->GetXaxis()->GetBinUpEdge(ix+1) > xmax) isMask=true;
      }

      if (ymin < ymax){
	if (fMaskND->GetYaxis()->GetBinLowEdge(iy+1) < ymin) isMask=true;
	if (fMaskND->GetYaxis()->GetBinUpEdge(iy+1) > ymax) isMask=true;
      }

      if(!isMask) fMask1D->SetBinContent(i+1, 1);
    }
    return fMask1D;
  }
}
