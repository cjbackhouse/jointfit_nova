#include "CAFAna/Core/SpectrumLoader.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Core/GenieWeightList.h"

#include "StandardRecord/StandardRecord.h"

#include <cassert>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

namespace ana
{
  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::string& wildcard, DataSource src)
    : SpectrumLoaderBase(wildcard, src)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::vector<std::string>& fnames,
                                 DataSource src)
    : SpectrumLoaderBase(fnames, src)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(DataSource src)
    : SpectrumLoaderBase(src)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader SpectrumLoader::FromSAMProject(const std::string& proj,
                                                DataSource src,
                                                int fileLimit)
  {
    SpectrumLoader ret;
    ret.fSource = src;
    ret.fWildcard = "project "+proj;
    ret.fFileSource = std::unique_ptr<IFileSource>(new SAMProjectSource(proj, fileLimit));
    return ret;
  }

  //----------------------------------------------------------------------
  SpectrumLoader::~SpectrumLoader()
  {
  }

  struct CompareByID
  {
    bool operator()(const Cut& a, const Cut& b)
    {
      return a.ID() < b.ID();
    }
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::Go()
  {
    if(fGone){
      std::cerr << "Error: can only call Go() once on a SpectrumLoader" << std::endl;
      abort();
    }
    fGone = true;

    // Find all the unique cuts
    std::set<Cut, CompareByID> cuts;
    for(auto& shiftdef: fHistDefs)
      for(auto& cutdef: shiftdef.second)
        cuts.insert(cutdef.first);
    for(const Cut& cut: cuts) fAllCuts.push_back(cut);

    fLivetimeByCut.resize(fAllCuts.size());
    fPOTByCut.resize(fAllCuts.size());


    const int Nfiles = NFiles();

    Progress* prog = 0;

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      if(Nfiles >= 0 && !prog) prog = new Progress(TString::Format("Filling %lu spectra from %d files matching '%s'", fHistDefs.TotalSize(), Nfiles, fWildcard.c_str()).Data());

      //std::cout<<"starting 1 HandleFile() "<<std::endl;
      HandleFile(f); //, Nfiles == 1 ? prog : 0);
      //std::cout<<"finished 1 HandleFile() "<<std::endl;

      if(Nfiles > 1 && prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end for fileIdx

    StoreExposures();

    if(prog){
      prog->Done();
      delete prog;
    }

    ReportExposures();

    fHistDefs.RemoveLoader(this);
    fHistDefs.Clear();
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());
    TTree* tr;
    //    if(f->GetListOfKeys()->Contains("cafmaker")){
    //      tr = (TTree*)f->Get("cafmaker/caf");
    //    }
    //    else{
    //      tr = (TTree*)f->Get("mvaselect/MVASelection");
    //    }
    //tr = (TTree*)f->Get("caf");
    tr = (TTree*)f->Get("mvaselect/MVASelection");
    assert(tr);

    // Surely no-one will generate 1000 universes?
    std::vector<std::array<double, 1000>> genie_tmp;
    const std::vector<std::string> genie_names = GetGenieWeightNames();
    genie_tmp.resize(genie_names.size());
    std::vector<int> genie_size_tmp;
    genie_size_tmp.resize(genie_names.size());

    FloatingExceptionOnNaN fpnan(false);

    caf::StandardRecord sr;

    tr->SetBranchAddress("Ev_reco", &sr.dune.Ev_reco);
    tr->SetBranchAddress("Ev_reco_nue", &sr.dune.Ev_reco_nue);
    tr->SetBranchAddress("Ev_reco_numu", &sr.dune.Ev_reco_numu);
    tr->SetBranchAddress("Elep_reco", &sr.dune.Elep_reco);
    tr->SetBranchAddress("mvaresult", &sr.dune.mvaresult);
    tr->SetBranchAddress("mvanue", &sr.dune.mvanue);
    tr->SetBranchAddress("mvanumu", &sr.dune.mvanumu);
    tr->SetBranchAddress("cvnnue", &sr.dune.cvnnue);
    tr->SetBranchAddress("cvnnumu", &sr.dune.cvnnumu);
    tr->SetBranchAddress("numu_pid", &sr.dune.numu_pid);
    tr->SetBranchAddress("nue_pid", &sr.dune.nue_pid);
    tr->SetBranchAddress("reco_q", &sr.dune.reco_q);
    tr->SetBranchAddress("vtx_x", &sr.dune.vtx_x);
    tr->SetBranchAddress("det_x", &sr.dune.det_x);

    // CW: add variables that Chris (M) wants for ND selections
    //tr->SetBranchAddress("muon_exit", &sr.dune.muon_exit);
    //tr->SetBranchAddress("Ehad_veto", &sr.dune.Ehad_veto);

    tr->SetBranchAddress("Ev", &sr.dune.Ev);

    tr->SetBranchAddress("Elep", &sr.dune.Elep);
    //    tr->SetBranchAddress("ccnc", &sr.dune.ccnc);
    tr->SetBranchAddress("isCC", &sr.dune.isCC);
    //    tr->SetBranchAddress("beamPdg", &sr.dune.beamPdg);
    //    tr->SetBranchAddress("neu", &sr.dune.neu);
    tr->SetBranchAddress("nuPDG", &sr.dune.nuPDG);
    tr->SetBranchAddress("nuPDGunosc", &sr.dune.nuPDGunosc);
    tr->SetBranchAddress("LepPDG", &sr.dune.LepPDG);
    tr->SetBranchAddress("mode", &sr.dune.mode);
    tr->SetBranchAddress("nipi0", &sr.dune.nipi0);
    tr->SetBranchAddress("nipip", &sr.dune.nipip);
    tr->SetBranchAddress("nipim", &sr.dune.nipim);
    tr->SetBranchAddress("Q2", &sr.dune.Q2);
    tr->SetBranchAddress("W", &sr.dune.W);
    tr->SetBranchAddress("Y", &sr.dune.Y);
    //    tr->SetBranchAddress("cc", &sr.dune.cc);

    tr->SetBranchAddress("nuvtxx_truth",  &sr.dune.nuvtxx_truth);
    tr->SetBranchAddress("nuvtxy_truth",  &sr.dune.nuvtxy_truth);
    tr->SetBranchAddress("nuvtxz_truth",  &sr.dune.nuvtxz_truth);


    tr->SetBranchAddress("run", &sr.dune.run);
    tr->SetBranchAddress("isFD", &sr.dune.isFD);
    tr->SetBranchAddress("isFHC", &sr.dune.isFHC);

    tr->SetBranchAddress("sigma_Ev_reco", &sr.dune.sigma_Ev_reco);
    tr->SetBranchAddress("sigma_Elep_reco", &sr.dune.sigma_Elep_reco);
    tr->SetBranchAddress("sigma_numu_pid", &sr.dune.sigma_numu_pid);
    tr->SetBranchAddress("sigma_nue_pid", &sr.dune.sigma_nue_pid);

    sr.dune.genie_wgt.resize(genie_names.size());

    for(unsigned int i = 0; i < genie_names.size(); ++i){
      tr->SetBranchAddress(("wgt_"+genie_names[i]).c_str(),
                           &genie_tmp[i]);
      tr->SetBranchAddress((genie_names[i]+"_nshifts").c_str(),
                           &genie_size_tmp[i]);
    }


    const int Nentries = tr->GetEntries();
    std::cout<<"number of entries "<<Nentries<<std::endl;
    for(int n = 0; n < Nentries; ++n){
      tr->GetEntry(n);

      // // Patch up isFD and isFHC which aren't in MVAResult files
      // if(sr.dune.run == 20000001 ||
      //    sr.dune.run == 20000002 ||
      //    sr.dune.run == 20000003){
      //   sr.dune.isFD = true;
      //   sr.dune.isFHC = true;
      // }
      // else if(sr.dune.run == 20000004 ||
      //         sr.dune.run == 20000005 ||
      //         sr.dune.run == 20000006){
      //   sr.dune.isFD = true;
      //   sr.dune.isFHC = false;
      // }
      // else if(sr.dune.run == 1){
      //   // ND all is already set
      // }
      // else{
      //   std::cout << "SpectrumLoader: Unrecognized run: "
      //             << sr.dune.run << std::endl;
      //   abort();
      // }

      // Reformat the genie systs
      for(unsigned int i = 0; i < genie_names.size(); ++i){
        const int Nuniv = genie_size_tmp[i];
        assert(Nuniv >= 0 && Nuniv <= int(genie_tmp[i].size()));
        sr.dune.genie_wgt[i].resize(Nuniv);
        for(int j = 0; j < Nuniv; ++j){
          sr.dune.genie_wgt[i][j] = genie_tmp[i][j];
        }
      }

      //std::cout<<"starting a HandleRecord() "<<std::endl;
      HandleRecord(&sr);
      //std::cout<<"finished a HandleRecord() "<<std::endl;

      if(prog && n%10000 == 0) prog->SetProgress(double(n)/Nentries);
    } // end for n
    //std::cout<<"finished a set of HandleRecord()"<<std::endl;
  }

  //----------------------------------------------------------------------
  /// Helper for \ref HandleRecord
  template<class T, class U> class CutVarCache
  {
  public:
    CutVarCache() : fVals(U::MaxID()+1), fValsSet(U::MaxID()+1, false) {}

    inline T Get(const U& var, const caf::StandardRecord* sr)
    {
      const unsigned int id = var.ID();

      if(fValsSet[id]){
        return fVals[id];
      }
      else{
        const T val = var(sr);
        fVals[id] = val;
        fValsSet[id] = true;
        return val;
      }
    }

  protected:
    // Seems to be faster to do this than [unordered_]map
    std::vector<T> fVals;
    std::vector<bool> fValsSet;
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleRecord(caf::StandardRecord* sr)
  {
    // Some shifts only adjust the weight, so they're effectively nominal, but
    // aren't grouped with the other nominal histograms. Keep track of the
    // results for nominals in these caches to speed those systs up.
    CutVarCache<bool, Cut> nomCutCache;
    CutVarCache<double, Var> nomWeiCache;
    CutVarCache<double, Var> nomVarCache;
    //std::cout<<"in HandleRecord() "<<std::endl;

    for(auto& shiftdef: fHistDefs){
      const SystShifts& shift = shiftdef.first;

      // Need to provide a clean slate for each new set of systematic shifts to
      // work from. Unfortunately, copying the whole StandardRecord is pretty
      // expensive. So we need to rely on this slightly dangerous "Restorer"
      // mechanism.

      // Spot checks to try and make sure no-one misses adding a variable to
      // Restorer
      static int iterationNo = 0;
      // Prime means we should get good coverage over all combinations
      const int kTestIterations = 9973;

      const TestVals* save = 0;
      if(++iterationNo % kTestIterations == 0)
        save = GetVals(sr, shiftdef.second);

      Restorer* restore = 0;
      double systWeight = 1;
      bool shifted = false;
      // Can special-case nominal to not pay cost of Shift() or Restorer
      if(!shift.IsNominal()){
        restore = new Restorer;
        shift.Shift(*restore, sr, systWeight);
        // Did the Shift actually modify the event at all?
        shifted = !restore->Empty();
      }

      for(auto& cutdef: shiftdef.second){
        const Cut& cut = cutdef.first;

        const bool pass = shifted ? cut(sr) : nomCutCache.Get(cut, sr);
        // Cut failed, skip all the histograms that depended on it
        if(!pass) continue;

        for(auto& weidef: cutdef.second){
          const Var& weivar = weidef.first;

          double wei = shifted ? weivar(sr) : nomWeiCache.Get(weivar, sr);

          wei *= systWeight;
          if(wei == 0) continue;

          for(auto& vardef: weidef.second){
            if(vardef.first.IsMulti()){
              for(double val: vardef.first.GetMultiVar()(sr)){
                for(Spectrum* s: vardef.second.spects)
                  s->Fill(val, wei);
              }
              continue;
            }

            const Var& var = vardef.first.GetVar();

            const double val = shifted ? var(sr) : nomVarCache.Get(var, sr);

            if(std::isnan(val) || std::isinf(val)){
              std::cerr << "Warning: Bad value: " << val
                        << " returned from a Var. The input variable(s) could "
                        << "be NaN in the CAF, or perhaps your "
                        << "Var code computed 0/0?";
              std::cout << " Not filling into this histogram for this slice." << std::endl;
              continue;
            }

            for(Spectrum* s: vardef.second.spects) s->Fill(val, wei);

            for(ReweightableSpectrum* rw: vardef.second.rwSpects){
              const double yval = rw->ReweightVar()(sr);

              if(std::isnan(yval) || std::isinf(yval)){
                std::cerr << "Warning: Bad value: " << yval
                          << " for reweighting Var";
                std::cout << ". Not filling into histogram." << std::endl;
                continue;
              }

              // TODO: ignoring events with no true neutrino etc
              if(yval != 0) rw->fHist->Fill(val, yval, wei);
            } // end for rw
          } // end for vardef
        } // end for weidef
      } // end for cutdef

      // Delete Restorer at this point and return StandardRecord to its
      // unshifted form ready for the next histogram.
      delete restore;

      // Make sure the record went back the way we found it
      if(save){
        CheckVals(save, sr, shift.ShortName(), shiftdef.second);
        delete save;
      }
    } // end for shiftdef
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::ReportExposures()
  {
    // The POT member variables we use here were filled as part of
    // SpectrumLoaderBase::GetNextFile() as we looped through the input files.

    // Let's just assume no-one is using the Cut::POT() function yet, so this
    // printout remains relevant...

    std::cout << fPOT << " POT" << std::endl;
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::AccumulateExposures(const caf::SRSpill* spill)
  {
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    for(auto& shiftdef: fHistDefs){
      for(auto& cutdef: shiftdef.second){
        for(auto& weidef: cutdef.second){
          for(auto& vardef: weidef.second){
            for(Spectrum* s: vardef.second.spects) s->fPOT += fPOT;
            for(ReweightableSpectrum* rw: vardef.second.rwSpects) rw->fPOT += fPOT;
          }
        }
      }
    }

    // std::map<int, double> livetime;
    // std::map<int, double> pot;

    // for(unsigned int i = 0; i < fAllCuts.size(); ++i){
    //   const int id = fAllCuts[i].ID();
    //   if(fLivetimeByCut[i] < 0){
    //     fLivetimeByCut[i] = 0;
    //     std::cout << "WARNING: no way to compute livetime for FD data spectrum. If you want a livetime you need to be applying one of the cuts from TimingCuts.h or similar. You probably should be anyway to remove bad data near the spill ends." << std::endl;
    //   }
    //   livetime.emplace(id, fLivetimeByCut[i]);
    //   pot.emplace(id, fPOTByCut[i]);
    // }

    // for(auto& shiftdef: fHistDefs){
    //   for(auto& cutdef: shiftdef.second){
    //     const Cut& cut = cutdef.first;
    //     const int id = cut.ID();

    //     for(auto& weidef: cutdef.second){
    //       for(auto& vardef: weidef.second){
    //         for(Spectrum* s: vardef.second.spects){
    //           s->fPOT += pot[id];
    //           s->fLivetime += livetime[id];
    //         }

    //         for(ReweightableSpectrum* rw: vardef.second.rwSpects){
    //           rw->fPOT += pot[id];
    //           rw->fLivetime += livetime[id];
    //         }
    //       }
    //     }
    //   }
    // }
  }

  //----------------------------------------------------------------------
  const SpectrumLoader::TestVals* SpectrumLoader::
  GetVals(const caf::StandardRecord* sr,
          IDMap<Cut, IDMap<Var, IDMap<VarOrMultiVar, SpectList>>>& hists) const
  {
    TestVals* ret = new TestVals;

    // Store values for all Vars and Cuts of interest
    for(auto& cutdef: hists){
      const bool cutval = cutdef.first(sr);
      ret->cuts.push_back(cutval);
      // Don't evaluate Vars when the Cut fails, might not be safe
      if(!cutval) continue;

      for(auto& weidef: cutdef.second){
        ret->weis.push_back(weidef.first(sr));

        for(auto& vardef: weidef.second){
          if(!vardef.first.IsMulti())
            ret->vars.push_back((vardef.first.GetVar())(sr));
        }
      }
    }

    return ret;
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::ValError(const std::string& type,
                                const std::string& shift,
                                const std::set<std::string>& /*req*/,
                                double orig, double now) const
  {
    // Try and print a comprehensive error message, I imagine this might be
    // hard to track down.

    std::cerr << std::endl;

    std::cerr << "Error. Value of " << type
              << " changed after it was shifted and then restored."
              << std::endl;

    std::cerr << "While applying shift " << shift;

    std::cerr << " initially had value " << orig
              << " now has " << now << std::endl;

    std::cerr << "Please check your use of Restorer very carefully"
              << std::endl;

    abort();
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::CheckVals(const TestVals* v,
                                 const caf::StandardRecord* sr,
                                 const std::string& shiftName,
                                 IDMap<Cut, IDMap<Var, IDMap<VarOrMultiVar, SpectList>>>& hists) const
  {
    unsigned int cutIdx = 0;
    unsigned int weiIdx = 0;
    unsigned int varIdx = 0;

    // Ensure everything is as TestVals says it should be

    for(auto& cutdef: hists){
      const bool cutval = cutdef.first(sr);

      if(cutval != v->cuts[cutIdx]){
        ValError("Cut", shiftName, {},
                 v->cuts[cutIdx], cutval);
      }
      ++cutIdx;

      // Don't evaluate Vars when the Cut fails, might not be safe
      if(!cutval) continue;

      for(auto& weidef: cutdef.second){
        const double weival = weidef.first(sr);
        if(!std::isnan(weival) && weival != v->weis[weiIdx]){
          ValError("Cut", shiftName, {},
                   v->weis[weiIdx], weival);
        }
        ++weiIdx;

        for(auto& vardef: weidef.second){
          if(vardef.first.IsMulti()) continue;
          const double varval = vardef.first.GetVar()(sr);
          if(!std::isnan(varval) && varval != v->vars[varIdx]){
            ValError("Var", shiftName, {},
                     v->vars[varIdx], varval);
          }
          ++varIdx;
        } // end for vardef
      } // end for weidef
    } // end for cutdef
  }
} // namespace
