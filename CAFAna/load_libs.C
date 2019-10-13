void load(std::string lib)
{
  std::cout << "." << std::flush;
  int ret = gSystem->Load(("lib"+lib).c_str());
  // In case of error, exit immediately with the error clearly showing, instead
  // of with a confusing secondary error about a page of output later.
  if(ret != 0){
    std::cout << std::endl << "gSystem->Load(lib"+lib+") failed with code " << ret << std::endl;
    exit(ret);
  }
}

void load_libs()
{
  // This magic incantation prevents ROOT doing slow cleanup work in
  // TList::RecursiveRemove() under ~TH1(). It also tends to lead to shutdown
  // crashes. This seems like a good compromise: go fast in batch mode
  // (probably fitting) and slow but not crashy in interactive move (probably
  // want to see the plots).
  if(gROOT->IsBatch()) gROOT->SetMustClean(false);

  // Colorize error messages. Would be better if we could just pick up the
  // flags novasoft uses, but they don't seem to be in any env var.
  gSystem->SetFlagsDebug(TString(gSystem->GetFlagsDebug())+" -fdiagnostics-color=auto");
  gSystem->SetFlagsOpt(TString(gSystem->GetFlagsOpt())+" -fdiagnostics-color=auto -UNDEBUG"); // match gcc's maxopt behaviour of retaining assert()

  // Include path
  TString includes = "-I$JOINTFIT_DIR -I$ROOTSYS/include";

  const std::vector<std::string> libs =
    {
      "OscLibFunc",
      "CAFAnaCore",
      "CAFAnaVars",
      "CAFAnaExperiment",
    };

  // Actually load the libraries
  std::cout << "Loading libraries";
  for(const std::string& lib: libs) load(lib);
  std::cout << std::endl;

  gSystem->SetIncludePath(includes);
}
