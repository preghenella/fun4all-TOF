R__ADD_INCLUDE_PATH($HOME/include)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libtofreco.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <tofreco/TOFReco.h>

void Fun4All_TOFReco(const char *fname)
{
  auto se = Fun4AllServer::instance();
  auto ana = new TOFReco("TOFReco");
  se->registerSubsystem(ana);

  auto in = new Fun4AllDstInputManager("in");
  in->fileopen(fname);
  se->registerInputManager(in);
  
  se->run();
  se->End();
}
