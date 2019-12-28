#ifndef TOFRECO_H__
#define TOFRECO_H__

#include <fun4all/SubsysReco.h>

class TOFReco : public SubsysReco
{

public:
  TOFReco(const std::string &name = "TOFReco");
  virtual ~TOFReco();

  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

 protected:

};

#endif
