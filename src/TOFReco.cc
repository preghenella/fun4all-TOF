#include "TOFReco.h"

#include <fun4all/Fun4AllReturnCodes.h>

TOFReco::TOFReco(const std::string& name)
  : SubsysReco(name)
{
}

TOFReco::~TOFReco()
{
}

int TOFReco::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TOFReco::process_event(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TOFReco::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
