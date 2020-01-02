#ifndef TOFRECO_H__
#define TOFRECO_H__

#include <fun4all/SubsysReco.h>

class Fun4AllHistoManager;
class PHG4HitContainer;
class SvtxTrackMap;
class PHG4TruthInfoContainer;
class SvtxTrack;

class TOFReco : public SubsysReco
{

public:
  
  TOFReco(const std::string &name = "TOFReco");
  virtual ~TOFReco();

  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

 protected:

  enum EDestination_t {PLANE, RADIUS, POINT};

  void histoInit();
  void process_hits();
  void process_tracks();
  bool is_track_in_acceptance(SvtxTrack *track);
  
  Fun4AllHistoManager *mHistoManager;
  PHG4HitContainer *mHits;
  SvtxTrackMap *mTracks;
  PHG4TruthInfoContainer *mTruth;
  
};

#endif
