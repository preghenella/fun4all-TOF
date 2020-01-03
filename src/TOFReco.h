#ifndef TOFRECO_H__
#define TOFRECO_H__

#include <fun4all/SubsysReco.h>
#include <string>

class Fun4AllHistoManager;
class PHG4HitContainer;
class SvtxTrackMap;
class PHG4TruthInfoContainer;
class SvtxTrack;

class TFile;
class TTree;

class TOFReco : public SubsysReco
{

public:
  
  TOFReco(const std::string &name = "TOFReco", const std::string &oname = "TOFReco.root");
  virtual ~TOFReco();

  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

 protected:

  enum EDestination_t {PLANE, RADIUS, POINT};

  void treeInit();
  void histoInit();
  void process_hits();
  void process_tracks();
  bool is_track_in_acceptance(SvtxTrack *track);
  
  /** input data members **/
  PHG4HitContainer *mHits;
  SvtxTrackMap *mTracks;
  PHG4TruthInfoContainer *mTruth;

  /** output data mambers **/
  std::string mFileName;
  TFile *mOutFile;
  TTree *mOutTree;
  static const int MAXTRACKS = 1024; 
  struct {
    int N;
    int    charge[MAXTRACKS];
    float   chisq[MAXTRACKS];
    int       ndf[MAXTRACKS];
    float       p[MAXTRACKS];
    float      pt[MAXTRACKS];
    float     eta[MAXTRACKS];
    float     phi[MAXTRACKS];
    float   dcaxy[MAXTRACKS];
    float    dcaz[MAXTRACKS];
    float    time[MAXTRACKS];
    float  length[MAXTRACKS];
    bool _primary[MAXTRACKS];
    int      _pid[MAXTRACKS];
  } mOutTracks;

  
};

#endif
