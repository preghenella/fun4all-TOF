#include "TOFReco.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>
#include <phgenfit/Fitter.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <GenFit/RKTrackRep.h>
#include <GenFit/MeasuredStateOnPlane.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

TOFReco::TOFReco(const std::string& name, const std::string &oname)
  : SubsysReco(name)
  , mFileName(oname)
{
}

TOFReco::~TOFReco()
{
}

int TOFReco::Init(PHCompositeNode *topNode)
{

#if 0
  auto geo = PHGeomUtility::GetTGeoManager(topNode);
  auto field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
  auto fitter = PHGenFit::Fitter::getInstance(geo, field, "DafRef", "RKTrackRep", false);
  if (!fitter) {
    std::cout << "Cannot find PHGenFit::Fitter" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
#endif

  treeInit();
  histoInit();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TOFReco::treeInit()
{
  mOutFile = TFile::Open(mFileName.c_str(), "RECREATE");
  mOutTree = new TTree("tofreco", "TOF reco tree");
  mOutTree->Branch("N"        , &mOutTracks.N        , "N/I");
  mOutTree->Branch("charge"   , &mOutTracks.charge   , "charge[N]/I");
  mOutTree->Branch("chisq"    , &mOutTracks.chisq    , "chisq[N]/F");
  mOutTree->Branch("ndf"      , &mOutTracks.ndf      , "ndf[N]/I");
  mOutTree->Branch("p"        , &mOutTracks.p        , "p[N]/F");
  mOutTree->Branch("pt"       , &mOutTracks.pt       , "pt[N]/F");
  mOutTree->Branch("eta"      , &mOutTracks.eta      , "eta[N]/F");
  mOutTree->Branch("phi"      , &mOutTracks.phi      , "phi[N]/F");
  mOutTree->Branch("dcaxy"    , &mOutTracks.dcaxy    , "dcaxy[N]/F");
  mOutTree->Branch("dcaz"     , &mOutTracks.dcaz     , "dcaz[N]/F");
  mOutTree->Branch("time"     , &mOutTracks.time     , "time[N]/F");
  mOutTree->Branch("length"   , &mOutTracks.length   , "length[N]/F");
  mOutTree->Branch("_primary" , &mOutTracks._primary , "_primary[N]/O");
  mOutTree->Branch("_pid"     , &mOutTracks._pid     , "_pid[N]/I");
}

void TOFReco::histoInit()
{
}

int TOFReco::process_event(PHCompositeNode *topNode)
{
  
  /** get hits **/
  mHits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_PSTOF");
  if (!mHits) {
    std::cout << "Cannot find G4HIT_PSTOF" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

   /** get track map **/
  mTracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!mTracks) {
    std::cout << "Cannot find SvtxTrackMap" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /** get truth info **/
  mTruth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!mTruth) {
    std::cout << "Cannot find G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /** reset **/
  mOutTracks.N = 0;

  process_hits();
  process_tracks();

  /** fill tree **/
  mOutTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

void TOFReco::process_hits()
{

  /** group hits by detector id **/
  std::map<int, std::vector<PHG4Hit *>> detid_hits;
  for (auto ihit = mHits->getHits().first; ihit != mHits->getHits().second; ++ihit) {
    auto hit = dynamic_cast<PHG4Hit *>(ihit->second);
    if (!hit) continue;
    detid_hits[hit->get_detid()].push_back(hit);
  }

  /** flag detector hits that are not separated enough in time **/
  for (auto &hits : detid_hits) {

    /** sort hits by time form first to last hit on detector id **/
    std::sort(std::begin(hits.second),
	      std::end(hits.second),
	      [](PHG4Hit *a, PHG4Hit *b) {return a->get_t(0) < b->get_t(0);});
    
    auto prev_t = -999.;
    for (const auto &hit : hits.second) {
      if (hit->get_t(0) - prev_t < 10.)
	hit->set_hit_type(0);
      prev_t = hit->get_t(0);
    }
  }  
}

void TOFReco::process_tracks()
{

  TVector3 pos_h;
  TVector3 pos_i, pos_r, pos_f;
  TVector3 mom_i, mom_r, mom_f;
  TMatrixDSym cov_i(6), cov_r(6), cov_f(6);
  std::unique_ptr<genfit::AbsTrackRep> rep(nullptr);
  std::unique_ptr<genfit::MeasuredStateOnPlane> sop(nullptr);

  /** loop over tracks **/
  for (auto itrack = mTracks->begin(); itrack != mTracks->end(); ++itrack) {
    auto track = dynamic_cast<SvtxTrack *>(itrack->second);
    if (!track) continue;
    auto particle = mTruth->GetParticle(track->get_truth_track_id());
    if (!particle) continue;
    if (!is_track_in_acceptance(track)) continue;

    /** add track info **/
    auto N = mOutTracks.N++;
    mOutTracks.charge[N] = track->get_charge();
    mOutTracks.chisq[N] = track->get_chisq();
    mOutTracks.ndf[N] = track->get_ndf();
    mOutTracks.p[N] = track->get_p();
    mOutTracks.pt[N] = track->get_pt();
    mOutTracks.eta[N] = track->get_eta();
    mOutTracks.phi[N] = track->get_phi();
    mOutTracks.dcaxy[N] = track->get_dca3d_xy();
    mOutTracks.dcaz[N] = track->get_dca3d_z();
    /** reset TOF info **/
    mOutTracks.time[N] = 0.;
    mOutTracks.length[N] = 0.;
    /** add particle info **/
    mOutTracks._pid[N] = particle->get_pid();
    mOutTracks._primary[N] = particle->get_track_id() > 0;
    
    /** get last track state **/
    auto istate = --track->end_states();
    auto state = istate->second; 
    if (!state) continue;
    
    /** get state info **/      
    auto pidguess = track->get_charge() * 211;
    auto length_i = state->get_pathlength();
    pos_i.SetXYZ(state->get_x(), state->get_y(), state->get_z());
    mom_i.SetXYZ(state->get_px(), state->get_py(), state->get_pz());
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
	cov_i[i][j] = state->get_error(i, j);

    /** 
     ** propagate to average TOF radius
     ** this is a dummy step to help rejecting
     ** hits that are too far from the track
     **/

    rep.reset(new genfit::RKTrackRep(pidguess));
    sop.reset(new genfit::MeasuredStateOnPlane(rep.get()));    
    sop->setPosMomCov(pos_i, mom_i, cov_i);
    try {
      rep->extrapolateToCylinder(*sop, 85.); // this is not working properly, do not understand why
    }
    catch(...) {
      continue;
    }
    sop->getPosMomCov(pos_r, mom_r, cov_r);

    /** 
     ** propagate to all TOF hit planes 
     ** (so far is to hit point, to be done with plane)
     ** rejecting hits that are too far away
     **/

    std::vector<std::tuple<PHG4Hit *, float, float>> matchable;

    /** loop over hits **/
    for (auto ihit = mHits->getHits().first; ihit != mHits->getHits().second; ++ihit) {
      auto hit = dynamic_cast<PHG4Hit *>(ihit->second);
      if (!hit) continue;
      if (hit->get_hit_type() == 0) continue; // not digitised, need to find a better way
      pos_h.SetXYZ(hit->get_x(0), hit->get_y(0), hit->get_z(0));

      //      if ((pos_h - pos_r).Mag() > 50.) continue; // this is not working, do not understand why

      rep.reset(new genfit::RKTrackRep(pidguess));
      sop.reset(new genfit::MeasuredStateOnPlane(rep.get()));

      /** extrapolate track to hit point **/
      sop->setPosMomCov(pos_i, mom_i, cov_i);
      auto length = length_i;
      try {
	length += rep->extrapolateToPoint(*sop, pos_h, false, false);
      }
      catch(...) {
	continue;
      }
      sop->getPosMomCov(pos_f, mom_f, cov_f);

      auto distance = (pos_f - pos_h).Mag();
      if (distance > 10.) continue;

      matchable.push_back(std::make_tuple(hit, length, distance));
      
    }

    /** match the closest hit **/

    /** sort matchable hits by distance **/
    std::sort(std::begin(matchable),
	      std::end(matchable),
	      [](std::tuple<PHG4Hit *, float, float> a,
		 std::tuple<PHG4Hit *, float, float> b)
	      {return std::get<2>(a) < std::get<2>(b);});

    if (matchable.size() <= 0)
      continue;

    auto matched = matchable[0];
    auto hit = std::get<0>(matched);
    auto length = std::get<1>(matched);
    
    mOutTracks.time[N] = hit->get_t(0);
    mOutTracks.length[N] = length;

  }
  
}

bool TOFReco::is_track_in_acceptance(SvtxTrack *track)
{
  /** minimal-loose checks on track eta-pt **/
  if (fabs(track->get_eta()) > 1.5) return false;
  return true;
}

int TOFReco::End(PHCompositeNode* topNode)
{
  mOutFile->cd();
  mOutTree->Write();
  mOutFile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}
