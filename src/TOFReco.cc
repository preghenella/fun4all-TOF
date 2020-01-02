#include "TOFReco.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>

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

#include "TH1F.h"
#include "TH2F.h"

TOFReco::TOFReco(const std::string& name)
  : SubsysReco(name)
{
}

TOFReco::~TOFReco()
{
}

int TOFReco::Init(PHCompositeNode *topNode)
{
  auto geo = PHGeomUtility::GetTGeoManager(topNode);
  auto field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
  auto fitter = PHGenFit::Fitter::getInstance(geo, field, "DafRef", "RKTrackRep", false);
  if (!fitter) {
    std::cout << "Cannot find PHGenFit::Fitter" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  histoInit();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TOFReco::histoInit()
{

  mHistoManager = new Fun4AllHistoManager(Name());
  mHistoManager->registerHisto( new TH1F("hTime_all"    , "hTime_all"    , 1000 , 0. , 100.) );
  mHistoManager->registerHisto( new TH1F("hTime_digit"  , "hTime_digit"  , 1000 , 0. , 100.) );
  mHistoManager->registerHisto( new TH1F("hTime_match"  , "hTime_match"  , 1000 , 0. , 100.) );
  mHistoManager->registerHisto( new TH2F("hBetaP_match" , "hBetaP_match" , 100, 0. , 10. , 100 , 0.1 , 1.1) );

  mHistoManager->registerHisto( new TH2F("hRecoTrack"   , "hRecoTrack"   , 30 , -1.5 , 1.5 , 100 , 0. , 10.) );
  mHistoManager->registerHisto( new TH2F("hMatchTrack"  , "hMatchTrack"  , 30 , -1.5 , 1.5 , 100 , 0. , 10.) );
  
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

  process_hits();
  process_tracks();

  return Fun4AllReturnCodes::EVENT_OK;
}

void TOFReco::process_hits()
{

  /** group hits by detector id **/
  std::map<int, std::vector<PHG4Hit *>> detid_hits;
  for (auto ihit = mHits->getHits().first; ihit != mHits->getHits().second; ++ihit) {
    auto hit = dynamic_cast<PHG4Hit *>(ihit->second);
    if (!hit) continue;

    std::cout << "Edep = " << hit->get_edep() << std::endl;
    auto trkid = hit->get_trkid();
    auto particle = mTruth->GetParticle(trkid);
    std::cout << " PID = " << particle->get_pid() << std::endl;

    if (hit->get_edep() <= 0.) {
      hit->set_hit_type(0);      
      continue;
    }



    detid_hits[hit->get_detid()].push_back(hit);
    dynamic_cast<TH1F *>(mHistoManager->getHisto("hTime_all"))->Fill(hit->get_t(0));
  }

  /** flag detector hits that are not separated enough in time **/
  for (auto &hits : detid_hits) {

    /** sort hits by time form first to last hit on detector id **/
    std::sort(std::begin(hits.second),
	      std::end(hits.second),
	      [](PHG4Hit *a, PHG4Hit *b) {return a->get_t(0) < b->get_t(0);});
    
    auto prev_t = -999.;
    for (const auto &hit : hits.second) {
      if (hit->get_t(0) - prev_t < 10.) {
	hit->set_hit_type(0);
      } else {
	dynamic_cast<TH1F *>(mHistoManager->getHisto("hTime_digit"))->Fill(hit->get_t(0));
      }
      prev_t = hit->get_t(0);
    }
  }
  
}

void TOFReco::process_tracks()
{

  /** loop over tracks **/
  for (auto itrack = mTracks->begin(); itrack != mTracks->end(); ++itrack) {
    auto track = dynamic_cast<SvtxTrack *>(itrack->second);
    if (!track) continue;
    auto particle = mTruth->GetParticle(track->get_truth_track_id());
    if (!particle) continue;
    if (!is_track_in_acceptance(track)) continue;

    dynamic_cast<TH2F *>(mHistoManager->getHisto("hRecoTrack"))->Fill(track->get_eta(), track->get_pt());
    
    /** get last track state **/
    auto istate = --track->end_states();
    auto state = istate->second; 
    if (!state) continue;
    
    /** 
     ** propagate close to TOF 
     **/

    /** get state info **/      
    TVector3 _pos(state->get_x(), state->get_y(), state->get_z());
    TVector3 _mom(state->get_px(), state->get_py(), state->get_pz());
    TMatrixDSym _cov(6);
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
	_cov[i][j] = state->get_error(i, j);
    auto _length = state->get_pathlength();

#if 0
    
    /** extrapolate track **/
    try {
      _length += rep->extrapolateToCylinder(*sop, 85.);
    }
    catch(...) {
      continue;
    }

#endif
    
    /** get state info **/
    TVector3 pos, mom;
    TMatrixDSym cov(6);
    //    sop->getPosMomCov(pos, mom, cov);

    std::vector<std::tuple<PHG4Hit *, float, float>> matchable;

    /** loop over hits **/
    for (auto ihit = mHits->getHits().first; ihit != mHits->getHits().second; ++ihit) {
      auto hit = dynamic_cast<PHG4Hit *>(ihit->second);
      if (!hit) continue;
      if (hit->get_hit_type() == 0) continue;
      
      /** get hit point **/
      TVector3 pnt(hit->get_x(0), hit->get_y(0), hit->get_z(0));
      //      if ((pnt - pos).Mag() > 50.) continue;

      int pidguess = -211;
      auto rep = std::unique_ptr<genfit::AbsTrackRep>(new genfit::RKTrackRep(pidguess));
      auto sop = std::unique_ptr<genfit::MeasuredStateOnPlane>(new genfit::MeasuredStateOnPlane(rep.get()));

      /** extrapolate track to hit point **/
      sop->setPosMomCov(_pos, _mom, _cov);
      auto length = _length;
      try {
	length += rep->extrapolateToPoint(*sop, pnt, false, false);
      }
      catch(...) {
	continue;
      }
      
      /** get state info **/
      sop->getPosMomCov(pos, mom, cov);

      auto distance = (pos - pnt).Mag();
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
    auto beta = length / (hit->get_t(0) * 29.979246);     
    
    dynamic_cast<TH2F *>(mHistoManager->getHisto("hMatchTrack"))->Fill(track->get_eta(), track->get_pt());
    dynamic_cast<TH1F *>(mHistoManager->getHisto("hTime_match"))->Fill(hit->get_t(0));
    dynamic_cast<TH2F *>(mHistoManager->getHisto("hBetaP_match"))->Fill(track->get_p(), beta);
  
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
  mHistoManager->dumpHistos("TOFReco.root", "RECREATE");
  return Fun4AllReturnCodes::EVENT_OK;
}
