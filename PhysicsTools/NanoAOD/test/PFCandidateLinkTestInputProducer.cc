// Input for test_catch2_PFCandidateLinkTableProducer.cc: PF candidates with references to tracks, GSF tracks, muons
// and PF block elements of clusters (which needs the ProductIDs of those, so everything is made in one module).
//   "tracks" 3, "gsfTracks" 2, "muons" 2 (muon 1 on track 2), "ecal" 3 and "hcal" 2 and "other" 1 PF clusters,
//   "blocks": one block of elements ecal 2, hcal 0, ecal 0, other 0, and an element with a null cluster, and
//   "otherBlocks" a copy of it that no candidate refers to;
//   "ticl": a photon of pt 50 and a charged hadron of pt 60 on track 1, "reversedTicl" the same in reverse order;
//   "otherMuons": 2 muons that no candidate refers to;
//   "" (the candidates), pt 10 to 60:
//     0  charged hadron, track 0, elements ecal 2 and hcal 0
//     1  electron, GSF track 1, elements ecal 0 and the null cluster
//     2  muon, track 2, muon 1, element other 0
//     3  photon, nothing
//     4, 5  copies of the "ticl" candidates

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <vector>

class PFCandidateLinkTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit PFCandidateLinkTestInputProducer(edm::ParameterSet const&)
      : tracksToken_(produces("tracks")),
        gsfTracksToken_(produces("gsfTracks")),
        muonsToken_(produces("muons")),
        otherMuonsToken_(produces("otherMuons")),
        ecalToken_(produces("ecal")),
        hcalToken_(produces("hcal")),
        otherToken_(produces("other")),
        blocksToken_(produces("blocks")),
        otherBlocksToken_(produces("otherBlocks")),
        ticlToken_(produces("ticl")),
        reversedTiclToken_(produces("reversedTicl")),
        candidatesToken_(produces()) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const tracks = event.emplace(tracksToken_, 3);
    auto const gsfTracks = event.emplace(gsfTracksToken_, 2);
    std::vector<reco::Muon> muonList(2);
    muonList[1].setInnerTrack(reco::TrackRef(tracks, 2));
    auto const muons = event.emplace(muonsToken_, std::move(muonList));
    event.emplace(otherMuonsToken_, 2);
    auto const ecal = event.emplace(ecalToken_, 3);
    auto const hcal = event.emplace(hcalToken_, 2);
    auto const other = event.emplace(otherToken_, 1);

    reco::PFBlock block;
    for (auto const& cluster : {reco::PFClusterRef(ecal, 2),
                                reco::PFClusterRef(hcal, 0),
                                reco::PFClusterRef(ecal, 0),
                                reco::PFClusterRef(other, 0),
                                reco::PFClusterRef()}) {
      reco::PFBlockElementCluster element(cluster, reco::PFBlockElement::ECAL);
      block.addElement(&element);
    }
    auto const blocks = event.emplace(blocksToken_, 1, block);
    event.emplace(otherBlocksToken_, 1, block);
    auto const inBlock = [&](reco::PFCandidate& candidate, unsigned int element) {
      candidate.addElementInBlock(reco::PFBlockRef(blocks, 0), element);
    };

    std::vector<reco::PFCandidate> ticl;
    ticl.emplace_back(0, p4(50), reco::PFCandidate::gamma);
    ticl.emplace_back(1, p4(60), reco::PFCandidate::h);
    ticl.back().setTrackRef(reco::TrackRef(tracks, 1));

    std::vector<reco::PFCandidate> candidates;
    candidates.emplace_back(1, p4(10), reco::PFCandidate::h);
    candidates.back().setTrackRef(reco::TrackRef(tracks, 0));
    inBlock(candidates.back(), 0);
    inBlock(candidates.back(), 1);
    candidates.emplace_back(-1, p4(20), reco::PFCandidate::e);
    candidates.back().setGsfTrackRef(reco::GsfTrackRef(gsfTracks, 1));
    inBlock(candidates.back(), 2);
    inBlock(candidates.back(), 4);
    candidates.emplace_back(1, p4(30), reco::PFCandidate::mu);
    candidates.back().setTrackRef(reco::TrackRef(tracks, 2));
    candidates.back().setMuonRef(reco::MuonRef(muons, 1));
    inBlock(candidates.back(), 3);
    candidates.emplace_back(0, p4(40), reco::PFCandidate::gamma);
    candidates.insert(candidates.end(), ticl.begin(), ticl.end());
    event.emplace(reversedTiclToken_, ticl.rbegin(), ticl.rend());
    event.emplace(ticlToken_, std::move(ticl));
    event.emplace(candidatesToken_, std::move(candidates));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  static reco::Particle::LorentzVector p4(double pt) { return {pt, 0, 0, pt}; }

  const edm::EDPutTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDPutTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  const edm::EDPutTokenT<reco::MuonCollection> muonsToken_, otherMuonsToken_;
  const edm::EDPutTokenT<reco::PFClusterCollection> ecalToken_, hcalToken_, otherToken_;
  const edm::EDPutTokenT<reco::PFBlockCollection> blocksToken_, otherBlocksToken_;
  const edm::EDPutTokenT<reco::PFCandidateCollection> ticlToken_, reversedTiclToken_, candidatesToken_;
};

DEFINE_FWK_MODULE(PFCandidateLinkTestInputProducer);
