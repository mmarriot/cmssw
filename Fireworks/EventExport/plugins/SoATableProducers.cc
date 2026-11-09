// SoATableProducer instantiations. Types that exist only in some development areas (e.g. TICL v6)
// are guarded by __has_include, so the package builds against the release and against those areas.

#include "FWCore/Framework/interface/MakerMacros.h"
#include "Fireworks/EventExport/interface/SoATableProducer.h"

#include "DataFormats/HGCalDigi/interface/HGCalDigiHost.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"

using HGCalDigiSoATableProducer = SoATableProducer<hgcaldigi::HGCalDigiHost>;
using HGCalRecHitSoATableProducer = SoATableProducer<HGCalSoARecHitsHostCollection>;
DEFINE_FWK_MODULE(HGCalDigiSoATableProducer);
DEFINE_FWK_MODULE(HGCalRecHitSoATableProducer);

#if __has_include("DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h")
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
using CaloClusterSoATableProducer = SoATableProducer<reco::CaloClusterHostCollection>;
DEFINE_FWK_MODULE(CaloClusterSoATableProducer);
#endif

#if __has_include("DataFormats/HGCalReco/interface/TracksterHostCollection.h")
#include "DataFormats/HGCalReco/interface/TracksterHostCollection.h"
using TracksterSoATableProducer = SoATableProducer<ticl::TracksterHostCollection>;
DEFINE_FWK_MODULE(TracksterSoATableProducer);
#endif
