// Checks HGCalGeometry::getCellOutline for all valid silicon cells of a detector:
//  - every outline is a counter-clockwise polygon with at least three corners
//  - in a full wafer, a cell whose outline is a regular hexagon is centred on getPosition (in both
//    endcaps) and, where getNewCorners uses the same cell size (v18 on), coincides with it
//  - the outlines of the cells of a wafer tile it: their areas add up to the area of the wafer
//    less the mouse bites at its corners (sampled once per partial type on a fine grid)
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/HGCalCommonData/interface/HGCalTypes.h"
#include "Geometry/HGCalCommonData/interface/HGCalWaferMask.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

class HGCalCellOutlineTester : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit HGCalCellOutlineTester(const edm::ParameterSet&);
  ~HGCalCellOutlineTester() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override {}
  void endRun(edm::Run const&, edm::EventSetup const&) override {}

private:
  const std::string nameDetector_;
  const double tolerance_;
  const edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomToken_;
};

HGCalCellOutlineTester::HGCalCellOutlineTester(const edm::ParameterSet& iC)
    : nameDetector_(iC.getParameter<std::string>("detectorName")),
      tolerance_(iC.getParameter<double>("tolerance")),
      geomToken_(esConsumes<HGCalGeometry, IdealGeometryRecord, edm::Transition::BeginRun>(
          edm::ESInputTag{"", nameDetector_})) {}

void HGCalCellOutlineTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detectorName", "HGCalEESensitive");
  desc.add<double>("tolerance", 1.e-3);
  descriptions.add("hgcalCellOutlineTester", desc);
}

void HGCalCellOutlineTester::beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
  const auto& geom = iSetup.getData(geomToken_);
  const auto& ddd = geom.topology().dddConstants();
  auto area = [](const std::vector<GlobalPoint>& p) {
    double a(0);
    for (unsigned int i = 0; i < p.size(); ++i)
      a += (p[i].x() * p[(i + 1) % p.size()].y() - p[(i + 1) % p.size()].x() * p[i].y());
    return 0.5 * a;
  };

  // (layer, zside, waferU, waferV) -> sum of the cell areas, number of cells
  std::map<std::tuple<int, int, int, int>, std::pair<double, int>> wafers;
  unsigned int nCell(0), nBad(0), nHex(0);
  double dHex(0), dCentre(0);
  for (const DetId& id : geom.getValidDetIds()) {
    HGCSiliconDetId sid(id);
    auto outline = geom.getCellOutline(id);
    double a = area(outline);
    ++nCell;
    if ((outline.size() < 3) || (a <= 0)) {
      if (++nBad <= 10)
        edm::LogVerbatim("HGCalGeom") << "No outline for " << sid << " (" << outline.size() << " corners, area " << a
                                      << ")";
      continue;
    }
    auto& w = wafers[{sid.layer(), sid.zside(), sid.waferU(), sid.waferV()}];
    w.first += a;
    w.second++;
    auto side = [&outline](unsigned int k) {
      const GlobalPoint &p = outline[k], &q = outline[(k + 1) % outline.size()];
      return std::hypot(q.x() - p.x(), q.y() - p.y());
    };
    bool regular = (outline.size() == 6);
    for (unsigned int k = 1; regular && k < 6; ++k)
      regular = (std::abs(side(k) - side(0)) < 1.e-4 * side(0));
    if (regular && (ddd.partialWaferType(sid.layer(), sid.waferU(), sid.waferV()) == HGCalTypes::WaferFull)) {
      ++nHex;
      double cx(0), cy(0);
      for (const auto& p : outline) {
        cx += p.x() / outline.size();
        cy += p.y() / outline.size();
      }
      auto pos = geom.getPosition(id, false, false);
      dCentre = std::max(dCentre, std::hypot(cx - pos.x(), cy - pos.y()));
      auto corners = geom.getNewCorners(id);
      for (const auto& p : outline) {
        double d(1.e9);
        for (unsigned int k = 0; k + 1 < corners.size(); ++k)
          d = std::min(d, static_cast<double>(std::hypot(p.x() - corners[k].x(), p.y() - corners[k].y())));
        dHex = std::max(dHex, d);
      }
    }
  }

  std::map<int, double> waferArea;
  std::map<int, std::pair<int, int>> tiled;  // partial type -> wafers not tiled, wafers
  double dArea(0);
  for (const auto& [key, w] : wafers) {
    auto [layer, zside, waferU, waferV] = key;
    int part = ddd.partialWaferType(layer, waferU, waferV);
    if (waferArea.find(part) == waferArea.end()) {
      auto wxy = HGCalWaferMask::waferXY(part, 0, ddd.waferSize(true), 0.0, 0.0, 0.0, ddd.v17OrLess());
      auto fxy = HGCalWaferMask::waferXY(HGCalTypes::WaferFull, 0, ddd.waferSize(true), 0.0, 0.0, 0.0, ddd.v17OrLess());
      double rCorner = ddd.waferSize(true) / std::sqrt(3.0), h = rCorner - ddd.mouseBite(true);
      auto inside = [&wxy](double x, double y) {
        bool in(false);
        for (unsigned int i = 0, j = wxy.size() - 2; i + 1 < wxy.size(); j = i++)
          if (((wxy[i].second > y) != (wxy[j].second > y)) &&
              (x <
               (wxy[j].first - wxy[i].first) * (y - wxy[i].second) / (wxy[j].second - wxy[i].second) + wxy[i].first))
            in = !in;
        return in;
      };
      const int nGrid(2000);
      const double step = 2 * rCorner / nGrid;
      long nIn(0);
      for (int ix = 0; ix < nGrid; ++ix) {
        for (int iy = 0; iy < nGrid; ++iy) {
          double x = -rCorner + (ix + 0.5) * step, y = -rCorner + (iy + 0.5) * step;
          bool bitten(false);
          for (unsigned int k = 0; k + 1 < fxy.size(); ++k)
            bitten |= ((x * fxy[k].first + y * fxy[k].second) / rCorner > h);
          if ((!bitten) && inside(x, y))
            ++nIn;
        }
      }
      waferArea[part] = nIn * step * step;
    }
    double rel = std::abs(w.first - waferArea[part]) / waferArea[part];
    dArea = std::max(dArea, rel);
    ++tiled[part].second;
    if (rel > tolerance_) {
      if (++tiled[part].first <= 3)
        edm::LogVerbatim("HGCalGeom") << "Wafer " << layer << ":" << zside << ":" << waferU << ":" << waferV
                                      << " partial type " << part << " with " << w.second << " cells of area "
                                      << w.first << " for a wafer area of " << waferArea[part];
    }
  }

  edm::LogVerbatim("HGCalGeom") << nameDetector_ << ": " << nCell << " cells, " << nBad << " without an outline; "
                                << nHex << " regular hexagons in full wafers, max distance to getNewCorners " << dHex
                                << " cm, of their centre to getPosition " << dCentre << " cm";
  for (const auto& [part, n] : tiled)
    edm::LogVerbatim("HGCalGeom") << nameDetector_ << ": partial type " << part << ": " << n.first << " of " << n.second
                                  << " wafers not tiled within " << tolerance_;
  edm::LogVerbatim("HGCalGeom") << nameDetector_ << ": max relative difference of cell and wafer areas " << dArea;
}

DEFINE_FWK_MODULE(HGCalCellOutlineTester);
