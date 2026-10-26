#include <numbers>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoTICL/Linking/interface/CornettoPairTest.h"

ticl::cornetto::PairTest ticl::cornetto::makePairTest(const edm::ParameterSet& conf) {
  const Cone cone{conf.getParameter<float>("maxBackwardDistance"),
                  conf.getParameter<float>("maxLongitudinalDistance"),
                  conf.getParameter<float>("transverseRadius0"),
                  conf.getParameter<float>("transverseSlope")};
  const float timeNSigma = conf.getParameter<float>("timeCompatibilityNSigma");
  const float timeFloor = conf.getParameter<float>("timeResolutionFloor");
  const PairTest test{conf.getParameter<float>("etaWindow"),
                      cone,
                      {cone.maxBackward,
                       cone.maxLongitudinal,
                       conf.getParameter<float>("emTransverseRadius0"),
                       conf.getParameter<float>("emTransverseSlope")},
                      timeNSigma > 0.f,
                      timeNSigma * timeNSigma,
                      timeFloor * timeFloor,
                      conf.getParameter<float>("typeVetoProbability")};
  const bool vetoValid =
      test.typeVetoProbability == 0.f || (test.typeVetoProbability >= 0.5f && test.typeVetoProbability <= 1.f);
  if (!(test.etaWindow > 0.f && test.etaWindow < std::numbers::pi_v<float>) || !coneValid(test.cone, test.emCone) ||
      !(timeFloor >= 0.f) || !vetoValid)
    throw cms::Exception("Configuration")
        << "Cornetto: invalid pair test parameters, see the parameter descriptions. etaWindow " << test.etaWindow
        << ", maxLongitudinalDistance " << test.cone.maxLongitudinal << ", maxBackwardDistance "
        << test.cone.maxBackward << ", transverseRadius0 " << test.cone.radius0 << ", transverseSlope "
        << test.cone.slope << ", emTransverseRadius0 " << test.emCone.radius0 << ", emTransverseSlope "
        << test.emCone.slope << ", timeResolutionFloor " << timeFloor << ", typeVetoProbability "
        << test.typeVetoProbability;
  return test;
}

void ticl::cornetto::fillPairTestDescription(edm::ParameterSetDescription& iDesc) {
  iDesc.add<float>("etaWindow", 0.3f)->setComment("Max |deta| and |dphi| between the barycenters. Must be in (0, pi).");
  iDesc.add<float>("maxLongitudinalDistance", 60.f)
      ->setComment("Max separation [cm] downstream of the anchor along its axis. Must be positive.");
  iDesc.add<float>("maxBackwardDistance", 20.f)
      ->setComment("Max separation [cm] upstream of the anchor along its axis. Must not be negative.");
  iDesc.add<float>("transverseRadius0", 4.f)
      ->setComment("Max transverse distance at zero separation [cm]. Must not be negative.");
  iDesc.add<float>("transverseSlope", 0.05f)
      ->setComment("Growth of the max transverse distance per cm of separation. Must not be negative.");
  iDesc.add<float>("emTransverseRadius0", 1.5f)
      ->setComment(
          "Max transverse distance at zero separation [cm] for an EM anchor. Must not be negative or above "
          "transverseRadius0.");
  iDesc.add<float>("emTransverseSlope", 0.05f)
      ->setComment(
          "Growth of the max transverse distance per cm of separation for an EM anchor. Must not be negative or "
          "above transverseSlope.");
  iDesc.add<float>("timeCompatibilityNSigma", 3.f)
      ->setComment(
          "Max |time difference at the anchor barycenter depth| in combined sigmas. Applies when both time "
          "errors are positive. Values <= 0 disable the gate.");
  iDesc.add<float>("timeResolutionFloor", 0.05f)
      ->setComment("Time resolution [ns] added in quadrature to the combined time error. Must not be negative.");
  iDesc.add<float>("typeVetoProbability", 0.6f)
      ->setComment(
          "EM/hadronic veto threshold: 0 disables the veto, else it must be in [0.5, 1]. A pair is refused if "
          "one EM probability is above this value and the other below 1 - this value. Applies only when both "
          "carry a PID.");
}
