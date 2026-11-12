###############################################################################
# Way to use this:
#   cmsRun testHGCalCellOutline_cfg.py geometry=D122 detector=HGCalEESensitive
#
#   Options for geometry D110, D120, D121, D122, D123, D124, ...
#   Options for detector HGCalEESensitive, HGCalHESiliconSensitive
#
###############################################################################
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('geometry',
                 "D122",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "geometry of operations: D110, D120, D121, D122, D123, D124, ...")
options.register('detector',
                 "HGCalEESensitive",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "detector: HGCalEESensitive, HGCalHESiliconSensitive")
options.parseArguments()
print(options)

geomName = "Run4" + options.geometry
geomFile = "Configuration.Geometry.GeometryExtended" + geomName + "Reco_cff"
import Configuration.Geometry.defaultPhase2ConditionsEra_cff as _settings
GLOBAL_TAG, ERA = _settings.get_era_and_conditions(geomName)
print("Geometry Name:   ", geomName)
print("Geom file Name:  ", geomFile)
print("Era Name:        ", ERA)

process = cms.Process("HGCalCellOutline",ERA)

process.load(geomFile)
process.load('Geometry.HGCalGeometry.hgcalCellOutlineTester_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')

if hasattr(process,'MessageLogger'):
    process.MessageLogger.HGCalGeom=dict()

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.hgcalCellOutlineTester.detectorName = options.detector
process.p1 = cms.Path(process.hgcalCellOutlineTester)
