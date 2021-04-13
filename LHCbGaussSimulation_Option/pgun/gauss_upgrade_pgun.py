from Gaudi.Configuration import *

importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu7.6-HorExtAngle.py")
##importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu11.4-HorExtAngle.py")
importOptions("$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py")
importOptions("11102201.py") #### Bd -> Kst gamma
#importOptions("12123002.py")  #### Bu -> K+  e+ e-
#importOptions("$LBPYTHIA8ROOT/options/Pythia8.py")
importOptions("$APPCONFIGOPTS/Gauss/Gauss-Upgrade-Baseline-20150522.py")
importOptions("$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py")
#importOptions("$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py")
importOptions("$APPCONFIGOPTS/Gauss/Upgrade.py")

from Gauss.Configuration import GenInit

GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1000
GaussGen.RunNumber = 1000

from Configurables import CondDB, Gauss ,DumpHepMCDecay
CondDB().Upgrade = True

from Configurables import LHCbApp
LHCbApp().DDDBtag    = "dddb-20171009"
LHCbApp().CondDBtag  = "sim-20170301-vc-md100" 


#LHCbApp().DDDBtag = 'upgrade/master'
#LHCbApp().DDDBtag = 'upgrade/upgrade/test-upgradeDB-yanxi'
#LHCbApp().CondDBtag ='upgrade/sim-20171127-vc-md100' 
LHCbApp().Simulation = True

LHCbApp().EvtMax = 50000


Gauss().DetectorGeo  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon', 'Magnet' ] }
Gauss().DetectorSim  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon', 'Magnet' ] }
Gauss().DetectorMoni = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Muon' ] }
Gauss().DetectorGeo  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Magnet' ] }
Gauss().DetectorSim  = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal', 'Magnet' ] }
Gauss().DetectorMoni = { "Detectors": ['VP', 'UT', 'FT', 'Rich1Pmt', 'Rich2Pmt', 'Ecal', 'Hcal' ] }
#updated
#Gauss.DetectorSim = {"Detectors": [ 'Magnet'] }
#Gauss.DetectorMoni ={"Detectors": [ 'Magnet'] }


## speed-up it
Gauss().SpilloverPaths = [] 
#Gauss().Phases = [ 'Generator' , 'Simulation' ] 
Gauss().Phases = ["Generator","GenToMCTree"]
Gauss().Production = 'PGUN'
#Gauss().OutputType = 'NONE' 


#smear primary vertices timing information
def CkvPVTimeSmearActivate():
    from Configurables import Generation
    from Configurables import BeamSpotMarkovChainSampleVertex
    Generation("Generation").addTool(BeamSpotMarkovChainSampleVertex, name="BeamSpotMarkovChainSampleVertex")
    Generation("Generation").VertexSmearingTool = "BeamSpotMarkovChainSampleVertex"
appendPostConfigAction(CkvPVTimeSmearActivate)

"""
generator = GaudiSequencer("generator")
generator.Members += ["DumpHepMCDecay"]
generator.DumpHepMCDecay.Depth = 99
generator.DumpHepMCDecay.PrintVertex = True
generator.DumpHepMCDecay.PrintEndVertex = True
generator.DumpHepMCDecay.Mode = 2
GenMonitor = GaudiSequencer( "GenMonitor" )
GenMonitor.Members += [ DumpHepMCDecay("DumpHepMCDecay")]
"""


from Configurables import GiGaGeo, GaussSensPlaneDet, TupleTool

#gg = GiGaGeo()
#gg.addTool ( GaussSensPlaneDet , name = 'PlaneDet' )
#pd = gg.PlaneDet 
#pd.OutputLevel = 1 
#pd.addTool ( TupleTool )
#
from Configurables import ApplicationMgr, NTupleSvc 
output_name = 'B2Kee_test'
output_name = 'B2KstGamma'
#output_name = 'test'
NTupleSvc().Output = [ 
        "FILE1 DATAFILE='%s.root' TYP='ROOT' OPT='NEW'"%(output_name,)
        ]
ApplicationMgr().ExtSvc += [ NTupleSvc() ] 


OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"%(output_name,)


#print Gauss() 
