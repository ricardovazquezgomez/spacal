from Gaudi.Configuration import *

importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu7.6-HorExtAngle.py")
##importOptions("$APPCONFIGOPTS/Gauss/Beam7000GeV-md100-nu11.4-HorExtAngle.py")
importOptions("$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py")
#importOptions("$DECFILESROOT/options/30000000.py")
importOptions("11102432.py")
importOptions("$LBPYTHIA8ROOT/options/Pythia8.py")
importOptions("$APPCONFIGOPTS/Gauss/Gauss-Upgrade-Baseline-20150522.py")
importOptions("$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py")
#importOptions("$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py")
importOptions("$APPCONFIGOPTS/Gauss/Upgrade.py")

from Gauss.Configuration import GenInit

GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 0
GaussGen.RunNumber = 1000

from Configurables import CondDB, Gauss 
CondDB().Upgrade = True

from Configurables import LHCbApp
LHCbApp().DDDBtag    = "dddb-20171009"
LHCbApp().CondDBtag  = "sim-20170301-vc-md100" 


#LHCbApp().DDDBtag = 'upgrade/master'
#LHCbApp().DDDBtag = 'upgrade/upgrade/test-upgradeDB-yanxi'
#LHCbApp().CondDBtag ='upgrade/sim-20171127-vc-md100' 
LHCbApp().Simulation = True

LHCbApp().EvtMax = 100


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
#Gauss().OutputType = 'NONE' 


#smear primary vertices timing information
def CkvPVTimeSmearActivate():
    from Configurables import Generation
    from Configurables import BeamSpotMarkovChainSampleVertex
    Generation("Generation").addTool(BeamSpotMarkovChainSampleVertex, name="BeamSpotMarkovChainSampleVertex")
    Generation("Generation").VertexSmearingTool = "BeamSpotMarkovChainSampleVertex"

appendPostConfigAction(CkvPVTimeSmearActivate)


from Configurables import GiGaGeo, GaussSensPlaneDet, TupleTool

#gg = GiGaGeo()
#gg.addTool ( GaussSensPlaneDet , name = 'PlaneDet' )
#pd = gg.PlaneDet 
#pd.OutputLevel = 1 
#pd.addTool ( TupleTool )
#
from Configurables import ApplicationMgr, NTupleSvc 
name = "nInter=1_Bd2KpPimPi0"
name = name+"_%ievt_Run%i_Evt%i"%(LHCbApp().EvtMax,GaussGen.RunNumber,GaussGen.FirstEventNumber)
NTupleSvc().Output = [ 
        "FILE1 DATAFILE='%s.root' TYP='ROOT' OPT='NEW'"%(name,)
        ]
ApplicationMgr().ExtSvc += [ NTupleSvc() ] 


OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"%(name,)


#print Gauss() 
