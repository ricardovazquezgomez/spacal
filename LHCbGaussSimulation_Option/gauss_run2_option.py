from Configurables import LHCbApp, CondDB, UpdateManagerSvc,DDDBConf, GenInit
from Gaudi.Configuration import * 


#from model Sim09g Model for 2018 - MD - 
importOptions("$APPCONFIGOPTS/Gauss/Beam6500GeV-md100-2018-nu1.6.py")
importOptions("$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py")
importOptions("$APPCONFIGOPTS/Gauss/DataType-2017.py")
importOptions("$APPCONFIGOPTS/Gauss/RICHRandomHits.py")
importOptions("$LBPYTHIA8ROOT/options/Pythia8.py")
importOptions("$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py")
importOptions("$DECFILESROOT/options/30000000.py")


DDDBConf().DbRoot = "/afs/cern.ch/user/z/zhangy/workdir/EcalUpgrade/GaussDev_v52r2/DDDB-runI-II/lhcb.xml" 
CondDB().Tags['DDDB'] = 'HEAD'

#--Set database tags 
from Configurables import LHCbApp


LHCbApp().DDDBtag   = "dddb-20170721-3"
LHCbApp().CondDBtag = "sim-20190430-vc-mu100"

LHCbApp().DDDBtag   = "dddb-20171030-3"
LHCbApp().CondDBtag = "sim-20170721-2-vc-md100"

from Configurables import Gauss
Gauss().DataType  = "2012"

####Gauss-Job.py
#--Generator phase, set random numbers
GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1
GaussGen.RunNumber        = 1082

#--Number of events
nEvts = 5
LHCbApp().EvtMax = nEvts

#Gauss().OutputType = 'NONE'
#Gauss().Histograms = 'NONE'
#--Set name of output files for given job (uncomment the lines)
#  Note that if you do not set it Gauss will make a name based on event type,
#  number of events and the date
#idFile = 'GaussTest'
#HistogramPersistencySvc().OutputFile = idFile+'-histos.root'
#
#OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.sim' TYP='POOL_ROOTTREE' OPT='RECREATE'"%idFile

#GenMonitor = GaudiSequencer( "GenMonitor" )
#SimMonitor = GaudiSequencer( "SimMonitor" )
#GenMonitor.Members += [ "GaussMonitor::CheckLifeTimeHepMC/HepMCLifeTime" ]
#SimMonitor.Members += [ "GaussMonitor::CheckLifeTimeMC/MCLifeTime" ]



Gauss().SpilloverPaths = [] 
Gauss().Phases = [ 'Generator' , 'Simulation' ] 


## configuration of detector
from Configurables import GiGaGeo, SpdPrsSensDet, TupleTool
gg = GiGaGeo()
gg.addTool ( SpdPrsSensDet, name = 'PlaneDet' )
pd = gg.PlaneDet 
pd.OutputLevel = 1 
Spd = "/dd/Structure/LHCb/DownstreamRegion/Spd"
Prs = "/dd/Structure/LHCb/DownstreamRegion/Prs"
Ecal = "/dd/Structure/LHCb/DownstreamRegion/Ecal"
Hcal = "/dd/Structure/LHCb/DownstreamRegion/Hcal"

pd.Detector = Spd

