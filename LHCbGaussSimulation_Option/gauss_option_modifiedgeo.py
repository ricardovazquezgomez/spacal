from Configurables import LHCbApp, CondDB, UpdateManagerSvc,DDDBConf
from Gaudi.Configuration import * 

#### $GAUSSOPTS/Gauss-2012.py \
from Gaudi.Configuration import *
importOptions("$APPCONFIGOPTS/Gauss/Sim08-Beam4000GeV-md100-2012-nu2.5.py")
importOptions("$APPCONFIGOPTS/Gauss/DataType-2012.py")
importOptions("$APPCONFIGOPTS/Gauss/RICHRandomHits.py")
importOptions("$APPCONFIGOPTS/Gauss/NoPacking.py")


DDDBConf().DbRoot = "/afs/cern.ch/user/z/zhangy/workdir/EcalUpgrade/GaussDev_v52r2/myDDDB-tutorial/lhcb.xml" 
CondDB().Tags['DDDB'] = 'HEAD'

#--Set database tags 
from Configurables import LHCbApp
LHCbApp().DDDBtag   = "dddb-20171030-3"
LHCbApp().CondDBtag = "sim-20170721-2-vc-md100"

#### $APPCONFIGOPTS/Gauss/Sim08-Beam4000GeV-mu100-2012-nu2.5.py
importOptions("$APPCONFIGOPTS/Gauss/Sim08-Beam4000GeV-mu100-2012-nu2.5.py")

#### $APPCONFIGOPTS/Gauss/DataType-2012.py
from Configurables import Gauss
Gauss().DataType  = "2012"

####Gauss-Job.py
from Gauss.Configuration import *

#--Generator phase, set random numbers
GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1
GaussGen.RunNumber        = 1082

#--Number of events
nEvts = 10
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

##### mine
importOptions("$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmNoCuts.py")
importOptions("$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py")
importOptions("$DECFILESROOT/options/30000000.py")
importOptions("$LBPYTHIA8ROOT/options/Pythia8.py")

## speed-up it
Gauss().SpilloverPaths = [] 
Gauss().Phases = [ 'Generator' , 'Simulation' ] 
