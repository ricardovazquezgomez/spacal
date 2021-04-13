from os import environ
from GaudiKernel.SystemOfUnits import *
from Gaudi.Configuration import *
from Configurables import GaudiSequencer, CombineParticles
from Configurables import DecayTreeTuple, EventTuple, TupleToolTrigger, TupleToolTISTOS,FilterDesktop, MCDecayTreeTuple,PrintMCTree
from Configurables import BackgroundCategory, TupleToolDecay, TupleToolVtxIsoln,TupleToolPid,EventCountHisto,TupleToolRecoStats
from Configurables import LoKi__Hybrid__TupleTool, TupleToolVeto
# Unit
SeqPhys = GaudiSequencer("SeqPhys")

mct = MCDecayTreeTuple('mct')
mct.Decay = "[beauty => K- pi+ ^(pi0 -> ^gamma ^gamma)]CC"
#mct.Branches = { }
#mct = MCDecayTreeTuple('mct')
#mct.Decay = "[Lambda_c+ => ^p+ ^K- ^pi+]CC"
#mct.Branches = {
#        "Lc" :"[Lambda_c+ => p+ K- pi+ ]CC" ,
#        "Lcp":"[Lambda_c+ => ^p+ K- pi+]CC" ,
#        "LcK":"[Lambda_c+ => p+ ^K- pi+]CC" ,
#        "LcH":"[Lambda_c+ => p+ K- ^pi+]CC" ,
#        }
#mct = MCDecayTreeTuple('mct')
#mct.Decay = "[D0 => ^K- ^pi+]CC"
#mct.Branches = {
#        "Lc" :"[D0 => K- pi+]CC" ,
#        "LcK":"[D0 => ^K- pi+]CC" ,
#        "LcH":"[D0 => K- ^pi+]CC" ,
#        }
def doIt():
    """
    specific post-config action for (x)GEN-files 
    """
    extension = "xgen"
    ext = extension.upper()

    from Configurables import DataOnDemandSvc
    dod  = DataOnDemandSvc ()
    from copy import deepcopy 
    algs = deepcopy ( dod.AlgMap ) 
    bad  = set() 
    for key in algs :
        if     0 <= key.find ( 'Rec'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'Raw'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'DAQ'     )                  : bad.add ( key )
        elif   0 <= key.find ( 'Trigger' )                  : bad.add ( key )
        elif   0 <= key.find ( 'Phys'    )                  : bad.add ( key )
        elif   0 <= key.find ( 'Prev/'   )                  : bad.add ( key )
        elif   0 <= key.find ( 'Next/'   )                  : bad.add ( key )
        elif   0 <= key.find ( '/MC/'    ) and 'GEN' == ext : bad.add ( key )
    for b in bad :
        del algs[b]

    dod.AlgMap = algs

    from Configurables import EventClockSvc, CondDB 
    EventClockSvc ( EventTimeDecoder = "FakeEventTime" )
    CondDB  ( IgnoreHeartBeat = True )

appendPostConfigAction( doIt)
mctl=[  'MCTupleToolHierarchy', 'MCTupleToolKinematic', 'MCTupleToolPrimaries']
mct.ToolList=mctl 

printMC = PrintMCTree()
printMC.ParticleNames = ["Lambda_c+","Lambda_c~-"]
printMC.ParticleNames = ["D0","D~0"]
printMC.ParticleNames = ["Xi_cc++","Xi_cc~--"]
printMC.ParticleNames = ["B_c+","B_c-"]



########################################################################
from Configurables import DaVinci
DaVinci().EvtMax = -1
DaVinci().PrintFreq = 1
DaVinci().SkipEvents = 4                       # Events to skip
DaVinci().DataType = "2016"
DaVinci().Simulation   =  True
#DaVinci().DDDBtag ="dddb-20150703"
#DaVinci().CondDBtag = "sim-20150703-vc-md100"
DaVinci().TupleFile = "Tuple.root"             # Ntuple
DaVinci().UserAlgorithms = [mct,printMC]  
#DaVinci().UserAlgorithms = [mct]  
#DaVinci().Input=["Gauss-26166050-100ev-20170208.xgen"]
#DaVinci().Input=["Lc.xgen"]
#DaVinci().Input=["Xicc.xgen"]
DaVinci().Input=["nInter=1_Bd2KpPimPi0_100evt_Run1000_Evt0.sim"]
#DaVinci().Input=["D0.xgen"]
