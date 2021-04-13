from os import environ
from GaudiKernel.SystemOfUnits import *
from Gaudi.Configuration import *
from Configurables import GaudiSequencer, CombineParticles
from Configurables import DecayTreeTuple, EventTuple, TupleToolTrigger, TupleToolTISTOS,FilterDesktop, MCDecayTreeTuple,PrintMCTree
from Configurables import BackgroundCategory, TupleToolDecay, TupleToolVtxIsoln,TupleToolPid,EventCountHisto,TupleToolRecoStats
from Configurables import LoKi__Hybrid__TupleTool, TupleToolVeto
# Unit
SeqPhys = GaudiSequencer("SeqPhys")

#mct = MCDecayTreeTuple('mct')
#mct.Decay = "gamma"
#mct.Decay = "[Beauty -> ^K*(892)0 ^gamma]CC"
#mct.Branches = {
#        "gamma" :"[Beauty -> K*(892)0 ^gamma]CC" ,
#        "Beauty" :"[Beauty -> K*(892)0 gamma]CC" ,
#        "Kst" :"[Beauty -> ^K*(892)0 gamma]CC" ,
#	}
mct = MCDecayTreeTuple('mct')
# B0 -> K+pi-pi0
#mct.Decay = "[Beauty -> ^K+ ^pi- ^(pi0 -> ^gamma ^gamma)]CC"
#mct.Branches = {
#        "gamma1" :"[Beauty -> K+ pi- (pi0 -> gamma ^gamma)]CC" ,
#        "gamma2" :"[Beauty -> K+ pi- (pi0 -> ^gamma gamma)]CC" ,
#        "B" :"[Beauty -> K+ pi- (pi0 -> gamma gamma)]CC" ,
#        "K" :"[Beauty -> ^K+ pi- (pi0 -> gamma gamma)]CC" ,
#        "Pi" :"[Beauty -> K+ ^pi- (pi0 -> gamma gamma)]CC" ,
#	}
# B+ -> K+pi0
mct.Decay = "[Beauty -> ^K+ ^(pi0 -> ^gamma ^gamma)]CC"
mct.Branches = {
        "gamma1" :"[Beauty -> K+ (pi0 -> gamma ^gamma)]CC" ,
        "gamma2" :"[Beauty -> K+ (pi0 -> ^gamma gamma)]CC" ,
        "B" :"[Beauty -> K+ (pi0 -> gamma gamma)]CC" ,
        "K" :"[Beauty -> ^K+ (pi0 -> gamma gamma)]CC" ,
	}
#mctB = MCDecayTreeTuple('mctB')
#mctB.Decay = "[B0]CC"
#mctB.Branches = {
#        "B0" :"[B0]CC" ,
#        }
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
#mctl=[  'MCTupleToolKinematic']
mct.ToolList=mctl 
#mctB.ToolList=mctl 

printMC = PrintMCTree()
printMC.ParticleNames = ["Lambda_c+","Lambda_c~-"]
printMC.ParticleNames = ["D0","D~0"]
printMC.ParticleNames = ["Xi_cc++","Xi_cc~--"]
printMC.ParticleNames = ["B0","B~0"]
printMC.ParticleNames = ["B+","B-"]




########################################################################
from Configurables import DaVinci
DaVinci().EvtMax = -1
#DaVinci().PrintFreq = 1
#DaVinci().SkipEvents = 4                       # Events to skip
DaVinci().DataType = "2018"
DaVinci().Simulation   =  True
#DaVinci().DDDBtag ="dddb-20150703"
#DaVinci().CondDBtag = "sim-20150703-vc-md100"
#name = "B2KstGamma"
DaVinci().TupleFile = "/eos/experiment/spacal/users/zhangy/GeneratorPhoton_B2KstGammaRun1000Evt1000.root"
DaVinci().TupleFile = "/eos/experiment/spacal/users/zhangy/GeneratorPhoton_B2KPiPi0Run1000Evt1000.root"
#DaVinci().UserAlgorithms = [printMC]  
DaVinci().UserAlgorithms = [mct]  
#DaVinci().Input=["Gauss-26166050-100ev-20170208.xgen"]
#DaVinci().Input=["Lc.xgen"]
DaVinci().Input=["/eos/experiment/spacal/users/zhangy/B2KstGammaGenerator_Run1000FirstEvent1000.sim"]
DaVinci().Input=["/eos/experiment/spacal/users/zhangy/B2KPiPi0.sim"]
#DaVinci().Input=["B2Kee.sim"]
#DaVinci().Input=["D0.xgen"]
# B0 -> K+pi-pi0
#DaVinci().Input=["B2KPiPi0_500k.sim"]
#DaVinci().TupleFile="B2KPiPi0_500k.root"
# B+ -> K+pi0
DaVinci().Input=["B2KPi0_500k.sim"]
DaVinci().TupleFile="B2KPi0_500k.root"

