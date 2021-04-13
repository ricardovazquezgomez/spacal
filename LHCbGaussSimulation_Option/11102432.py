# file /build/jenkins-tests/workspace/nightly-builds/checkout/tmp/checkout/DBASE/Gen/DecFiles/v30r41/options/11102432.py generated: Mon, 16 Dec 2019 12:50:06
#
# Event Type: 11102432
#
# ASCII decay Descriptor: {[[B0]nos -> K+ pi- (pi0 -> gamma gamma)]cc, [[B0]os -> K- pi+ (pi0 -> gamma gamma)]cc}
#
from Configurables import Generation
Generation().EventType = 11102432
Generation().SampleGenerationTool = "SignalRepeatedHadronization"
from Configurables import SignalRepeatedHadronization
Generation().addTool( SignalRepeatedHadronization )
Generation().SignalRepeatedHadronization.ProductionTool = "PythiaProduction"
from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "$DECFILESROOT/dkfiles/Bd_K+pi-pi0=DecProdCut,sqDalitz.dec"
Generation().SignalRepeatedHadronization.CutTool = "DaughtersInLHCb"
Generation().SignalRepeatedHadronization.SignalPIDList = [ 511,-511 ]

Generation().PileUpTool = "FixedNInteractions"
