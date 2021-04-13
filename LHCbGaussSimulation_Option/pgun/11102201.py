# file /afs/cern.ch/user/z/zhangy/GaussDev_v49r12/Gen/DecFiles/options/11102201.py generated: Wed, 29 May 2019 17:51:14
#
# Event Type: 11102201
#
# ASCII decay Descriptor: {[[B0]nos -> (K*(892)0 -> K+ pi-) gamma]cc, [[B0]os -> (K*(892)~0 -> K- pi+) gamma]cc}
#
#from Configurables import Generation
#Generation().EventType = 11102201
#Generation().SampleGenerationTool = "SignalRepeatedHadronization"
#from Configurables import SignalRepeatedHadronization
#Generation().addTool( SignalRepeatedHadronization )
#Generation().SignalRepeatedHadronization.ProductionTool = "PythiaProduction"
#Generation().SignalRepeatedHadronization.CutTool = "DaughtersInLHCb"
#Generation().SignalRepeatedHadronization.SignalPIDList = [ 511,-511 ]



# Ad-hoc particle gun code

from Configurables import ParticleGun
pgun = ParticleGun("ParticleGun")
pgun.SignalPdgCode = 511


from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "$DECFILESROOT/dkfiles/Bd_Kstgamma=DecProdCut.dec"

pgun.DecayTool = "EvtGenDecay"
pgun.GenCutTool = "DaughtersInLHCb"

from Configurables import FlatNParticles
pgun.NumberOfParticlesTool = "FlatNParticles"
pgun.addTool( FlatNParticles , name = "FlatNParticles" )

from Configurables import MomentumSpectrum
pgun.ParticleGunTool = "MomentumSpectrum"
pgun.addTool( MomentumSpectrum , name = "MomentumSpectrum" )
pgun.MomentumSpectrum.PdgCodes = [ 511,-511 ]
pgun.MomentumSpectrum.InputFile = "$PGUNSDATAROOT/data/Ebeam6500GeV/MomentumSpectrum_521.root"
pgun.MomentumSpectrum.BinningVariables = "pteta"
pgun.MomentumSpectrum.HistogramPath = "h_pteta"

from Configurables import BeamSpotSmearVertex
pgun.addTool(BeamSpotSmearVertex, name="BeamSpotSmearVertex")
pgun.VertexSmearingTool = "BeamSpotSmearVertex"
pgun.EventType = 11102201
