# file /afs/cern.ch/user/z/zhangy/GaussDev_v49r12/Gen/DecFiles/options/12123002.py generated: Wed, 29 May 2019 17:51:41
#
# Event Type: 12123002
#
# ASCII decay Descriptor: [B+ -> K+ e+ e-]cc
#
#from Configurables import Generation
#Generation().EventType = 12123002
#Generation().SampleGenerationTool = "SignalRepeatedHadronization"
#from Configurables import SignalRepeatedHadronization
#Generation().addTool( SignalRepeatedHadronization )
#Generation().SignalRepeatedHadronization.ProductionTool = "PythiaProduction"
#Generation().SignalRepeatedHadronization.CutTool = "DaughtersInLHCb"
#Generation().SignalRepeatedHadronization.SignalPIDList = [ 521,-521 ]



# Ad-hoc particle gun code

from Configurables import ParticleGun
pgun = ParticleGun("ParticleGun")
pgun.SignalPdgCode = 521


from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "$DECFILESROOT/dkfiles/Bu_Kee=MS,DecProdCut.dec"

pgun.DecayTool = "EvtGenDecay"
pgun.GenCutTool = "DaughtersInLHCb"

from Configurables import FlatNParticles
pgun.NumberOfParticlesTool = "FlatNParticles"
pgun.addTool( FlatNParticles , name = "FlatNParticles" )

from Configurables import MomentumSpectrum
pgun.ParticleGunTool = "MomentumSpectrum"
pgun.addTool( MomentumSpectrum , name = "MomentumSpectrum" )
pgun.MomentumSpectrum.PdgCodes = [ 521,-521 ]
pgun.MomentumSpectrum.InputFile = "$PGUNSDATAROOT/data/Ebeam6500GeV/MomentumSpectrum_521.root"
pgun.MomentumSpectrum.BinningVariables = "pteta"
pgun.MomentumSpectrum.HistogramPath = "h_pteta"

from Configurables import BeamSpotSmearVertex
pgun.addTool(BeamSpotSmearVertex, name="BeamSpotSmearVertex")
pgun.VertexSmearingTool = "BeamSpotSmearVertex"
pgun.EventType = 12123002
