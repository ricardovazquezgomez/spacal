# file /build/jenkins-tests/workspace/nightly-builds/checkout/tmp/checkout/DBASE/Gen/DecFiles/v30r41/options/11102432.py generated: Mon, 16 Dec 2019 12:50:06
#
# Event Type: 11102432
#
# ASCII decay Descriptor: {[[B0]nos -> K+ pi- (pi0 -> gamma gamma)]cc, [[B0]os -> K- pi+ (pi0 -> gamma gamma)]cc}
#


from Configurables import ParticleGun
pgun = ParticleGun("ParticleGun")
pgun.SignalPdgCode = 511


from Configurables import ToolSvc
from Configurables import EvtGenDecay
ToolSvc().addTool( EvtGenDecay )
ToolSvc().EvtGenDecay.UserDecayFile = "$DECFILESROOT/dkfiles/Bd_K+pi-pi0=DecProdCut,sqDalitz.dec"

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
pgun.EventType = 11102432
