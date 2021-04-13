from ROOT import *
from array import array

import sys
print "Ecal: 12520 mm"
print "Neutron shielding: 12110 mm"
print "Distance: 12520-12110 = 410 m"

gROOT.ProcessLine(".x ~/lhcbStyle.C")


#sys.exit()
# B0 -> K+pi-pi0
#fi=TFile("B2KPiPi0_500k.root")
# B+ -> K+pi0
fi=TFile("B2KPi0_500k.root")
oldt = fi.Get("mct/MCDecayTree")

# B0 -> K+pi-pi0
#fs=TFile("Flux_B2KPiPi0_500k.root","recreate")
# B+ -> K+pi0
fs=TFile("Flux_B2KPi0_500k.root","recreate")
nt=TTree("tree","")

prod_vertex_x = array("d",[0.])
prod_vertex_y = array("d",[0.])
prod_vertex_z = array("d",[0.])

entry_x = array("d",[0.])
entry_y = array("d",[0.])
entry_z = array("d",[0.])

px = array("d",[0.])
py = array("d",[0.])
pz = array("d",[0.])

G4_index = array("i",[0])
pdgID = array("i",[0])

mother_index=array ("i",[0])

charge = array("d",[0])
timing = array("d",[0])
eKinetic = array("d",[0])
et = array("d",[0])
eTot = array("d",[0])

runNumber = array("i",[0])
evtNumber = array("i",[0])
nt.Branch("runNumber",runNumber,"runNumber/I")
nt.Branch("evtNumber",evtNumber,"evtNumber/I")
evtIndex = array("i",[0])
nt.Branch("evtIndex",evtIndex,"evtIndex/I")

nt.Branch("prod_vertex_x",prod_vertex_x,"prod_vertex_x/D")
nt.Branch("prod_vertex_y",prod_vertex_y,"prod_vertex_y/D")
nt.Branch("prod_vertex_z",prod_vertex_z,"prod_vertex_z/D")

nt.Branch("entry_x",entry_x,"entry_x/D")
nt.Branch("entry_y",entry_y,"entry_y/D")
nt.Branch("entry_z",entry_z,"entry_z/D")


nt.Branch("px",px,"px/D")
nt.Branch("py",py,"py/D")
nt.Branch("pz",pz,"pz/D")

nt.Branch("G4index",G4_index,"G4index/I")
nt.Branch("pdgID",pdgID,"pdgID/I")
nt.Branch("mother_index",mother_index,"mother_index/I")

nt.Branch("charge",charge,"charge/D")
nt.Branch("timing",timing,"timing/D")
nt.Branch("eKinetic",eKinetic,"eKinetic/D")
nt.Branch("eTot",eTot,"eTot/D")

currentID = -9999

evtIndex[0] = -1
runN,evtN=0,0
from math import sqrt,atan
evtIndex[0] = -1
for evt in oldt:
    if not abs(evt.pi0_MC_MOTHER_ID) == 511 and not abs(evt.pi0_MC_MOTHER_ID) == 521:continue
    if evt.gamma1_TRUEP_Z<0.1 or evt.gamma1_TRUEP_Z<0.1:continue
    if atan(evt.gamma1_TRUEPT/evt.gamma1_TRUEP_Z)>0.25 or atan(evt.gamma1_TRUEPT/evt.gamma1_TRUEP_Z)>0.25 :continue

    prod_vertex_x[0] = evt.gamma1_TRUEORIGINVERTEX_X
    prod_vertex_y[0] = evt.gamma1_TRUEORIGINVERTEX_Y
    prod_vertex_z[0] = evt.gamma1_TRUEORIGINVERTEX_Z
    px[0] = evt.gamma1_TRUEP_X/1000.
    py[0] = evt.gamma1_TRUEP_Y/1000.
    pz[0] = evt.gamma1_TRUEP_Z/1000.
    #### propagated to z=12620
    entry_z[0] = 12620.
    entry_x[0] = px[0]/pz[0]*(12620.-prod_vertex_z[0]) + prod_vertex_x[0]
    entry_y[0] = py[0]/pz[0]*(12620.-prod_vertex_z[0]) + prod_vertex_y[0]

    eTot[0] = evt.gamma1_TRUEP_E/1000.
    eKinetic[0] = eTot[0]
    et[0] = evt.gamma1_TRUEPT/1000.



    G4_index[0] = 0
    pdgID[0] = 22
    mother_index[0] = 0
    timing[0] = (12620.-prod_vertex_z[0])*1.E+9/(300000.*1000.*1000.) * (eTot[0]/pz[0]) + evt.gamma1_TRUEORIGINVERTEX_E
    charge[0] = 0

    runNumber[0] = 0
    evtNumber[0] = 0
    evtIndex[0]  += 1

    nt.Fill()

    prod_vertex_x[0] = evt.gamma2_TRUEORIGINVERTEX_X
    prod_vertex_y[0] = evt.gamma2_TRUEORIGINVERTEX_Y
    prod_vertex_z[0] = evt.gamma2_TRUEORIGINVERTEX_Z
    px[0] = evt.gamma2_TRUEP_X/1000.
    py[0] = evt.gamma2_TRUEP_Y/1000.
    pz[0] = evt.gamma2_TRUEP_Z/1000.
    #### propagated to z=12620
    entry_z[0] = 12620.
    entry_x[0] = px[0]/pz[0]*(12620.-prod_vertex_z[0]) + prod_vertex_x[0]
    entry_y[0] = py[0]/pz[0]*(12620.-prod_vertex_z[0]) + prod_vertex_y[0]

    eTot[0] = evt.gamma2_TRUEP_E/1000.
    eKinetic[0] = eTot[0]
    et[0] = evt.gamma2_TRUEPT/1000.



    G4_index[0] = 0
    pdgID[0] = 22
    mother_index[0] = 0
    timing[0] = (12620.-prod_vertex_z[0])*1.E+9/(300000.*1000.*1000.) * (eTot[0]/pz[0]) + evt.gamma2_TRUEORIGINVERTEX_E
    charge[0] = 0

    runNumber[0] = 0
    evtNumber[0] = 0
    #evtIndex[0]  += 1   #not ht next event yet

    nt.Fill()
        
nt.Write()
print "Output Entires: ",nt.GetEntries()
