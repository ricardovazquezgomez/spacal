from ROOT import *
from array import array

import sys
print "Ecal: 12520 mm"
print "Neutron shielding: 12110 mm"
print "Distance: 12520-12110 = 410 m"

gROOT.ProcessLine(".x ~/lhcbStyle.C")

input_file = "/eos/lhcb/user/z/zhangy/EcalUpgrade/mu7.6_3000evt_Run1000_Evt0.root"

if len(sys.argv)>1:
    input_file = sys.argv[-1]


fo=TFile(input_file)
gDirectory.cd("GaussGeo.Ecal.Tuple")

tr = gDirectory.Get("Hits")

print "Input Entries: ",tr.GetEntries()

#sys.exit()
fs=TFile("Flux_mu7.6_Run1000_Evt0.root","recreate")
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
nt.Branch("et",et,"et/D")   #transverse energy
nt.Branch("eTot",eTot,"eTot/D")

currentID = -9999

#txt_output = open("source.dat","w")
evtIndex[0] = -1
runN,evtN=0,0
hNColl = TH1F("hNColl","",20,-0.5,19.5)
from math import sqrt
total = 0
filled = -1
for evt in tr:
    if evt.ID == filled: continue
    filled  = evt.ID
    if evt.pz<1.E-6:continue
    prod_vertex_x[0] = evt.prod_vertex_x
    prod_vertex_y[0] = evt.prod_vertex_y
    prod_vertex_z[0] = evt.prod_vertex_z

    px[0] = evt.px
    py[0] = evt.py
    pz[0] = evt.pz

    #### z between  12600 - 12629, propagated to 12620
    dz = 12620.-evt.entry_z
    #print evt.pz,evt.px,evt.py,evt.PDGID,evt.eKinetic,evt.eTot
    dt = dz/1000./(evt.pz/evt.eTot*299792458.)*1.E+9 #ns
    dx = evt.px/evt.pz*dz
    dy = evt.py/evt.pz*dz

    entry_x[0] = evt.entry_x + dx
    entry_y[0] = evt.entry_y + dy
    entry_z[0] = evt.entry_z + dz

    timing[0] = evt.timing + dt

    #total  += 1; if total >100: break
    #print evt.entry_z,12620,dx,dy,dz,dt


    G4_index[0] = evt.ID
    pdgID[0] = evt.PDGID
    mother_index[0] = evt.mother_ID
    charge[0] = evt.charge
    eKinetic[0] = evt.eKinetic/1000.  #to GeV
    et[0] = sqrt(pow(evt.entry_x,2.)+pow(evt.entry_y,2.))/evt.entry_z*evt.eKinetic
    eTot[0] = evt.eTot/1000.

    runNumber[0] = evt.runNumber
    evtNumber[0] = evt.evtNumber
    if (runN != evt.runNumber) or ( evtN != evt.evtNumber):
        evtIndex[0] += 1
        hNColl.Fill(evt.nColl)
        runN = evt.runNumber
        evtN = evt.evtNumber

    if G4_index[0] != currentID:
        nt.Fill()
        currentID = G4_index[0]
        #txt_output.write(
        #        #("%i "*3)%(runNumber[0], evtNumber[0], evtIndex[0])+
        #        ("%i "*1)%(evtIndex[0])+
        #        ("%f "*3)%(prod_vertex_x[0], prod_vertex_y[0], prod_vertex_z[0])+
        #        ("%f "*3)%(entry_x[0], entry_y[0], entry_z[0])+
        #        ("%f "*3)%(px[0], py[0], pz[0])+
        #        ("%f "*4)%(eKinetic[0],eTot[0],et[0],timing[0])+
        #        (("%i "*3)[:-1])%(G4_index[0], pdgID[0],mother_index[0])+
        #        "\n"
        #        )
        
nt.Write()
print "Output Entires: ",nt.GetEntries()

can = TCanvas("can","",600,500)
can.cd()
hNColl.Draw("e")
hNColl.SetXTitle("No. of collisions")
can.SaveAs("NCollisions.pdf")
