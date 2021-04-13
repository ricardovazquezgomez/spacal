from ROOT import *
from array import array

import sys
print "Ecal: 12520 mm"
print "Neutron shielding: 12110 mm"
print "Distance: 12520-12110 = 410 m"

gROOT.ProcessLine(".x ~/lhcbStyle.C")



#sys.exit()
fixed = True
name = "input_fakeGamma.root"
if fixed:
    name = "test.root"
fs=TFile(name,"recreate")
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

evtIndex[0] = -1
runN,evtN=0,0
from math import sqrt,atan,cos,sin
rnd = TRandom3()
evtIndex[0] = -1
for ii in range(1):


    prod_vertex_x[0] = 0.
    prod_vertex_y[0] = 0.
    prod_vertex_z[0] = 0.

    if fixed:
        eTot[0] = 1.+1.*(ii/100) #in GeV
        entry_x[0] = rnd.Uniform(-1.7,1.7)
        entry_y[0] = rnd.Uniform(-1.7,1.7)
        entry_z[0] = 12520.

        eKinetic[0] = eTot[0]

        ct = rnd.Uniform(-3.1416,3.1416)
        #ct = cos(45./180.*TMath.Pi())
        px[0] = eTot[0]*sin(4./180.*TMath.Pi()) * cos(ct)
        py[0] = eTot[0]*sin(4./180.*TMath.Pi()) * sin(ct)
        if rnd.Uniform(0.,1.)>0.5:
            py[0] = -py[0]
        pz[0] = eTot[0]*cos(4./180.*TMath.Pi())

        et[0] = 0.
    if not fixed:
        entry_x[0] = rnd.Uniform(-250.,250)  #in mm
        entry_y[0] = rnd.Uniform(-250.,250) 
        entry_z[0] = 12520.

        eTot[0] = rnd.Uniform(10.,20.) #inGeV
        eKinetic[0] = eTot[0]

        px[0] = eTot[0]*entry_x[0]/sqrt(entry_x[0]*entry_x[0] + entry_y[0]*entry_y[0] +entry_z[0]*entry_z[0])
        py[0] = eTot[0]*entry_y[0]/sqrt(entry_x[0]*entry_x[0] + entry_y[0]*entry_y[0] +entry_z[0]*entry_z[0])
        pz[0] = eTot[0]*entry_z[0]/sqrt(entry_x[0]*entry_x[0] + entry_y[0]*entry_y[0] +entry_z[0]*entry_z[0])

        et[0] = sqrt(px[0]*px[0]+py[0]*py[0])

    if atan(et[0]/pz[0])>0.4:continue


    G4_index[0] = 0
    pdgID[0] = 22
    mother_index[0] = 0
    timing[0] = 12.*1.E+9/(300000.*1000.)
    charge[0] = 0

    runNumber[0] = 0
    evtNumber[0] = 0
    evtIndex[0]  += 1

    nt.Fill()
        
nt.Write()
print "Output Entires: ",nt.GetEntries()

