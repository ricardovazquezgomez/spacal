// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch
//

#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "MyMaterials.hh"
#include "CreateTree.hh"
#include "PrimaryGeneratorAction.hh"
#include "TRandom3.h"


#include <vector>

using namespace CLHEP;



EventAction::EventAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


EventAction::~EventAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::BeginOfEventAction (const G4Event* evt)
{
  G4int evtNb = evt->GetEventID () ;

   // G4cout << "---> Begin of Event: " << evtNb << G4endl ;
  // if ( evtNb%printModulo == 0 )
  // {
  //   G4cout << "---> Begin of Event: " << evtNb << G4endl ;
  // }


  CreateTree::Instance ()->Clear () ;




  G4PrimaryVertex * vertex = evt->GetPrimaryVertex () ;
  G4double x = vertex->GetX0 () ;
  G4double y = vertex->GetY0 () ;
  G4double z = vertex->GetZ0 () ;
  G4PrimaryParticle * particle = vertex->GetPrimary () ;
  G4double InitEnergy = particle->GetKineticEnergy () ;
  // G4cout << x << " " << y << " " << z << " " <<  InitEnergy << G4endl ;
  G4double px = particle->GetMomentumDirection().getX();
  G4double py = particle->GetMomentumDirection().getY();
  G4double pz = particle->GetMomentumDirection().getZ();
  G4int pdg = particle->GetPDGcode();

  // for the first event of the run
  if(evtNb == 0)
  {
    if(Parameters::Instance()->main_config.simulationType == 3)
    {
      CreateTree::Instance()->optX.push_back(x);
      CreateTree::Instance()->optY.push_back(y);
      CreateTree::Instance()->optZ.push_back(z);
      CreateTree::Instance()->optE.push_back(InitEnergy/CLHEP::eV);
    }
  }

  // --------------------- STORE INFO FOR X_0 / R_M ----------------------------- //
  int Radial_nSteps       = 5000;
  int Longitudinal_nSteps = 5000;
  CreateTree::Instance() -> Radial_stepLength       = 2500. / Radial_nSteps;       // in mm
  CreateTree::Instance() -> Longitudinal_stepLength = 2500. / Longitudinal_nSteps; // in mm

  // INSTANCE RUN/EVENT IN TREE
  CreateTree::Instance ()->Event = evt->GetEventID () ;
  CreateTree::Instance ()->inputMomentum->at (0) = px;
  CreateTree::Instance ()->inputMomentum->at (1) = py;
  CreateTree::Instance ()->inputMomentum->at (2) = pz;
  CreateTree::Instance ()->inputMomentum->at (3) = InitEnergy/MeV;
  CreateTree::Instance ()->inputInitialPosition->at (0) = x/mm;
  CreateTree::Instance ()->inputInitialPosition->at (1) = y/mm;
  CreateTree::Instance ()->inputInitialPosition->at (2) = z/mm;
  CreateTree::Instance ()->primary_type = pdg;

  // CreateTree::Instance()->events_per_run++;

}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void EventAction::EndOfEventAction(const G4Event* evt)
{
  evt -> GetEventID();

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  {
    CreateTree::Instance ()->EnergyPerModuleFront->Fill();
    CreateTree::Instance ()->EnergyPerModuleRear->Fill();

  }

  if(CreateTree::Instance ()->GetTree ()) CreateTree::Instance ()->Fill ();
  if(CreateTree::Instance()->summary) CreateTree::Instance()->summary->Fill();


  if(Parameters::Instance()->main_config.saveSimple)
  {
    // G4cout << ">>>>>>>>>>>>>>>>>>>>>>> " <<  CreateTree::Instance()->listOfTimestamps_front_back_0.size() << std::endl;
    // G4cout << ">>>>>>>>>>>>>>>>>>>>>>> " <<  CreateTree::Instance()->listOfTimestamps_front_back_1.size() << std::endl;
    //quick timestamp calculation


    //apply PDE
    std::vector<double> tAtfetPDE_0;
    std::vector<double> tAtfetPDE_1;
    double pde = CreateTree::Instance()->pde;
    std::cout << ">>>>><<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>< " << pde << std::endl;
    TRandom3 *rand = new TRandom3(0);
    int numb_of_phot_for_time_average = 5;
    for(unsigned int i = 0 ; i < CreateTree::Instance()->listOfTimestamps_front_back_0.size(); i++)
    {
      double numb = rand->Uniform(1.0);
      if(numb < pde)
      {
        tAtfetPDE_0.push_back(CreateTree::Instance()->listOfTimestamps_front_back_0.at(i));
      }
    }
    for(unsigned int i = 0 ; i < CreateTree::Instance()->listOfTimestamps_front_back_1.size(); i++)
    {
      double numb = rand->Uniform(1.0);
      if(numb < pde)
      {
        tAtfetPDE_1.push_back(CreateTree::Instance()->listOfTimestamps_front_back_1.at(i));
      }
    }
    // then, order the listOfTimestamps
    std::sort(tAtfetPDE_0.begin(),tAtfetPDE_0.end());
    std::sort(tAtfetPDE_1.begin(),tAtfetPDE_1.end());


    // set the k number, checking that it is not greater than the number of timestamps available
    // 0
    int effectiveN = numb_of_phot_for_time_average;
    if(numb_of_phot_for_time_average > tAtfetPDE_0.size())
    {
      effectiveN = tAtfetPDE_0.size();
    }
    for(int j = 0 ; j < effectiveN; j++)
    {
      CreateTree::Instance()->timestamp_front_back_0 +=   tAtfetPDE_0.at(j) / effectiveN; // avg
    }
    // 1
    if(numb_of_phot_for_time_average > tAtfetPDE_1.size())
    {
      effectiveN = tAtfetPDE_1.size();
    }
    for(int j = 0 ; j < effectiveN; j++)
    {
      CreateTree::Instance()->timestamp_front_back_1 +=   tAtfetPDE_1.at(j) / effectiveN; // avg
    }
    CreateTree::Instance()->simple_time->Fill();

    // G4cout << ">>>>>>>>>>>>>>>>>>>>>>> " <<  CreateTree::Instance()->timestamp_front_back_0 << std::endl;
    // G4cout << ">>>>>>>>>>>>>>>>>>>>>>> " <<  CreateTree::Instance()->timestamp_front_back_1 << std::endl;
  }

}
