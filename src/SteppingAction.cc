// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#include "SteppingAction.hh"

using namespace std;
using namespace CLHEP;

int to_int (string name)
{
  int Result ;             // int which will contain the result
  stringstream convert (name) ;
  string dummy ;
  convert >> dummy ;
  convert >> Result ;
  return Result ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


SteppingAction::SteppingAction (PrimaryGeneratorAction* gen_action):
fPrimaryGeneratorAction(gen_action)
{
  saveAll               = Parameters::Instance()->main_config.saveAll            ;
  saveTree              = Parameters::Instance()->main_config.saveTree           ;
  saveShower            = Parameters::Instance()->main_config.saveShower         ;
  saveStructure         = Parameters::Instance()->main_config.saveStructure      ;
  savePrimaries         = Parameters::Instance()->main_config.savePrimaries      ;
  savePhotons           = Parameters::Instance()->main_config.savePhotons        ;
  savePhotonAbsPoint    = Parameters::Instance()->main_config.savePhotonAbsPoint ;
  saveSummary           = Parameters::Instance()->main_config.saveSummary        ;
  saveLAPPD             = Parameters::Instance()->main_config.saveLAPPD          ;
  logging_volume        = Parameters::Instance()->main_config.logging_volume     ;
  pre_volume            = Parameters::Instance()->main_config.pre_volume         ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
SteppingAction::SteppingAction (const string& configFileName)
{
  ConfigFile config (configFileName) ;
}
SteppingAction::~SteppingAction ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

float SteppingAction::propagate(float zstart, float zend, float refractive){
    //G4cout<<zstart<<"\t"<<zend<<"\t"<<abs(zend-zstart)/(CLHEP::c_light/refractive)<<"\t"<<CLHEP::c_light<<G4endl;
    //use default unit: mm, ns, c_light=299...mm/ns
    return abs(zend-zstart)/(CLHEP::c_light/refractive);
}



void SteppingAction::UserSteppingAction (const G4Step * theStep)
{
  int photonDetected = 0;


  //---------------------------//
  // STEP INFO                 //
  //---------------------------//
  G4Track* theTrack = theStep->GetTrack () ;
  G4int trackID = theTrack->GetTrackID();
  
  // kill particle if the user has set an event time threshold
  if(Parameters::Instance()->main_config.event_time_cut > 0)
  {
    if(theTrack->GetGlobalTime()/CLHEP::ns > Parameters::Instance()->main_config.event_time_cut)
    { 
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      return;
    }
  }
  

  G4double energy = theStep->GetTotalEnergyDeposit()/MeV;
  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());
  G4ParticleDefinition* particleType = theTrack->GetDefinition () ;

  G4StepPoint * thePrePoint  = theStep->GetPreStepPoint () ;
  G4StepPoint * thePostPoint = theStep->GetPostStepPoint () ;
  // const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
  G4String thePrePVName  = "None" ; if ( thePrePV )  thePrePVName  = thePrePV  -> GetName () ;
  G4String thePostPVName = "None" ; if ( thePostPV ) thePostPVName = thePostPV -> GetName () ;

  G4TouchableHandle theTouchable = thePrePoint->GetTouchableHandle();
  G4ThreeVector origin(0.,0.,0.);
  G4ThreeVector globalOrigin = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
  // G4cout << thePrePV->GetName() << " "
  //        << globalOrigin.getX() << " "
  //        << globalOrigin.getY() << " "
  //        << globalOrigin.getZ() << " "
  //        << G4endl;

  // local coordinates

  // G4cout << "theTrackInfo->GetPrimaryID() = " << theTrackInfo->GetPrimaryID() << G4endl;
  // G4int nStep = theTrack -> GetCurrentStepNumber();
  G4String materialName = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName();
  G4String materialNamePre = ""; if(thePrePV) materialNamePre = thePrePV->GetLogicalVolume()->GetMaterial()->GetName();
  G4String materialNamePost = ""; if(thePostPV) materialNamePost = thePostPV->GetLogicalVolume()->GetMaterial()->GetName();

  //-------------------
  // get local position
  // G4double global_x = thePrePosition.x()/mm;
  // G4double global_y = thePrePosition.y()/mm;
  // G4double global_z = thePrePosition.z()/mm;
  G4ThreeVector OnDetectorPosition = thePostPoint->GetPosition();
  G4ThreeVector PreOnDetectorMomentum = thePrePoint->GetMomentumDirection();
  G4ThreeVector PostOnDetectorMomentum = thePostPoint->GetMomentumDirection();
  G4ThreeVector positionVector = thePrePoint->GetPosition(); //get the position vector
  G4ThreeVector momentumVector = thePrePoint->GetMomentum(); //get the position vector

  G4ThreeVector localPosition;
  localPosition.setX(positionVector.getX() - globalOrigin.getX());
  localPosition.setY(positionVector.getY() - globalOrigin.getY());
  localPosition.setZ(positionVector.getZ() - globalOrigin.getZ());


  int crystalID = -1;
  int cellID    = -1;
  int moduleID  = -1;
  Float_t module_x = 0;
  Float_t module_y = 0;
  Float_t module_z = 0;
  Int_t module_type = 0;
  int showerIsInCrystal;
  // G4cout<<"PV name, x: "<<thePrePVName<<"\t "<<positionVector.getX()<<G4endl;


  // G4int ModuleID_x=-1;
  // G4int ModuleID_y=-1;
  G4String volume_name="";
  // G4int volume_cpID=0;
  G4TouchableHandle preHistory =  thePrePoint->GetTouchableHandle();
  // G4bool frontCell = true;
  // G4int module_index=-1;
  // G4VPVParameterisation *module_parametrization;


  // if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  // {
  //
  //
  //   if(showerIsInCrystal==1)
  //   {
  //
  //     if(true)
  //     {
  //       // std::cout << "------------" << std::endl;
  //       int offsetX = 0;
  //       int offsetY = 0;
  //
  //       for(int ii=0;ii<preHistory->GetHistoryDepth();ii++)
  //       {
  //         volume_name = preHistory->GetVolume(ii)->GetName();
  //         if(volume_name.contains("Module"))
  //         {
  //           G4ThreeVector modulePosition = preHistory->GetVolume(ii)->GetTranslation();
  //           // preHistory->GetVolume(ii)->GetParameterisation();
  //           G4cout << modulePosition.getX() << " "
  //                  << modulePosition.getY() << " "
  //                  << modulePosition.getZ() << " "
  //                  << G4endl;
  //           // G4cout << ii << " " << volume_name << " " << preHistory->GetReplicaNumber(ii) << G4endl;
  //         }
  //       }
  //
  //       // if(CreateTree::Instance()->pipe == 0)
  //       // {
  //       //   for(int ii=0;ii<preHistory->GetHistoryDepth();ii++)
  //       //   {
  //       //     // no offset
  //       //     volume_name = preHistory->GetVolume(ii)->GetName();
  //       //
  //       //     // std::cout << ii << " " << volume_name << " " << preHistory->GetReplicaNumber(ii) << std::endl;
  //       //     if(volume_name.contains("modules_replica")) ModuleID_x = preHistory->GetReplicaNumber(ii);
  //       //     if(volume_name.contains("calo_full_PV")) ModuleID_y = preHistory->GetReplicaNumber(ii);
  //       //     // The material between two front and rear cell is not considered
  //       //
  //       //     // std::cout << " --------> " << ModuleID_x << " " << ModuleID_y << std::endl;
  //       //   }
  //       // }
  //       // else
  //       // {
  //       //   for(int ii=0;ii<preHistory->GetHistoryDepth();ii++)
  //       //   {
  //       //     volume_name = preHistory->GetVolume(ii)->GetName();
  //       //     // x
  //       //     if(volume_name.contains("module_PV")) ModuleID_x = preHistory->GetReplicaNumber(ii);
  //       //     if(volume_name.contains("calo_right_PV")) // right calo
  //       //     {
  //       //       // so x index is not replica but replica + offset
  //       //       offsetX = CreateTree::Instance()->pipe_modules_nx + (CreateTree::Instance()->module_array_x - CreateTree::Instance()->pipe_modules_nx)/2;
  //       //     }
  //       //
  //       //     // y
  //       //     if(volume_name.contains("calo_full_bottom_PV")) ModuleID_y = preHistory->GetReplicaNumber(ii);
  //       //     if(volume_name.contains("calo_short_left_row_PV"))
  //       //     {
  //       //       ModuleID_y = preHistory->GetReplicaNumber(ii);
  //       //       offsetY    = (CreateTree::Instance()->module_array_y - CreateTree::Instance()->pipe_modules_ny)/2;
  //       //     }
  //       //     if(volume_name.contains("calo_short_right_row_PV"))
  //       //     {
  //       //       ModuleID_y = preHistory->GetReplicaNumber(ii);
  //       //       offsetY    = (CreateTree::Instance()->module_array_y - CreateTree::Instance()->pipe_modules_ny)/2;
  //       //     }
  //       //     if(volume_name.contains("calo_full_top_PV"))
  //       //     {
  //       //       ModuleID_y = preHistory->GetReplicaNumber(ii);
  //       //       offsetY    = CreateTree::Instance()->pipe_modules_ny + (CreateTree::Instance()->module_array_y - CreateTree::Instance()->pipe_modules_ny)/2;
  //       //     }
  //       //   }
  //       // }
  //       ModuleID_x = ModuleID_x + offsetX;
  //       ModuleID_y = ModuleID_y + offsetY;
  //       // std::cout << " --------> " << ModuleID_x << " " << ModuleID_y << std::endl;
  //
  //       if(positionVector.getZ() > CreateTree::Instance()->zSeparationAbsolutePosition) frontCell=false;
  //
  //
  //       if(ModuleID_x!=-1 && ModuleID_y!=-1)
  //       {
  //         //It seems  that root doesn't support 2D array well, so convert into 1D array: index = ix+nx*iy
  //         if(ModuleID_x>=CreateTree::Instance()->module_array_x || ModuleID_y>=CreateTree::Instance()->module_array_y)
  //         G4cout<<"Number of modules in the calorimeter exceeds inputs !!!"<<G4endl;
  //         module_index = ModuleID_x + ModuleID_y*CreateTree::Instance()->module_array_x;
  //         if(frontCell)
  //         {
  //           CreateTree::Instance()->depositedEnergyByModuleFront[module_index] += energy;
  //           //need propagation
  //           CreateTree::Instance()->depositionTimingByModuleFront[module_index] += energy*(theTrack->GetGlobalTime() + propagate(positionVector.getZ(),CreateTree::Instance()->zSeparationAbsolutePosition -40.)); // this actually "propagates" the photons to front surface
  //         }
  //         if(!frontCell)
  //         {
  //           CreateTree::Instance()->depositedEnergyByModuleRear[module_index] += energy;
  //           //need propagation
  //           CreateTree::Instance()->depositionTimingByModuleRear[module_index] += energy*(theTrack->GetGlobalTime() + propagate(positionVector.getZ(),CreateTree::Instance()->zSeparationAbsolutePosition+100.)); //in ns*MeV // this actually "propagates" the photons to end surface.
  //         }
  //         //G4cout<<" Module_index "<<module_index<<"  : energy "<<energy/GeV<<" : Front? "<<frontCell<<" : "<<CreateTree::Instance()->depositedEnergyByModuleFront[module_index]<<" : "<<CreateTree::Instance()->depositedEnergyByModuleRear[module_index]<<G4endl;
  //       }
  //     }
  //
  //   }
  // }


  // ---------------------------//
  // OPTICAL PHOTONS            //
  // ---------------------------//
  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  {
    G4String processName ;
    // force scintillation as creator process if the optical photon
    // was created as a primary, check for the creator otherwise.
    // this is done because if the optical is a primary, the search for
    // creator process is impossible, and it would crash
    // the simulation

    // check if the primary was an optical photon
    // G4ParticleDefinition* particleDefinition = fPrimaryGeneratorAction->GetParticleDefinition();
    // std::string primaryParticleDefinitionString = Parameters::Instance()->particleDefinition;

    G4ParticleDefinition* particleDefinition = fPrimaryGeneratorAction->GetParticleDefinition();
    if(particleDefinition->GetParticleName() == "opticalphoton") // special case: optical created manually somewhere. will do some special stuff
    {
      processName = "Manual";
    }
    else
    {
      processName = theTrack->GetCreatorProcess()->GetProcessName();
    }

    // G4cout << ">>>>>>>>>>>>>>>>>> qui" << G4endl;

    // boundary process
    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

    //find the boundary process only once
    if(!boundary)
    {
      G4ProcessManager* pm = theTrack->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int i;
      for( i=0;i<nprocesses;i++)
      {
        if((*pv)[i]->GetProcessName()=="OpBoundary")
        {
          boundary = (G4OpBoundaryProcess*)(*pv)[i];
          break;
        }
      }
    }
    boundaryStatus = boundary->GetStatus();



    // more complex boundary check
    bool saveThePhoton = false;
    if((boundaryStatus == FresnelRefraction) && (thePostPVName.contains(logging_volume.c_str())))
    {
      if(pre_volume == "") // if it is not requested to check for pre volume name..
      {
        saveThePhoton = true;
        // photonDetected = 1;
      }
      else // check also pre volume name
      {
        if( thePrePVName.contains(pre_volume.c_str()))
        {
          saveThePhoton = true;
          // photonDetected = 1;
        }
      }
    }

    // //if optical calibration run, save the photon even if it's lost and not detected, but write down that it's not detected
    // if(opticalCalibrationRun)
    // {
    //   if((theTrack->GetTrackStatus() == fStopAndKill) || (theTrack->GetTrackStatus() == fKillTrackAndSecondaries))
    //   {
    //     saveThePhoton = true;
    //     photonDetected = 0;
    //   }
    // }



    if(saveThePhoton)
    {
      // get pmt

      // get mother logical volume name
      // set as detected

      // showerIsInCrystal = 1;
      // crystal number and cell number
      // crystal module name has this structure
      // crystal_N_cell_M_absorber_name-module
      // we wanto to extract N and M, the crystal and cell number
      // take the full volume name

      // //mod to test time
      // std::string namePV = thePrePVName;
      // // crystal_N_cell_M_absorber_name-module
      // // find position of first underscore
      // std::size_t cryBegin = namePV.find_first_of("_");
      // // take a substring that starts with the crystal number
      // std::string cryFullString = namePV.substr(cryBegin+1,namePV.size());
      // // N_cell_M_absorber_name-module
      // // now find the first underscore in this new string
      // std::size_t cryEnd = cryFullString.find_first_of("_");
      // // take a string with just the crystal number
      // std::string cryString = cryFullString.substr(0,cryEnd);
      // // N
      // // now take a string that starts with "cell"
      // std::string cellBaseString = cryFullString.substr(cryEnd+1,cryFullString.size());
      // // cell_M_absorber_name-module
      // // find first underscore in this string
      // std::size_t cellBegin = cellBaseString.find_first_of("_");
      // // extract a substring that starts with the cell number
      // std::string cellFullString = cellBaseString.substr(cellBegin+1,cellBaseString.size());
      // // M_absorber_name-module
      // std::size_t cellEnd = cellFullString.find_first_of("_");
      // // take a string with just the cell number
      // std::string cellString = cellFullString.substr(0,cellEnd);
      //
      // // G4cout << "---------------------" << G4endl;
      // // G4cout << namePV << " "
      // //        << cryFullString << " "
      // //        << cryString << " "
      // //        << cellBaseString << " "
      // //        << cellFullString << " "
      // //        << cellString << " "
      // //        << G4endl;
      //
      // std::istringstream ( cryString )  >> crystalID;
      // std::istringstream ( cellString ) >> cellID;
      //
      for(int ii=0;ii<preHistory->GetHistoryDepth();ii++)
      {
        volume_name = preHistory->GetVolume(ii)->GetName();
        if(volume_name.contains("Module"))
        {
          // find module type
          std::size_t modBegin = volume_name.find_first_of("_");
          std::string modString = volume_name.substr(modBegin+1,volume_name.size());
          std::istringstream ( modString ) >> module_type;
        }
      }
      // end of mod to test time

      G4String readoutName = thePostPV->GetMotherLogical()->GetName();

      // G4int front_back       = -1;
      // if(readoutName.contains("Negative"))
      // {
      //   front_back = 0;
      // }
      // if(readoutName.contains("Positive"))
      // {
      //   front_back = 1;
      // }

      // find pmt and module number from name
      // structure of name:
      // PMT_XX_module_YY
      // so
      // std::string nameOfPMT =  (std::string) thePostPVName;
      // // YY
      // std::string str_module_number = nameOfPMT.substr(nameOfPMT.find_last_of("_")+1,100);
      // int module_number = atoi(str_module_number.c_str());
      // // PMT_XX_module
      // std::string cutString = nameOfPMT.substr(0,nameOfPMT.find_last_of("_"));
      //
      // std::string str_pmt_number = cutString.substr(cutString.find_first_of("_") + 1,cutString.find_last_of("_") - cutString.find_first_of("_") - 1);
      // int pmt_number = atoi(str_pmt_number.c_str());


      if(saveAll || savePhotons)
      {

        // CreateTree::Instance()->module_number = module_number;
        // CreateTree::Instance()->front_back    = front_back;
        // CreateTree::Instance()->pmt_number    = pmt_number;
        CreateTree::Instance()->vertX         = theTrack->GetVertexPosition().x();
        CreateTree::Instance()->vertY         = theTrack->GetVertexPosition().y();
        CreateTree::Instance()->vertZ         = theTrack->GetVertexPosition().z();
        CreateTree::Instance()->vertMomentumX = theTrack->GetVertexMomentumDirection().x();
        CreateTree::Instance()->vertMomentumY = theTrack->GetVertexMomentumDirection().y();
        CreateTree::Instance()->vertMomentumZ = theTrack->GetVertexMomentumDirection().z();
        CreateTree::Instance()->PositionX     = OnDetectorPosition.getX();
        CreateTree::Instance()->PositionY     = OnDetectorPosition.getY();
        CreateTree::Instance()->PositionZ     = OnDetectorPosition.getZ();
        CreateTree::Instance()->PreMomentumX  = PreOnDetectorMomentum.getX();
        CreateTree::Instance()->PreMomentumY  = PreOnDetectorMomentum.getY();
        CreateTree::Instance()->PreMomentumZ  = PreOnDetectorMomentum.getZ();
        CreateTree::Instance()->PostMomentumX = PostOnDetectorMomentum.getX();
        CreateTree::Instance()->PostMomentumY = PostOnDetectorMomentum.getY();
        CreateTree::Instance()->PostMomentumZ = PostOnDetectorMomentum.getZ();
        CreateTree::Instance()->globalTime    = theTrack->GetGlobalTime()/CLHEP::ns;
        CreateTree::Instance()->localTime     = theTrack->GetLocalTime()/CLHEP::ns;
        CreateTree::Instance()->PhotonEnergy  = theTrack->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV;
        CreateTree::Instance()->processName    = processName;
        CreateTree::Instance()->photonDetected  = photonDetected;

        // CreateTree::Instance()->photon_cellID = (Int_t) cellID;
        // CreateTree::Instance()->photon_moduleID = (Int_t) moduleID;
        CreateTree::Instance()->photon_moduleType = (Int_t) module_type;
        // CreateTree::Instance()->photon_module_x = module_x;
        // CreateTree::Instance()->photon_module_y = module_y;
        // CreateTree::Instance()->photon_module_z = module_z;
        CreateTree::Instance()->photons->Fill();
      }

      // if(front_back == 0)
      // {
      //   CreateTree::Instance()->listOfTimestamps_front_back_0.push_back(theTrack->GetGlobalTime()/CLHEP::ns);
      // }
      // else
      // {
      //   CreateTree::Instance()->listOfTimestamps_front_back_1.push_back(theTrack->GetGlobalTime()/CLHEP::ns);
      // }

      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }

    // should be now done by stacking action, useless
    // if( (theTrack->GetLogicalVolumeAtVertex()->GetName().contains("crystal")) && (nStep == 1) && (processName == "Scintillation") )
    // {
    //   if( !propagateScintillation ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    // }
    //
    // if( (theTrack->GetLogicalVolumeAtVertex()->GetName().contains("crystal")) && (nStep == 1) && (processName == "Cerenkov") )
    // {
    //   if( !propagateCerenkov ) theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    // }


    if((theTrack->GetTrackStatus() == fStopAndKill) || (theTrack->GetTrackStatus() == fKillTrackAndSecondaries))
    {
      if(saveAll || savePhotonAbsPoint)
      {
        // save death point in space
        CreateTree::Instance()->abs_x             = OnDetectorPosition.getX();
        CreateTree::Instance()->abs_y             = OnDetectorPosition.getY();
        CreateTree::Instance()->abs_z             = OnDetectorPosition.getZ();
        CreateTree::Instance()->abs_globalTime    = theTrack->GetGlobalTime()/CLHEP::ns;
        CreateTree::Instance()->abs_PhotonEnergy  = theTrack->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV;
        CreateTree::Instance()->photonsAbsPoint->Fill();
      }
    }
  } // optical photon
  else// non optical photon
  {
    if(theTrack->GetTrackID() == 1) // primary particle (not optical photon)
    {
      if(saveAll || savePrimaries)
      {
        //check if it entered the absborber already (in that case it's been already saved)
        if(  !CreateTree::Instance()->primaryFirstEntranceFound) // not yet found
        {
          if(thePrePVName.contains("abs") || thePrePVName.contains("Abs") ) // enteting now
          {
            CreateTree::Instance()->primaryFirstEntranceFound = true;
            //save primary vertex in TTree
            G4ThreeVector particleVertexPosition = theTrack->GetVertexPosition ();
            G4ThreeVector particleVertexMomentum = theTrack->GetVertexMomentumDirection ();
            G4double      particleVertexEnergy   = theTrack->GetVertexKineticEnergy ();

            G4ThreeVector particlePosition = thePrePoint->GetPosition();
            G4ThreeVector particleMomentum = thePrePoint->GetMomentumDirection();
            G4double      particleEnergy   = theTrack->GetTotalEnergy ();
            CreateTree::Instance()->primaryPositionAtVertexX = particleVertexPosition.getX();
            CreateTree::Instance()->primaryPositionAtVertexY = particleVertexPosition.getY();
            CreateTree::Instance()->primaryPositionAtVertexZ = particleVertexPosition.getZ();
            CreateTree::Instance()->primaryMomentumAtVertexX = particleVertexMomentum.getX();
            CreateTree::Instance()->primaryMomentumAtVertexY = particleVertexMomentum.getY();
            CreateTree::Instance()->primaryMomentumAtVertexZ = particleVertexMomentum.getZ();
            CreateTree::Instance()->primaryEnergyAtVertex    = particleVertexEnergy;

            CreateTree::Instance()->primaryPositionOnAbsorberX = particlePosition.getX();
            CreateTree::Instance()->primaryPositionOnAbsorberY = particlePosition.getY();
            CreateTree::Instance()->primaryPositionOnAbsorberZ = particlePosition.getZ();
            CreateTree::Instance()->primaryPositionOnAbsorberT = theTrack->GetGlobalTime()/CLHEP::ns;
            CreateTree::Instance()->primaryMomentumOnAbsorberX = particleMomentum.getX();
            CreateTree::Instance()->primaryMomentumOnAbsorberY = particleMomentum.getY();
            CreateTree::Instance()->primaryMomentumOnAbsorberZ = particleMomentum.getZ();
            CreateTree::Instance()->primaryEnergyOnAbsorber    = particleEnergy;
            CreateTree::Instance()->primaries->Fill();

            // CreateTree::Instance()->primary_type = theTrack->GetDynamicParticle()->GetPDGcode();
          }
        }
      }
    }
    // G4cout << ">>> begin non optical photon" << G4endl;

    // find process
    // boundary process
    // G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    // static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

    //find the boundary process only once
    // if(!boundary)
    // {
    // G4double energy = theStep->GetTotalEnergyDeposit()/MeV;
    G4double IonEnergy = (theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit())/MeV;
    G4double NonIonEnergy = theStep->GetNonIonizingEnergyDeposit()/MeV;

	G4String processName = "none"; 
    if ( thePrePoint->GetProcessDefinedStep() ) processName = thePrePoint->GetProcessDefinedStep ()->GetProcessName();

    // G4cout << "Event = " << CreateTree::Instance ()->Event << G4endl;
    //if ( IonEnergy == 0. ) return ;
    // G4cout << "IonEnergy = " << IonEnergy << G4endl;

    // search for material number. bad but faster than looping on the std::map...
    int materialNumber = 0 ;
    if      ( materialName == "Quartz" )     materialNumber = 1;
    else if ( materialName == "SiO2Ce" )     materialNumber = 2;
    else if ( materialName == "DSB" )        materialNumber = 3;
    else if ( materialName == "LuAG_Ce" )    materialNumber = 4;
    else if ( materialName == "YAG_Ce" )     materialNumber = 5;
    else if ( materialName == "GAGG_Ce_Mg" ) materialNumber = 6;
    else if ( materialName == "Water" )      materialNumber = 7;
    else if ( materialName == "GFAG" )           materialNumber = 8;
    else if ( materialName == "GAGG_very_fast" ) materialNumber = 9;
    else if ( materialName == "GYAGG" )          materialNumber = 10;
    else if ( materialName == "GAGG_slow" )      materialNumber = 11;
    else if ( materialName == "Polystyrene" )      materialNumber = 12;
    else if ( materialName == "Lead_crystal" )      materialNumber = 13;
    else if ( materialName == "PureTungsten2_crystal" )      materialNumber = 14;
    else if ( materialName == "Al_crystal" )      materialNumber = 15;
    else if ( materialName == "Cu_crystal" )      materialNumber = 16;
    else if ( materialName == "GAGG_ILM" )        materialNumber = 17;

    CreateTree::Instance ()->depositedEnergyTotal += energy/GeV;

    if(thePrePVName.contains("crystal"))
    {
        showerIsInCrystal = 1;
    }
    else
    {
        showerIsInCrystal = 0;
    }

    // if(showerIsInCrystal)
    // {
      for(int ii=0;ii<preHistory->GetHistoryDepth();ii++)
      {
        volume_name = preHistory->GetVolume(ii)->GetName();
        if(volume_name.contains("Module"))
        {
          // find module type
          modBegin = volume_name.find_first_of("_");
          modString = volume_name.substr(modBegin+1,volume_name.size());
          std::istringstream ( modString ) >> module_type;
        }
      }
      // end of time test
      CreateTree::Instance()->total_energy_deposited += energy;
    // }


    if(saveAll || saveShower)
    {
      CreateTree::Instance()->showerIsInCrystal      = (Int_t) showerIsInCrystal;
      // CreateTree::Instance()->crystalID              = (Int_t) crystalID;
      // CreateTree::Instance()->shower_cellID          = (Int_t) cellID;
      // CreateTree::Instance()->shower_moduleID        = (Int_t) moduleID;
      CreateTree::Instance()->shower_moduleType      = (Int_t) module_type;
      // CreateTree::Instance()->shower_module_x        = module_x;
      // CreateTree::Instance()->shower_module_y        = module_y;
      // CreateTree::Instance()->shower_module_z        = module_z;
      CreateTree::Instance()->pdgID                  = (Int_t) theTrack->GetParticleDefinition()->GetPDGEncoding();
      CreateTree::Instance()->trackID                  = (Int_t) theTrack->GetTrackID();
      CreateTree::Instance()->primaryID              = theTrackInfo->GetPrimaryID();
      CreateTree::Instance()->primaryPDGID           = theTrackInfo->GetPrimaryPDGID();
      CreateTree::Instance()->primaryEnergy          = theTrackInfo->GetPrimaryEnergy();
      CreateTree::Instance()->showerX                = positionVector.getX();
      CreateTree::Instance()->showerY                = positionVector.getY();
      CreateTree::Instance()->showerZ                = positionVector.getZ();
      CreateTree::Instance()->showerPx               = momentumVector.getX();
      CreateTree::Instance()->showerPy               = momentumVector.getY();
      CreateTree::Instance()->showerPz               = momentumVector.getZ();
      CreateTree::Instance()->showerT                = theTrack->GetGlobalTime()/CLHEP::ns;;
      CreateTree::Instance()->showerTotalEnDep       = energy;
      CreateTree::Instance()->showerIonizingEnDep    = IonEnergy;
      CreateTree::Instance()->showerNonIonizingEnDep = NonIonEnergy;
      CreateTree::Instance()->showerProcessName      = processName;
      CreateTree::Instance()->showerMaterialName     = materialName;
      CreateTree::Instance()->showerMaterialNamePre  = materialNamePre;
      CreateTree::Instance()->showerMaterialNamePost = materialNamePost;
      CreateTree::Instance()->showerMaterialNumber   = materialNumber;
      CreateTree::Instance()->showerLocalX           = localPosition.getX();
      CreateTree::Instance()->showerLocalY           = localPosition.getY();
      CreateTree::Instance()->showerLocalZ           = localPosition.getZ();
      CreateTree::Instance()->shower->Fill();
    }


    //cout << "Process: " << processName.data() << ", material: " << materialName.data() << " particle ID: " << theTrack->GetParticleDefinition()->GetPDGEncoding() << " z = " << positionVector.getZ() << endl;

    if(strcmp(processName.data(),"Transportation")!=0) return;
    if(strcmp(materialName.data(),"LAPPD_Window")!=0) return;

    //cout << "Process: " << processName.data() << ", material: " << materialName.data() << " particle ID: " << theTrack->GetParticleDefinition()->GetPDGEncoding() << " z = " << positionVector.getZ() << endl;

    //if(abs(positionVector.getZ() - Parameters::Instance()->main_config.cell_separator_position + 0.5*Parameters::Instance()->main_config.separation_thickness)>0.01) return;

    if(abs(theTrack->GetParticleDefinition()->GetPDGEncoding())!=11) return;

    //cout << "Got entry, saving..." << endl;

    if(saveAll || saveLAPPD)
    {
           
      CreateTree::Instance()->pdgID                  = theTrack->GetParticleDefinition()->GetPDGEncoding();
      CreateTree::Instance()->trackID                = (Int_t) theTrack->GetTrackID();
      CreateTree::Instance()->showerX                = positionVector.getX();
      CreateTree::Instance()->showerY                = positionVector.getY();
      CreateTree::Instance()->showerZ                = positionVector.getZ();
      CreateTree::Instance()->showerPx               = momentumVector.getX();
      CreateTree::Instance()->showerPy               = momentumVector.getY();
      CreateTree::Instance()->showerPz               = momentumVector.getZ();
      CreateTree::Instance()->showerT                = theTrack->GetGlobalTime()/CLHEP::ns;
      CreateTree::Instance()->showerProcessName      = processName;
      CreateTree::Instance()->LAPPD->Fill();
      //cout << "Entry saved!" << endl;

    }


  } // non optical photon


  return ;
}
