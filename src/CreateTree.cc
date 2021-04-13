// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#include "CreateTree.hh"
#include <algorithm>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if ( fInstance )
  {
    return ;
  }

  this->fInstance = this ;
  this->fname     = name ;

  // saveAll            = pass_saveAll           ;
  // saveTree           = pass_saveTree          ;
  // saveShower         = pass_saveShower        ;
  // saveStructure      = pass_saveStructure     ;
  // savePrimaries      = pass_savePrimaries     ;
  // savePhotons        = pass_savePhotons       ;
  // savePhotonAbsPoint = pass_savePhotonAbsPoint;
  // saveSummary        = pass_saveSummary;
  // saveSimple         = pass_saveSimple;
  // savePhotonGen      = pass_savePhotonGen;
  // saveEnergyPerModule = pass_saveEnergyPerModule;
  // opticalCalibrationRun = pass_opticalCalibrationRun;

  // std::cout << "----------------------------------------" << std::endl;
  // std::cout << "----------------------------------------" << std::endl;
  // std::cout << saveAll              << " "
  //           << saveTree             << " "
  //           << saveShower           << " "
  //           << saveStructure        << " "
  //           << savePrimaries        << " "
  //           << savePhotons          << " "
  //           << savePhotonAbsPoint   << " "
  //           << saveSummary          << " "
  //           << saveSimple           << " "
  //           << savePhotonGen << " "
  //           << std::endl;
  // std::cout << "----------------------------------------" << std::endl;
  // std::cout << "----------------------------------------" << std::endl;
  //----------------------------------//
  // STANDARD TREE                    //
  //----------------------------------//
  //legacy from old sim, will be removed later
  this->ftree = NULL;
  energyPerCrystal = new vector<float> (1800, 0.) ;
  inputMomentum = new vector<double> (4, 0.) ;
  inputInitialPosition = new vector<double> (3, 0.) ;
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveTree)
  {
    this->ftree     = new TTree (name,name) ;

    this->GetTree ()->Branch ("Event", &this->Event, "Event/I") ;

    this->GetTree ()->Branch ("inputMomentum","vector<double>",&inputMomentum) ;
    this->GetTree ()->Branch ("inputInitialPosition","vector<double>",&inputInitialPosition) ;
    this->GetTree ()->Branch ("depositedEnergyTotal",        &this->depositedEnergyTotal,                "depositedEnergyTotal/F") ;
    this->GetTree ()->Branch ("Radial_stepLength",               &Radial_stepLength,                                     "Radial_stepLength/F");
    this->GetTree ()->Branch ("Longitudinal_stepLength",         &Longitudinal_stepLength,                         "Longitudinal_stepLength/F");
    this->GetTree ()->Branch ("Radial_ion_energy_absorber",       Radial_ion_energy_absorber,             "Radial_ion_energy_absorber[5000]/F");
    this->GetTree ()->Branch ("Longitudinal_ion_energy_absorber", Longitudinal_ion_energy_absorber, "Longitudinal_ion_energy_absorber[5000]/F");
    this->GetTree()->Branch("PrimaryParticleX",PrimaryParticleX,"PrimaryParticleX[1000]/F");
    this->GetTree()->Branch("PrimaryParticleY",PrimaryParticleY,"PrimaryParticleY[1000]/F");
    this->GetTree()->Branch("PrimaryParticleZ",PrimaryParticleZ,"PrimaryParticleZ[1000]/F");
    this->GetTree()->Branch("PrimaryParticleE",PrimaryParticleE,"PrimaryParticleE[1000]/F");
  }
  // end of STANDARD TREE
  //------------------------------------


  //----------------------------------//
  // STRUCTURE TREE                   //
  //----------------------------------//
  // saving the simulation structure into output files
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
  {
    fibres = new TTree ("fibres", "fibres") ;
    fibres->Branch("ID"          , &this->fibreID           , "ID/I");
    fibres->Branch("type"        , &this->fibreType         , "type/I");
    fibres->Branch("material"    , &this->fibreMaterial     , "material/I");
    fibres->Branch("x"           , &this->fibreX            , "x/F");
    fibres->Branch("y"           , &this->fibreY            , "y/F");
    fibres->Branch("z"           , &this->fibreZ            , "z/F");
    fibres->Branch("dx"          , &this->fibreDX           , "dx/F");
    fibres->Branch("dy"          , &this->fibreDY           , "dy/F");
    fibres->Branch("dz"          , &this->fibreDZ           , "dz/F");
    fibres->Branch("separation_z", &this->fibreSeparationZ  , "separation_z/F");
    fibres->Branch("sections"    , &this->fibreSections     , "sections/I");

    holes = new TTree ("holes", "holes") ;
    holes->Branch("ID"          , &this->holeID           , "ID/I");
    holes->Branch("type"        , &this->holeType         , "type/I");
    holes->Branch("material"    , &this->holeMaterial     , "material/I");
    holes->Branch("x"           , &this->holeX            , "x/F");
    holes->Branch("y"           , &this->holeY            , "y/F");
    holes->Branch("z"           , &this->holeZ            , "z/F");
    holes->Branch("dx"          , &this->holeDX           , "dx/F");
    holes->Branch("dy"          , &this->holeDY           , "dy/F");
    holes->Branch("dz"          , &this->holeDZ           , "dz/F");
    holes->Branch("separation_z", &this->holeSeparationZ  , "separation_z/F");
    holes->Branch("sections"    , &this->holeSections     , "sections/I");


    // cells->Branch("moduleID" , &this->moduleID     , "moduleID/I");
    cells = new TTree ("cells", "cells") ;
    cells->Branch("ID"          , &this->cellID           , "ID/I");
    cells->Branch("type"        , &this->cellType         , "type/I");
    cells->Branch("material"    , &this->cellMaterial     , "material/I");
    cells->Branch("x"           , &this->cellX            , "x/F");
    cells->Branch("y"           , &this->cellY            , "y/F");
    cells->Branch("z"           , &this->cellZ            , "z/F");
    cells->Branch("dx"          , &this->cellDX           , "dx/F");
    cells->Branch("dy"          , &this->cellDY           , "dy/F");
    cells->Branch("dz"          , &this->cellDZ           , "dz/F");
    cells->Branch("separation_z", &this->cellSeparationZ  , "separation_z/F");
    cells->Branch("sections"    , &this->cellSections     , "sections/I");

    absorbers = new TTree ("absorbers", "absorbers") ;
    absorbers->Branch("ID"          , &this->absorberID           , "ID/I");
    absorbers->Branch("type"        , &this->absorberType         , "type/I");
    absorbers->Branch("material"    , &this->absorberMaterial     , "material/I");
    absorbers->Branch("x"           , &this->absorberX            , "x/F");
    absorbers->Branch("y"           , &this->absorberY            , "y/F");
    absorbers->Branch("z"           , &this->absorberZ            , "z/F");
    absorbers->Branch("dx"          , &this->absorberDX           , "dx/F");
    absorbers->Branch("dy"          , &this->absorberDY           , "dy/F");
    absorbers->Branch("dz"          , &this->absorberDZ           , "dz/F");
    absorbers->Branch("separation_z", &this->absorberSeparationZ  , "separation_z/F");
    absorbers->Branch("sections"    , &this->absorberSections     , "sections/I");

    modules = new TTree ("modules", "modules") ;
    modules->Branch("ID"          , &this->moduleID           , "ID/I");
    modules->Branch("type"        , &this->moduleType         , "type/I");
    modules->Branch("material"    , &this->moduleMaterial     , "material/I");
    modules->Branch("x"           , &this->moduleX            , "x/F");
    modules->Branch("y"           , &this->moduleY            , "y/F");
    modules->Branch("z"           , &this->moduleZ            , "z/F");
    modules->Branch("dx"          , &this->moduleDX           , "dx/F");
    modules->Branch("dy"          , &this->moduleDY           , "dy/F");
    modules->Branch("dz"          , &this->moduleDZ           , "dz/F");
    modules->Branch("separation_z", &this->moduleSeparationZ  , "separation_z/F");
    modules->Branch("sections"    , &this->moduleSections     , "sections/I");

    calorimeter = new TTree ("calorimeter", "calorimeter") ;
    calorimeter->Branch("ID"          ,&this->calorimeterID          , "ID/I");
    calorimeter->Branch("type"        ,&this->calorimeterType        , "type/I");
    calorimeter->Branch("material"    ,&this->calorimeterMaterial    , "material/I");
    calorimeter->Branch("x"           ,&this->calorimeterX           , "x/F");
    calorimeter->Branch("y"           ,&this->calorimeterY           , "y/F");
    calorimeter->Branch("z"           ,&this->calorimeterZ           , "z/F");
    calorimeter->Branch("dx"          ,&this->calorimeterDX          , "dx/F");
    calorimeter->Branch("dy"          ,&this->calorimeterDY          , "dy/F");
    calorimeter->Branch("dz"          ,&this->calorimeterDZ          , "dz/F");
    calorimeter->Branch("separation_z",&this->calorimeterSeparationZ , "separation_z/F");
    calorimeter->Branch("sections"    ,&this->calorimeterSections    , "sections/I");
  }
  // end of STRUCTURE TREE
  //------------------------------------


  //----------------------------------//
  // SHOWER TREE                      //
  //----------------------------------//
  // saving the all info of energy depositions in the output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveShower)
  {
    shower = new TTree ("shower" , "shower") ;
    shower->Branch("run"                , &this->Run                    , "run/I");
    shower->Branch("event"              , &this->Event                  , "event/I");
    shower->Branch("isInCrystal"        , &this->showerIsInCrystal      , "isInCrystal/I");
    // shower->Branch("crystalID"          , &this->crystalID              , "crystalID/I");
    // shower->Branch("cellID"             , &this->shower_cellID , "cellID/I");
    // shower->Branch("sectionID"             , &this->shower_moduleID , "sectionID/I");
    // shower->Branch("moduleID"             , &this->shower_moduleID , "moduleID/I");
    shower->Branch("moduleType"             , &this->shower_moduleType , "moduleType/I");
    // shower->Branch("module_x"             , &this->shower_module_x , "module_x/F");
    // shower->Branch("module_y"             , &this->shower_module_y , "module_y/F");
    // shower->Branch("module_z"             , &this->shower_module_z , "module_z/F");
    shower->Branch("pdgID"        , &this->pdgID, "pdgID/I");
    shower->Branch("trackID"        , &this->trackID, "trackID/I");
    shower->Branch("primaryID"          , &this->primaryID              , "primaryID/I");
    shower->Branch("primaryPDGID"          , &this->primaryPDGID              , "primaryPDGID/I");
    shower->Branch("primaryEnergy"      , &this->primaryEnergy          , "primaryEnergy/F");
    shower->Branch("x"                  , &this->showerX                , "x/F");
    shower->Branch("y"                  , &this->showerY                , "y/F");
    shower->Branch("z"                  , &this->showerZ                , "z/F");
    shower->Branch("t"                  , &this->showerT                , "t/F");
    shower->Branch("px"                 , &this->showerPx               , "px/F");
    shower->Branch("py"                 , &this->showerPy               , "py/F");
    shower->Branch("pz"                 , &this->showerPz               , "pz/F");
    shower->Branch("totalEnDep"         , &this->showerTotalEnDep       , "totalEnDep/F");
    shower->Branch("ionizingEnDep"      , &this->showerIonizingEnDep    , "ionizingEnDep/F");
    shower->Branch("nonIonizingEnDep"   , &this->showerNonIonizingEnDep , "nonIonizingEnDep/F");
    shower->Branch("processName"        , &this->showerProcessName);
    shower->Branch("materialName"       , &this->showerMaterialName);
    shower->Branch("materialNamePre"    , &this->showerMaterialNamePre);
    shower->Branch("materialNamePost"   , &this->showerMaterialNamePost);
    shower->Branch("materialNumber"     , &this->showerMaterialNumber,"materialNumber/I");
    shower->Branch("localX"                  , &this->showerLocalX                , "localX/F");
    shower->Branch("localY"                  , &this->showerLocalY                , "localY/F");
    shower->Branch("localZ"                  , &this->showerLocalZ                , "localZ/F");
    shower->Branch("primary_PositionOnAbsorberX", &this->primaryPositionOnAbsorberX,"primary_PositionOnAbsorberX/F");
    shower->Branch("primary_PositionOnAbsorberY", &this->primaryPositionOnAbsorberY,"primary_PositionOnAbsorberY/F");
    shower->Branch("primary_PositionOnAbsorberZ", &this->primaryPositionOnAbsorberZ,"primary_PositionOnAbsorberZ/F");
    shower->Branch("primary_MomentumOnAbsorberX", &this->primaryMomentumOnAbsorberX,"primary_MomentumOnAbsorberX/F");
    shower->Branch("primary_MomentumOnAbsorberY", &this->primaryMomentumOnAbsorberY,"primary_MomentumOnAbsorberY/F");
    shower->Branch("primary_MomentumOnAbsorberZ", &this->primaryMomentumOnAbsorberZ,"primary_MomentumOnAbsorberZ/F");
    shower->Branch("primary_EnergyOnAbsorber"   , &this->primaryEnergyOnAbsorber   ,"primary_EnergyOnAbsorber/F");
  }
  // end of SHOWER TREE
  //------------------------------------

  //----------------------------------//
  // LAPPD TREE                      //
  //----------------------------------//
  // saving the all info of particle position and momentum in the output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveLAPPD)
  {
    LAPPD = new TTree ("LAPPD" , "LAPPD") ;
    LAPPD->Branch("pdgID"              , &this->pdgID                  , "pdgID/I");
    LAPPD->Branch("trackID"            , &this->trackID                , "trackID/I");
    LAPPD->Branch("x"                  , &this->showerX                , "x/F");
    LAPPD->Branch("y"                  , &this->showerY                , "y/F");
    LAPPD->Branch("z"                  , &this->showerZ                , "z/F");
    LAPPD->Branch("t"                  , &this->showerT                , "t/F");
    LAPPD->Branch("px"                 , &this->showerPx               , "px/F");
    LAPPD->Branch("py"                 , &this->showerPy               , "py/F");
    LAPPD->Branch("pz"                 , &this->showerPz               , "pz/F");
    LAPPD->Branch("primary_PositionOnAbsorberX", &this->primaryPositionOnAbsorberX,"primary_PositionOnAbsorberX/F");
    LAPPD->Branch("primary_PositionOnAbsorberY", &this->primaryPositionOnAbsorberY,"primary_PositionOnAbsorberY/F");
    LAPPD->Branch("primary_PositionOnAbsorberZ", &this->primaryPositionOnAbsorberZ,"primary_PositionOnAbsorberZ/F");
    LAPPD->Branch("primary_PositionOnAbsorberT", &this->primaryPositionOnAbsorberT,"primary_PositionOnAbsorberT/F");
    LAPPD->Branch("processName"        , &this->showerProcessName);
  }
  // end of LAPPD TREE
  //------------------------------------


  //----------------------------------//
  // PRIMARIES TREE                   //
  //----------------------------------//
  // saving the all info of primary particles in output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePrimaries || Parameters::Instance()->main_config.saveSummary)
  {
    primaries = new TTree ("primaries" , "primaries") ;
    primaries->Branch("run"                , &this->Run                       , "run/I");
    primaries->Branch("event"              , &this->Event                     , "event/I");
    primaries->Branch("PositionAtVertexX"  , &this->primaryPositionAtVertexX  ,"PositionAtVertexX/F");
    primaries->Branch("PositionAtVertexY"  , &this->primaryPositionAtVertexY  ,"PositionAtVertexY/F");
    primaries->Branch("PositionAtVertexZ"  , &this->primaryPositionAtVertexZ  ,"PositionAtVertexZ/F");
    primaries->Branch("MomentumAtVertexX"  , &this->primaryMomentumAtVertexX  ,"MomentumAtVertexX/F");
    primaries->Branch("MomentumAtVertexY"  , &this->primaryMomentumAtVertexY  ,"MomentumAtVertexY/F");
    primaries->Branch("MomentumAtVertexZ"  , &this->primaryMomentumAtVertexZ  ,"MomentumAtVertexZ/F");
    primaries->Branch("EnergyAtVertex"     , &this->primaryEnergyAtVertex     ,"EnergyAtVertex/F");
    primaries->Branch("PositionOnAbsorberX", &this->primaryPositionOnAbsorberX,"PositionOnAbsorberX/F");
    primaries->Branch("PositionOnAbsorberY", &this->primaryPositionOnAbsorberY,"PositionOnAbsorberY/F");
    primaries->Branch("PositionOnAbsorberZ", &this->primaryPositionOnAbsorberZ,"PositionOnAbsorberZ/F");
    primaries->Branch("MomentumOnAbsorberX", &this->primaryMomentumOnAbsorberX,"MomentumOnAbsorberX/F");
    primaries->Branch("MomentumOnAbsorberY", &this->primaryMomentumOnAbsorberY,"MomentumOnAbsorberY/F");
    primaries->Branch("MomentumOnAbsorberZ", &this->primaryMomentumOnAbsorberZ,"MomentumOnAbsorberZ/F");
    primaries->Branch("EnergyOnAbsorber"   , &this->primaryEnergyOnAbsorber   ,"EnergyOnAbsorber/F");
  }
  // end of PRIMARIES TREE
  //------------------------------------

  // if(opticalCalibrationRun)
  // {
  //   optCaliInfo = new TTree("optCaliInfo","optCaliInfo");
  //   optCaliInfo->Branch("run"           , &this->Run           , "run/I");
  //   optCaliInfo->Branch("events_per_run"         , &this->events_per_run         , "events_per_run/I");
  // }
  //----------------------------------//
  // PHOTONS TREE                     //
  //----------------------------------//
  // saving the all info of photon detection in output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotons)
  {
    // branch only what's needed
    photons = new TTree ("photons", "photons");
    photons->Branch("run"           , &this->Run           , "run/I");
    photons->Branch("event"         , &this->Event         , "event/I");
    // if(Parameters::Instance()->main_config.simulationType != 3)
    // {
    //   photons->Branch("front_back"    , &this->front_back    , "front_back/I");
    //   photons->Branch("pmt_number"    , &this->pmt_number    , "pmt_number/I");
    //   photons->Branch("module_number" , &this->module_number , "module_number/I");
    // }
    photons->Branch("vertX"         , &this->vertX         , "vertX/F");
    photons->Branch("vertY"         , &this->vertY         , "vertY/F");
    photons->Branch("vertZ"         , &this->vertZ         , "vertZ/F");
    if(Parameters::Instance()->main_config.simulationType != 3)
    {
      photons->Branch("vertMomentumX" , &this->vertMomentumX , "vertMomentumX/F");
      photons->Branch("vertMomentumY" , &this->vertMomentumY , "vertMomentumY/F");
      photons->Branch("vertMomentumZ" , &this->vertMomentumZ , "vertMomentumZ/F");
    }
    photons->Branch("PositionX"     , &this->PositionX     , "PositionX/F");
    photons->Branch("PositionY"     , &this->PositionY     , "PositionY/F");
    photons->Branch("PositionZ"     , &this->PositionZ     , "PositionZ/F");
    if(Parameters::Instance()->main_config.simulationType != 3)
    {
      photons->Branch("PreMomentumX"  , &this->PreMomentumX  , "PreMomentumX/F");
      photons->Branch("PreMomentumY"  , &this->PreMomentumY  , "PreMomentumY/F");
      photons->Branch("PreMomentumZ"  , &this->PreMomentumZ  , "PreMomentumZ/F");
      photons->Branch("PostMomentumX" , &this->PostMomentumX , "PostMomentumX/F");
      photons->Branch("PostMomentumY" , &this->PostMomentumY , "PostMomentumY/F");
      photons->Branch("PostMomentumZ" , &this->PostMomentumZ , "PostMomentumZ/F");
    }
    photons->Branch("globalTime"    , &this->globalTime    , "globalTime/F");
    if(Parameters::Instance()->main_config.simulationType != 3)
    {
      photons->Branch("localTime"     , &this->localTime     , "localTime/F");
    }
    photons->Branch("PhotonEnergy"  , &this->PhotonEnergy  ,      "PhotonEnergy/F");
    // photons->Branch("cellID"        , &this->photon_cellID ,     "cellID/I");
    // photons->Branch("moduleID"      , &this->photon_moduleID ,   "moduleID/I");
    photons->Branch("moduleType"    , &this->photon_moduleType , "moduleType/I");
    // photons->Branch("module_x"      , &this->photon_module_x ,   "module_x/F");
    // photons->Branch("module_y"      , &this->photon_module_y ,   "module_y/F");
    // photons->Branch("module_z"      , &this->photon_module_z ,   "module_z/F");
    if(Parameters::Instance()->main_config.simulationType != 3)
    {
      photons->Branch("processName"   ,  &this->processName );
      // photons->Branch("detected"      , &this->photonDetected           , "photonDetected/I");
      photons->Branch("primary_PositionOnAbsorberX", &this->primaryPositionOnAbsorberX,"primary_PositionOnAbsorberX/F");
      photons->Branch("primary_PositionOnAbsorberY", &this->primaryPositionOnAbsorberY,"primary_PositionOnAbsorberY/F");
      photons->Branch("primary_PositionOnAbsorberZ", &this->primaryPositionOnAbsorberZ,"primary_PositionOnAbsorberZ/F");
      photons->Branch("primary_MomentumOnAbsorberX", &this->primaryMomentumOnAbsorberX,"primary_MomentumOnAbsorberX/F");
      photons->Branch("primary_MomentumOnAbsorberY", &this->primaryMomentumOnAbsorberY,"primary_MomentumOnAbsorberY/F");
      photons->Branch("primary_MomentumOnAbsorberZ", &this->primaryMomentumOnAbsorberZ,"primary_MomentumOnAbsorberZ/F");
      photons->Branch("primary_EnergyOnAbsorber"   , &this->primaryEnergyOnAbsorber   ,"primary_EnergyOnAbsorber/F");
    }
  }
  // end of PHOTONS TREE
  //------------------------------------


  //----------------------------------//
  // PHOTONS ABS TREE                 //
  //----------------------------------//
  // saving the all info of photon absorption point in output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotonAbsPoint)
  {
    photonsAbsPoint = new TTree ("photonsAbsPoint", "photonsAbsPoint");
    photonsAbsPoint->Branch("run"           , &this->Run          , "run/I");
    photonsAbsPoint->Branch("event"         , &this->Event   , "event/I");
    photonsAbsPoint->Branch("x"             , &this->abs_x             , "x/F");
    photonsAbsPoint->Branch("y"             , &this->abs_y             , "y/F");
    photonsAbsPoint->Branch("z"             , &this->abs_z             , "z/F");
    photonsAbsPoint->Branch("globalTime"    , &this->abs_globalTime    , "globalTime/F");
    photonsAbsPoint->Branch("PhotonEnergy"  , &this->abs_PhotonEnergy  , "PhotonEnergy/F");
  }
  // end of PHOTONS ABS TREE
  //------------------------------------


  //----------------------------------//
  // SUMMARY OUTPUT TREE              //
  //----------------------------------//
  // saving summary output in output file
  summary = NULL;
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveSummary)
  {
    summary = new TTree ("summary", "summary");
    summary->Branch("run"                    , &this->Run                       ,"run/I");
    summary->Branch("event"                  , &this->Event                     ,"event/I");
    summary->Branch("PositionAtVertexX"      , &this->inputInitialPosition->at(0)  ,"PositionAtVertexX/F");
    summary->Branch("PositionAtVertexY"      , &this->inputInitialPosition->at(1)  ,"PositionAtVertexY/F");
    summary->Branch("PositionAtVertexZ"      , &this->inputInitialPosition->at(2)  ,"PositionAtVertexZ/F");
    summary->Branch("MomentumAtVertexX"      , &this->inputMomentum->at(0)  ,"MomentumAtVertexX/F");
    summary->Branch("MomentumAtVertexY"      , &this->inputMomentum->at(1)  ,"MomentumAtVertexY/F");
    summary->Branch("MomentumAtVertexZ"      , &this->inputMomentum->at(2)  ,"MomentumAtVertexZ/F");
    summary->Branch("EnergyAtVertex"         , &this->inputMomentum->at(3)     ,"EnergyAtVertex/F");
    summary->Branch("primaryType"            , &this->primary_type              ,"primaryType/I");
    summary->Branch("totalEnergyDeposited"   , &this->total_energy_deposited    ,"totalEnergyDeposited/F");
  }
  // end of SUMMARY OUTPUT TREE
  //------------------------------------

  // listOfTimestamps_front_back_0 = new std::vector<double>;
  // listOfTimestamps_front_back_1 = new std::vector<double>;

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveSimple)
  {
    simple_time = new TTree ("simple_time", "simple_time");
    simple_time->Branch("run"                    , &this->Run                       ,"run/I");
    simple_time->Branch("event"                  , &this->Event                     ,"event/I");
    simple_time->Branch("timestamp_front_back_0" , &this->timestamp_front_back_0    ,"timestamp_front_back_0/F");
    simple_time->Branch("timestamp_front_back_1" , &this->timestamp_front_back_1    ,"timestamp_front_back_1/F");
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotonGen)
  {
    opticalPhotonGen = new TTree ("opticalPhotonGen", "opticalPhotonGen");
    opticalPhotonGen->Branch("run"           , &this->Run            ,"run/I");
    opticalPhotonGen->Branch("event"         , &this->Event          ,"event/I");
    opticalPhotonGen->Branch("x"             , &this->gen_x          ,"x/F");
    opticalPhotonGen->Branch("y"             , &this->gen_y          ,"y/F");
    opticalPhotonGen->Branch("z"             , &this->gen_z          ,"z/F");
    opticalPhotonGen->Branch("MomentumX"     , &this->gen_ux         ,"MomwntumX/F");
    opticalPhotonGen->Branch("MomentumY"     , &this->gen_uy         ,"MomwntumY/F");
    opticalPhotonGen->Branch("MomentumZ"     , &this->gen_uz         ,"MomwntumZ/F");
    opticalPhotonGen->Branch("PhotonEnergy"  , &this->gen_energy     ,"PhotonEnergy/F");
    opticalPhotonGen->Branch("globalTime"    , &this->gen_globalTime ,"globalTime/F");
    opticalPhotonGen->Branch("localTime"     , &this->gen_localTime  ,"localTime/F");
    opticalPhotonGen->Branch("processName"  ,  &this->gen_process );
    opticalPhotonGen->Branch("primary_PositionOnAbsorberX", &this->primaryPositionOnAbsorberX,"primary_PositionOnAbsorberX/F");
    opticalPhotonGen->Branch("primary_PositionOnAbsorberY", &this->primaryPositionOnAbsorberY,"primary_PositionOnAbsorberY/F");
    opticalPhotonGen->Branch("primary_PositionOnAbsorberZ", &this->primaryPositionOnAbsorberZ,"primary_PositionOnAbsorberZ/F");
    opticalPhotonGen->Branch("primary_MomentumOnAbsorberX", &this->primaryMomentumOnAbsorberX,"primary_MomentumOnAbsorberX/F");
    opticalPhotonGen->Branch("primary_MomentumOnAbsorberY", &this->primaryMomentumOnAbsorberY,"primary_MomentumOnAbsorberY/F");
    opticalPhotonGen->Branch("primary_MomentumOnAbsorberZ", &this->primaryMomentumOnAbsorberZ,"primary_MomentumOnAbsorberZ/F");
    opticalPhotonGen->Branch("primary_EnergyOnAbsorber"   , &this->primaryEnergyOnAbsorber   ,"primary_EnergyOnAbsorber/F");

  }

  SaveOpticalRunVector = false;

  // YANXI MOD
  module_array_x        = Parameters::Instance()->main_config.modules_nx;
  module_array_y        = Parameters::Instance()->main_config.modules_ny;
  pipe_modules_nx       = Parameters::Instance()->main_config.pipe_modules_nx;
  pipe_modules_ny       = Parameters::Instance()->main_config.pipe_modules_ny;
  // CreateTree::Instance()->module_array_x = modules_nx;
  // CreateTree::Instance()->module_array_y = modules_ny;
  // CreateTree::Instance()->pipe_modules_nx = pipe_modules_nx;
  // CreateTree::Instance()->pipe_modules_ny = pipe_modules_ny;

  module_total = module_array_x*module_array_y;
  G4cout<<"Number of modules "<< module_array_x<<"x"<<module_array_y<<G4endl;
  // these numbers define the 4 or 1 regions of the calorimeter. We need to save this info
  // as usual, we use CreateTree as storing place
  if((pipe_modules_nx == 0) && (pipe_modules_ny == 0))
  {
    pipe = 0;
  }
  else
  {
    pipe = 1;
  }
  depositedEnergyByModuleFront  = new float[module_total];
  depositedEnergyByModuleRear   = new float[module_total];
  depositionTimingByModuleFront = new float[module_total];
  depositionTimingByModuleRear  = new float[module_total];
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  {
    CreateTreeEnergyDepositByModule();
  }

  // END OF YANXI MOD



  // inizialize by cleaning
  this->Clear () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::~CreateTree ()
{}



int CreateTree::Fill ()
{
  return this->GetTree ()->Fill () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


bool CreateTree::Write (TFile * outfile)
{
  outfile->cd () ;

  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.saveTree)           ftree->Write () ;
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.saveShower)         shower->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.savePrimaries)      primaries->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.savePhotons)        photons->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.savePhotonAbsPoint) photonsAbsPoint->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.saveLAPPD)          LAPPD->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.saveSummary)        summary->Write();
  if(Parameters::Instance()->main_config.saveSimple  || Parameters::Instance()->main_config.saveAll)            simple_time->Write();
  if(Parameters::Instance()->main_config.saveAll     || Parameters::Instance()->main_config.savePhotonGen)      opticalPhotonGen->Write();
  if(Parameters::Instance()->main_config.simulationType == 3)
  {
    gDirectory->WriteObject(&optN    , "optN");
    gDirectory->WriteObject(&optX    , "optX");
    gDirectory->WriteObject(&optY    , "optY");
    gDirectory->WriteObject(&optZ    , "optZ");
    gDirectory->WriteObject(&optE    , "optE");
  }
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  {
    EnergyPerModuleFront->Write();
    EnergyPerModuleRear->Write();
  }
  // FIXME for now always on
  // int opticalMaterial = 0;

  for(int i = 0 ; i < Parameters::Instance()->enHisto.size(); i++)
  {
    Parameters::Instance()->enHisto[i]->Write();
  }
  for(int i = 0 ; i < Parameters::Instance()->tHisto.size(); i++)
  {
    Parameters::Instance()->tHisto[i]->Write();
  }
  for(int i = 0 ; i < Parameters::Instance()->absHisto.size(); i++)
  {
    Parameters::Instance()->absHisto[i]->Write();
  }
  // if(abs_GaGG_Ce_Mg)
  // {
  //   abs_GaGG_Ce_Mg->Write();
  // }
  // if(abs_GaGG_Ce_Mg_old)
  // {
  //   abs_GaGG_Ce_Mg_old->Write();
  // }
  // if(abs_YAG_Ce)
  // {
  //   abs_YAG_Ce->Write();
  // }
  // also write a std::vector of crystal numbers and of crystal light yields
  std::vector<float> v_crystalLightYield;
  std::vector<float> v_crystalResolutionScale;
  // CrystalMaterialList
  // std::cout << "Found " << CrystalMaterialList.size() << " crystal materials "<< std::endl;
  for(unsigned int i = 0 ; i < Parameters::Instance()->CrystalMaterialList.size(); i++)
  {
    CrystalMaterialList.push_back (Parameters::Instance()->CrystalMaterialList[i] );
    v_crystalLightYield.push_back( (Parameters::Instance()->GetCrystalMaterial(Parameters::Instance()->CrystalMaterialList[i]))->GetMaterialPropertiesTable()->GetConstProperty("SCINTILLATIONYIELD") );
    v_crystalResolutionScale.push_back( (Parameters::Instance()->GetCrystalMaterial(Parameters::Instance()->CrystalMaterialList[i]))->GetMaterialPropertiesTable()->GetConstProperty("RESOLUTIONSCALE") );
  }
  gDirectory->WriteObject(&CrystalMaterialList, "crystalMaterialList");
  gDirectory->WriteObject(&v_crystalLightYield, "crystalLightYield");
  gDirectory->WriteObject(&v_crystalResolutionScale, "crystalResolutionScale");

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
  {
    calorimeter->Write();
    modules->Write();
    absorbers->Write();
    cells->Write();
    holes->Write();
    fibres->Write () ;
  }

  // if(opticalMaterial != 2) // so if it's NOT a instantaneous monocromatic source
  // {
  //   if(enHisto)
  //   {
  //     enHisto->Write();
  //   }
  //   if(tHisto)
  //   {
  //     tHisto->Write();
  //   }
  //   if(abs_GaGG_Ce_Mg)
  //   {
  //     abs_GaGG_Ce_Mg->Write();
  //   }
  //   if(abs_YAG_Ce)
  //   {
  //     abs_YAG_Ce->Write();
  //   }
  // }
  return true ;
}

// bool CreateTree::WriteStructure (TFile * outfile)
// {
//   outfile->cd () ;
//   modules->Write();
//   absorbers->Write();
//   cells->Write();
//   holes->Write();
//   fibres->Write () ;
// }
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::CreateTreeEnergyDepositByModule(){
    //Save energy deposit by module
    EnergyPerModuleFront = new TTree("EnergyPerModuleFront","EnergyPerModuleFront");
    EnergyPerModuleFront->Branch("module_total",&module_total,"module_total/I");
    EnergyPerModuleFront->Branch("de",this->depositedEnergyByModuleFront,"de[module_total]/F");
    EnergyPerModuleFront->Branch("te",this->depositionTimingByModuleFront,"te[module_total]/F");
    EnergyPerModuleRear = new TTree("EnergyPerModuleRear","EnergyPerModuleRear");
    EnergyPerModuleRear->Branch("module_total",&module_total,"module_total/I");
    EnergyPerModuleRear->Branch("de",this->depositedEnergyByModuleRear,"de[module_total]/F");
    EnergyPerModuleRear->Branch("te",this->depositionTimingByModuleRear,"te[module_total]/F");
}

// void CreateTree::AddCrystalMaterial(int num, G4Material* aMaterial)
// {
//   CrystalMaterialList.push_back(num);
//   CrystalMaterialMap.insert(std::make_pair(num,aMaterial));
// }
//
// G4Material* CreateTree::GetCrystalMaterial(int num)
// {
//   return CrystalMaterialMap[num];
// }



void CreateTree::Clear ()
{
  Event	= 0 ;

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  {
    for(int ii=0;ii<module_total;ii++)
    {
        depositedEnergyByModuleFront[ii] = 0.;
        depositedEnergyByModuleRear[ii] = 0.;
        depositionTimingByModuleFront[ii] = 0.;
        depositionTimingByModuleRear[ii] = 0.;
    }
  }




  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveTree)
  {
    depositedEnergyTotal = 0. ;
    for (int i = 0 ; i < 4 ; ++i)
    {
      inputMomentum->at (i) = 0. ;
    }
    for (int i = 0 ; i < 3 ; ++i)
    {
      inputInitialPosition->at (i) = 0. ;
    }

    Radial_stepLength = 0.;
    Longitudinal_stepLength = 0.;
    for(int i = 0; i < 5000; ++i)
    {
      Radial_ion_energy_absorber[i] = 0.;
      Longitudinal_ion_energy_absorber[i] = 0.;
    }

    for(int i = 0; i < 1000; ++i)
    {
      PrimaryParticleX[i] = 0.;
      PrimaryParticleY[i] = 0.;
      PrimaryParticleZ[i] = 0.;
      PrimaryParticleE[i] = 0.;
    }
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveShower)
  {
    showerIsInCrystal      = 0 ;
    // crystalID              = 0 ;
    // shower_cellID          = 0 ;
    // shower_moduleID        = 0 ;
    shower_moduleType      = 0 ;
    // shower_module_x        = 0 ;
    // shower_module_y        = 0 ;
    // shower_module_z        = 0 ;
    showerX                = 0 ;
    showerY                = 0 ;
    showerZ                = 0 ;
    showerT                = 0 ;
    showerPx               = 0 ;
    showerPy               = 0 ;
    showerPz               = 0 ;
    showerTotalEnDep       = 0 ;
    showerIonizingEnDep    = 0 ;
    showerNonIonizingEnDep = 0 ;
    showerProcessName      = "";
    showerMaterialName     = "";
    showerMaterialNamePre     = "";
    showerMaterialNamePost     = "";
    showerMaterialNumber   = 0;
    showerLocalX           = 0 ;
    showerLocalY           = 0 ;
    showerLocalZ           = 0 ;
    primaryPositionOnAbsorberX = 0  ;
    primaryPositionOnAbsorberY = 0  ;
    primaryPositionOnAbsorberZ = 0  ;
    primaryMomentumOnAbsorberX = 0  ;
    primaryMomentumOnAbsorberY = 0  ;
    primaryMomentumOnAbsorberZ = 0  ;
    primaryEnergyOnAbsorber    = 0  ;
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveLAPPD)
  {
    showerX                = 0 ;
    showerY                = 0 ;
    showerZ                = 0 ;
    showerT                = 0 ;
    showerPx               = 0 ;
    showerPy               = 0 ;
    showerPz               = 0 ;
    showerProcessName      = "";
    primaryPositionOnAbsorberX = 0  ;
    primaryPositionOnAbsorberY = 0  ;
    primaryPositionOnAbsorberZ = 0  ;
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotons)
  {
    // Nprimaries = 0;
    // front_back    = 0 ;
    // pmt_number    = 0 ;
    // module_number = 0 ;
    vertX         = 0 ;
    vertY         = 0 ;
    vertZ         = 0 ;
    vertMomentumX = 0 ;
    vertMomentumY = 0 ;
    vertMomentumZ = 0 ;
    PositionX     = 0 ;
    PositionY     = 0 ;
    PositionZ     = 0 ;
    PreMomentumX  = 0 ;
    PreMomentumY  = 0 ;
    PreMomentumZ  = 0 ;
    PostMomentumX = 0 ;
    PostMomentumY = 0 ;
    PostMomentumZ = 0 ;
    globalTime    = 0 ;
    localTime     = 0 ;
    PhotonEnergy  = 0 ;
    processName = "";
    // photon_cellID = 0;
    // photon_moduleID = 0;
    photon_moduleType = 0;
    // photon_module_x = 0;
    // photon_module_y = 0;
    // photon_module_z = 0;
    primaryPositionOnAbsorberX = 0  ;
    primaryPositionOnAbsorberY = 0  ;
    primaryPositionOnAbsorberZ = 0  ;
    primaryMomentumOnAbsorberX = 0  ;
    primaryMomentumOnAbsorberY = 0  ;
    primaryMomentumOnAbsorberZ = 0  ;
    primaryEnergyOnAbsorber    = 0  ;
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotonAbsPoint)
  {
    abs_x            = 0;
    abs_y            = 0;
    abs_z            = 0;
    abs_globalTime   = 0;
    abs_PhotonEnergy = 0;
  }


  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePrimaries)
  {
    primaryPositionAtVertexX   = 0  ;
    primaryPositionAtVertexY   = 0  ;
    primaryPositionAtVertexZ   = 0  ;
    primaryMomentumAtVertexX   = 0  ;
    primaryMomentumAtVertexY   = 0  ;
    primaryMomentumAtVertexZ   = 0  ;
    primaryEnergyAtVertex      = 0  ;
    primaryPositionOnAbsorberX = 0  ;
    primaryPositionOnAbsorberY = 0  ;
    primaryPositionOnAbsorberZ = 0  ;
    primaryMomentumOnAbsorberX = 0  ;
    primaryMomentumOnAbsorberY = 0  ;
    primaryMomentumOnAbsorberZ = 0  ;
    primaryEnergyOnAbsorber    = 0  ;
    primaryFirstEntranceFound = false;
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.savePhotonGen)
  {
    gen_x            = 0 ;
    gen_y            = 0 ;
    gen_z            = 0 ;
    gen_ux           = 0 ;
    gen_uy           = 0 ;
    gen_uz           = 0 ;
    gen_energy       = 0 ;
    gen_globalTime   = 0 ;
    gen_localTime    = 0 ;
    gen_process      = "";
    primaryPositionOnAbsorberX = 0  ;
    primaryPositionOnAbsorberY = 0  ;
    primaryPositionOnAbsorberZ = 0  ;
    primaryMomentumOnAbsorberX = 0  ;
    primaryMomentumOnAbsorberY = 0  ;
    primaryMomentumOnAbsorberZ = 0  ;
    primaryEnergyOnAbsorber    = 0  ;
  }

  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveSummary)
  {
    for (int i = 0 ; i < 4 ; ++i)
    {
      inputMomentum->at (i) = 0. ;
    }
    for (int i = 0 ; i < 3 ; ++i)
    {
      inputInitialPosition->at (i) = 0. ;
    }
    total_energy_deposited = 0;
    primary_type = 0;
  }

  if(Parameters::Instance()->main_config.saveSimple || Parameters::Instance()->main_config.saveAll)
  {
    listOfTimestamps_front_back_0.clear();
    listOfTimestamps_front_back_1.clear();
    timestamp_front_back_0 = 0;
    timestamp_front_back_1 = 0;
  }


}
