// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"

// #include "DetectorParameterisation.hh"
#include "CreateTree.hh"


#include <algorithm>
#include <string>
#include <sstream>
#include <limits>

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4Cons.hh>
#include <G4UnionSolid.hh>
#include <Cell.hh>
#include <Absorber.hh>
#include <ModuleParametrization.hh>
#include <utility>
#include "recallCAD.hh"


using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //-----------------------------------------------------
  //------------- Define colors --------------
  //-----------------------------------------------------
  white      = G4Colour(1.00, 1.00, 1.00) ;  // white
  grey       = G4Colour(0.50, 0.50, 0.50) ;  // grey
  black      = G4Colour(0.00, 0.00, 0.00) ;  // black
  red        = G4Colour(1.00, 0.00, 0.00) ;  // red
  green      = G4Colour(0.00, 1.00, 0.00) ;  // green
  palegreen  = G4Colour(130.0/255.0, 1.00, 130.0/255.0) ;  // pale green
  blue       = G4Colour(0.00, 0.00, 1.00) ;  // blue
  cyan       = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  air        = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  magenta    = G4Colour(1.00, 0.00, 1.00) ;  // magenta
  yellow     = G4Colour(1.00, 1.00, 0.00) ;  // yellow
  brass      = G4Colour(0.80, 0.60, 0.40) ;  // brass
  brown      = G4Colour(0.70, 0.40, 0.10) ;  // brown
  orange     = G4Colour(1.00, 0.33, 0.00) ;  // pale orange
  paleorange = G4Colour(1.00, 170.0/255.0, 0.00) ;  // orange

  color_vector.push_back(blue);
  color_vector.push_back(cyan);
  color_vector.push_back(green);
  color_vector.push_back(palegreen);
  color_vector.push_back(paleorange);
  color_vector.push_back(red);
  color_vector.push_back(black);
  color_vector.push_back(orange);
  color_vector.push_back(grey);
  color_vector.push_back(brown);


  G4cout << ">>>>>> TEST <<<<<<" << G4endl ;
  Parameters::Instance()->PrintConfig(Parameters::Instance()->main_config);
  for(int iC = 0 ; iC < Parameters::Instance()->module_config.size() ; iC++ )
  {
    Parameters::Instance()->PrintConfig(Parameters::Instance()->module_config[iC]);
  }
  G4cout << ">>>>>> END OF TEST <<<<<<" << G4endl ;

  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  // ConfigFile config (configFileName) ;

  saveAll            = Parameters::Instance()->main_config.saveAll;
  saveTree           = Parameters::Instance()->main_config.saveTree;
  saveShower         = Parameters::Instance()->main_config.saveShower;
  saveStructure      = Parameters::Instance()->main_config.saveStructure;
  savePrimaries      = Parameters::Instance()->main_config.savePrimaries;
  savePhotons        = Parameters::Instance()->main_config.savePhotons;
  savePhotonAbsPoint = Parameters::Instance()->main_config.savePhotons;
  saveLAPPD          = Parameters::Instance()->main_config.saveLAPPD;

  checkOverlaps      = Parameters::Instance()->main_config.checkOverlaps;
  world_material     = Parameters::Instance()->main_config.world_material;
  B_field_intensity  = Parameters::Instance()->main_config.B_field_intensity;

  //---------------------------------------//
  // CALORIMETER                           //
  //---------------------------------------//
  modules_nx           = Parameters::Instance()->main_config.modules_nx;
  modules_ny           = Parameters::Instance()->main_config.modules_ny;
  calorimeter_position = Parameters::Instance()->main_config.calorimeter_position;
  //visibility

  worldVisibility     = Parameters::Instance()->main_config.worldVisibility;
  caloVisibility      = Parameters::Instance()->main_config.caloVisibility;

  scaleFactor         = Parameters::Instance()->main_config.plex_abs_length_scale_factor;



  if(Parameters::Instance()->calorimeter_method == 1)
  {
    // create the one and only module here
    Spacal *Spacal_module = new Spacal(Parameters::Instance()->main_config);
    // store it in the det const modules list
    modules_list.push_back(Spacal_module);

    // we need to use it's dimensions to create the calo and the world
    calorimeter_size_x = modules_nx*Spacal_module->GetModuleSizeX();
    calorimeter_size_y = modules_ny*Spacal_module->GetModuleSizeY();
    calorimeter_size_z = Spacal_module->GetModuleSizeZ() + 2.0*fabs(Spacal_module->GetModuleZshift());

    module_size_x = Spacal_module->GetModuleSizeX();
    module_size_y = Spacal_module->GetModuleSizeY();
    module_size_z = Spacal_module->GetModuleSizeZ();

    cell_separator_position = Parameters::Instance()->main_config.cell_separator_position;
    pipe_modules_nx = Parameters::Instance()->main_config.pipe_modules_nx;
    pipe_modules_ny = Parameters::Instance()->main_config.pipe_modules_ny;

    // check if module_n - pipe_n is > 0
    if((modules_nx - pipe_modules_nx) <= 0 )
    {
      G4cerr << "ERROR: (modules_nx - pipe_modules_nx) cannot be 0 or less! Aborting..." << std::endl;
      exit(-1);
    }
    if((modules_ny - pipe_modules_ny) <= 0 )
    {
      G4cerr << "ERROR: (modules_ny - pipe_modules_ny) cannot be 0 or less! Aborting..." << std::endl;
      exit(-1);
    }
    // check if module_n - pipe_n is even
    if(pipe_modules_nx != 0)
    {
      if(((modules_nx - pipe_modules_nx) % 2) != 0)
      {
        G4cerr << "ERROR: The value of (modules_nx - pipe_modules_nx) cannot be odd! Aborting..." << std::endl;
        exit(-1);
      }
    }
    if(pipe_modules_ny != 0)
    {
      if(((modules_ny - pipe_modules_ny) % 2) != 0)
      {
        G4cerr << "ERROR: The value of (modules_ny - pipe_modules_ny) cannot be odd! Aborting..." << std::endl;
        exit(-1);
      }
    }
  }
  else
  {
    calorimeter_size_x = Parameters::Instance()->main_config.calorimeter_size_x;
    calorimeter_size_y = Parameters::Instance()->main_config.calorimeter_size_y;
    calorimeter_size_z = Parameters::Instance()->main_config.calorimeter_size_z;

    // and create as many module prototypes as needed
    for(int iC = 0 ; iC < Parameters::Instance()->module_config.size() ; iC++ )
    {
      // create the module
      Spacal *Spacal_module = new Spacal(Parameters::Instance()->module_config[iC]);
      // store it in the det const modules list
      modules_list.push_back(Spacal_module);
    }
  }

  expHall_x = 4.0*calorimeter_size_x;
  expHall_y = 4.0*calorimeter_size_y;
  expHall_z = 4.0*calorimeter_size_z;

  B_field_IsInitialized = false ;
  // set calorimeter position
  calorimeter_pos_x = 0;
  calorimeter_pos_y = 0;
  if(Parameters::Instance()->main_config.calorimeter_position == 0)
  {
    calorimeter_pos_z = 0;
  }
  else
  {
    calorimeter_pos_z = calorimeter_size_z/2.0;
  }
  // store absolute position of separation front/back, to be used in energy accumulation
  // worth saving in output too, to be used by propagateHybrid!
  CreateTree::Instance()->zSeparationAbsolutePosition = calorimeter_pos_z + cell_separator_position;
  std::cout << "zSeparationAbsolutePosition " << " " << CreateTree::Instance()->zSeparationAbsolutePosition << std::endl;

  containerVisibility = Parameters::Instance()->main_config.containerVisibility;





  initializeMaterials () ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


DetectorConstruction::~DetectorConstruction ()
{}


  //---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;




  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4cout << ">>>>>>>>>>>>> Build Experimental Hall "<< G4endl;
  G4VSolid * worldS = new G4Box ("worldS", 0.5 * expHall_x, 0.5 * expHall_y, 0.5 * expHall_z) ;
  G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
  G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (0.,0.,0.), worldLV, "World", 0, false, 0, checkOverlaps) ;
  G4VisAttributes* VisAttWorld = new G4VisAttributes (black) ;
  VisAttWorld->SetVisibility (Parameters::Instance()->main_config.worldVisibility) ;
  VisAttWorld->SetForceWireframe (true) ;
  worldLV->SetVisAttributes (VisAttWorld) ;

  // The calorimeter
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // the calorimeter is made of 1 or more identical modules. They are arranged in a grid, with
  // the possibility to "remove" some module (to make space for a beam pipe)
  G4cout << ">>>>>>>>>>>>> Build calorimeter "<< G4endl;
  G4VSolid * calorimeterS = new G4Box ("calorimeterS", 0.5 * calorimeter_size_x, 0.5 * calorimeter_size_y, 0.5 * calorimeter_size_z ) ;
  G4LogicalVolume * calorimeterLV = new G4LogicalVolume (calorimeterS, WoMaterial, "calorimeterLV") ;
  new G4PVPlacement (0, G4ThreeVector (calorimeter_pos_x, calorimeter_pos_y, calorimeter_pos_z), calorimeterLV, "calorimeterPV", worldLV, false, 0, checkOverlaps) ;
  G4VisAttributes* Calo_VisAtt = new G4VisAttributes(orange);
  Calo_VisAtt->SetVisibility(Parameters::Instance()->main_config.caloVisibility);
  Calo_VisAtt->SetForceWireframe(Parameters::Instance()->main_config.wireFrame);
  calorimeterLV->SetVisAttributes(Calo_VisAtt);

  // save calorimeter
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
  {
    CreateTree::Instance()->calorimeterID                = 0; // only 1 calo per simulation, please...
    CreateTree::Instance()->calorimeterMaterial          = 0; // irrelevant
    CreateTree::Instance()->calorimeterType              = 0; // irrelevant
    CreateTree::Instance()->calorimeterX                 = calorimeter_pos_x;
    CreateTree::Instance()->calorimeterY                 = calorimeter_pos_y;
    CreateTree::Instance()->calorimeterZ                 = calorimeter_pos_z;
    CreateTree::Instance()->calorimeterDX                = calorimeter_size_x;
    CreateTree::Instance()->calorimeterDY                = calorimeter_size_y;
    CreateTree::Instance()->calorimeterDZ                = calorimeter_size_z;
    CreateTree::Instance()->calorimeterSeparationZ       = 0; // irrelevant
    CreateTree::Instance()->calorimeterSections = 1;
    CreateTree::Instance()->calorimeter->Fill();
  }

  //-----------------------------------------//
  // Modules
  //-----------------------------------------//

  // int iForbid = 1;
  // int jForbid = 1;

  int moduleID = 0;
  int absorberID = 0;
  // int cellID = 0;
  // int crystalID = 0;
  G4cout << ">>>>>>>>>>>>> Build modules "<< G4endl;

  if(Parameters::Instance()->calorimeter_method == 1)
  {

    G4LogicalVolume *moduleLV  = modules_list[0]->GetModuleLV();
    G4VisAttributes* Module_VisAtt = new G4VisAttributes(black);  // color
    Module_VisAtt->SetVisibility(Parameters::Instance()->main_config.moduleVisibility);
    Module_VisAtt->SetForceWireframe(Parameters::Instance()->main_config.wireFrame);
    moduleLV->SetVisAttributes(Module_VisAtt);

    // change paradigm, pass to G4VPVParameterisation
    // makes code more symmetric, and simplifies saving of module positions
    // all we need to do is create a map of the module centers
    std::vector<G4double> x;
    std::vector<G4double> y;
    std::vector<G4double> z;
    // loop on modules x and y
    for(int iX = 0 ; iX < modules_nx; iX++)
    {
      for(int iY = 0 ; iY < modules_ny; iY++)
      {
        if((pipe_modules_nx + pipe_modules_ny) != 0)
        {
          // skip pipe modules
          if( (iX > ((modules_nx - pipe_modules_nx)/2 -1) ) && ( iX < (modules_nx - (modules_nx - pipe_modules_nx)/2 ) ) ) // if
          {
            if( (iY > ((modules_ny - pipe_modules_ny)/2 -1) ) && ( iY < (modules_ny - (modules_ny - pipe_modules_ny)/2 ) ) )
            {
              continue;
            }
          }
        }

        // calc module coordinates, assuming no gap among modules
        double module_center_x = module_size_x*(iX - (modules_nx-1)/(2.0)); // distance of module from x.y center of calo
        double module_center_y = module_size_y*(iY - (modules_ny-1)/(2.0)); // distance of module from x.y center of calo
        x.push_back(module_center_x * mm);
        y.push_back(module_center_y * mm);
        z.push_back(modules_list[0]->GetModuleZshift()); // always in the center of calo in z
        // std::cout << iX << " " << iY << std::endl;
        // std::cout << module_center_x << " " << module_center_y << std::endl;
      }
    }
    // G4cout << "points n " << x.size() << G4endl;
    // place modules
    G4VPVParameterisation* module_parametrization =
                           new ModuleParametrization (x.size(),
                                                      x,y,z,
                                                      G4ThreeVector( modules_list[0]->GetModuleSizeX(), modules_list[0]->GetModuleSizeY(), modules_list[0]->GetModuleSizeZ() ),
                                                      modules_list[0]->GetECALposition(),
                                                      modules_list[0]->GetSeparationPositionZ(),
                                                      modules_list[0]->GetNofSections() );
    //
    // we use as name the ecal identifier, it will be useful in shower saving and hybrid propagation
    std::stringstream moduleNameStream;
    moduleNameStream << "Module_" << modules_list[0]->GetECALposition();
    new G4PVParameterised(moduleNameStream.str().c_str(),       // their name
                          moduleLV,   // their logical volume
                          calorimeterLV,       // Mother logical volume
                          kUndefined,          // Are placed along this axis
                          x.size(),    // Number of chambers
                          module_parametrization,    // The parametrisation
                          false); // checking overlaps
    // PositionModulesWithReplica(moduleLV,calorimeterLV);

  }
  else
  {
    // for each spacal type,
    // 1. import CAD volume
    // 2. place modules
    for(int i = 0 ; i < modules_list.size() ; i++ )
    {
      // G4cout << "module n " << i << G4endl;

      // bool container_visibility = true;
      // bool container_wireframe = true;
      std::stringstream sname;
      sname << "volumeLV " << modules_list[i]->GetECALposition();
      recallCAD CADvolume;
      G4LogicalVolume* vol_logical = CADvolume.returnMeshedLV(modules_list[i]->GetContainerVolumeFileName().c_str()
                                                        , WoMaterial
                                                        , sname.str().c_str());
      // cow_mesh->SetScale(500);   

      sname.str("");
      sname << "volumePV " << modules_list[i]->GetECALposition();

      G4Colour container_color;
      if(i < color_vector.size()-1)
      {
        container_color = color_vector[i];
      }
      else
      {
        container_color = grey;
      }

      G4VisAttributes* vol_VisAtt = new G4VisAttributes(container_color);  // color
      vol_VisAtt->SetVisibility(containerVisibility);
      vol_VisAtt->SetForceWireframe(true);
      // vol_VisAtt->SetForceWireframe(false);
      vol_logical->SetVisAttributes(vol_VisAtt);

      new G4PVPlacement( 0
                       , G4ThreeVector(0, 0, 0)
                       , vol_logical
                       , sname.str().c_str()
                       , calorimeterLV
                       , false, 0
      );
      sname.str("");


      // place the module
      // get the module
      G4LogicalVolume *moduleLV  = modules_list[i]->GetModuleLV();
      //
      // get the points
      std::vector<Parameters::Point> coordinates = modules_list[i]->GetCoordinates();
        // G4cout << coordinates.size() << G4endl;
      std::vector<G4double> x;
      std::vector<G4double> y;
      std::vector<G4double> z;
      for(int iP = 0 ; iP < coordinates.size() ; iP++)
      {
        x.push_back(coordinates[iP].x * mm);
        y.push_back(coordinates[iP].y * mm);
        z.push_back(modules_list[i]->GetModuleZshift());
      }
      // place with parametrisation
      // G4cout << "points n " << x.size() << G4endl;
      G4VPVParameterisation* module_parametrization =
                             new ModuleParametrization (x.size(),
                                                        x,y,z,
                                                        G4ThreeVector( modules_list[i]->GetModuleSizeX(), modules_list[i]->GetModuleSizeY(), modules_list[i]->GetModuleSizeZ() ),
                                                        modules_list[i]->GetECALposition(),
                                                        modules_list[i]->GetSeparationPositionZ(),
                                                        modules_list[i]->GetNofSections() );
      //
      G4VisAttributes* Module_VisAtt = new G4VisAttributes(container_color);  // color
      Module_VisAtt->SetVisibility(modules_list[i]->GetModuleVisibility());
      Module_VisAtt->SetForceWireframe(modules_list[i]->GetModuleWireFrame());
      moduleLV->SetVisAttributes(Module_VisAtt);
      std::stringstream moduleNameStream;
      moduleNameStream << "Module_" << modules_list[i]->GetECALposition();
      new G4PVParameterised(moduleNameStream.str().c_str(),       // their name
                            moduleLV,   // their logical volume
                            vol_logical,       // Mother logical volume
                            kUndefined,          // Are placed along this axis
                            x.size(),    // Number of chambers
                            module_parametrization,    // The parametrisation
                            false); // checking overlaps
    }
  }

  //PG call the magnetic field initialisation
  if (B_field_intensity > 0.1 * tesla) ConstructField () ;

  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}

// void DetectorConstruction::PositionModulesWithReplica(G4LogicalVolume * moduleLV,G4LogicalVolume * calorimeterLV) // method not used anymore
// {
//
//   // int iMod = 0;
//   // int jMod = 0;
//   // int module_number = 0;
//   std::stringstream sname;
//   // build the container volumes
//   G4VSolid        *calo_bottom_S;
//   G4LogicalVolume *calo_bottom_LV;
//   G4VSolid        *calo_top_S;
//   G4LogicalVolume *calo_top_LV;
//   G4VSolid        *calo_left_S;
//   G4LogicalVolume *calo_left_LV;
//   G4VSolid        *calo_right_S;
//   G4LogicalVolume *calo_right_LV;
//
//   // full row box, will be used anyway
//   G4VSolid        *calo_full_row_S;
//   G4LogicalVolume *calo_full_row_LV;
//   calo_full_row_S = new G4Box ("calo_full_row_S", calorimeter_size_x/2.,module_size_y/2.,module_size_z/2.);
//   calo_full_row_LV = new G4LogicalVolume (calo_full_row_S, WoMaterial,"calo_full_row_LV") ;
//   G4VisAttributes* calo_full_row_VisAtt = new G4VisAttributes(black);  // color
//   calo_full_row_VisAtt->SetVisibility(false);
//   calo_full_row_VisAtt->SetForceWireframe(true);
//   calo_full_row_LV->SetVisAttributes(calo_full_row_VisAtt);
//   // short row box
//   G4VSolid        *calo_short_row_S;
//   G4LogicalVolume *calo_short_row_LV;
//   // sname << moduleName << "_PV";
//   // sname.str("")
//   // decide on replicas chain
//   if((pipe_modules_nx == 0) && (pipe_modules_ny == 0))// no pipe
//   {
//
//     new G4PVReplica("calo_full_PV",
//           calo_full_row_LV,
//           calorimeterLV,
//           kYAxis,
//           modules_ny,
//           module_size_y
//         );
//     //
//     new G4PVReplica("modules_replica",
//           moduleLV,
//           calo_full_row_LV,
//           kXAxis,
//           modules_nx,
//           module_size_x
//           );
//   }
//   else // with pipe
//   {
//
//     calo_bottom_S  = new G4Box ("calo_bottom_S", calorimeter_size_x/2.,(modules_ny - pipe_modules_ny)*module_size_y/4.,module_size_z/2.);
//     calo_top_S     = new G4Box ("calo_top_S"   , calorimeter_size_x/2.,(modules_ny - pipe_modules_ny)*module_size_y/4.,module_size_z/2.);
//     calo_left_S    = new G4Box ("calo_left_S"  , (modules_nx -pipe_modules_nx)*module_size_x/4., pipe_modules_ny * module_size_y/2.,module_size_z/2.);
//     calo_right_S   = new G4Box ("calo_right_S" , (modules_nx -pipe_modules_nx)*module_size_x/4., pipe_modules_ny * module_size_y/2.,module_size_z/2.);
//     calo_bottom_LV = new G4LogicalVolume (calo_bottom_S , WoMaterial,"calo_bottom_LV") ;
//     calo_top_LV    = new G4LogicalVolume (calo_top_S    , WoMaterial,"calo_top_LV") ;
//     calo_left_LV   = new G4LogicalVolume (calo_left_S   , WoMaterial,"calo_left_LV") ;
//     calo_right_LV  = new G4LogicalVolume (calo_right_S  , WoMaterial,"calo_right_LV") ;
//
//     new G4PVPlacement (0, G4ThreeVector (0,-calorimeter_size_y/2. + (modules_ny -pipe_modules_ny)*module_size_y/4.,
//                                          0),
//                                          calo_bottom_LV, "calo_bottom_PV",
//                                          calorimeterLV, false, 0, checkOverlaps) ;
//     //
//     new G4PVPlacement (0, G4ThreeVector (0,
//                                          calorimeter_size_y/2. - (modules_ny - pipe_modules_ny)*module_size_y/4.,
//                                          0),
//                                          calo_top_LV, "calo_top_PV",
//                                          calorimeterLV, false, 0, checkOverlaps) ;
//     //
//     new G4PVPlacement (0, G4ThreeVector (-calorimeter_size_x/2. + (modules_nx - pipe_modules_nx)*module_size_x/4.,
//                                          0,
//                                          0),
//                                          calo_left_LV, "calo_left_PV",
//                                          calorimeterLV, false, 0, checkOverlaps) ;
//     //
//     new G4PVPlacement (0, G4ThreeVector (+calorimeter_size_x/2. - (modules_nx - pipe_modules_nx)*module_size_x/4.,
//                                          0,
//                                          0),
//                                          calo_right_LV, "calo_right_PV",
//                                          calorimeterLV, false, 0, checkOverlaps) ;
//     //
//
//     G4VisAttributes* calo_elements_VisAtt = new G4VisAttributes(black);  // color
//     calo_elements_VisAtt->SetVisibility(false);
//     calo_elements_VisAtt->SetForceWireframe(true);
//     calo_bottom_LV->SetVisAttributes(calo_elements_VisAtt);
//     calo_top_LV   ->SetVisAttributes(calo_elements_VisAtt);
//     calo_left_LV  ->SetVisAttributes(calo_elements_VisAtt);
//     calo_right_LV ->SetVisAttributes(calo_elements_VisAtt);
//
//     // sname << ;
//     calo_short_row_S = new G4Box ("calo_short_row_S",(modules_nx - pipe_modules_nx)*module_size_x/4.,module_size_y/2.,module_size_z/2.);
//     calo_short_row_LV = new G4LogicalVolume (calo_short_row_S, WoMaterial,"calo_short_row_LV") ;
//
//     //
//     G4VisAttributes* calo_short_row_VisAtt = new G4VisAttributes(black);  // color
//     calo_short_row_VisAtt->SetVisibility(false);
//     calo_short_row_VisAtt->SetForceWireframe(true);
//     calo_short_row_LV->SetVisAttributes(calo_short_row_VisAtt);
//
//
//     new G4PVReplica("calo_full_bottom_PV",
//           calo_full_row_LV,
//           calo_bottom_LV,
//           kYAxis,
//           (modules_ny - pipe_modules_ny)/2.0,
//           module_size_y
//         );
//     //
//     new G4PVReplica("calo_full_top_PV",
//           calo_full_row_LV,
//           calo_top_LV,
//           kYAxis,
//           (modules_ny - pipe_modules_ny)/2.0,
//           module_size_y
//         );
//
//     new G4PVReplica("calo_short_left_row_PV",
//           calo_short_row_LV,
//           calo_left_LV,
//           kYAxis,
//           pipe_modules_ny,
//           module_size_y
//           );
//     //
//     new G4PVReplica("calo_short_right_row_PV",
//           calo_short_row_LV,
//           calo_right_LV,
//           kYAxis,
//           pipe_modules_ny,
//           module_size_y
//           );
//
//     new G4PVReplica("module_PV",
//           moduleLV,
//           calo_full_row_LV,
//           kXAxis,
//           modules_nx,
//           module_size_x
//           );
//     //
//     new G4PVReplica("module_PV",
//           moduleLV,
//           calo_short_row_LV,
//           kXAxis,
//           (modules_nx-pipe_modules_nx)/2.0,
//           module_size_x
//           );
//   }
//
//
//
//
//   sname.str("");
// }




void DetectorConstruction::initializeMaterials ()
{
  // define materials
  AirMaterial = MyMaterials::Air();
  AirKillerMaterial = MyMaterials::AirKiller();
  //LAPPD_average = MyMaterials::LAPPD_average();

  WoMaterial = NULL ;
  if      ( world_material == 1 ) WoMaterial = MyMaterials::Air () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "Wo. material: "<< WoMaterial << G4endl ;

  InterfaceContainerMaterial = NULL;
  if(readoutType == 0 || readoutType == 3 || readoutType == 4)
  {
    InterfaceContainerMaterial = AirMaterial;
  }
  else
  {
    InterfaceContainerMaterial = AirKillerMaterial;
  }



  InterfaceMaterial  = NULL;
  if      ( cone_material == 0 ) InterfaceMaterial = MyMaterials::Air () ;
  else if ( cone_material == 1 ) InterfaceMaterial = MyMaterials::PLEX (scaleFactor) ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid cone material specifier " << cone_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "Interface material: "<< InterfaceMaterial << G4endl ;

  GapAbsToInterfaceMaterial  = NULL;
  if      ( gap_abs_interface_material == 0 ) GapAbsToInterfaceMaterial = MyMaterials::Air () ;
  else if ( gap_abs_interface_material == 1 ) GapAbsToInterfaceMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapAbsToInterfaceMaterial material specifier " << gap_abs_interface_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "GapAbsToInterfaceMaterial material: "<< GapAbsToInterfaceMaterial << G4endl ;

  GapInterfaceToReadoutMaterial  = NULL;
  if      ( gap_interface_readout_material == 0 ) GapInterfaceToReadoutMaterial = MyMaterials::Air () ;
  else if ( gap_interface_readout_material == 1 ) GapInterfaceToReadoutMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapInterfaceToReadoutMaterial material specifier " << gap_interface_readout_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "GapInterfaceToReadoutMaterial material: "<< GapInterfaceToReadoutMaterial << G4endl ;


  PLEXMaterial = NULL ;
  PLEXMaterial = MyMaterials::PLEX (scaleFactor) ;
  //
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "PLEX. material: "<< PLEXMaterial << G4endl ;
  //
  PVCMaterial = NULL ;
  PVCMaterial = MyMaterials::PVC () ;
  //
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "PVC. material: "<< PVCMaterial << G4endl ;


}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void DetectorConstruction::ConstructField ()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl ;
  static G4TransportationManager * trMgr = G4TransportationManager::GetTransportationManager () ;

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager * globalFieldMgr = trMgr->GetFieldManager () ;

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector (
      0.0522 * B_field_intensity,
      0.0522 * B_field_intensity,
      0.9973 * B_field_intensity
    ) ;

    B_field = new G4UniformMagField (fieldVector) ;
    globalFieldMgr->SetDetectorField (B_field) ;
    globalFieldMgr->CreateChordFinder (B_field) ;
    globalFieldMgr->GetChordFinder ()->SetDeltaChord (0.005 * mm) ;
    B_field_IsInitialized = true ;
  }
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl ;
  return ;
}


