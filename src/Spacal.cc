#include "Parameters.hh"
#include "Spacal.hh"
#include "CreateTree.hh"
#include "ConfigFile.hh"
#include "HoleParametrization.hh"

#include <algorithm>
#include <string>
#include <sstream>

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
#include "recallCAD.hh"


using namespace CLHEP;

Spacal::Spacal()
{

}


Spacal::~Spacal(){}

Spacal::Spacal(Data_t input_parameters)
{

  moduleName = "module";

  // prepare variables for Spacal module, by reading config file
  parameters = input_parameters;

  // start filling the module data saving structure
  module.fileName   = parameters.fileName;
  ecal_position     = parameters.ecal_position;
  volume_file       = parameters.volume_file;
  module_visibility = parameters.moduleVisibility;
  module_wireframe  = parameters.wireFrame;
  // sections          = parameters.NofSpacalSegments;
  separation_z_position = parameters.cell_separator_position;
  moduleZshift = parameters.moduleZshift;
  if(Parameters::Instance()->calorimeter_method == 2)
  {
    ScanMap(Parameters::Instance()->mapEcal);
    if(ecal_position != 0)
    {
      // run on all positions, save only the ones of this module
      for(int iP = 0 ; iP < Parameters::Instance()->coordinates.size() ; iP++ )
      {
        if(Parameters::Instance()->coordinates[iP].z == ecal_position)
        {
          coordinates.push_back(Parameters::Instance()->coordinates[iP]);
        }
      }
    }
  }
  // find xy min and max on the map



  //-----------------------------------------------------
  //------------- Define colors --------------
  //-----------------------------------------------------
  white   = G4Colour(1.00, 1.00, 1.00) ;  // white
  grey    = G4Colour(0.50, 0.50, 0.50) ;  // grey
  black   = G4Colour(0.00, 0.00, 0.00) ;  // black
  red     = G4Colour(1.00, 0.00, 0.00) ;  // red
  green   = G4Colour(0.00, 1.00, 0.00) ;  // green
  blue    = G4Colour(0.00, 0.00, 1.00) ;  // blue
  cyan    = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  air     = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  magenta = G4Colour(1.00, 0.00, 1.00) ;  // magenta
  yellow  = G4Colour(1.00, 1.00, 0.00) ;  // yellow
  brass   = G4Colour(0.80, 0.60, 0.40) ;  // brass
  brown   = G4Colour(0.70, 0.40, 0.10) ;  // brown
  orange  = G4Colour(1.00, 0.33, 0.00) ;  // orange

  // make volume
  ComputeVolumes();

  // materials
  InitializeMaterials() ;


  //declear a stringstream for composition of names
  std::stringstream sname;
  std::stringstream smodule;
  std::string moduleName = "module";

  sname << moduleName << "_S";
  moduleS = new G4Box (sname.str().c_str(), 0.5 * module_size_x, 0.5 * module_size_y, 0.5 * module_size_z ) ;
  sname.str("");
  sname << moduleName << "_LV";
  moduleLV = new G4LogicalVolume (moduleS, WoMaterial, sname.str().c_str()) ;
  sname.str("");

  G4VisAttributes* Module_VisAtt = new G4VisAttributes(black);  // color
  Module_VisAtt->SetVisibility(parameters.moduleVisibility);
  Module_VisAtt->SetForceWireframe(true);
  moduleLV->SetVisAttributes(Module_VisAtt);

  // build the module
  BuildModule(moduleLV);

}

void Spacal::ComputeVolumes()
{
  ComputeAbsorber();
  ComputeCells();
  ComputeInterface();
  ComputeReadout();
  ComputeModule();
}

// routine to create the absorber object
void Spacal::ComputeAbsorber()
{
  // rescale AbsorberSizeZ, adding separation_thickness if it's there
  CreateTree::Instance()->absorberLength = parameters.AbsSizeZ;
  AbsSizeX = parameters.AbsSizeX;
  AbsSizeY = parameters.AbsSizeY;
  AbsSizeZ = parameters.AbsSizeZ + parameters.separation_thickness;
  absorber.SetName(parameters.AbsName);
  absorber.SetDimensions( parameters.AbsSizeX,
                          parameters.AbsSizeY,
                          AbsSizeZ);
  absorber.SetPosition(parameters.AbsPositionX,
                       parameters.AbsPositionY,
                       parameters.AbsPositionZ);
  absorber.SetMaterial(parameters.AbsMaterial,parameters.W_fraction);
  absorber.SetSingleSeparationLayer(parameters.SingleSeparationLayer);
  absorber.SetSeparationMaterial(parameters.separation_material,parameters.W_fraction);
  absorber.SetLAPPDMaterials(parameters.LAPPD_materials,parameters.W_fraction);
  absorber.SetOpticalSurface(parameters.absorber_reflectivity,
                             parameters.absorber_specular_lobe,
                             parameters.absorber_specular_spike,
                             parameters.absorber_backscatter,
                             parameters.absorber_sigma_alpha);
  absorber.SetEsrThickness(parameters.separation_thickness);
}



// routine to create the cell objects
void Spacal::ComputeCells()
{

  for(unsigned int iCell = 0 ; iCell < parameters.CellNames.size(); iCell++)
  {
    Cell cell;
    cell.SetID(iCell);
    cell.SetName(parameters.CellNames[iCell]);
    G4double modPositionZ;
    if(parameters.CellPositionZ[iCell] > 0)
    {
      modPositionZ = parameters.CellPositionZ[iCell] + (absorber.GetEsrThickness())/2.0;
    }
    if(parameters.CellPositionZ[iCell] < 0)
    {
      modPositionZ = parameters.CellPositionZ[iCell] - (absorber.GetEsrThickness())/2.0;
    }
    if(parameters.CellPositionZ[iCell] == 0)
    {
      modPositionZ = 0;
    }
    cell.SetPosition          (parameters.CellPositionX[iCell],
                               parameters.CellPositionY[iCell],
                               modPositionZ);
    cell.SetCrystalMaterial   (parameters.CellMaterial[iCell],
                               parameters.user_lightyield,
                               parameters.abs_length_scale_factor,
                               parameters.cladding_abs_length_scale_factor);
    cell.SetGapSize(parameters.gapSize);
    // add crystal material to list of crystal materials for this module and to the global one, only if not already there
    bool isThere = false;
    for(unsigned int i = 0 ; i < parameters.crystalMaterialList.size();i++)
    {
      if(parameters.CellMaterial[iCell] == parameters.crystalMaterialList[i])
      {
        isThere = true;
      }
    }
    if(!isThere)
    {
      //add to this module
      parameters.crystalMaterialList.push_back(parameters.CellMaterial[iCell]);
      // and add to global calorimeters
      Parameters::Instance()->AddCrystalMaterial(parameters.CellMaterial[iCell],cell.GetCrystalMaterial());
    }

    cell.SetXelements                (parameters.CellXelements[iCell]);
    cell.SetYelements                (parameters.CellYelements[iCell]);
    cell.SetNominalCrystalDimensions (parameters.CellCrystalSizeX[iCell],
                                      parameters.CellCrystalSizeY[iCell],
                                      parameters.CellCrystalSizeZ[iCell]);
    cell.SetCrystalPitch             (parameters.CellCrystalPitchX[iCell],
                                      parameters.CellCrystalPitchY[iCell]);
    cell.SetAirLayer                 (parameters.CellAirLayer[iCell]);
    // cell.SetExtGapMaterial    (CellExtGapMaterial[iCell]);
    cell.SetIntGapMaterial           (parameters.CellIntGapMaterial[iCell]);
    cell.SetStaggering               (parameters.CellStaggering[iCell]);
    cell.SetStaggeringAxis           (parameters.CellStaggeringAxis[iCell]);
    cell.SetStaggeringSize           (parameters.CellStaggeringSize[iCell]);
    cell.SetStaggeringParity         (parameters.CellStaggeringParity[iCell]);
    cell.SetStaggeringRemove         (parameters.CellStaggeringRemove[iCell]);
    cell.SetHoleShape                (parameters.CellHoleShape[iCell]);
    cell.SetCrystalShape             (parameters.CellCrystalShape[iCell]);
    cell.SetCrystalCladding          (parameters.CellCrystalCladding[iCell]);

    // not very elegant since it's a module parameter, but ok
    cell.SetCrystalInnerCladdingFraction(parameters.crystal_inner_cladding_fraction);
    cell.SetCrystalOuterCladdingFraction(parameters.crystal_outer_cladding_fraction);


    cell.MakeCellStruture();
    cell.CalculateCellDimensions();

    // // find if the absorber is made of one or two sections
    // for(int iCell = 0 ; iCell < data.CellPositionZ.size(); iCell++)
    // {
    //   // double meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2])/2.0;
    double positionZ = cell.GetPositionZ();
    double dimZ = cell.GetSizeZ();
    bool zCellIsAlreadyThere = false;
    for(unsigned int iz = 0; iz < absorber.zPos.size(); iz++)
    {
      if(positionZ == absorber.zPos[iz])
      {
        zCellIsAlreadyThere = true;
      }
    }
    if(!zCellIsAlreadyThere)
    {
      absorber.zPos.push_back(positionZ);
      absorber.zDim.push_back(dimZ);
    }





    absorber.AddCell(cell);
  }

  absorber.NofSpacalSegments = absorber.zPos.size();
  sections                   = absorber.NofSpacalSegments;

  G4cout << "Absorber " << absorber.GetName()
        <<  " cells = " << absorber.GetNumberOfCells() << G4endl;

}


void Spacal::ComputeInterface()
{
  InterfaceSizeX = parameters.AbsSizeX;
  InterfaceSizeY = parameters.AbsSizeY;
  InterfaceSizeZ = parameters.InterfaceSizeZ;
  if(parameters.gapSize > 0) // rescale if there is a gap
  {
    InterfaceSizeZ = InterfaceSizeZ+2.0*parameters.gapSize;
  }
  if(parameters.abs_interface_extraGap > 0) // rescale if there is an extra gap only between the absorber and the interface
  {
    InterfaceSizeZ = InterfaceSizeZ+parameters.abs_interface_extraGap;
  }
}


void Spacal::ComputeReadout()
{
  ReadoutSizeX = parameters.AbsSizeX;
  ReadoutSizeY = parameters.AbsSizeY;
  ReadoutSizeZ = parameters.ReadoutSizeZ;

  if(parameters.readoutType == -1)
  {
    G4cout << "No valid readout choice given, setting to default, no readout volume" << G4endl ;
    parameters.readoutType = 5;
  }

  if(parameters.readoutType == 0 || parameters.readoutType == 3)// just 1 cm air gap
  {
    LGrelevantDimension = 10.0 *mm;
  }
  if(parameters.readoutType == 4)
  {
    LGrelevantDimension = 0.1 *mm;
  }
  if(parameters.readoutType == 1)// tb 2018
  {
    LGrelevantDimension = 20.0 *mm;
  }
  if (parameters.readoutType == 2 | parameters.readoutType == 7 ) // tb 2019
  {
    LGrelevantDimension = 15.3 *mm;
  }

}

void Spacal::ComputeModule()
{
  module_size_x = (AbsSizeX);
  module_size_y = (AbsSizeY);
  module_size_z = (AbsSizeZ + 2.0*InterfaceSizeZ + 2.0*ReadoutSizeZ);

  // modules will be replicas, no need to generate a complex
  // array of positions. this will be just a "fake" module always in 0
  module_pos_x = 0;
  module_pos_y = 0;
  module_pos_z = 0;
}



void Spacal::BuildModule(G4LogicalVolume * moduleLV)
{
  int crystalID = 0;
  int cellID = 0;


  //-----------------------------------//
  //---------------------------------- //
  // The interface volumes
  //-----------------------------------//
  //-----------------------------------//

  //-----------------------------------//
  // positive z inteface
  //-----------------------------------//

  sname << "Interface_Positive_Z_"
        << moduleName
        << "_S";
  Interface_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * (InterfaceSizeZ)) ;
  sname.str("");
  sname << "Interface_Positive_Z_"
        << moduleName
        << "_LV";
  Interface_Positive_Z_LV = new G4LogicalVolume (Interface_Positive_Z_S, InterfaceContainerMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Interface_Positive_Z_"
        << moduleName
        << "_PV";
  Interface_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 +(AbsSizeZ+InterfaceSizeZ)/2.0),
                                                                 Interface_Positive_Z_LV, sname.str().c_str(),
                                                                 moduleLV, false, 0, parameters.checkOverlaps) ;
  sname.str("");
  cout << "Z POSITION OF Interface_Positive_Z_" << moduleName << "_PV "<< AbsSizeZ << " " << InterfaceSizeZ << endl;

  //-----------------------------------//
  // elements of positive z interface
  //-----------------------------------//
  // place an air/grease layer at the beginning ad at the end of the interfaces
  // prepare pointers

  sname << "Gap_Abs_Interface_Positive_Z_"
        << moduleName
        << "_S";
  Gap_Abs_Interface_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * (parameters.gapSize + parameters.abs_interface_extraGap)) ;
  sname.str("");
  sname << "Gap_Abs_Interface_Positive_Z_"
        << moduleName
        << "_LV";
  Gap_Abs_Interface_Positive_Z_LV = new G4LogicalVolume (Gap_Abs_Interface_Positive_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Gap_Abs_Interface_Positive_Z_"
        << moduleName
        << "_PV";
  Gap_Abs_Interface_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 -(InterfaceSizeZ-(parameters.gapSize + parameters.abs_interface_extraGap))/2.0),
                                                                 Gap_Abs_Interface_Positive_Z_LV, sname.str().c_str(),
                                                                 Interface_Positive_Z_LV, false, 0, parameters.checkOverlaps) ;
  sname.str("");

  if(parameters.esr_on_positive_exit)
  {
    sname << "reflector_single_readout_positive";
    G4OpticalSurface* reflector_single_readout_positive = new G4OpticalSurface(sname.str().c_str());
    reflector_single_readout_positive->SetType(dielectric_metal);
    reflector_single_readout_positive->SetFinish(polished);
    reflector_single_readout_positive->SetModel(unified);
    reflector_single_readout_positive->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
    G4LogicalBorderSurface* l_reflector_single_readout_positive = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                   Gap_Abs_Interface_Positive_Z_PV,
                                                                   Interface_Positive_Z_PV,
                                                                   reflector_single_readout_positive);
    sname.str("");
  }

  sname << "Gap_Interface_Readout_Positive_Z_"
        << moduleName
        << "_S";
  Gap_Interface_Readout_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * parameters.gapSize) ;
  sname.str("");
  sname << "Gap_Interface_Readout_Positive_Z_"
        << moduleName
        << "_LV";
  Gap_Interface_Readout_Positive_Z_LV = new G4LogicalVolume (Gap_Interface_Readout_Positive_Z_S, GapInterfaceToReadoutMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Gap_Interface_Readout_Positive_Z_"
        << moduleName
        << "_PV";
  Gap_Interface_Readout_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 +(InterfaceSizeZ-parameters.gapSize)/2.0),
                                                                 Gap_Interface_Readout_Positive_Z_LV, sname.str().c_str(),
                                                                 Interface_Positive_Z_LV, false, 0, parameters.checkOverlaps) ;
  sname.str("");
  // end of positive z
  //-----------------------------------//


  //-----------------------------------//
  // negative z
  //-----------------------------------//

  sname << "Interface_Negative_Z_"
        << moduleName
        << "_S";
  Interface_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * InterfaceSizeZ ) ;
  sname.str("");

  sname << "Interface_Negative_Z_"
        << moduleName
        << "_LV";
  Interface_Negative_Z_LV = new G4LogicalVolume (Interface_Negative_Z_S, InterfaceContainerMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Interface_Negative_Z_"
        << moduleName
        << "_PV";
  Interface_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 -(AbsSizeZ+InterfaceSizeZ)/2.0),
                                                                 Interface_Negative_Z_LV, sname.str().c_str(),
                                                                 moduleLV, false, 0, parameters.checkOverlaps) ;
  sname.str("");
  //


  sname << "Gap_Abs_Interface_Negative_Z_"
        << moduleName
        << "_S";
  Gap_Abs_Interface_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * (parameters.gapSize + parameters.abs_interface_extraGap)) ;
  sname.str("");
  sname << "Gap_Abs_Interface_Negative_Z_"
        << moduleName
        << "_LV";
  Gap_Abs_Interface_Negative_Z_LV = new G4LogicalVolume (Gap_Abs_Interface_Negative_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Gap_Abs_Interface_Negative_Z_"
        << moduleName
        << "_PV";
  Gap_Abs_Interface_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 +(InterfaceSizeZ- (parameters.gapSize + parameters.abs_interface_extraGap))/2.0),
                                                                 Gap_Abs_Interface_Negative_Z_LV, sname.str().c_str(),
                                                                 Interface_Negative_Z_LV, false, 0, parameters.checkOverlaps) ;
  sname.str("");

  // set here already a reflector at the exit of fibres, if required
  
  if(parameters.esr_on_negative_exit)
  {
    sname << "reflector_single_readout_negative";
    G4OpticalSurface* reflector_single_readout_negative = new G4OpticalSurface(sname.str().c_str());
    reflector_single_readout_negative->SetType(dielectric_metal);
    reflector_single_readout_negative->SetFinish(polished);
    reflector_single_readout_negative->SetModel(unified);
    reflector_single_readout_negative->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
    G4LogicalBorderSurface* l_reflector_single_readout_negative = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                   Gap_Abs_Interface_Negative_Z_PV,
                                                                   Interface_Negative_Z_PV,
                                                                   reflector_single_readout_negative);
    sname.str("");
  }
  

  sname << "Gap_Interface_Readout_Negative_Z_"
        << moduleName
        << "_S";
  Gap_Interface_Readout_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * parameters.gapSize) ;
  sname.str("");
  sname << "Gap_Interface_Readout_Negative_Z_"
        << moduleName
        << "_LV";
  Gap_Interface_Readout_Negative_Z_LV = new G4LogicalVolume (Gap_Interface_Readout_Negative_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Gap_Interface_Readout_Negative_Z_"
        << moduleName
        << "_PV";
  Gap_Interface_Readout_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 -(InterfaceSizeZ-parameters.gapSize)/2.0),
                                                                 Gap_Interface_Readout_Negative_Z_LV, sname.str().c_str(),
                                                                 Interface_Negative_Z_LV, false, 0, parameters.checkOverlaps) ;
  sname.str("");
  // end of negative z
  //-----------------------------------//

  // vis att.
  G4VisAttributes* Interface_VisAtt = new G4VisAttributes(red);  // color
  Interface_VisAtt->SetVisibility(parameters.interfaceVisibility);
  Interface_VisAtt->SetForceWireframe(true);
  Interface_Positive_Z_LV->SetVisAttributes(Interface_VisAtt);
  Interface_Negative_Z_LV->SetVisAttributes(Interface_VisAtt);

  G4VisAttributes* Gap_VisAtt = new G4VisAttributes(cyan);  // color
  Gap_VisAtt->SetVisibility(parameters.gapsVisibility);
  Gap_VisAtt->SetForceWireframe(true);
  Gap_Abs_Interface_Positive_Z_LV->SetVisAttributes(Gap_VisAtt);
  Gap_Interface_Readout_Positive_Z_LV->SetVisAttributes(Gap_VisAtt);
  Gap_Abs_Interface_Negative_Z_LV->SetVisAttributes(Gap_VisAtt);
  Gap_Interface_Readout_Negative_Z_LV->SetVisAttributes(Gap_VisAtt);
  // end of the interfaces
  //--------------------------------------- //




  //--------------------------------------- //
  // The readout volumes
  // prepare pointers

  sname << "Readout_Positive_Z_"
        << moduleName
        << "_S";
  Readout_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * ReadoutSizeX, 0.5 * ReadoutSizeY, 0.5 * ReadoutSizeZ ) ;
  sname.str("");
  sname << "Readout_Positive_Z_"
        << moduleName
        << "_LV";
  Readout_Positive_Z_LV = new G4LogicalVolume (Readout_Positive_Z_S, AirKillerMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Readout_Positive_Z_"
        << moduleName
        << "_PV";
  Readout_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 +(AbsSizeZ+2.0*InterfaceSizeZ+ReadoutSizeZ)/2.0),
                                                                 Readout_Positive_Z_LV, sname.str().c_str(),
                                                                 moduleLV, false, 0, parameters.checkOverlaps) ;
  sname.str("");


  sname << "Readout_Negative_Z_"
        << moduleName
        << "_S";
  Readout_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * ReadoutSizeX, 0.5 * ReadoutSizeY, 0.5 * ReadoutSizeZ ) ;
  sname.str("");
  sname << "Readout_Negative_Z_"
        << moduleName
        << "_LV";
  Readout_Negative_Z_LV = new G4LogicalVolume (Readout_Negative_Z_S, AirKillerMaterial, sname.str().c_str()) ;
  sname.str("");

  sname << "Readout_Negative_Z_"
        << moduleName
        << "_PV";
  Readout_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                 0,
                                                                 -(AbsSizeZ+2.0*InterfaceSizeZ+ReadoutSizeZ)/2.0),
                                                                 Readout_Negative_Z_LV, sname.str().c_str(),
                                                                 moduleLV, false, 0, parameters.checkOverlaps) ;
  sname.str("");
  G4VisAttributes* Readout_VisAtt = new G4VisAttributes(green);  // color
  Readout_VisAtt->SetVisibility(parameters.readoutVisibility);
  Readout_VisAtt->SetForceWireframe(true);
  Readout_Positive_Z_LV->SetVisAttributes(Readout_VisAtt);
  Readout_Negative_Z_LV->SetVisAttributes(Readout_VisAtt);
  // end of the readout
  //--------------------------------------- //



  //--------------------------------------- //
  // Specific optical readout part
  // this is intentionally hardcoded (but the meaningful job is done in different routines)
  // because there is no way to modularize every possible readout
  // the user will have to take care of writing another function to create the volumes of
  // light concentrators and photo-detectors, and modify the CreateTree and SteppingAction classes
  // to have a personalized readout of the light transport

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
  rotationMatrix->rotateY(180.*deg);
  G4double shiftCAD = - InterfaceSizeZ/2.0 + (parameters.gapSize + parameters.abs_interface_extraGap);


  if(parameters.readoutType == 0 || parameters.readoutType == 3 || parameters.readoutType == 4 || parameters.readoutType == 5|| parameters.readoutType == 6 || parameters.readoutType == 8)
  {
    // no need for any additional volume
  }
  else
  {
    ConstructLightGuides(Interface_Negative_Z_LV,Interface_Negative_Z_PV,moduleName,LGrelevantDimension,rotationMatrix,-shiftCAD);
    // // on the other side, we need to rotate the light guides...
    ConstructLightGuides(Interface_Positive_Z_LV,Interface_Positive_Z_PV,moduleName,LGrelevantDimension,0, +shiftCAD);
  }

  // // and pmts
  // // no rotation needed
  if(parameters.readoutType == 5 || parameters.readoutType == 8)
  {
    // no volumes needed
  }
  else
  {
    if(parameters.readoutType == 3 || parameters.readoutType == 4)
    {
      ConstructOneLargePMT(Readout_Positive_Z_LV,moduleName,InterfaceSizeX/2.0,0);
      ConstructOneLargePMT(Readout_Negative_Z_LV,moduleName,InterfaceSizeX/2.0,0);
    }
    else
    {
      if(parameters.readoutType == 6)
      {
        ConstructClearFiber(Interface_Positive_Z_LV,Interface_Positive_Z_PV,Gap_Abs_Interface_Positive_Z_PV,Gap_Interface_Readout_Positive_Z_PV,moduleName);
        ConstructClearFiber(Interface_Negative_Z_LV,Interface_Negative_Z_PV,Gap_Abs_Interface_Negative_Z_PV,Gap_Interface_Readout_Negative_Z_PV,moduleName);
         //with this configuration, we also always want to have a realistic readout:
        ConstructPMTs(Readout_Positive_Z_LV,moduleName,LGrelevantDimension,0);
        ConstructPMTs(Readout_Negative_Z_LV,moduleName,LGrelevantDimension,0);
      }
      else
      {
        ConstructPMTs(Readout_Positive_Z_LV,moduleName,LGrelevantDimension,0);
        ConstructPMTs(Readout_Negative_Z_LV,moduleName,LGrelevantDimension,0);
      }
    }
  }






  //--------------------------------------- //
  // The absorber
  //--------------------------------------- //
  // the absorber is 1, but if there is a cell longitudinal separation
  // it means that it is logically divided in two
  // now, if this separation is of type 0 or 1 (just air, or aluminization)
  // nothing need to be done here
  // otherwise, if the separation is type 2, with a foil of esr, the absorber need to be effectively
  // separated in two
  // to do this, it's created larger by a factor equal to the esr thickness
  // then a air volume is created in the middle

  // prepare pointers

  // create the shape
  sname << absorber.GetName() << "_S"
        << "_" << moduleName ;
  Absorber_S    = new G4Box (sname.str().c_str(),
                         0.5 * absorber.GetSizeX(),
                         0.5 * absorber.GetSizeY(),
                         0.5 * absorber.GetSizeZ());
  sname.str("");
  // create the logical volume
  sname << absorber.GetName() << "_LV"
        << "_" << moduleName;
  Absorber_LV   = new G4LogicalVolume (Absorber_S,
                                     absorber.GetMaterial(),
                                     sname.str().c_str()) ;
  sname.str("");
  // create the physical volume (i.e. place it in space)
  sname << absorber.GetName() << "_PV"
        << "_" << moduleName;
  Absorber_PV  =  new G4PVPlacement (0,
                                 G4ThreeVector (absorber.GetPositionX(),
                                                absorber.GetPositionY(),
                                                absorber.GetPositionZ()),
                                 Absorber_LV,
                                 sname.str().c_str(),
                                 moduleLV,
                                 false, 0, parameters.checkOverlaps);
  sname.str("");

  // G4cout << "abs x " << std::setprecision(std::numeric_limits<G4double>::digits10 + 1) << absorber.GetPositionX() << G4endl;
  // G4cout << "abs y " << absorber.GetPositionY() << G4endl;
  // G4cout << "abs z " << absorber.GetPositionZ() << G4endl;

  G4VisAttributes* Absorber_VisAtt = new G4VisAttributes(cyan);  // color
  Absorber_VisAtt->SetVisibility(parameters.absorberVisibility);
  Absorber_VisAtt->SetForceWireframe(true);
  Absorber_LV->SetVisAttributes(Absorber_VisAtt);




  // save the absorber in the TTree in output file
  if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
  {
    CreateTree::Instance()->absorberID       = 0;   // only one abs per module anyway...
    CreateTree::Instance()->absorberMaterial = absorber.GetMaterialIdentifier();
    CreateTree::Instance()->absorberType     = ecal_position;
    CreateTree::Instance()->absorberX        = module_pos_x + absorber.GetPositionX();
    CreateTree::Instance()->absorberY        = module_pos_y + absorber.GetPositionY();
    CreateTree::Instance()->absorberZ        = module_pos_z + absorber.GetPositionZ();
    CreateTree::Instance()->absorberDX       = absorber.GetSizeX();
    CreateTree::Instance()->absorberDY       = absorber.GetSizeY();
    CreateTree::Instance()->absorberDZ       = absorber.GetSizeZ();
    CreateTree::Instance()->absorberSeparationZ = separation_z_position;
    CreateTree::Instance()->absorberSections = sections;
    CreateTree::Instance()->absorbers->Fill();
  }

  // G4cout << "abs x " << module_pos_x[iMod][jMod] + absorber.GetPositionX() << G4endl;
  // G4cout << "abs y " << module_pos_y[iMod][jMod] + absorber.GetPositionY() << G4endl;
  // G4cout << "abs z " << module_pos_z[iMod][jMod] + absorber.GetPositionZ() << G4endl;


  // now do sub absorber volumes, to contain holes. also do the esr volume, if needed.
  // base the choice on the number of sections in the absorber, and on the separation_thickness



  abs_section_S  = new G4VSolid        *[sections];
  abs_section_LV = new G4LogicalVolume *[sections];
  abs_section_PV = new G4PVPlacement   *[sections];

  for(int iSec = 0; iSec < sections; iSec++)
  {
    sname.str("");
    sname << "Abs_section_S_" << iSec;
    abs_section_S[iSec]    = new G4Box (sname.str().c_str(),
                             0.5 * absorber.GetSizeX(),
                             0.5 * absorber.GetSizeY(),
                             0.5 * absorber.zDim[iSec]);
    sname.str("");
    sname << "Abs_section_LV_" << iSec;
    abs_section_LV[iSec]   = new G4LogicalVolume (abs_section_S[iSec],
                                         absorber.GetMaterial(),
                                         sname.str().c_str()) ;
    sname.str("");
    sname << "Abs_section_PV_" << iSec;
    abs_section_PV[iSec]  =  new G4PVPlacement (0,
                                   G4ThreeVector (absorber.GetPositionX(),
                                                  absorber.GetPositionY(),
                                                  absorber.zPos[iSec]),
                                   abs_section_LV[iSec],
                                   sname.str().c_str(),
                                   Absorber_LV,
                                   false, 0, parameters.checkOverlaps);
    sname.str("");

    // // the absorber does not have an index of refraction, but it will be touched by optical photons
    // // if no surface is defined, it will kill them. So we define a G4LogicalSkinSurface
    // sname << "SkinSurf_Abs_section_" << iSec << "_" << absorber.GetName()
    //       << "_" << moduleName;
    // G4OpticalSurface* absSurface = new G4OpticalSurface(sname.str().c_str());
    // sname.str("");
    // absSurface->SetType(dielectric_metal);
    // absSurface->SetFinish(ground);
    // absSurface->SetModel(unified);
    // absSurface->SetSigmaAlpha(absorber.GetSigmaAlpha());
    // absSurface->SetMaterialPropertiesTable(absorber.GetOpticalSurface());
    // sname << "LSS_"
    //       << absorber.GetName() << "_LSS"
    //       << "_" << moduleName;
    //
    // G4LogicalSkinSurface *abs_LSS;
    // // create a skin surface, unless the absorber is made of air (special case to study individual fibers in air)
    // if(parameters.AbsMaterial != 9) // 9) is Air
    // {
    //   abs_LSS = new G4LogicalSkinSurface(sname.str().c_str(),abs_section_LV[iSec],absSurface);
    // }
    // sname.str("");

    G4VisAttributes* Absorber_Section_VisAtt = new G4VisAttributes(blue);  // color
    Absorber_Section_VisAtt->SetVisibility(parameters.absorberVisibility);      // absorbers are always visible
    Absorber_Section_VisAtt->SetForceWireframe(parameters.wireFrame);
    abs_section_LV[iSec]->SetVisAttributes(Absorber_VisAtt);
  }


  //-------------------------------//
  // Layer diving the absorber, on which reflectors are
  // Must be set to LAPPD
  //-------------------------------//

  

  if((parameters.separation_thickness > 0) && (sections > 1)) // only if reflector thickness is > 0, remember that it is set to 0 if the separation is made by aluminization
  {

    // do a container volume anyway,
    // just be careful of the material

    sname << "esr_S_"
        << absorber.GetName()
        << "_" << moduleName;
    esr_S    = new G4Box (sname.str().c_str(),
                            0.5 * absorber.GetSizeX(),
                            0.5 * absorber.GetSizeY(),
                            0.5 * (parameters.separation_thickness));
    sname.str("");
    // create the logical volume
    sname << "esr_LV_" << absorber.GetName()
          << "_" << moduleName;
    // choose material: 
    // if single separation 
    // G4cout << "DEBUG: " << absorber.GetSeparationMaterial() << G4endl;
    esr_LV   = new G4LogicalVolume (esr_S,
                                        MyMaterials::Air(),
                                        sname.str().c_str()) ;
    sname.str("");
    // create the physical volume (i.e. place it in space)
    sname << "esr_PV_"
          << absorber.GetName()
          << "_" << moduleName;;
    esr_PV  =  new G4PVPlacement (0,
                                    G4ThreeVector (0,
                                                  0,
                                                  parameters.cell_separator_position),
                                    esr_LV,
                                    sname.str().c_str(),
                                    Absorber_LV,
                                    false, 0, parameters.checkOverlaps);
    sname.str("");
    G4VisAttributes* esr_VisAtt = new G4VisAttributes(green);  // color
    esr_VisAtt->SetVisibility(parameters.esrVisibility);      // absorbers are always visible
    esr_VisAtt->SetForceWireframe(parameters.wireFrame);
    esr_LV->SetVisAttributes(esr_VisAtt);    

    if(!(parameters.SingleSeparationLayer))
    { //LAPPD HERE

      double cumulative_thickness = 0.;
      double layer_pos = 0.;

      LAPPD_materials_vec = absorber.GetLAPPDMaterials();
      if (LAPPD_materials_vec.size() != parameters.LAPPD_layers.size()) cout << "ERROR: VECTOR SIZES DO NOT MATCH" << endl;

      // Set VisibilityAttributes for the LAPPD
      G4VisAttributes* lappd_VisAtt = new G4VisAttributes(brown);  // color
      lappd_VisAtt->SetVisibility(parameters.lappdVisibility);      // absorbers are always visible
      lappd_VisAtt->SetForceWireframe(parameters.wireFrame);

      for (int layer_num=0; layer_num<parameters.LAPPD_layers.size(); layer_num++)
      {
        // layer_pos = parameters.cell_separator_position - 0.5*parameters.separation_thickness + cumulative_thickness + 0.5*parameters.LAPPD_layers.at(layer_num);
        layer_pos = - 0.5*parameters.separation_thickness + cumulative_thickness + 0.5*parameters.LAPPD_layers.at(layer_num);
        
        sname.str("");
        sname << "LAPPD_S_LAYER_"
              << to_string(layer_num)
              << "_"
              << absorber.GetName()
              << "_" << moduleName;
        LAPPD_S.push_back( new G4Box (sname.str().c_str(),
                                      0.5 * absorber.GetSizeX(),
                                      0.5 * absorber.GetSizeY(),
                                      0.5 * parameters.LAPPD_layers.at(layer_num)));
        sname.str("");
        // create the logical volume
        sname << "LAPPD_LV_LAYER_"
              << to_string(layer_num)
              << "_" 
              << absorber.GetName()
              << "_" << moduleName;     
        if(Parameters::Instance()->main_config.verbosity > 0) G4cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << G4endl << LAPPD_materials_vec.at(layer_num) << G4endl;
        LAPPD_LV.push_back(new G4LogicalVolume (LAPPD_S.at(layer_num),
                                                LAPPD_materials_vec.at(layer_num),
                                                sname.str().c_str()) );
        // Set VisibilityAttributes for the LAPPD
        LAPPD_LV[layer_num]->SetVisAttributes(lappd_VisAtt);

        sname.str("");
        // create the physical volume (i.e. place it in space)
        sname << "LAPPD_PV_LAYER_"
              << to_string(layer_num)
              << "_"          
              << absorber.GetName()
              << "_" << moduleName;;
        LAPPD_PV.push_back(new G4PVPlacement (0,
                                              G4ThreeVector (0, 0, layer_pos),
                                              LAPPD_LV.at(layer_num),
                                              sname.str().c_str(),
                                              esr_LV,
                                              false, 
                                              0, 
                                              parameters.checkOverlaps));
        sname.str("");      


        
        cout << ">>> INFO: added LAPPD layer " << layer_num << " with thickness " << parameters.LAPPD_layers.at(layer_num) << " (cumulative: " << cumulative_thickness << ")" << endl;
        cout << ">>> INFO: its position is at " << layer_pos << endl;    

        cumulative_thickness += parameters.LAPPD_layers.at(layer_num);
      }
    }
  } 
  //--------------------------------------- //



  // build the hole parametrization
  // this is where the cell concept loses meaning, we are using it
  // to separate front and back, and to arrage crystals, but no more than that
  // position vectors. how to distinguish front and back cells?
  // we know how many sections
  // we can check cellz position vs abs separation

  std::vector<Cell> CellCollection[2];
  std::vector<G4double> hole_x[2];
  std::vector<G4double> hole_y[2];
  std::vector<G4double> hole_z[2];
  // dimensions NEED to be the same for elements in the same cell, so...
  std::vector<G4double> hole_dx[2];
  std::vector<G4double> hole_dy[2];
  std::vector<G4double> hole_dz[2];

  std::vector<G4double> crystal_x[2];
  std::vector<G4double> crystal_y[2];
  std::vector<G4double> crystal_z[2];
  // dimensions NEED to be the same for elements in the same cell, so...
  std::vector<G4double> crystal_dx[2];
  std::vector<G4double> crystal_dy[2];
  std::vector<G4double> crystal_dz[2];

  std::cout << "sections = " << sections << std::endl;

  // run on all cells in the absorber
  for(int iCell = 0; iCell < absorber.GetNumberOfCells() ; iCell++)
  {
    Cell cell = absorber.GetCell(iCell);
    if(sections == 1)
    {
      // all cells in the same collection, use front
      CellCollection[0].push_back(cell);
    }
    else
    {
      // distinguish front and back cells
      if((module_pos_z + absorber.GetPositionZ() + cell.GetPositionZ()) > separation_z_position)
      {
        // z is positive -> back cell
        CellCollection[1].push_back(cell);
      }
      else
      {
        // z is negative -> front cell
        CellCollection[0].push_back(cell);
      }
    }
  }

  for(int iSec = 0; iSec < sections; iSec++) // run on sections
  {
    // dimensions NEED to be the same for elements in the same cell, so
    // we take only the first
    hole_dx[iSec].push_back(CellCollection[iSec][0].GetLayer(0).hole_dimension_x[0]);
    hole_dy[iSec].push_back(CellCollection[iSec][0].GetLayer(0).hole_dimension_y[0]);
    hole_dz[iSec].push_back(CellCollection[iSec][0].GetLayer(0).hole_dimension_z[0]);
    crystal_dx[iSec].push_back(CellCollection[iSec][0].GetLayer(0).crystal_dimension_x[0]);
    crystal_dy[iSec].push_back(CellCollection[iSec][0].GetLayer(0).crystal_dimension_y[0]);
    crystal_dz[iSec].push_back(CellCollection[iSec][0].GetLayer(0).crystal_dimension_z[0]);

    for(int iCell = 0; iCell < CellCollection[iSec].size() ; iCell++)
    {
      // now run on all cell elements
      Cell cell = CellCollection[iSec][iCell];
      G4int NumberOfLayers = cell.GetNumberOfLayers();
      for(int iLay = 0; iLay < NumberOfLayers ; iLay++)
      {
        //get the row
        row_t layer = cell.GetLayer(iLay);
        // loop on elements of the row, aka the columns, aka the individual holes (elements)
        int NumberOfLayerElements = layer.elements;
        for(int iEl = 0; iEl < NumberOfLayerElements ; iEl++)
        {
          hole_x[iSec].push_back(cell.GetPositionX() + layer.hole_pos_x[iEl] );
          hole_y[iSec].push_back(cell.GetPositionY() + layer.hole_pos_y[iEl]);
          hole_z[iSec].push_back(0);
          crystal_x[iSec].push_back(cell.GetPositionX() + layer.hole_pos_x[iEl]);
          crystal_y[iSec].push_back(cell.GetPositionY() + layer.hole_pos_y[iEl]);
          crystal_z[iSec].push_back(cell.GetPositionZ() + layer.crystal_pos_z[iEl]);

        }
      }
    }
  }

  // // debug
  // for(int iSec = 0; iSec < sections; iSec++) // run on sections
  // {
  //   for(int iP = 0 ; iP < x[iSec].size(); iP++)
  //   {
  //     std::cout << iSec << " " << iP << " "
  //               << x[iSec][iP] << " "
  //               << y[iSec][iP] << " "
  //               << z[iSec][iP] << " "
  //               << std::endl;
  //   }
  // }
  // std::cout << "---------------------------------" << std::endl;
  // std::cout << "---------------------------------" << std::endl;
  // for(int iSec = 0; iSec < sections; iSec++) // run on sections
  // {
  //   for(int iP = 0 ; iP < dx[iSec].size(); iP++)
  //   {
  //     std::cout << iSec << " " << iP << " "
  //               << dx[iSec][iP] << " "
  //               << dy[iSec][iP] << " "
  //               << dz[iSec][iP] << " "
  //               << std::endl;
  //   }
  // }


  // now save cells (needed for next steps)
  // run on cells of the absorber
  for(int iCell = 0; iCell < absorber.GetNumberOfCells() ; iCell++)
  {
    Cell cell = absorber.GetCell(iCell);

    // save the cell in the cell TTree in output file
    if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
    {
      CreateTree::Instance()->cellID       = cellID;
      CreateTree::Instance()->cellMaterial = 0;
      CreateTree::Instance()->cellType     = ecal_position;
      CreateTree::Instance()->cellX        = module_pos_x + absorber.GetPositionX() + cell.GetPositionX();
      CreateTree::Instance()->cellY        = module_pos_y + absorber.GetPositionY() + cell.GetPositionY();
      CreateTree::Instance()->cellZ        = module_pos_z + absorber.GetPositionZ() + cell.GetPositionZ();
      CreateTree::Instance()->cellDX       = cell.GetSizeX();
      CreateTree::Instance()->cellDY       = cell.GetSizeY();
      CreateTree::Instance()->cellDZ       = cell.GetSizeZ();
      CreateTree::Instance()->absorberSeparationZ = separation_z_position;
      CreateTree::Instance()->cellSections = sections;
      CreateTree::Instance()->cells->Fill();
    }
  }



  // create the hole volumes (with all the other volumes inside)
  // needs to be done separately for front and back section, if they exist


  // place them with parametrisation
  for(int iSec = 0; iSec < sections; iSec++) // run on sections
  {
    // create the hole volume
    G4LogicalVolume *hole_LV = MakeHole(CellCollection[iSec][0]);


    // G4PVPlacement *hole_PV;
    // create the physical volume (i.e. place the hole in the absorber)
    // sname << "hole_PV";
          // << cell.GetName() << "_"
          // << absorber.GetName() << "_"
          // << moduleName;


    G4VPVParameterisation* hole_parametrization =
                                 new HoleParametrization (
                                                            hole_x[iSec],hole_y[iSec],hole_z[iSec],
                                                            G4ThreeVector( hole_dx[iSec][0], hole_dy[iSec][0], hole_dz[iSec][0] ),
                                                            crystal_x[iSec],crystal_y[iSec],crystal_z[iSec],
                                                            G4ThreeVector( crystal_dx[iSec][0], crystal_dy[iSec][0], crystal_dz[iSec][0] ),
                                                            ecal_position,CellCollection[iSec][0].GetCrystalMaterialIdentifier());
    // hole_PV  =  new G4PVPlacement (0,
    //                                G4ThreeVector (0,
    //                                               0,
    //                                               0),
    //                                hole_LV,
    //                                sname.str().c_str(),
    //                                abs_section_LV[iSec],
    //                                false, 0, parameters.checkOverlaps);
    G4PVParameterised* hole_PV = new G4PVParameterised("hole",       // their name
                          hole_LV,   // their logical volume
                          abs_section_LV[iSec],       // Mother logical volume
                          kUndefined,          // Are placed along this axis
                          hole_x[iSec].size(),    // Number of chambers
                          hole_parametrization,    // The parametrisation
                          false); // checking overlaps
    sname.str("");



    // add an optical surface between hole and absorber
    // create a skin surface, unless the absorber is made of air (special case to study individual fibers in air)
    sname << "surf_HoleToAbs" ;
          // << cell.GetName() << "_"
          // << absorber.GetName() << "_"
          // << moduleName;
    G4OpticalSurface* surf_HoleToAbs = new G4OpticalSurface(sname.str().c_str());
    surf_HoleToAbs->SetType(dielectric_metal);
    surf_HoleToAbs->SetFinish(ground);
    surf_HoleToAbs->SetModel(unified);
    surf_HoleToAbs->SetSigmaAlpha(absorber.GetSigmaAlpha()); // irrelevant
    surf_HoleToAbs->SetMaterialPropertiesTable(absorber.GetOpticalSurface());
    G4LogicalBorderSurface* l_surf_HoleToAbs = NULL;

    if(absorber.GetMaterialIdentifier() != 9) // 9) is Air
    {
      l_surf_HoleToAbs = new G4LogicalBorderSurface(sname.str().c_str(),
                                                    hole_PV,
                                                    abs_section_PV[iSec],
                                                    surf_HoleToAbs);
    }
    sname.str("");


  }



  // // run on cells of the absorber
  // for(int iCell = 0; iCell < absorber.GetNumberOfCells() ; iCell++)
  // {
  //   // cells are a collection of holes in an absorber
  //   // first, get the cell specifications
  //   Cell cell = absorber.GetCell(iCell);
  //
  //   // save the cell in the cell TTree in output file
  //   // if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
  //   // {
  //   //   CreateTree::Instance()->cellID       = cellID;
  //   //   CreateTree::Instance()->cellMaterial = 0;
  //   //   CreateTree::Instance()->cellType     = ecal_position;
  //   //   CreateTree::Instance()->cellX        = module_pos_x + absorber.GetPositionX() + cell.GetPositionX();
  //   //   CreateTree::Instance()->cellY        = module_pos_y + absorber.GetPositionY() + cell.GetPositionY();
  //   //   CreateTree::Instance()->cellZ        = module_pos_z + absorber.GetPositionZ() + cell.GetPositionZ();
  //   //   CreateTree::Instance()->cellDX       = cell.GetSizeX();
  //   //   CreateTree::Instance()->cellDY       = cell.GetSizeY();
  //   //   CreateTree::Instance()->cellDZ       = cell.GetSizeZ();
  //   //   CreateTree::Instance()->absorberSeparationZ = separation_z_position;
  //   //   CreateTree::Instance()->cellSections = sections;
  //   //   CreateTree::Instance()->cells->Fill();
  //   // }
  //
  //
  //
  //   // G4cout << "cell x " << module_pos_x[iMod][jMod] + absorber.GetPositionX() + cell.GetPositionX() << G4endl;
  //   // G4cout << "cell y " << module_pos_y[iMod][jMod] + absorber.GetPositionY() + cell.GetPositionY() << G4endl;
  //   // G4cout << "cell z " << module_pos_z[iMod][jMod] + absorber.GetPositionZ() + cell.GetPositionZ() << G4endl;
  //
  //
  //
  //
  //
  //   // now for each cell run on rows
  //   // get the number of rows of this cell
  //   G4int NumberOfLayers = cell.GetNumberOfLayers();
  //   for(int iLay = 0; iLay < NumberOfLayers ; iLay++)
  //   {
  //     //get the row
  //     row_t layer = cell.GetLayer(iLay);
  //
  //     // loop on elements of the row, aka the columns, aka the individual holes (elements)
  //     int NumberOfLayerElements = layer.elements;
  //     for(int iEl = 0; iEl < NumberOfLayerElements ; iEl++)
  //     {
  //
  //
  //     }// end of loop on columns
  //   }//end loop on rows
  //   cellID++;
  // }// end loop on cells
  // moduleID++;
  // absorberID++;
  // }//end loop on modules y
  // }//end loop on modules x

  // CreateTree::Instance()->WriteStructure();

}


G4LogicalVolume* Spacal::MakeHole(Cell cell )
{

  // now, each element is made by a hole filled by a crystal.
  // however, we need also to specify the air gaps
  // the thickness of air gaps is defined by the user (airLayer key in cfg file)
  // the hole transvers dimensions are already made larger than the crystal, by a
  // factor 2*airLayer, to fill the crystal plus the air gap. But the hole cannot
  // be longer than the absorber (and hopefully the user didn't define a crystal, or the 3 cry, shorter
  // than a hole).
  // Here, we also need to create physical volumes for the air gaps, because we will need
  // optical surfaces (and also to put optical grease/glue if needed).

  // NEW FEATURE, 1/9/2020: possibility to have round fibers and/or round holes.
  // In principle the solution is simple, it just requires to use other constructors
  // on the basis of some configuration keys. However, 4 combinations now become possible:
  // 1) Square hole, square fiber
  // 2) Round hole, round fiber
  // 3) Square hole, round fiber
  // 4) Round hole, square fiber
  // The first 2 options are the most straighforward and also probably the only ones that will be used
  // Option 3 could be taken into consideration in some application, and option 4 it's
  // indeed quite crazy but hey, let's implement it anyway, you never know what your boss
  // might ask you to test...
  // So in fact option 4 introduces a problem: we need an air_hole to be able to properly
  // define the lateral optical surface crystal to air gap. The dimension of this air_hole
  // was defined just as the average between crystal and hole (x,y) dims -->.
  // air_hole = 0.5(cry+hole)
  // Now, this might not work for option 4 (simple drawings will show you that). 
  // Furthermore, the possible options for hole type and crystals (and especially, since 
  // the crystal could have cladding), are implemented

  row_t layer = cell.GetLayer(0);

  //--------------------------------------- //
  // the hole
  G4VSolid        *hole_S;
  G4LogicalVolume *hole_LV;
  G4PVPlacement   *hole_PV;
  // create the hole shape
  sname << "hole_S";
  int holeShape = cell.GetHoleShape();
  int crystalShape = cell.GetCrystalShape();
  bool crystalCladding = cell.GetCrystalCladding();
  if(holeShape == 0) // squared holes
  {
    hole_S = new G4Box(sname.str().c_str(),
                       0.5 * layer.hole_dimension_x[0],
                       0.5 * layer.hole_dimension_y[0],
                       0.5 * layer.hole_dimension_z[0]);
  }
  else if (holeShape == 1) // round holes
  {
    if(crystalShape == 0)
    {
      G4double diagonal = sqrt(pow(layer.crystal_dimension_x[0],2) + pow(layer.crystal_dimension_y[0],2));
      hole_S = new G4Tubs(sname.str().c_str(),
                            0.,
                            0.5 * (layer.hole_dimension_x[0] + diagonal) / 2.0,
                            0.5 * layer.hole_dimension_z[0],
                            0. * deg, 360. * deg);

    }
    else 
    {
      hole_S = new G4Tubs(sname.str().c_str(),
                        0.,
                        0.5 * layer.hole_dimension_x[0],
                        0.5 * layer.hole_dimension_z[0],
                        0. * deg, 360. * deg);
    }
  }
  else 
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid hole shape specifier " << holeShape << G4endl;
    exit(-1);
  }
  
  sname.str("");
  // create the logical volume
  sname << "hole_LV";
  hole_LV   = new G4LogicalVolume (hole_S,
                                   AirMaterial,     // always air
                                   sname.str().c_str()) ;
  sname.str("");
  // set the Visualization properties of the hole
  G4VisAttributes* Hole_VisAtt = new G4VisAttributes(air);
  Hole_VisAtt->SetVisibility(parameters.holeVisibility);
  Hole_VisAtt->SetForceWireframe(true);
  hole_LV->SetVisAttributes(Hole_VisAtt);

  // do another volume just like the hole, it's needed for the PV, for the crystal to hole surface
  // is it ok? and what about the exit?
  G4VSolid        *air_hole_S;
  G4LogicalVolume *air_hole_LV;
  G4PVPlacement   *air_hole_PV;
  // create the hole shape
  sname << "air_hole_S";
  if(holeShape == 0) // squared holes
  {
    air_hole_S    = new G4Box (sname.str().c_str(),
                             0.5 * (layer.hole_dimension_x[0] + layer.crystal_dimension_x[0])/2.0,
                             0.5 * (layer.hole_dimension_y[0] + layer.crystal_dimension_y[0])/2.0,
                             0.5 * layer.crystal_dimension_z[0]);
  }
  else if (holeShape == 1)
  {
    // special case if hole is round and crystal is rectangular 
    if(crystalShape == 0)
    {
      G4double diagonal = sqrt(pow(layer.crystal_dimension_x[0],2) + pow(layer.crystal_dimension_y[0],2));
      air_hole_S = new G4Tubs(sname.str().c_str(),
                            0.,
                            0.5 * (layer.hole_dimension_x[0] + diagonal) / 2.0,
                            0.5 * layer.hole_dimension_z[0],
                            0. * deg, 360. * deg);

    }
    else 
    {
      air_hole_S = new G4Tubs(sname.str().c_str(),
                            0.,
                            0.5 * (layer.hole_dimension_x[0] + layer.crystal_dimension_x[0]) / 2.0,
                            0.5 * layer.hole_dimension_z[0],
                            0. * deg, 360. * deg);
    }
  }
  else
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid hole shape specifier " << holeShape << G4endl;
    exit(-1);
  }

  sname.str("");
  // create the logical volume
  sname << "air_hole_LV";
  air_hole_LV   = new G4LogicalVolume (air_hole_S,
                                       AirMaterial,     // always air
                                       sname.str().c_str()) ;
  sname.str("");
  // create the physical volume (i.e. place the hole in the absorber)
  sname << "air_hole_PV";
  air_hole_PV  =  new G4PVPlacement (0,
                                 G4ThreeVector (layer.crystal_pos_x[0],
                                                layer.crystal_pos_y[0],
                                                layer.crystal_pos_z[0]),
                                 air_hole_LV,
                                 sname.str().c_str(),
                                 hole_LV,
                                 false, 0, parameters.checkOverlaps);
  //
  G4VisAttributes* air_Hole_VisAtt = new G4VisAttributes(blue);
  air_Hole_VisAtt->SetVisibility(false); // always hide this volume, it's just instrumental to the optical surfaces, not needed for other purpose
  air_Hole_VisAtt->SetForceWireframe(true);
  air_hole_LV->SetVisAttributes(air_Hole_VisAtt);

  // end of the hole
  //--------------------------------------- //



  //--------------------------------------- //
  // create the crystal
  G4VSolid *crystal_S;
  G4LogicalVolume *crystal_LV;
  G4PVPlacement *crystal_PV;
  G4VisAttributes *Cry_VisAtt;

  // create the shape
  // needs to distinguish the case of cladding 
  // in this case the crystal volume becomes the outer cladding volume,
  // then two more volumes are placed inside: the inner cladding and the core,
  // with the latter that will get the name "crystal" (necessary for stepping action)

  G4VSolid        *inner_cladding_S;
  G4LogicalVolume *inner_cladding_LV;
  G4PVPlacement   *inner_cladding_PV;

  G4VSolid        *core_S;
  G4LogicalVolume *core_LV;
  G4PVPlacement   *core_PV;

  // if(crystalCladding && crystalShape == 1) // cladding only for round fibers!
  if(crystalCladding) // cladding 
  {
    sname << "outer_cladding";
  }
  else 
  {
    sname << "crystal";
  }
  
  if(crystalShape == 0) // squared 
  {
    crystal_S = new G4Box(sname.str().c_str(),
                        0.5 * layer.crystal_dimension_x[0],
                        0.5 * layer.crystal_dimension_y[0],
                        0.5 * layer.crystal_dimension_z[0]);
  }
  else if (crystalShape == 1) // round 
  {
    crystal_S = new G4Tubs(sname.str().c_str(),
                        0.,
                        0.5 * layer.crystal_dimension_x[0],
                        0.5 * layer.crystal_dimension_z[0],
                        0. * deg, 360. * deg);
  }
  else 
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid crystal shape specifier " << crystalShape << G4endl;
    exit(-1);
  }
  sname.str("");

  
  
  if(crystalCladding) // cladding 
  {
    sname << "outer_cladding";
    crystal_LV = new G4LogicalVolume(crystal_S,
                                    cell.GetOuterCladdingMaterial(),
                                    sname.str().c_str());
    sname.str("");
    
    if(crystalShape == 0) // for squared 
    {
      sname << "inner_cladding";
      // make the two other volumes 
      inner_cladding_S = new G4Box(sname.str().c_str(),
                          0.5 * cell.GetInnerCladdingDiameter(),
                          0.5 * cell.GetInnerCladdingDiameter(),
                          0.5 * layer.crystal_dimension_z[0]);
      inner_cladding_LV = new G4LogicalVolume(inner_cladding_S,
                                      cell.GetInnerCladdingMaterial(),
                                      sname.str().c_str());
      inner_cladding_PV = new G4PVPlacement(0,
                                    G4ThreeVector(0, 
                                                  0,
                                                  0),
                                    inner_cladding_LV,
                                    sname.str().c_str(),
                                    crystal_LV, //mother volume is crystal_LV
                                    false, 0, parameters.checkOverlaps);
      sname.str("");
      sname << "crystal";
      core_S = new G4Box(sname.str().c_str(),
                          0.5 * cell.GetCoreDiameter(),
                          0.5 * cell.GetCoreDiameter(),
                          0.5 * layer.crystal_dimension_z[0]);
      core_LV = new G4LogicalVolume(core_S,
                                      cell.GetCrystalMaterial(),
                                      sname.str().c_str());
      core_PV = new G4PVPlacement(0,
                                    G4ThreeVector(0, 
                                                  0,
                                                  0),
                                    core_LV,
                                    sname.str().c_str(),
                                    inner_cladding_LV, //mother volume is crystal_LV
                                    false, 0, parameters.checkOverlaps);
      sname.str("");

    }
    else if(crystalShape == 1) // for round fibers
    {
      
      sname << "inner_cladding";
      // make the two other volumes 
      inner_cladding_S = new G4Tubs(sname.str().c_str(),
                          0.,
                          0.5 * cell.GetInnerCladdingDiameter(),
                          0.5 * layer.crystal_dimension_z[0],
                          0. * deg, 360. * deg);      
      inner_cladding_LV = new G4LogicalVolume(inner_cladding_S,
                                      cell.GetInnerCladdingMaterial(),
                                      sname.str().c_str());
      inner_cladding_PV = new G4PVPlacement(0,
                                    G4ThreeVector(0, 
                                                  0,
                                                  0),
                                    inner_cladding_LV,
                                    sname.str().c_str(),
                                    crystal_LV, //mother volume is crystal_LV
                                    false, 0, parameters.checkOverlaps);
      sname.str("");
      sname << "crystal";
      core_S = new G4Tubs(sname.str().c_str(),
                          0.,
                          0.5 * cell.GetCoreDiameter(),
                          0.5 * layer.crystal_dimension_z[0],
                          0. * deg, 360. * deg);
      core_LV = new G4LogicalVolume(core_S,
                                      cell.GetCrystalMaterial(),
                                      sname.str().c_str());
      core_PV = new G4PVPlacement(0,
                                    G4ThreeVector(0, 
                                                  0,
                                                  0),
                                    core_LV,
                                    sname.str().c_str(),
                                    inner_cladding_LV, //mother volume is crystal_LV
                                    false, 0, parameters.checkOverlaps);
      sname.str("");
    }
    else // should never happen, already checked. code here just for simmetry 
    {
      G4cerr << "<Spacal::InitializeMaterials()>: Invalid crystal shape specifier " << crystalShape << G4endl;
      exit(-1);
    }

    // set vis attriubutes
    G4VisAttributes *inCl_VisAtt = new G4VisAttributes(cell.GetInnerCladdingColor());
    inCl_VisAtt->SetVisibility(parameters.crystalsVisibility);
    inCl_VisAtt->SetForceWireframe(parameters.wireFrame);
    inner_cladding_LV->SetVisAttributes(inCl_VisAtt);
    
    G4VisAttributes *coreCl_VisAtt = new G4VisAttributes(cell.GetCrystalColor());
    coreCl_VisAtt->SetVisibility(parameters.crystalsVisibility);
    coreCl_VisAtt->SetForceWireframe(parameters.wireFrame);
    core_LV->SetVisAttributes(coreCl_VisAtt);

    Cry_VisAtt = new G4VisAttributes(cell.GetOuterCladdingColor());

    // also define internal optical surfaces
    // one (back anf forth, so in fact 2) can be defined here:
    // surface between core and inner cladding 
    // the other one needs to be defined in the mother method, 
    // because the constructor needs the physical volume of "crystal" and 
    // this will be placed only afterwards  
    if(parameters.crystal_lateral_depolishing > 0)
    {
      sname << "latSurf_core_to_inner_1";
      G4OpticalSurface* opt_surf_core_to_inner_1 = new G4OpticalSurface(sname.str().c_str());
      opt_surf_core_to_inner_1->SetType(dielectric_dielectric);
      opt_surf_core_to_inner_1->SetFinish(ground);
      opt_surf_core_to_inner_1->SetModel(unified);
      opt_surf_core_to_inner_1->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
      opt_surf_core_to_inner_1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_border_core_to_inner_1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    core_PV,
                                                                    inner_cladding_PV,
                                                                    opt_surf_core_to_inner_1);
      sname.str("");
      //
      sname << "latSurf_core_to_inner_2";
      G4OpticalSurface* opt_surf_core_to_inner_2 = new G4OpticalSurface(sname.str().c_str());
      opt_surf_core_to_inner_2->SetType(dielectric_dielectric);
      opt_surf_core_to_inner_2->SetFinish(ground);
      opt_surf_core_to_inner_2->SetModel(unified);
      opt_surf_core_to_inner_2->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
      opt_surf_core_to_inner_2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_border_core_to_inner_2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    inner_cladding_PV,
                                                                    core_PV,
                                                                    opt_surf_core_to_inner_2);
      sname.str("");
      //
    }
  }
  else 
  {
    sname << "crystal";
    crystal_LV = new G4LogicalVolume(crystal_S,
                                   cell.GetCrystalMaterial(),
                                   sname.str().c_str());
    sname.str("");
    Cry_VisAtt = new G4VisAttributes(cell.GetCrystalColor());
  }
  
  sname.str("");
   
  Cry_VisAtt->SetVisibility(parameters.crystalsVisibility);
  Cry_VisAtt->SetForceWireframe(parameters.wireFrame);
  crystal_LV->SetVisAttributes(Cry_VisAtt);
  // end of the crystal
  //--------------------------------------- //

  // if(crystalCladding && crystalShape == 1) // cladding only for round fibers!
  if(crystalCladding) // cladding
  {
    sname << "outer_cladding";
  }
  else 
  {
    sname << "crystal";
  }
  crystal_PV = new G4PVPlacement(0,
                                 G4ThreeVector(0, // crystal always in the center of air_hole!
                                               0,
                                               0),
                                 crystal_LV,
                                 sname.str().c_str(),
                                 air_hole_LV, //mother volume is air_hole_LV
                                 false, 0, parameters.checkOverlaps);
  sname.str("");

  // only if cladding and round fibers 
  // optical surface from inner to outer cladding 
  // if(crystalCladding && crystalShape == 1)
  if(crystalCladding)
  {
    if(parameters.crystal_lateral_depolishing > 0)
    {
      sname << "latSurf_inner_to_outer_1";
      G4OpticalSurface* opt_surf_inner_to_outer_1 = new G4OpticalSurface(sname.str().c_str());
      opt_surf_inner_to_outer_1->SetType(dielectric_dielectric);
      opt_surf_inner_to_outer_1->SetFinish(ground);
      opt_surf_inner_to_outer_1->SetModel(unified);
      opt_surf_inner_to_outer_1->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
      opt_surf_inner_to_outer_1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_border_inner_to_outer_1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    inner_cladding_PV,
                                                                    crystal_PV,
                                                                    opt_surf_inner_to_outer_1);
      sname.str("");
      //
      sname << "latSurf_inner_to_outer_2";
      G4OpticalSurface* opt_surf_inner_to_outer_2 = new G4OpticalSurface(sname.str().c_str());
      opt_surf_inner_to_outer_2->SetType(dielectric_dielectric);
      opt_surf_inner_to_outer_2->SetFinish(ground);
      opt_surf_inner_to_outer_2->SetModel(unified);
      opt_surf_inner_to_outer_2->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
      opt_surf_inner_to_outer_2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    crystal_PV,
                                                                    inner_cladding_PV,
                                                                    opt_surf_inner_to_outer_2);
      sname.str("");
      //
    }
  }
    

  // opt surf if depolishing of exit faces
  if(parameters.crystal_exit_depolishing > 0)
  {
    //
    // surf to cry air gap positive z interface
    sname << "surf_CryToPositiveInterface1";
    G4OpticalSurface* surf_CryToPositiveInterface1 = new G4OpticalSurface(sname.str().c_str());
    surf_CryToPositiveInterface1->SetType(dielectric_dielectric);
    surf_CryToPositiveInterface1->SetFinish(ground);
    surf_CryToPositiveInterface1->SetModel(unified);
    surf_CryToPositiveInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
    surf_CryToPositiveInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  crystal_PV,
                                                                  Gap_Abs_Interface_Positive_Z_PV,
                                                                  surf_CryToPositiveInterface1);
    sname.str("");
    //
    sname << "surf_CryToPositiveInterface2";
    G4OpticalSurface* surf_CryToPositiveInterface2 = new G4OpticalSurface(sname.str().c_str());
    surf_CryToPositiveInterface2->SetType(dielectric_dielectric);
    surf_CryToPositiveInterface2->SetFinish(ground);
    surf_CryToPositiveInterface2->SetModel(unified);
    surf_CryToPositiveInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
    surf_CryToPositiveInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  Gap_Abs_Interface_Positive_Z_PV,
                                                                  crystal_PV,
                                                                  surf_CryToPositiveInterface2);
    sname.str("");



    // surf to cry air gap negative z interface
    sname << "surf_CryToNegativeInterface1";
    G4OpticalSurface* surf_CryToNegativeInterface1 = new G4OpticalSurface(sname.str().c_str());
    surf_CryToNegativeInterface1->SetType(dielectric_dielectric);
    surf_CryToNegativeInterface1->SetFinish(ground);
    surf_CryToNegativeInterface1->SetModel(unified);
    surf_CryToNegativeInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
    surf_CryToNegativeInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  crystal_PV,
                                                                  Gap_Abs_Interface_Negative_Z_PV,
                                                                  surf_CryToNegativeInterface1);
    sname.str("");
    //
    sname << "surf_CryToNegativeInterface2";
    G4OpticalSurface* surf_CryToNegativeInterface2 = new G4OpticalSurface(sname.str().c_str());
    surf_CryToNegativeInterface2->SetType(dielectric_dielectric);
    surf_CryToNegativeInterface2->SetFinish(ground);
    surf_CryToNegativeInterface2->SetModel(unified);
    surf_CryToNegativeInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
    surf_CryToNegativeInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  Gap_Abs_Interface_Negative_Z_PV,
                                                                  crystal_PV,
                                                                  surf_CryToNegativeInterface2);
    sname.str("");

    // if(crystalCladding && crystalShape == 1) // cladding only for round fibers!
    if(crystalCladding ) // cladding 
    {
      // INNER CLADDING TO INTERFACES
      //
      // surf inner cladding to air gap positive z interface
      sname << "surf_InnerCladdingToPositiveInterface1";
      G4OpticalSurface* surf_InnerCladdingToPositiveInterface1 = new G4OpticalSurface(sname.str().c_str());
      surf_InnerCladdingToPositiveInterface1->SetType(dielectric_dielectric);
      surf_InnerCladdingToPositiveInterface1->SetFinish(ground);
      surf_InnerCladdingToPositiveInterface1->SetModel(unified);
      surf_InnerCladdingToPositiveInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_InnerCladdingToPositiveInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_InnerCladdingborder1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    inner_cladding_PV,
                                                                    Gap_Abs_Interface_Positive_Z_PV,
                                                                    surf_InnerCladdingToPositiveInterface1);
      sname.str("");
      //
      sname << "surf_InnerCladdingToPositiveInterface2";
      G4OpticalSurface* surf_InnerCladdingToPositiveInterface2 = new G4OpticalSurface(sname.str().c_str());
      surf_InnerCladdingToPositiveInterface2->SetType(dielectric_dielectric);
      surf_InnerCladdingToPositiveInterface2->SetFinish(ground);
      surf_InnerCladdingToPositiveInterface2->SetModel(unified);
      surf_InnerCladdingToPositiveInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_InnerCladdingToPositiveInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_InnerCladdingborder2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    Gap_Abs_Interface_Positive_Z_PV,
                                                                    inner_cladding_PV,
                                                                    surf_InnerCladdingToPositiveInterface2);
      sname.str("");
  
  
  
      // surf inner cladding to air gap negative z interface
      sname << "surf_InnerCladdingToNegativeInterface1";
      G4OpticalSurface* surf_InnerCladdingToNegativeInterface1 = new G4OpticalSurface(sname.str().c_str());
      surf_InnerCladdingToNegativeInterface1->SetType(dielectric_dielectric);
      surf_InnerCladdingToNegativeInterface1->SetFinish(ground);
      surf_InnerCladdingToNegativeInterface1->SetModel(unified);
      surf_InnerCladdingToNegativeInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_InnerCladdingToNegativeInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_InnerCladdingborder3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    inner_cladding_PV,
                                                                    Gap_Abs_Interface_Negative_Z_PV,
                                                                    surf_InnerCladdingToNegativeInterface1);
      sname.str("");
      //
      sname << "surf_InnerCladdingToNegativeInterface2";
      G4OpticalSurface* surf_InnerCladdingToNegativeInterface2 = new G4OpticalSurface(sname.str().c_str());
      surf_InnerCladdingToNegativeInterface2->SetType(dielectric_dielectric);
      surf_InnerCladdingToNegativeInterface2->SetFinish(ground);
      surf_InnerCladdingToNegativeInterface2->SetModel(unified);
      surf_InnerCladdingToNegativeInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_InnerCladdingToNegativeInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_InnerCladdingborder4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    Gap_Abs_Interface_Negative_Z_PV,
                                                                    inner_cladding_PV,
                                                                    surf_InnerCladdingToNegativeInterface2);
      

      // CORE TO INTERFACES
      //
      // surf core to air gap positive z interface
      sname << "surf_CoreToPositiveInterface1";
      G4OpticalSurface* surf_CoreToPositiveInterface1 = new G4OpticalSurface(sname.str().c_str());
      surf_CoreToPositiveInterface1->SetType(dielectric_dielectric);
      surf_CoreToPositiveInterface1->SetFinish(ground);
      surf_CoreToPositiveInterface1->SetModel(unified);
      surf_CoreToPositiveInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_CoreToPositiveInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_Coreborder1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    core_PV,
                                                                    Gap_Abs_Interface_Positive_Z_PV,
                                                                    surf_CoreToPositiveInterface1);
      sname.str("");
      //
      sname << "surf_CoreToPositiveInterface2";
      G4OpticalSurface* surf_CoreToPositiveInterface2 = new G4OpticalSurface(sname.str().c_str());
      surf_CoreToPositiveInterface2->SetType(dielectric_dielectric);
      surf_CoreToPositiveInterface2->SetFinish(ground);
      surf_CoreToPositiveInterface2->SetModel(unified);
      surf_CoreToPositiveInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_CoreToPositiveInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_Coreborder2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    Gap_Abs_Interface_Positive_Z_PV,
                                                                    core_PV,
                                                                    surf_CoreToPositiveInterface2);
      sname.str("");
  
  
  
      // surf outer cladding to air gap negative z interface
      sname << "surf_CoreToNegativeInterface1";
      G4OpticalSurface* surf_CoreToNegativeInterface1 = new G4OpticalSurface(sname.str().c_str());
      surf_CoreToNegativeInterface1->SetType(dielectric_dielectric);
      surf_CoreToNegativeInterface1->SetFinish(ground);
      surf_CoreToNegativeInterface1->SetModel(unified);
      surf_CoreToNegativeInterface1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_CoreToNegativeInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_Coreborder3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    core_PV,
                                                                    Gap_Abs_Interface_Negative_Z_PV,
                                                                    surf_CoreToNegativeInterface1);
      sname.str("");
      //
      sname << "surf_CoreToNegativeInterface2";
      G4OpticalSurface* surf_CoreToNegativeInterface2 = new G4OpticalSurface(sname.str().c_str());
      surf_CoreToNegativeInterface2->SetType(dielectric_dielectric);
      surf_CoreToNegativeInterface2->SetFinish(ground);
      surf_CoreToNegativeInterface2->SetModel(unified);
      surf_CoreToNegativeInterface2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
      surf_CoreToNegativeInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
      G4LogicalBorderSurface* l_Coreborder4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                    Gap_Abs_Interface_Negative_Z_PV,
                                                                    core_PV,
                                                                    surf_CoreToNegativeInterface2);
    }

  }








  // add air gap on internal side, only if there is an internal side...
  if(cell.GetAirLayer() > 0)
  {
    G4double intGapPosition;
    if(cell.GetPositionZ() != 0)
    {
      if(cell.GetPositionZ() >  0) // then the exit must be in positive position
      {
        // extGapPosition = + (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
        intGapPosition = - (cell.GetHoleDimensionZ()/2.0 - cell.GetGapSize() / 2.0 );
      }
      else  // then the exit must be in negative position
      {
        // extGapPosition = - (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
        intGapPosition = + (cell.GetHoleDimensionZ()/2.0 - cell.GetGapSize() / 2.0 );
      }


      //--------------------------------------- //
      // the internal air gap
      G4VSolid        *air_gap_int_S;
      G4LogicalVolume *air_gap_int_LV;
      G4PVPlacement   *air_gap_int_PV;

      // create the shape
      // section as the hole, thickness as the airLayer
      sname << "air_gap_int";
      
      if (holeShape == 0) // squared holes
      {
        air_gap_int_S = new G4Box(sname.str().c_str(),
                                  0.5 * layer.hole_dimension_x[0],
                                  0.5 * layer.hole_dimension_y[0],
                                  0.5 * cell.GetGapSize());
      }
      else if (holeShape == 1)
      {
        air_gap_int_S = new G4Tubs(sname.str().c_str(),
                                   0.,
                                   0.5 * layer.hole_dimension_x[0],
                                   0.5 * cell.GetGapSize(),
                                   0. * deg, 360. * deg);
      }
      else
      {
        G4cerr << "<Spacal::InitializeMaterials()>: Invalid hole shape specifier " << holeShape << G4endl;
        exit(-1);
      }
      sname.str("");
      // create the logical volume
      sname << "air_gap_int";
      air_gap_int_LV  = new G4LogicalVolume (air_gap_int_S,
                                              cell.GetIntGapMaterial(),
                                              sname.str().c_str()) ;
      sname.str("");
      // create the physical volume (i.e. place it in space)
      sname << "air_gap_int";
      air_gap_int_PV = new G4PVPlacement  (0,
                                            G4ThreeVector (0,0,
                                            intGapPosition),
                                            air_gap_int_LV,
                                            sname.str().c_str(),
                                            hole_LV,    //mother is hole
                                            false, 0, parameters.checkOverlaps);
      sname.str("");
      G4VisAttributes* air_gap_int_VisAtt = new G4VisAttributes(cell.GetIntGapMaterialColor());
      air_gap_int_VisAtt->SetVisibility(parameters.crystalsVisibility);
      air_gap_int_VisAtt->SetForceWireframe(parameters.wireFrame);
      air_gap_int_LV->SetVisAttributes(air_gap_int_VisAtt);
      // end of air gap "back"
      //--------------------------------------- //

      // optical surfaces
      // define the optical surf that separates cells longitudinally
      // sepation can be
      // 0) nothing = air
      // 1) aluminization
      // 2) esr


      if(parameters.cell_separation_type == 0)  // air gap
      {
      // CASE 0)
      // if exit surface is polished ----> nothing to do
      // if exit surface depolished  ----> set a depolished surface between the crystal and the air layer just added
        if(parameters.crystal_exit_depolishing > 0)
        {
          //
          // surf to int air gap
          sname << "intSurf1";
          G4OpticalSurface* int_opt_surf1 = new G4OpticalSurface(sname.str().c_str());
          int_opt_surf1->SetType(dielectric_dielectric);
          int_opt_surf1->SetFinish(ground);
          int_opt_surf1->SetModel(unified);
          int_opt_surf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
          int_opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
          G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                        crystal_PV,
                                                                        air_gap_int_PV,
                                                                        int_opt_surf1);
          sname.str("");
          //
          sname << "intSurf2";
          G4OpticalSurface* int_opt_surf2 = new G4OpticalSurface(sname.str().c_str());
          int_opt_surf2->SetType(dielectric_dielectric);
          int_opt_surf2->SetFinish(ground);
          int_opt_surf2->SetModel(unified);
          int_opt_surf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
          int_opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
          G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                        air_gap_int_PV,
                                                                        crystal_PV,
                                                                        int_opt_surf2);
          sname.str("");

          if(crystalCladding) // cladding 
          // if(crystalCladding && crystalShape == 1) // cladding only for round fibers!
          {
            //
            // inner cladding surf to int air gap
            sname << "innerSurf1";
            G4OpticalSurface* innerSurf1 = new G4OpticalSurface(sname.str().c_str());
            innerSurf1->SetType(dielectric_dielectric);
            innerSurf1->SetFinish(ground);
            innerSurf1->SetModel(unified);
            innerSurf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            innerSurf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_Innerborder1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          inner_cladding_PV,
                                                                          air_gap_int_PV,
                                                                          innerSurf1);
            sname.str("");
            //
            sname << "innerSurf2";
            G4OpticalSurface* innerSurf2 = new G4OpticalSurface(sname.str().c_str());
            innerSurf2->SetType(dielectric_dielectric);
            innerSurf2->SetFinish(ground);
            innerSurf2->SetModel(unified);
            innerSurf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            innerSurf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_Innerborder2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          air_gap_int_PV,
                                                                          inner_cladding_PV,
                                                                          innerSurf2);
            sname.str("");
  
            //
            // core surf to int air gap
            sname << "coreSurf1";
            G4OpticalSurface* coreSurf1 = new G4OpticalSurface(sname.str().c_str());
            coreSurf1->SetType(dielectric_dielectric);
            coreSurf1->SetFinish(ground);
            coreSurf1->SetModel(unified);
            coreSurf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            coreSurf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_core1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          core_PV,
                                                                          air_gap_int_PV,
                                                                          coreSurf1);
            sname.str("");
            //
            sname << "coreSurf2";
            G4OpticalSurface* coreSurf2 = new G4OpticalSurface(sname.str().c_str());
            coreSurf2->SetType(dielectric_dielectric);
            coreSurf2->SetFinish(ground);
            coreSurf2->SetModel(unified);
            coreSurf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            coreSurf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_core2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          air_gap_int_PV,
                                                                          crystal_PV,
                                                                          coreSurf2);
            sname.str("");
          }
          
        }
      }

      // CASE 1)
      // aluminization
      if(parameters.cell_separation_type == 1)  // aluminization
      {
        // aluminization means polished front painted
        //
        sname << "extSurf1";
        G4OpticalSurface* opt_surf1 = new G4OpticalSurface(sname.str().c_str());
        opt_surf1->SetType(dielectric_dielectric);
        opt_surf1->SetFinish(polishedfrontpainted);
        opt_surf1->SetModel(unified);
        // opt_surf1->SetSigmaAlpha(crystal_exit_depolishing);
        opt_surf1->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
        G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                       crystal_PV,
                                                                       air_gap_int_PV,
                                                                       opt_surf1);
        sname.str("");
        //
        // if(crystalCladding && crystalShape == 1) // cladding only for round fibers!
        if(crystalCladding ) // cladding
        {
          //
          sname << "alumInnerCladding";
          G4OpticalSurface* alumInnerCladding = new G4OpticalSurface(sname.str().c_str());
          alumInnerCladding->SetType(dielectric_dielectric);
          alumInnerCladding->SetFinish(polishedfrontpainted);
          alumInnerCladding->SetModel(unified);
          alumInnerCladding->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
          G4LogicalBorderSurface* l_alumInnerCladdingborder1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                         inner_cladding_PV,
                                                                         air_gap_int_PV,
                                                                         alumInnerCladding);
          sname.str("");

          //
          sname << "alumCoreCladding";
          G4OpticalSurface* alumCoreCladding = new G4OpticalSurface(sname.str().c_str());
          alumCoreCladding->SetType(dielectric_dielectric);
          alumCoreCladding->SetFinish(polishedfrontpainted);
          alumCoreCladding->SetModel(unified);
          alumCoreCladding->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
          G4LogicalBorderSurface* l_alumCoreCladdingborder1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                         core_PV,
                                                                         air_gap_int_PV,
                                                                         alumCoreCladding);
          sname.str("");

        }
      }



      if(parameters.cell_separation_type == 2)  // esr
      {
        // CASE 2)
        // esr
        // if this is the case, there MUST be a esr volume, and there must be an air gap
        // first, the depolishing, if it's there
        if(parameters.crystal_exit_depolishing > 0)
        {
          //
          // surf to int air gap
          sname << "intSurf1";
          G4OpticalSurface* int_opt_surf1 = new G4OpticalSurface(sname.str().c_str());
          int_opt_surf1->SetType(dielectric_dielectric);
          int_opt_surf1->SetFinish(ground);
          int_opt_surf1->SetModel(unified);
          int_opt_surf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
          int_opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
          G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                        crystal_PV,
                                                                        air_gap_int_PV,
                                                                        int_opt_surf1);
          sname.str("");
          //
          sname << "intSurf2";
          G4OpticalSurface* int_opt_surf2 = new G4OpticalSurface(sname.str().c_str());
          int_opt_surf2->SetType(dielectric_dielectric);
          int_opt_surf2->SetFinish(ground);
          int_opt_surf2->SetModel(unified);
          int_opt_surf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
          int_opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
          G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                        air_gap_int_PV,
                                                                        crystal_PV,
                                                                        int_opt_surf2);
          sname.str("");

          if(crystalCladding) // cladding 
          {
            // inner cladding surf to int air gap
            sname << "InnerCladdingintSurf1";
            G4OpticalSurface* InnerCladdingintSurf1 = new G4OpticalSurface(sname.str().c_str());
            InnerCladdingintSurf1->SetType(dielectric_dielectric);
            InnerCladdingintSurf1->SetFinish(ground);
            InnerCladdingintSurf1->SetModel(unified);
            InnerCladdingintSurf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            InnerCladdingintSurf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_InnerCladding_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          inner_cladding_PV,
                                                                          air_gap_int_PV,
                                                                          InnerCladdingintSurf1);
            sname.str("");
            //
            sname << "InnerCladdingintSurf2";
            G4OpticalSurface* InnerCladdingintSurf2 = new G4OpticalSurface(sname.str().c_str());
            InnerCladdingintSurf2->SetType(dielectric_dielectric);
            InnerCladdingintSurf2->SetFinish(ground);
            InnerCladdingintSurf2->SetModel(unified);
            InnerCladdingintSurf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            InnerCladdingintSurf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_InnerCladding_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          air_gap_int_PV,
                                                                          inner_cladding_PV,
                                                                          InnerCladdingintSurf2);
            sname.str("");

            // core surf to int air gap
            sname << "CoreintSurf1";
            G4OpticalSurface* CoreintSurf1 = new G4OpticalSurface(sname.str().c_str());
            CoreintSurf1->SetType(dielectric_dielectric);
            CoreintSurf1->SetFinish(ground);
            CoreintSurf1->SetModel(unified);
            CoreintSurf1->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            CoreintSurf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_core_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          core_PV,
                                                                          air_gap_int_PV,
                                                                          CoreintSurf1);
            sname.str("");
            //
            sname << "CoreintSurf2";
            G4OpticalSurface* CoreintSurf2 = new G4OpticalSurface(sname.str().c_str());
            CoreintSurf2->SetType(dielectric_dielectric);
            CoreintSurf2->SetFinish(ground);
            CoreintSurf2->SetModel(unified);
            CoreintSurf2->SetSigmaAlpha(parameters.crystal_exit_depolishing);
            CoreintSurf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
            G4LogicalBorderSurface* l_core_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                          air_gap_int_PV,
                                                                          core_PV,
                                                                          CoreintSurf2);
            sname.str("");
          }


        }

        // then, the esr surface between air_gap_int_PV and esr volume (just in one direction)
        // here two cases: 
        // 1) the separation is single type, so optical surface between air_gap_int_PV and esr_PV
        // 2) the separation is multiple type, so optical surface is between air_gap_int_PV and the 
        //    adjacent volume in the separation stack
        // BUT given some limitation in doubleing point numbers, reflecting on G4 volume position precision
        // there could be misalignement so it's better to always define the air_gap_int_PV - esr_PV 
        // surface, then also define the other two

        sname << "esrSurf1";
        G4OpticalSurface* esrSurf1 = new G4OpticalSurface(sname.str().c_str());
        esrSurf1->SetType(dielectric_metal);
        esrSurf1->SetFinish(polished);
        esrSurf1->SetModel(unified);
        esrSurf1->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
        G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                       air_gap_int_PV,
                                                                       esr_PV,
                                                                       esrSurf1);
        sname.str("");
        
        if(!(parameters.SingleSeparationLayer)) 
        {
          if(cell.GetPositionZ() >  0) // then the LAPPD adjacent module is the last one in the stack (positive position)
          {
            sname << "esrSurf2";
            G4OpticalSurface* esrSurf2 = new G4OpticalSurface(sname.str().c_str());
            esrSurf2->SetType(dielectric_metal);
            esrSurf2->SetFinish(polished);
            esrSurf2->SetModel(unified);
            esrSurf2->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
            G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                           air_gap_int_PV,
                                                                           LAPPD_PV[LAPPD_PV.size()-1],
                                                                           esrSurf2);
            sname.str("");
          }
          else  // then the LAPPD adjacent module is the first one in the stack (negative position)
          {
            sname << "esrSurf3";
            G4OpticalSurface* esrSurf3 = new G4OpticalSurface(sname.str().c_str());
            esrSurf3->SetType(dielectric_metal);
            esrSurf3->SetFinish(polished);
            esrSurf3->SetModel(unified);
            esrSurf3->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
            G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                           air_gap_int_PV,
                                                                           LAPPD_PV[0],
                                                                           esrSurf3);
            sname.str("");
          }
        }
        
      }
    }
  } // end if airLayer


  // optical surfaces
  // optical photons are going to meet, in the cell
  // 1) surface crystal - hole when they exit from the side  -> surface could be
  //    depolished, so need for surface dielectric_dielectric with depolishing,
  //    if user requires it - back and forth!
  // 2) surface  crystal - ext gap material                  -> surface could be
  //    depolished, so need for surface dielectric_dielectric with depolishing,
  //    if user requires it - back and forth!
  // 3) surface  crystal - int gap material                  -> surface could be
  //    depolished, so need for surface dielectric_dielectric with depolishing,
  //    if user requires it - back and forth!
  // 4) surface int - reflector   -> if there is a reflector between modules
  // 5) surface int               -> if there is aluminization


  //---------------------------------------- //
  // 1) surface crystal - hole (crystal lateral surf)
  if(parameters.crystal_lateral_depolishing > 0)
  {
    sname << "latSurf1";
    G4OpticalSurface* opt_surf1 = new G4OpticalSurface(sname.str().c_str());
    opt_surf1->SetType(dielectric_dielectric);
    opt_surf1->SetFinish(ground);
    opt_surf1->SetModel(unified);
    opt_surf1->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
    opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  crystal_PV,
                                                                  air_hole_PV,
                                                                  opt_surf1);
    sname.str("");
    //
    sname << "latSurf2";
    G4OpticalSurface* opt_surf2 = new G4OpticalSurface(sname.str().c_str());
    opt_surf2->SetType(dielectric_dielectric);
    opt_surf2->SetFinish(ground);
    opt_surf2->SetModel(unified);
    opt_surf2->SetSigmaAlpha(parameters.crystal_lateral_depolishing);
    opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
    G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                  air_hole_PV,
                                                                  crystal_PV,
                                                                  opt_surf2);
    sname.str("");
    //
  }
  // end of 1)
  //---------------------------------------- //
  return hole_LV;
}

void Spacal::InitializeMaterials ()
{
  // define materials
  AirMaterial = MyMaterials::Air();
  AirKillerMaterial = MyMaterials::AirKiller();
  WoMaterial = MyMaterials::Air () ;

  InterfaceContainerMaterial = NULL;
  if(parameters.readoutType == 0 || parameters.readoutType == 3  || parameters.readoutType == 7 || parameters.readoutType == 5)
  {
    InterfaceContainerMaterial = AirMaterial;
  }
  else
  {
    InterfaceContainerMaterial = AirKillerMaterial;
  }
  
  if(parameters.readoutType == 4)
  {
    // set material to the same as GapAbsToInterfaceMaterial and GapInterfaceToReadoutMaterial
    if ( parameters.gap_abs_interface_material == 1 || parameters.gap_interface_readout_material == 1 )
    {
      InterfaceContainerMaterial = MyMaterials::OpticalGrease () ;
    }
    else 
    {
      InterfaceContainerMaterial = MyMaterials::Air () ;
    }
  }

  if(parameters.readoutType == 8)
  {
    // set interface container to stainless
    InterfaceContainerMaterial = MyMaterials::StainlessSteel();

  }

  InterfaceMaterial  = NULL;
  if      ( parameters.cone_material == 0 ) InterfaceMaterial = MyMaterials::Air () ;
  else if ( parameters.cone_material == 1 ) InterfaceMaterial = MyMaterials::PLEX (parameters.plex_abs_length_scale_factor) ;
  else
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid cone material specifier " << parameters.cone_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "Interface material: "<< InterfaceMaterial << G4endl ;

  GapAbsToInterfaceMaterial  = NULL;
  if      ( parameters.gap_abs_interface_material == 0 ) GapAbsToInterfaceMaterial = MyMaterials::Air () ;
  else if ( parameters.gap_abs_interface_material == 1 ) GapAbsToInterfaceMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid GapAbsToInterfaceMaterial material specifier " << parameters.gap_abs_interface_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "GapAbsToInterfaceMaterial material: "<< GapAbsToInterfaceMaterial << G4endl ;

  GapInterfaceToReadoutMaterial  = NULL;
  if      ( parameters.gap_interface_readout_material == 0 ) GapInterfaceToReadoutMaterial = MyMaterials::Air () ;
  else if ( parameters.gap_interface_readout_material == 1 ) GapInterfaceToReadoutMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<Spacal::InitializeMaterials()>: Invalid GapInterfaceToReadoutMaterial material specifier " << parameters.gap_interface_readout_material << G4endl ;
    exit (-1) ;
  }
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "GapInterfaceToReadoutMaterial material: "<< GapInterfaceToReadoutMaterial << G4endl ;

  PLEXMaterial = NULL ;
  PLEXMaterial = MyMaterials::PLEX (parameters.plex_abs_length_scale_factor) ; 
  //
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "PLEX. material: "<< PLEXMaterial << G4endl ;
  //
  PVCMaterial = NULL ;
  PVCMaterial = MyMaterials::PVC () ;
  //
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "PVC. material: "<< PVCMaterial << G4endl ;

}


// spacal 2019 test beam readout and 2020
void Spacal::ConstructLightGuides(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4double lguide_edge,G4RotationMatrix *rot, G4double shiftCAD)
{
  // now the light guides
  // Light guide structure
  // Truncated of cone, cut on the sides
  // goes from 14.4 x14.4 mm2 (on the fibers) to 5 mm radius on the PMT
  // in reality, the 14.4 mm is reduced to make space for a small air gap
  // first, do a simple truncated cone
  G4LogicalVolume* LguideLV;
  G4PVPlacement *LightGuide_PV;
  recallCAD cone3cm;
  std::stringstream sname;
  // G4double lguide_edge = 14.4 *mm;
  G4double airGap = 0.05*mm; // 50 microns on each side
  G4double PLEX_depth = InterfaceSizeZ - 2.0*parameters.gapSize;
  G4double PMT_radius = 5 *mm;
  G4VSolid* coneS = new G4Cons("aCone",0,PMT_radius,0, 0.5*sqrt(2)*(lguide_edge-2.0*airGap) ,0.5*PLEX_depth,0.0 * deg,360.0* deg);
  // then prepare a hollow box, section of the hole 20x20 mm2
  if(parameters.readoutType != 7)
  {
    G4Box* innerBoxS = new G4Box ("innerBoxS", 0.5*(lguide_edge-2.0*airGap), 0.5*(lguide_edge-2.0*airGap), 0.5*PLEX_depth) ;
    G4Box* outerBoxS = new G4Box ("outerBoxS", 0.5*sqrt(2)*(lguide_edge-2.0*airGap), 0.5*sqrt(2)*(lguide_edge-2.0*airGap), 0.5*PLEX_depth) ;
    G4VSolid* subtract = new G4SubtractionSolid("Hollow-Box", outerBoxS, innerBoxS,0,G4ThreeVector(0.,0.,0.));
    // subtract hollow box from truncated cone
    G4VSolid* coneSolid = new G4SubtractionSolid("coneSolid", coneS, subtract,0,  G4ThreeVector(0.,0.,0.));
    LguideLV = new G4LogicalVolume (coneSolid,InterfaceMaterial,"LguideLV");
  }
  else // tb 2020, 3cm LG fish tail
  {
    LguideLV = cone3cm.returnMeshedLV(parameters.light_guide_file, PLEXMaterial, "LguideLV"); // "LguideLV"
  }
 

  int nLGx = 3;
  int nLGy = 3;
  G4double lg_pitch = lguide_edge;
  G4double lg_x[3] = {-lg_pitch,0,lg_pitch};
  G4double lg_y[3] = {-lg_pitch,0,lg_pitch};
  for(int iLG = 0; iLG < nLGx; iLG++)
  {
    for(int jLG = 0; jLG < nLGy; jLG++)
    {
      sname << "LightGuide_"
            << iLG << "_" << jLG << "_"
            << moduleName
            << "_PV";
      if(parameters.readoutType == 7) //tb 2020, we use the variable lguide_edge to shift the CAD back in place
      {
        LightGuide_PV = new G4PVPlacement (rot, G4ThreeVector (lg_x[iLG],lg_y[jLG],shiftCAD), LguideLV, sname.str().c_str(), Interface_LV, false, 0, parameters.checkOverlaps);
      }
      else
      {
        LightGuide_PV = new G4PVPlacement (rot, G4ThreeVector (lg_x[iLG],lg_y[jLG],0), LguideLV, sname.str().c_str(), Interface_LV, false, 0, parameters.checkOverlaps);
      }
      
      sname.str("");

      G4VisAttributes* LG_VisAtt = new G4VisAttributes(green);  // color
      LG_VisAtt->SetVisibility(parameters.lgVisibility);
      LG_VisAtt->SetForceWireframe(parameters.wireFrame);
      LguideLV->SetVisAttributes(LG_VisAtt);

      if(parameters.esr_on_cones || parameters.teflon_on_cones)
      {
        sname << "esrSurf_LightGuide_"
              << iLG << "_" << jLG << "_"
              << moduleName;
        G4OpticalSurface* esrSurf = new G4OpticalSurface(sname.str().c_str());
        esrSurf->SetType(dielectric_metal);
        
        if(parameters.esr_on_cones && parameters.teflon_on_cones)
        {
          G4cout << "WARNING: both esr_on_cones = 1 and teflon_on_cones = 1 are given! Setting to just teflon_on_cones = 1" << G4endl ;
          // then there is no need to change the parameter value, given the ordering of 
          // next if statements the teflon one will prevail, if both are true
        }
        if(parameters.esr_on_cones)
        {
          esrSurf->SetFinish(polished);
          esrSurf->SetMaterialPropertiesTable(MyMaterials::ESR(parameters.esrTransmittance));
        }
        if(parameters.teflon_on_cones)
        {
          esrSurf->SetFinish(ground);
          esrSurf->SetMaterialPropertiesTable(MyMaterials::Teflon());
        }
        
        
        esrSurf->SetModel(unified);
        
        G4LogicalBorderSurface* l_border = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      LightGuide_PV,
                                                                      Interface_PV,
                                                                      esrSurf);
        sname.str("");

      }
    }
  }
}

// void DetectorConstruction::JustAirGap(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4double lguide_edge,G4RotationMatrix *rot)
// {
//   std::stringstream sname;
//
// }



void Spacal::ConstructPMTs(G4LogicalVolume *Readout_LV, std::string moduleName,G4double lg_pitch,G4RotationMatrix *rot)
{
  std::stringstream sname;
  G4double PMT_radius = 5*mm;
  G4double PMT_length = ReadoutSizeZ;
  G4VSolid * TubeS = new G4Tubs ("TubeS", 0., PMT_radius, 0.5*PMT_length, 0.*deg, 360.*deg);
  G4LogicalVolume * PMTLV = new G4LogicalVolume (TubeS,PLEXMaterial, "PMTLV");

  int nLGx = 3;
  int nLGy = 3;
  // G4double lg_pitch = 14.4 * mm;
  G4double lg_x[3] = {-lg_pitch,0,lg_pitch};
  G4double lg_y[3] = {-lg_pitch,0,lg_pitch};

  for(int iLG = 0; iLG < nLGx; iLG++)
  {
    for(int jLG = 0; jLG < nLGy; jLG++)
    {
      int pmt_number = iLG*nLGy + jLG;
      sname << "PMT_"
            << pmt_number
            << "_"
            << moduleName;
      new G4PVPlacement (rot, G4ThreeVector (lg_x[iLG],lg_y[jLG],0), PMTLV, sname.str().c_str(), Readout_LV, false, 0, parameters.checkOverlaps);
      sname.str("");
      G4VisAttributes* LG_VisAtt = new G4VisAttributes(grey);  // color
      LG_VisAtt->SetVisibility(parameters.lgVisibility);
      LG_VisAtt->SetForceWireframe(parameters.wireFrame);
      PMTLV->SetVisAttributes(LG_VisAtt);
    }
  }

}


void Spacal::ConstructOneLargePMT(G4LogicalVolume *Readout_LV, std::string moduleName,G4double PMT_radius,G4RotationMatrix *rot)
{
  std::stringstream sname;
  // G4double PMT_radius = 5*mm;
  G4double PMT_length = ReadoutSizeZ;
  G4VSolid * pmtBox = new G4Box ("pmtBox", PMT_radius, PMT_radius, 0.5*PMT_length) ;
  // G4VSolid * TubeS = new G4Tubs ("TubeS", 0., PMT_radius, 0.5*PMT_length, 0.*deg, 360.*deg);
  G4LogicalVolume * PMTLV = new G4LogicalVolume (pmtBox,PLEXMaterial, "PMTLV");
  sname << "PMT_"
        << 0
        << "_"
        << moduleName;
  new G4PVPlacement (rot, G4ThreeVector (00,0), PMTLV, sname.str().c_str(), Readout_LV, false, 0, parameters.checkOverlaps);
  sname.str("");
  G4VisAttributes* LG_VisAtt = new G4VisAttributes(grey);  // color
  LG_VisAtt->SetVisibility(parameters.lgVisibility);
  LG_VisAtt->SetForceWireframe(parameters.wireFrame);
  PMTLV->SetVisAttributes(LG_VisAtt);
}



void Spacal::ConstructClearFiber(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV,  G4PVPlacement *Gap_Interface_PV,  G4PVPlacement *Gap_Interface_Readout_PV, std::string moduleName)
{
  G4double fiber_length = InterfaceSizeZ - 2.0*parameters.gapSize;
  //using the same construction method used for the crystals and the absorbers
  for(int iCell = 0; iCell < absorber.GetNumberOfCells() ; iCell++)
  {
     Cell cell = absorber.GetCell(iCell);
     G4int NumberOfLayers = cell.GetNumberOfLayers();
     for(int iLay = 0; iLay < NumberOfLayers ; iLay++)
     {
      std::stringstream sname;
      row_t layer = cell.GetLayer(iLay);
      int NumberOfLayerElements = layer.elements;
      for(int i = 0; i < NumberOfLayerElements ; i++)
      {
        G4double x_dim = layer.crystal_dimension_x[i]; // square geometry only!
        // in the x-y position, we have to include all the information of the absorber-cell-layer structure
        G4double x_pos = absorber.GetPositionX()+cell.GetPositionX()+layer.hole_pos_x[i]+layer.crystal_pos_x[i];
        G4double y_pos = absorber.GetPositionY()+cell.GetPositionY()+layer.hole_pos_y[i]+layer.crystal_pos_y[i];
        G4double SecondCladdingSize = 0.03 * x_dim; // 3 % of the fiber diameter, which is the same as AbsSizeX
        G4double FirstCladdingSize = 0.06 * x_dim;
        // volumes
        // The definition of the Solid and logical are left intentionally inside the loops
        // allows to change the dimension of the crystal hole by hole, if selected in the config.
        G4VSolid *Interface_Z_2nd_Cladding_S;
        G4LogicalVolume *Interface_Z_2nd_Cladding_LV;
        G4PVPlacement   *Interface_Z_2nd_Cladding_PV;
        G4VSolid *Interface_Z_1st_Cladding_S;
        G4LogicalVolume *Interface_Z_1st_Cladding_LV;
        G4PVPlacement   *Interface_Z_1st_Cladding_PV;
        G4VSolid *Interface_Z_Core_S;
        G4LogicalVolume *Interface_Z_Core_LV;
        G4PVPlacement   *Interface_Z_Core_PV;

        // 2nd cladding of the optical fiber (external one)
        Interface_Z_2nd_Cladding_S = new G4Tubs (sname.str().c_str(), 0., x_dim / 2., 0.5 * fiber_length, 0., 2 * M_PI) ;
        sname.str("");
        sname << "Interface_Z_2nd_Cladding_"
              << moduleName
              << "_LV";
        Interface_Z_2nd_Cladding_LV = new G4LogicalVolume (Interface_Z_2nd_Cladding_S, PLEXMaterial, sname.str().c_str()) ;
        sname.str("");

        sname << "Interface_Z_2nd_Cladding_"
              << moduleName
              << "_PV";
        Interface_Z_2nd_Cladding_PV = new G4PVPlacement (0, G4ThreeVector (x_pos,
                                                                      y_pos,
                                                                      0),
                                                                      Interface_Z_2nd_Cladding_LV, sname.str().c_str(),
                                                                      Interface_LV, false, 0, parameters.checkOverlaps) ;
        sname.str("");
        // 1st cladding of the optical fiber (intermediate one)
        Interface_Z_1st_Cladding_S = new G4Tubs (sname.str().c_str(), 0., x_dim / 2. - SecondCladdingSize, 0.5 * fiber_length, 0., 2 * M_PI) ;
        sname.str("");
        sname << "Interface_Z_1st_Cladding_"
              << moduleName
              << "_LV";
        Interface_Z_1st_Cladding_LV = new G4LogicalVolume (Interface_Z_1st_Cladding_S, PLEXMaterial, sname.str().c_str()) ;
        sname.str("");

        sname << "Interface_Z_1st_Cladding_"
              << moduleName
              << "_PV";
        Interface_Z_1st_Cladding_PV = new G4PVPlacement (0, G4ThreeVector (x_pos,
                                                                      y_pos,
                                                                      0),
                                                                      Interface_Z_1st_Cladding_LV, sname.str().c_str(),
                                                                      Interface_Z_2nd_Cladding_LV, false, 0, parameters.checkOverlaps) ;
        sname.str("");
        // core of the optical fiber (intermediate one)
        Interface_Z_Core_S = new G4Tubs (sname.str().c_str(), 0., x_dim / 2. - FirstCladdingSize, 0.5 * fiber_length, 0., 2 * M_PI) ;
        sname.str("");
        sname << "Interface_Z_Core_"
              << moduleName
              << "_LV";
        Interface_Z_Core_LV = new G4LogicalVolume (Interface_Z_Core_S, MyMaterials::Polystyrene(parameters.user_lightyield,parameters.abs_length_scale_factor), sname.str().c_str()) ;
        sname.str("");

        sname << "Interface_Z_Core_"
              << moduleName
              << "_PV";
        Interface_Z_Core_PV = new G4PVPlacement (0, G4ThreeVector (x_pos,
                                                                      y_pos,
                                                                      0),
                                                                      Interface_Z_Core_LV, sname.str().c_str(),
                                                                      Interface_Z_1st_Cladding_LV, false, 0, parameters.checkOverlaps) ;
        // optical surfaces. We need to define:
        // gap to fiber, fiber to gap for all the fibers components and on both sides of the Clear Fiber
        sname.str("");

        sname << "fiberSurf"
              << i << "_"
              << cell.GetName() << "_"
              << absorber.GetName() << "_"
              << moduleName;
        G4OpticalSurface* surf_GapToClearFiber = new G4OpticalSurface(sname.str().c_str());
        surf_GapToClearFiber->SetType(dielectric_dielectric);
        surf_GapToClearFiber->SetFinish(polished);
        surf_GapToClearFiber->SetModel(unified);
        surf_GapToClearFiber->SetMaterialPropertiesTable(MyMaterials::clear_fiber_optical());
        // core to first cladding
        // incoming
        G4LogicalBorderSurface* l_borderGapCFiber1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_Core_PV,
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // outgoing
        G4LogicalBorderSurface* l_borderGapCFiber2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      Interface_Z_Core_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // first cladding to second
        // incoming
        G4LogicalBorderSurface* l_borderGapCFiber3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // outgoing
        G4LogicalBorderSurface* l_borderGapCFiber4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // second cladding to interface
        // incoming
        G4LogicalBorderSurface* l_borderGapCFiber5 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      Interface_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // outgoing
        G4LogicalBorderSurface* l_borderGapCFiber6 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_PV,
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // gap (toward absorber) to fiber
        // incoming
        // second cladding (external, FP)
        G4LogicalBorderSurface* l_borderGapCFiber7 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_PV,
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // first cladding (PMMA)
        G4LogicalBorderSurface* l_borderGapCFiber8 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_PV,
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // fiber core (PS)
        G4LogicalBorderSurface* l_borderGapCFiber9 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_PV,
                                                                      Interface_Z_Core_PV,
                                                                      surf_GapToClearFiber);
        // outgoing
        // second cladding (external, FP)
        G4LogicalBorderSurface* l_borderGapCFiber10 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      Gap_Interface_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // first cladding (PMMA)
        G4LogicalBorderSurface* l_borderGapCFiber11 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      Gap_Interface_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // fiber core (PS)
        G4LogicalBorderSurface* l_borderGapCFiber12 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_Core_PV,
                                                                      Gap_Interface_PV,
                                                                      surf_GapToClearFiber);
        // fiber to gap (toward readout)
        // incoming
        // second cladding (external, FP)
        G4LogicalBorderSurface* l_borderGapCFiber13 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      Gap_Interface_Readout_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // first cladding (PMMA)
        G4LogicalBorderSurface* l_borderGapCFiber14 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      Gap_Interface_Readout_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // fiber core (PS)
        G4LogicalBorderSurface* l_borderGapCFiber15 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Interface_Z_Core_PV,
                                                                      Gap_Interface_Readout_PV,
                                                                      surf_GapToClearFiber);
        // outgoing
        // second cladding (external, FP)
        G4LogicalBorderSurface* l_borderGapCFiber16 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_Readout_PV,
                                                                      Interface_Z_2nd_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // first cladding (PMMA)
        G4LogicalBorderSurface* l_borderGapCFiber17 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_Readout_PV,
                                                                      Interface_Z_1st_Cladding_PV,
                                                                      surf_GapToClearFiber);
        sname.str("");
        // fiber core (PS)
        G4LogicalBorderSurface* l_borderGapCFiber18 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                      Gap_Interface_Readout_PV,
                                                                      Interface_Z_Core_PV,
                                                                      surf_GapToClearFiber);
      }
    }
  }
}


void Spacal::ScanMap(TH2D *mapEcal)
{
  // set module block sizes
  G4double module_block_size_x = mapEcal->GetXaxis()->GetBinWidth(1) * cm;
  G4double module_block_size_y = mapEcal->GetYaxis()->GetBinWidth(1) * cm;

  xMapMin = +INFINITY;
  xMapMax = -INFINITY;
  yMapMin = +INFINITY;
  yMapMax = -INFINITY;

  // get nbisx and y
  int nBinsX = mapEcal->GetXaxis()->GetNbins();
  int nBinsY = mapEcal->GetYaxis()->GetNbins();
  // run on all bins, find the vector of types

  for(int i = 1 ; i < nBinsX+1 ; i++)
  {
    for(int j = 1 ; j < nBinsY+1 ; j++)
    {

      G4double x   = mapEcal->GetXaxis()->GetBinCenter(i) * cm;
      G4double y   = mapEcal->GetYaxis()->GetBinCenter(j) * cm;
      int number   = (int) mapEcal->GetBinContent(i,j);

      if(number == ecal_position)
      {
        if(( x - module_block_size_x/2.0) < xMapMin )
        {
          xMapMin = ( x - module_block_size_x/2.0);
        }
        if(( x + module_block_size_x/2.0) > xMapMax )
        {
          xMapMax = ( x + module_block_size_x/2.0);
        }
        if(( y - module_block_size_y/2.0) < yMapMin )
        {
          yMapMin = ( y - module_block_size_y/2.0);
        }
        if(( y + module_block_size_y/2.0) > yMapMax )
        {
          yMapMax = ( y + module_block_size_y/2.0);
        }
      }
    }
  }

  G4cout << ecal_position << " "
         << xMapMin << " "
         << xMapMax << " "
         << yMapMin << " "
         << yMapMax << " "
         << G4endl;
}
