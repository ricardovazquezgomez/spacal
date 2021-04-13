// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#include "Parameters.hh"

#include "CreateTree.hh"
#include "ConfigFile.hh"

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

using namespace CLHEP;

Parameters* Parameters::fInstance = NULL ;

Parameters::Parameters()
{
  if(fInstance)
  {
    return ;
  }
  this->fInstance = this ;
}

Parameters::~Parameters(){}


void Parameters::ReadMainConfig(std::string fileName)
{
  // this routines triggers the read of the main config file

  // invoke the routine FillStructure to actually read the file
  main_config = FillStructure(fileName);
  main_config.ecal_position = 0;
}


void Parameters::ReadModuleConfig(std::string fileName,int ecal_type,std::string container_volume)
{
  // this routines triggers the read a module config file

  // invoke the routine FillStructure to actually read the file
  Data_t data = FillStructure(fileName);
  data.ecal_position = ecal_type;
  data.volume_file = container_volume;
  // add structure to list of module config files
  module_config.push_back(data);
}

struct Data_t Parameters::FillStructure(std::string fileName)
{
  // prepare a Data_t structure
  Data_t data;
  // open input config file
  ConfigFile config(fileName.c_str());

  //---------------------------------------------------------//
  // IMPORT INDIVUDUAL KEYS
  //---------------------------------------------------------//
  // save config file name
  data.fileName               = fileName;
  config.readIntoVect(data.module_config_list,"modules_config");
  config.readIntoVect(data.ecal_type,"ecal_type");
  config.readIntoVect(data.container_volume,"container_volume");
  data.ecal_map_file          = config.read<std::string>("ecal_map_file","");
  data.ecal_map_histo         = config.read<std::string>("ecal_map_histo","");
  // read all keys (see Parameters.hh for key definitions)
  data.seed                   = config.read<long int>("seed",-1);
  data.checkOverlaps          = config.read<bool>("checkOverlaps",0) ;
  data.B_field_intensity      = config.read<double>("B_field_intensity",0) * tesla ;
  data.use_gps                = config.read<bool>("useGPS",false);
  data.input_filename         = config.read<TString>("input_file",TString("LHCbGaussSimulation_Option/input_fakeGamma_fixedE50GeV.root"));
  data.input_treename         = config.read<TString>("input_tree",TString("DecayTree"));
  data.skipEvents             = config.read<int>("skipEvents",0);
  data.Zshift                 = config.read<double>("deltaZ",-12500.); //shifted from LHCb coordinate to module coordinate
  data.gps_instructions_file  = config.read<std::string>("gps_instructions_file");
  data.defaultCut             = config.read<double> ("defaultCut",0.025) * mm;
  data.simulationType         = config.read<int>("simulationType",-1); // default to -1
  data.switchOnScintillation  = config.read<int>("switchOnScintillation",1);
  data.propagateScintillation = config.read<int>("propagateScintillation",1);
  data.switchOnCerenkov       = config.read<int>("switchOnCerenkov",1);
  data.propagateCerenkov      = config.read<int>("propagateCerenkov",1);
  data.saveAll                = config.read<bool>("saveAll",0);
  data.saveTree               = config.read<bool>("saveTree",0);
  data.saveShower             = config.read<bool>("saveShower",0);
  data.savePrimaries          = config.read<bool>("savePrimaries",0);
  data.savePhotons            = config.read<bool>("savePhotons",0);
  data.savePhotonAbsPoint     = config.read<bool>("savePhotonAbsPoint",0);
  data.saveStructure          = config.read<bool>("saveStructure",0);
  data.saveSummary            = config.read<bool>("saveSummary",0);
  data.saveSimple             = config.read<bool>("saveSimple",0);
  data.savePhotonGen          = config.read<bool>("savePhotonGen",0);
  data.saveEnergyPerModule    = config.read<bool>("saveEnergyPerModule",0);
  data.saveLAPPD              = config.read<bool>("saveLAPPD",1);
  data.primaries              = config.read<int>("primaries",1);
  data.opticalMaterial        = config.read<int>("opticalMaterial",-1);
  data.world_material         = config.read<int>("world_material",-1);
  data.logging_volume         = config.read<std::string>("logging_volume","PMT");
  data.pre_volume             = config.read<std::string>("pre_volume","");
  data.gapSize                = config.read<double>("gap_size",0.01); // default air gap is 10 micron
  data.abs_interface_extraGap = config.read<double>("abs_interface_extraGap",0); 
  data.worldVisibility        = config.read<int>("worldVisibility",0);
  data.caloVisibility         = config.read<int>("caloVisibility",0);
  data.moduleVisibility       = config.read<int>("moduleVisibility",0);
  data.holeVisibility         = config.read<int>("holeVisibility",1);
  data.crystalsVisibility     = config.read<int>("crystalsVisibility",1);
  data.interfaceVisibility    = config.read<int>("interfaceVisibility",1);
  data.gapsVisibility         = config.read<int>("gapsVisibility",1);
  data.readoutVisibility      = config.read<int>("readoutVisibility",1);
  data.absorberVisibility     = config.read<int>("absorberVisibility",0);
  data.esrVisibility          = config.read<int>("esrVisibility",1);
  data.lgVisibility           = config.read<int>("lgVisibility",1);
  data.lappdVisibility        = config.read<int>("lappdVisibility",0);
  data.containerVisibility    = config.read<bool>("containerVisibility",1);
  data.wireFrame              = config.read<int>("wireFrame",0);
  data.calorimeter_position   = config.read<int>("calorimeter_position",0);
  data.calorimeter_size_x     = config.read<double>("calorimeter_size_x",0) * mm;
  data.calorimeter_size_y     = config.read<double>("calorimeter_size_y",0) * mm;
  data.calorimeter_size_z     = config.read<double>("calorimeter_size_z",0) * mm;
  data.modules_nx             = config.read<int>("modules_nx",0);
  data.modules_ny             = config.read<int>("modules_ny",0);
  data.pipe_modules_nx        = config.read<int>("pipe_modules_nx",0);
  data.pipe_modules_ny        = config.read<int>("pipe_modules_ny",0);
  data.InterfaceSizeZ         = config.read<double>("interface_length",30);
  data.gap_abs_interface_material = config.read<int>("gap_abs_interface_material",0);
  data.gap_interface_readout_material = config.read<int>("gap_interface_readout_material",0);
  data.ReadoutSizeZ           = config.read<double>("readout_length",50);
  data.AbsName                = config.read<std::string>("absorber_name","");
  data.AbsSizeX               = config.read<G4double>("absorber_size_x",0);
  data.AbsSizeY               = config.read<G4double>("absorber_size_y",0);
  data.AbsSizeZ               = config.read<G4double>("absorber_size_z",0);
  data.AbsPositionX           = config.read<G4double>("absorber_pos_x",0);
  data.AbsPositionY           = config.read<G4double>("absorber_pos_y",0);
  data.AbsPositionZ           = config.read<G4double>("absorber_pos_z",0);
  data.AbsMaterial            = config.read<G4double>("absorber_material",0);
  data.absorber_reflectivity  = config.read<G4double>("absorber_reflectivity",0);
  data.absorber_specular_lobe = config.read<G4double>("absorber_specular_lobe",0);
  data.absorber_specular_spike= config.read<G4double>("absorber_specular_spike",0);
  data.absorber_backscatter   = config.read<G4double>("absorber_backscatter",0);
  data.absorber_sigma_alpha   = config.read<G4double>("absorber_sigma_alpha",0);
  data.W_fraction             = config.read<G4double>("W_fraction",0.722);
  data.cone_material          = config.read<int>("cone_material",0);
  data.esr_on_cones            = config.read<bool>("esr_on_cones",0);
  data.teflon_on_cones            = config.read<bool>("teflon_on_cones",0);
  data.esrTransmittance        = config.read<double>("esr_transmittance",0);
  data.cell_separation_type    = config.read<int>("cell_separation_type",0);
  data.cell_separator_position = config.read<G4double>("cell_separator_position",0);
  //data.separation_thickness = config.read<double>("separation_thickness",0);
  // data.separation_material  = config.read<int>("separation_material",0);
  data.readoutType = config.read<int>("readoutType",-1);
  data.abs_length_scale_factor = config.read<double>("abs_length_scale_factor",1);
  data.cladding_abs_length_scale_factor = config.read<double>("cladding_abs_length_scale_factor",0);
  data.plex_abs_length_scale_factor = config.read<double>("plex_abs_length_scale_factor",0);
  data.user_lightyield = config.read<double>("user_lightyield",-1);
  data.crystal_lateral_depolishing = config.read<double>("crystal_lateral_depolishing",0);
  data.crystal_exit_depolishing    = config.read<double>("crystal_exit_depolishing",0);
  data.moduleZshift    = config.read<double>("moduleZshift",0);
  data.crystal_inner_cladding_fraction    = config.read<double>("crystal_inner_cladding_fraction",0.02);
  data.crystal_outer_cladding_fraction    = config.read<double>("crystal_outer_cladding_fraction",0.02);
  data.event_time_cut    = config.read<double>("event_time_cut_ns",0);
  data.light_guide_file  = config.read<std::string>("light_guide_file","");
  data.verbosity         = config.read<int>("verbosity",0);
  data.esr_on_positive_exit = config.read<bool>("esr_on_positive_exit",0);
  data.esr_on_negative_exit = config.read<bool>("esr_on_negative_exit",0);
  // data.optical_calibration_module_type = config.read<int>("optical_calibration_module_type",-1);
  
  if(data.cladding_abs_length_scale_factor == 0)
  {
    G4cout << "No cladding abs length rescale factor specified, setting it equal to abs_length_scale_factor = " << data.abs_length_scale_factor << G4endl;
    data.cladding_abs_length_scale_factor = data.abs_length_scale_factor;
  }
  if(data.cell_separation_type != 2)
  {
    data.separation_thickness = 0;
  }
  if(data.user_lightyield < 0)
  {
    G4cout << "Using default light yields for " << data.fileName << G4endl;
  }
  else
  {
    G4cout << "Using user light yield = " << data.user_lightyield << " for " << data.fileName << G4endl;
  }
  if(data.abs_length_scale_factor == 1)
  {
    G4cout << "Using standard abs length for" << data.fileName << G4endl;
  }
  else
  {
    G4cout << "Using abs length scale factor = " << data.abs_length_scale_factor << " for" << data.fileName << G4endl;
  }
  //---------------------------------------------------------//



  //---------------------------------------------------------//
  // IMPORT VECTORS
  //---------------------------------------------------------//

  config.readIntoVect(data.CellNames,        "cell_name");
  config.readIntoVect(data.CellPositionX,    "cell_pos_x");
  config.readIntoVect(data.CellPositionY,    "cell_pos_y");
  config.readIntoVect(data.CellPositionZ,    "cell_pos_z");
  config.readIntoVect(data.CellXelements,    "cell_x_elements");
  config.readIntoVect(data.CellYelements,    "cell_y_elements");
  config.readIntoVect(data.CellCrystalSizeX ,"cell_crystal_size_x");
  config.readIntoVect(data.CellCrystalSizeY ,"cell_crystal_size_y");
  config.readIntoVect(data.CellCrystalSizeZ ,"cell_crystal_size_z");
  config.readIntoVect(data.CellCrystalPitchX,"cell_crystal_pitch_x");
  config.readIntoVect(data.CellCrystalPitchY,"cell_crystal_pitch_y");
  config.readIntoVect(data.CellMaterial,     "cell_crystal_material");
  config.readIntoVect(data.LAPPD_layers,     "LAPPD_layers_thickness");
  config.readIntoVect(data.LAPPD_materials,  "LAPPD_layers_materials");
  


  if (data.LAPPD_layers.size()>0) {
    int LAPPD_thickness = 0.;
    for (auto & element : data.LAPPD_layers) {
      int element_in_nm = (int) roundf(element * 1000000);
      LAPPD_thickness += element_in_nm ;
      // std::cout << setprecision(20);
      // std::cout << " ------------ --------- " << element_in_mu << " " << element << " " << data.LAPPD_thickness << std::endl;
    }
    data.separation_thickness = LAPPD_thickness / 1000000.0; // * micrometer;
    std::cout << setprecision(20) << "data.separation_thickness " << data.separation_thickness << std::endl;
    data.separation_material = 10; // set air as separation material, but in fact useless
    data.SingleSeparationLayer = false;
  } else {
    data.separation_thickness = config.read<double>("separation_thickness",0);
    data.separation_material  = config.read<int>("separation_material",0);
    data.SingleSeparationLayer = true;
  }

  // optional keys
  // cell shape 
  //
  bool holeShapeGiven = config.read<bool>("cell_hole_shape", 0);
  bool crystalShapeGiven = config.read<bool>("cell_crystal_shape", 0);
  bool crystalCladdingGiven = config.read<bool>("cell_crystal_cladding", 0);
  if(holeShapeGiven)
  {
    config.readIntoVect(data.CellHoleShape, "cell_hole_shape");
  }
  else
  {
    for(int i = 0; i < data.CellNames.size(); i++)
    {
      data.CellHoleShape.push_back(0);
    }
  }
  if(crystalShapeGiven)
  {
    config.readIntoVect(data.CellCrystalShape, "cell_crystal_shape");
  }
  else
  {
    for(int i = 0; i < data.CellNames.size(); i++)
    {
      data.CellCrystalShape.push_back(0);
    }
  }
  if(crystalCladdingGiven)
  {
    config.readIntoVect(data.CellCrystalCladding, "cell_crystal_cladding");
  }
  else
  {
    for(int i = 0; i < data.CellNames.size(); i++)
    {
      data.CellCrystalCladding.push_back(0);
    }
  }


  // air layer. if nothing is given, no air layer (so build 3 meaningless arrays)
  bool airLayerGiven = config.read<bool>("cell_air_layer",0);
  if(airLayerGiven)
  {
    config.readIntoVect(data.CellAirLayer,     "cell_air_layer");
    //internal air layer? if it's not given, it means both layers are external
    bool intAirLayerGiven = config.read<bool>("cell_int_gap_material",0);
    if(intAirLayerGiven) // /internal air layer? if it's not given, it means both layers are external
    {
      config.readIntoVect(data.CellIntGapMaterial,  "cell_int_gap_material");
    }
    else //if it's not given, it means both layers are external - FIXME: is this needed now?
    {
      config.readIntoVect(data.CellIntGapMaterial,  "cell_ext_gap_material");
    }
  }
  else
  {
    for(int i = 0; i < data.CellNames.size(); i++)
    {
      data.CellAirLayer      .push_back(0);
      data.CellIntGapMaterial.push_back(1);
    }
  }
  // staggering. if nothing is given, no staggering
  bool stagGiven = config.read<bool>("cell_staggering",0);
  if(stagGiven)
  {
    config.readIntoVect(data.CellStaggering,         "cell_staggering");
    config.readIntoVect(data.CellStaggeringAxis,     "cell_staggering_axis");
    config.readIntoVect(data.CellStaggeringSize  ,   "cell_staggering_size");
    config.readIntoVect(data.CellStaggeringParity,   "cell_staggering_parity");
    config.readIntoVect(data.CellStaggeringRemove,   "cell_staggering_remove");
  }
  else
  {
    for(int i = 0; i < data.CellNames.size(); i++)
    {
      data.CellStaggering.push_back(0);
      data.CellStaggeringAxis.push_back(0);
      data.CellStaggeringSize.push_back(0);
      data.CellStaggeringParity.push_back(0);
      data.CellStaggeringRemove.push_back(0);
    }
  }



  //---------------------------------------------------------//
  // CHECK VECTORS
  //---------------------------------------------------------//
  // check that all vectors have the same length
  // pedantic, embarrassing implementation, but saves a lot of time when debugging config files, so...
  bool cell_allEquals = true;
  if(data.CellNames.size() != data.CellPositionX.size())        cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellPositionX size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellPositionY.size())        cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellPositionY size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellPositionZ.size())        cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellPositionZ size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellMaterial.size())         cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellMaterial size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellXelements.size())        cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellXelements size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellYelements.size())        cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellYelements size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellAirLayer.size())         cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellAirLayer size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellIntGapMaterial.size())   cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellIntGapMaterial size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellStaggering.size())       cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellStaggering size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellStaggeringAxis.size())   cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellStaggeringAxis size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellStaggeringSize.size())   cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellStaggeringSize size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellStaggeringParity.size()) cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellStaggeringParity size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellStaggeringRemove.size()) cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellStaggeringRemove size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalSizeX .size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalSizeX size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalSizeY .size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalSizeY size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalSizeZ .size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalSizeZ size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalPitchX.size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalPitchX size!" << G4endl ;
    exit (-1) ;
  }




  if(data.CellNames.size() != data.CellHoleShape.size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellHoleShape size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalShape.size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalShape size!" << G4endl ;
    exit (-1) ;
  }

  if(data.CellNames.size() != data.CellCrystalCladding.size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match in " << data.fileName << G4endl ;
    G4cerr << "Wrong CellCrystalCladding size!" << G4endl ;
    exit (-1) ;
  }
  


  // check light guide file
  if(data.readoutType == 7 && data.light_guide_file == "" )
  {
    G4cerr << "<DetectorConstruction>: You chose readoutType = 7 but you didn't provide an OBJ file for the light guide! " << G4endl ;
    G4cerr << "Aborting!" << G4endl ;
    exit (-1) ;
  }


  // // find if the module is made of one or two sections
  // // std::vector<double> zCell;
  // for(int iCell = 0 ; iCell < data.CellPositionZ.size(); iCell++)
  // {
  //   // double meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2])/2.0;
  //   bool zCellIsAlreadyThere = false;
  //   for(unsigned int iz = 0; iz < data.zCell.size(); iz++)
  //   {
  //     if(data.CellPositionZ[iCell] == data.zCell[iz])
  //     {
  //       zCellIsAlreadyThere = true;
  //     }
  //   }
  //   if(!zCellIsAlreadyThere)
  //   {
  //     data.zCell.push_back(data.CellPositionZ[iCell]);
  //     data.zDim.push_back(data.CellCrystalSizeZ[iCell]);
  //   }
  // }
  //
  // data.NofSpacalSegments = data.zCell.size();

  // // find if the dimensions of the 1 or 2 sections
  // // std::vector<double> zCell;
  // for(int iCell = 0 ; iCell < data.CellCrystalSizeZ.size(); iCell++)
  // {
  //   // double meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2])/2.0;
  //   bool zCellIsAlreadyThere = false;
  //   for(unsigned int iz = 0; iz < data.zDim.size(); iz++)
  //   {
  //     if(data.CellCrystalSizeZ[iCell] == data.zDim[iz])
  //     {
  //       zCellIsAlreadyThere = true;
  //     }
  //   }
  //   if(!zCellIsAlreadyThere)
  //   {
  //     data.zDim.push_back(data.CellCrystalSizeZ[iCell]);
  //   }
  // }

  // data.NofSpacalSegments = data.zDim.size();


  return data;
}


void Parameters::ScanMap(TH2D *histo2d)
{
  // set module block sizes
  module_block_size_x = histo2d->GetXaxis()->GetBinWidth(1) * cm;
  module_block_size_y = histo2d->GetYaxis()->GetBinWidth(1) * cm;

  xMapMin = +INFINITY;
  xMapMax = -INFINITY;
  yMapMin = +INFINITY;
  yMapMax = -INFINITY;

  // get nbisx and y
  int nBinsX = histo2d->GetXaxis()->GetNbins();
  int nBinsY = histo2d->GetYaxis()->GetNbins();
  // run on all bins, find the vector of types

  for(int i = 1 ; i < nBinsX+1 ; i++)
  {
    for(int j = 1 ; j < nBinsY+1 ; j++)
    {

      G4double x   = histo2d->GetXaxis()->GetBinCenter(i) * cm;
      G4double y   = histo2d->GetYaxis()->GetBinCenter(j) * cm;

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




      int number   = (int) histo2d->GetBinContent(i,j);
      Point p(x,y,number);
      coordinates.push_back(p);

      bool alreadyThere = false;
      for(int k = 0; k < ecal_type_number.size(); k++)
      {
        if(number == ecal_type_number[k])
        {
          alreadyThere = true;
        }
      }
      if(!alreadyThere)
      {
        ecal_type_number.push_back(number);
        ecal_type_xMapMin.push_back(+INFINITY);
        ecal_type_xMapMax.push_back(-INFINITY);
        ecal_type_yMapMin.push_back(+INFINITY);
        ecal_type_yMapMax.push_back(-INFINITY);
      }
    }
  }

  G4cout << xMapMin << " "
         << xMapMax << " "
         << yMapMin << " "
         << yMapMax << " "
         << G4endl;

  // now find x,y min and max, but for all the types
  for(int k = 0; k < ecal_type_number.size(); k++) // for each type
  {
    for(int i = 1 ; i < nBinsX+1 ; i++)
    {
      for(int j = 1 ; j < nBinsY+1 ; j++)
      {

        G4double x   = histo2d->GetXaxis()->GetBinCenter(i) * cm;
        G4double y   = histo2d->GetYaxis()->GetBinCenter(j) * cm;
        int number   = (int) histo2d->GetBinContent(i,j);

        if(number == ecal_type_number[k])
        {
          if(( x - module_block_size_x/2.0) < ecal_type_xMapMin[k] )
          {
            ecal_type_xMapMin[k] = ( x - module_block_size_x/2.0);
          }
          if(( x + module_block_size_x/2.0) > ecal_type_xMapMax[k] )
          {
            ecal_type_xMapMax[k] = ( x + module_block_size_x/2.0);
          }
          if(( y - module_block_size_y/2.0) < ecal_type_yMapMin[k] )
          {
            ecal_type_yMapMin[k] = ( y - module_block_size_y/2.0);
          }
          if(( y + module_block_size_y/2.0) > ecal_type_yMapMax[k] )
          {
            ecal_type_yMapMax[k] = ( y + module_block_size_y/2.0);
          }
        }
      }
    }
  }

  for(int k = 0; k < ecal_type_number.size(); k++) // for each type
  {
    G4cout << ecal_type_number[k]  << " "
           << ecal_type_xMapMin[k] << " "
           << ecal_type_xMapMax[k] << " "
           << ecal_type_yMapMin[k] << " "
           << ecal_type_yMapMax[k] << " "
           << G4endl;
  }


  // G4cout << "<<<<<<<<<<<<<<<<<<<<<<<<< coordinates "<<coordinates.size() << "\n";
  // for(int k = 0; k < ecal_type_number.size(); k++)
  // {
  //   G4cout << ecal_type_number[k] << " ";
  // }
  // G4cout << G4endl;

}

void Parameters::AddCrystalMaterial(int num, G4Material* aMaterial)
{
  // first check if the material is not already there!
  bool alreadyThere = false;
  for(int k = 0; k < CrystalMaterialList.size(); k++)
  {
    if(num == CrystalMaterialList[k])
    {
      alreadyThere = true;
    }
  }
  if(!alreadyThere)
  {
    CrystalMaterialList.push_back(num);
    CrystalMaterialMap.insert(std::make_pair(num,aMaterial));
  }
}

G4Material* Parameters::GetCrystalMaterial(int num)
{
  return CrystalMaterialMap[num];
}


void Parameters::PrintConfig(Data_t data)
{
  G4cout << "**********************************************"           << G4endl;
  G4cout << "* CONFIG FILE " << data.fileName                          << G4endl;
  G4cout << "**********************************************"           << G4endl;
  G4cout << "fileName                         = " << data.fileName                       << G4endl;
  G4cout << "seed                             = " << data.seed                           << G4endl;
  G4cout << "checkOverlaps                    = " << data.checkOverlaps                  << G4endl;
  G4cout << "B_field_intensity                = " << data.B_field_intensity              << G4endl;
  G4cout << "use_gps                          = " << data.use_gps                        << G4endl;
  G4cout << "input_filename                   = " << data.input_filename                 << G4endl;
  G4cout << "input_treename                   = " << data.input_treename                 << G4endl;
  G4cout << "skipEvents                       = " << data.skipEvents                     << G4endl;
  G4cout << "Zshift                           = " << data.Zshift                         << G4endl;
  G4cout << "gps_instructions_file            = " << data.gps_instructions_file          << G4endl;
  G4cout << "defaultCut                       = " << data.defaultCut                     << G4endl;
  G4cout << "simulationType                   = " << data.simulationType                 << G4endl;
  G4cout << "switchOnScintillation            = " << data.switchOnScintillation          << G4endl;
  G4cout << "propagateScintillation           = " << data.propagateScintillation         << G4endl;
  G4cout << "switchOnCerenkov                 = " << data.switchOnCerenkov               << G4endl;
  G4cout << "propagateCerenkov                = " << data.propagateCerenkov              << G4endl;
  G4cout << "saveAll                          = " << data.saveAll                        << G4endl;
  G4cout << "saveTree                         = " << data.saveTree                       << G4endl;
  G4cout << "saveShower                       = " << data.saveShower                     << G4endl;
  G4cout << "saveLAPPD                        = " << data.saveLAPPD                      << G4endl;
  G4cout << "savePrimaries                    = " << data.savePrimaries                  << G4endl;
  G4cout << "savePhotons                      = " << data.savePhotons                    << G4endl;
  G4cout << "savePhotonAbsPoint               = " << data.savePhotonAbsPoint             << G4endl;
  G4cout << "saveStructure                    = " << data.saveStructure                  << G4endl;
  G4cout << "saveSummary                      = " << data.saveSummary                    << G4endl;
  G4cout << "saveSimple                       = " << data.saveSimple                     << G4endl;
  G4cout << "savePhotonGen                    = " << data.savePhotonGen                  << G4endl;
  G4cout << "saveEnergyPerModule              = " << data.saveEnergyPerModule            << G4endl;
  G4cout << "primaries                        = " << data.primaries                      << G4endl;
  G4cout << "opticalMaterial                  = " << data.opticalMaterial                << G4endl;
  G4cout << "logging_volume                   = " << data.logging_volume                 << G4endl;
  G4cout << "pre_volume                       = " << data.pre_volume                     << G4endl;
  G4cout << "world_material                   = " << data.world_material                 << G4endl;
  G4cout << "gapSize                          = " << data.gapSize                        << G4endl;
  G4cout << "worldVisibility                  = " << data.worldVisibility                << G4endl;
  G4cout << "caloVisibility                   = " << data.caloVisibility                 << G4endl;
  G4cout << "moduleVisibility                 = " << data.moduleVisibility               << G4endl;
  G4cout << "holeVisibility                   = " << data.holeVisibility                 << G4endl;
  G4cout << "crystalsVisibility               = " << data.crystalsVisibility             << G4endl;
  G4cout << "interfaceVisibility              = " << data.interfaceVisibility            << G4endl;
  G4cout << "gapsVisibility                   = " << data.gapsVisibility                 << G4endl;
  G4cout << "readoutVisibility                = " << data.readoutVisibility              << G4endl;
  G4cout << "absorberVisibility               = " << data.absorberVisibility             << G4endl;
  G4cout << "esrVisibility                    = " << data.esrVisibility                  << G4endl;
  G4cout << "lgVisibility                     = " << data.lgVisibility                   << G4endl;
  G4cout << "containerVisibility              = " << data.containerVisibility           << G4endl;
  G4cout << "wireFrame                        = " << data.wireFrame                      << G4endl;
  G4cout << "calorimeter_position             = " << data.calorimeter_position           << G4endl;
  G4cout << "calorimeter_size_x               = " << data.calorimeter_size_x             << G4endl;
  G4cout << "calorimeter_size_y               = " << data.calorimeter_size_y             << G4endl;
  G4cout << "calorimeter_size_z               = " << data.calorimeter_size_z             << G4endl;
  G4cout << "modules_nx                       = " << data.modules_nx                     << G4endl;
  G4cout << "modules_ny                       = " << data.modules_ny                     << G4endl;
  G4cout << "pipe_modules_nx                  = " << data.pipe_modules_nx                << G4endl;
  G4cout << "pipe_modules_ny                  = " << data.pipe_modules_ny                << G4endl;
  G4cout << "InterfaceSizeZ                   = " << data.InterfaceSizeZ                 << G4endl;
  G4cout << "gap_abs_interface_material       = " << data.gap_abs_interface_material     << G4endl;
  G4cout << "gap_interface_readout_material   = " << data.gap_interface_readout_material << G4endl;
  G4cout << "ReadoutSizeZ                     = " << data.ReadoutSizeZ                   << G4endl;
  G4cout << "AbsName                          = " << data.AbsName                        << G4endl;
  G4cout << "AbsSizeX                         = " << data.AbsSizeX                       << G4endl;
  G4cout << "AbsSizeY                         = " << data.AbsSizeY                       << G4endl;
  G4cout << "AbsSizeZ                         = " << data.AbsSizeZ                       << G4endl;
  G4cout << "AbsPositionX                     = " << data.AbsPositionX                   << G4endl;
  G4cout << "AbsPositionY                     = " << data.AbsPositionY                   << G4endl;
  G4cout << "AbsPositionZ                     = " << data.AbsPositionZ                   << G4endl;
  G4cout << "AbsMaterial                      = " << data.AbsMaterial                    << G4endl;
  G4cout << "absorber_reflectivity            = " << data.absorber_reflectivity          << G4endl;
  G4cout << "absorber_specular_lobe           = " << data.absorber_specular_lobe         << G4endl;
  G4cout << "absorber_specular_spike          = " << data.absorber_specular_spike        << G4endl;
  G4cout << "absorber_backscatter             = " << data.absorber_backscatter           << G4endl;
  G4cout << "absorber_sigma_alpha             = " << data.absorber_sigma_alpha           << G4endl;
  G4cout << "W_fraction                       = " << data.W_fraction                     << G4endl;
  G4cout << "cone_material                    = " << data.cone_material                  << G4endl;
  G4cout << "esr_on_cones                     = " << data.esr_on_cones                   << G4endl;
  G4cout << "teflon_on_cones                  = " << data.teflon_on_cones                << G4endl;
  G4cout << "esrTransmittance                 = " << data.esrTransmittance               << G4endl;
  G4cout << "cell_separation_type             = " << data.cell_separation_type           << G4endl;
  G4cout << "cell_separator_position          = " << data.cell_separator_position        << G4endl;
  G4cout << "separation_thickness             = " << data.separation_thickness           << G4endl;
  G4cout << "separation_material              = " << data.separation_material            << G4endl;
  G4cout << "readoutType                      = " << data.readoutType                    << G4endl;
  G4cout << "abs_length_scale_factor          = " << data.abs_length_scale_factor        << G4endl;
  G4cout << "cladding_abs_length_scale_factor = " << data.cladding_abs_length_scale_factor << G4endl;
  G4cout << "plex_abs_length_scale_factor     = " << data.plex_abs_length_scale_factor << G4endl;
  G4cout << "user_lightyield                  = " << data.user_lightyield                << G4endl;
  G4cout << "crystal_lateral_depolishing      = " << data.crystal_lateral_depolishing    << G4endl;
  G4cout << "crystal_exit_depolishing         = " << data.crystal_exit_depolishing       << G4endl;
  // G4cout << "optical_calibration_module_type  = " << data.optical_calibration_module_type<< G4endl;
  G4cout << "ecal_position                    = " << data.ecal_position                  << G4endl;
  G4cout << "volume_file                      = " << data.volume_file                    << G4endl;
  G4cout << "moduleZshift                     = " << data.moduleZshift                   << G4endl;
  G4cout << "crystal_inner_cladding_fraction  = " << data.crystal_inner_cladding_fraction<< G4endl;
  G4cout << "crystal_outer_cladding_fraction  = " << data.crystal_outer_cladding_fraction<< G4endl;
  G4cout << "event_time_cut_ns                = " << data.event_time_cut<< G4endl;
  G4cout << "light_guide_file                 = " << data.light_guide_file << G4endl;
  G4cout << "esr_on_positive_exit             = " << data.esr_on_positive_exit << G4endl;
  G4cout << "esr_on_negative_exit             = " << data.esr_on_negative_exit << G4endl;
  // G4cout << "module_visibility                = " << data.module_visibility              << G4endl;
  G4cout << "**********************************************"           << G4endl;
  G4cout << "* End of CONFIG FILE " << data.fileName                   << G4endl;
  G4cout << "**********************************************"           << G4endl;

}


void Parameters::WriteParameters(TFile *outfile)
{
  CreateParamTree();
  // write main config parameters (always)
  // TDirectory *parametersDir = outfile->mkdir("Parameters");
  outfile->cd();
  // TDirectory *configDir = gDirectory->mkdir("ConfigFile");
  // configDir->cd();
  doWriteParameters(outfile,main_config);
  // configDir->cd("../..");
  // write all the module config parameters, if present
  for(unsigned int i = 0 ; i < module_config.size(); i++)
  {
    doWriteParameters(outfile,module_config[i]);
  }
  paramTTree->Write();
}

void Parameters::doWriteParameters(TFile *outfile, Data_t parameters)
{
  // easy way, use a TTree to write the struc in the output file
  // do a copy 
  paramStruct = parameters;
  // fill TTree
  paramTTree->Fill();
  // write to file
  // paramTTree->Write();

}

void Parameters::CreateParamTree()
{
  // create the ttree
  paramTTree = new TTree("parameters","parameters");
  //link the TTree to the struct  
  paramTTree->Branch("moduleZshift",&paramStruct.moduleZshift,"moduleZshift/D");
  paramTTree->Branch("verbosity",&paramStruct.verbosity,"verbosity/I");
  paramTTree->Branch("fileName",&paramStruct.fileName);
  paramTTree->Branch("ecal_position",&paramStruct.ecal_position,"ecal_position/I");
  paramTTree->Branch("volume_file",&paramStruct.volume_file);
  paramTTree->Branch("light_guide_file",&paramStruct.light_guide_file);
  paramTTree->Branch("ecal_map_file",&paramStruct.ecal_map_file);
  paramTTree->Branch("ecal_map_histo",&paramStruct.ecal_map_histo);
  // paramTTree->Branch("seed",&paramStruct.seed,"seed/G");
  paramTTree->Branch("checkOverlaps",&paramStruct.checkOverlaps,"checkOverlaps/O");
  paramTTree->Branch("B_field_intensity",&paramStruct.B_field_intensity,"B_field_intensity/D");
  paramTTree->Branch("use_gps",&paramStruct.use_gps,"use_gps/O");
  paramTTree->Branch("input_filename",&paramStruct.input_filename);
  paramTTree->Branch("input_treename",&paramStruct.input_treename);
  paramTTree->Branch("skipEvents",&paramStruct.skipEvents,"skipEvents/I");
  paramTTree->Branch("Zshift",&paramStruct.Zshift,"Zshift/D");
  paramTTree->Branch("gps_instructions_file",&paramStruct.gps_instructions_file);
  paramTTree->Branch("defaultCut",&paramStruct.defaultCut,"defaultCut/D");
  paramTTree->Branch("simulationType",&paramStruct.simulationType,"simulationType/I");
  paramTTree->Branch("switchOnScintillation",&paramStruct.switchOnScintillation,"switchOnScintillation/I");
  paramTTree->Branch("propagateScintillation",&paramStruct.propagateScintillation,"propagateScintillation/I");
  paramTTree->Branch("switchOnCerenkov",&paramStruct.switchOnCerenkov,"switchOnCerenkov/I");
  paramTTree->Branch("propagateCerenkov",&paramStruct.propagateCerenkov,"propagateCerenkov/I");
  paramTTree->Branch("saveAll",&paramStruct.saveAll,"saveAll/O");
  paramTTree->Branch("saveTree",&paramStruct.saveTree,"saveTree/O");
  paramTTree->Branch("saveShower",&paramStruct.saveShower,"saveShower/O");
  paramTTree->Branch("savePrimaries",&paramStruct.savePrimaries,"savePrimaries/O");
  paramTTree->Branch("savePhotons",&paramStruct.savePhotons,"savePhotons/O");
  paramTTree->Branch("savePhotonAbsPoint",&paramStruct.savePhotonAbsPoint,"savePhotonAbsPoint/O");
  paramTTree->Branch("saveStructure",&paramStruct.saveStructure,"saveStructure/O");
  paramTTree->Branch("saveSummary",&paramStruct.saveSummary,"saveSummary/O");
  paramTTree->Branch("saveSimple",&paramStruct.saveSimple,"saveSimple/O");
  paramTTree->Branch("savePhotonGen",&paramStruct.savePhotonGen,"savePhotonGen/O");
  paramTTree->Branch("saveEnergyPerModule",&paramStruct.saveEnergyPerModule,"saveEnergyPerModule/O");
  paramTTree->Branch("saveLAPPD",&paramStruct.saveLAPPD,"saveLAPPD/O");
  paramTTree->Branch("primaries",&paramStruct.primaries,"primaries/I");
  paramTTree->Branch("opticalMaterial",&paramStruct.opticalMaterial,"opticalMaterial/I");
  paramTTree->Branch("logging_volume",&paramStruct.logging_volume);
  paramTTree->Branch("pre_volume",&paramStruct.pre_volume);
  paramTTree->Branch("world_material",&paramStruct.world_material,"world_material/I");
  paramTTree->Branch("gapSize",&paramStruct.gapSize,"gapSize/D");
  paramTTree->Branch("abs_interface_extraGap",&paramStruct.abs_interface_extraGap,"abs_interface_extraGap/D");
  paramTTree->Branch("worldVisibility",&paramStruct.worldVisibility,"worldVisibility/I");
  paramTTree->Branch("caloVisibility",&paramStruct.caloVisibility,"caloVisibility/I");
  paramTTree->Branch("moduleVisibility",&paramStruct.moduleVisibility,"moduleVisibility/I");
  paramTTree->Branch("holeVisibility",&paramStruct.holeVisibility,"holeVisibility/I");
  paramTTree->Branch("crystalsVisibility",&paramStruct.crystalsVisibility,"crystalsVisibility/I");
  paramTTree->Branch("interfaceVisibility",&paramStruct.interfaceVisibility,"interfaceVisibility/I");
  paramTTree->Branch("gapsVisibility",&paramStruct.gapsVisibility,"gapsVisibility/I");
  paramTTree->Branch("readoutVisibility",&paramStruct.readoutVisibility,"readoutVisibility/I");
  paramTTree->Branch("absorberVisibility",&paramStruct.absorberVisibility,"absorberVisibility/I");
  paramTTree->Branch("esrVisibility",&paramStruct.esrVisibility,"esrVisibility/I");
  paramTTree->Branch("lgVisibility",&paramStruct.lgVisibility,"lgVisibility/I");
  paramTTree->Branch("lappdVisibility",&paramStruct.lappdVisibility,"lappdVisibility/I");
  paramTTree->Branch("containerVisibility",&paramStruct.containerVisibility,"containerVisibility/O");
  paramTTree->Branch("wireFrame",&paramStruct.wireFrame,"wireFrame/I");
  paramTTree->Branch("calorimeter_position",&paramStruct.calorimeter_position,"calorimeter_position/I");
  paramTTree->Branch("calorimeter_size_x",&paramStruct.calorimeter_size_x,"calorimeter_size_x/D");
  paramTTree->Branch("calorimeter_size_y",&paramStruct.calorimeter_size_y,"calorimeter_size_y/D");
  paramTTree->Branch("calorimeter_size_z",&paramStruct.calorimeter_size_z,"calorimeter_size_z/D");
  paramTTree->Branch("modules_nx",&paramStruct.modules_nx,"modules_nx/I");
  paramTTree->Branch("modules_ny",&paramStruct.modules_ny,"modules_ny/I");
  paramTTree->Branch("pipe_modules_nx",&paramStruct.pipe_modules_nx,"pipe_modules_nx/I");
  paramTTree->Branch("pipe_modules_ny",&paramStruct.pipe_modules_ny,"pipe_modules_ny/I");
  paramTTree->Branch("InterfaceSizeZ",&paramStruct.InterfaceSizeZ,"InterfaceSizeZ/D");
  paramTTree->Branch("gap_abs_interface_material",&paramStruct.gap_abs_interface_material,"gap_abs_interface_material/I");
  paramTTree->Branch("gap_interface_readout_material",&paramStruct.gap_interface_readout_material,"gap_interface_readout_material/I");
  paramTTree->Branch("ReadoutSizeZ",&paramStruct.ReadoutSizeZ,"ReadoutSizeZ/D");
  paramTTree->Branch("AbsName",&paramStruct.AbsName);
  paramTTree->Branch("AbsSizeX",&paramStruct.AbsSizeX,"AbsSizeX/D");
  paramTTree->Branch("AbsSizeY",&paramStruct.AbsSizeY,"AbsSizeY/D");
  paramTTree->Branch("AbsSizeZ",&paramStruct.AbsSizeZ,"AbsSizeZ/D");
  paramTTree->Branch("AbsPositionX",&paramStruct.AbsPositionX,"AbsPositionX/D");
  paramTTree->Branch("AbsPositionY",&paramStruct.AbsPositionY,"AbsPositionY/D");
  paramTTree->Branch("AbsPositionZ",&paramStruct.AbsPositionZ,"AbsPositionZ/D");
  paramTTree->Branch("AbsMaterial",&paramStruct.AbsMaterial,"AbsMaterial/I");
  paramTTree->Branch("absorber_reflectivity",&paramStruct.absorber_reflectivity,"absorber_reflectivity/D");
  paramTTree->Branch("absorber_specular_lobe",&paramStruct.absorber_specular_lobe,"absorber_specular_lobe/D");
  paramTTree->Branch("absorber_specular_spike",&paramStruct.absorber_specular_spike,"absorber_specular_spike/D");
  paramTTree->Branch("absorber_backscatter",&paramStruct.absorber_backscatter,"absorber_backscatter/D");
  paramTTree->Branch("absorber_sigma_alpha",&paramStruct.absorber_sigma_alpha,"absorber_sigma_alpha/D");
  paramTTree->Branch("W_fraction",&paramStruct.W_fraction,"W_fraction/D");
  paramTTree->Branch("cone_material",&paramStruct.cone_material,"cone_material/I");
  paramTTree->Branch("esr_on_cones",&paramStruct.esr_on_cones,"esr_on_cones/O");
  paramTTree->Branch("teflon_on_cones",&paramStruct.teflon_on_cones,"teflon_on_cones/O");
  // double esrTransmittance;
  // G4int cell_separation_type;
  // G4double cell_separator_position;
  // double separation_thickness;
  // int separation_material;
  paramTTree->Branch("esrTransmittance",&paramStruct.esrTransmittance,"esrTransmittance/D");
  paramTTree->Branch("cell_separation_type",&paramStruct.cell_separation_type,"cell_separation_type/I");
  paramTTree->Branch("cell_separator_position",&paramStruct.cell_separator_position,"cell_separator_position/D");
  paramTTree->Branch("separation_thickness",&paramStruct.separation_thickness,"separation_thickness/D");
  paramTTree->Branch("separation_material",&paramStruct.separation_material,"separation_material/I");
  paramTTree->Branch("readoutType",&paramStruct.readoutType,"readoutType/I");
  paramTTree->Branch("abs_length_scale_factor",&paramStruct.abs_length_scale_factor,"abs_length_scale_factor/D");
  paramTTree->Branch("cladding_abs_length_scale_factor",&paramStruct.cladding_abs_length_scale_factor,"cladding_abs_length_scale_factor/D");
  paramTTree->Branch("plex_abs_length_scale_factor",&paramStruct.plex_abs_length_scale_factor,"plex_abs_length_scale_factor/D");
  paramTTree->Branch("user_lightyield",&paramStruct.user_lightyield,"user_lightyield/D");
  paramTTree->Branch("crystal_lateral_depolishing",&paramStruct.crystal_lateral_depolishing,"crystal_lateral_depolishing/D");
  paramTTree->Branch("crystal_exit_depolishing",&paramStruct.crystal_exit_depolishing,"crystal_exit_depolishing/D");
  paramTTree->Branch("crystal_inner_cladding_fraction",&paramStruct.crystal_inner_cladding_fraction,"crystal_inner_cladding_fraction/D");
  paramTTree->Branch("crystal_outer_cladding_fraction",&paramStruct.crystal_outer_cladding_fraction,"crystal_outer_cladding_fraction/D");
  paramTTree->Branch("event_time_cut",&paramStruct.event_time_cut,"event_time_cut/D");
  paramTTree->Branch("LAPPD_thickness",&paramStruct.LAPPD_thickness,"LAPPD_thickness/D");
  paramTTree->Branch("SingleSeparationLayer",&paramStruct.SingleSeparationLayer,"SingleSeparationLayer/O");
  paramTTree->Branch("esr_on_positive_exit",&paramStruct.esr_on_positive_exit,"esr_on_positive_exit/O");
  paramTTree->Branch("esr_on_negative_exit",&paramStruct.esr_on_negative_exit,"esr_on_negative_exit/O");
  // t->Branch("",&);
  paramTTree->Branch("CellNames",&paramStruct.CellNames);
  paramTTree->Branch("CellPositionX",&paramStruct.CellPositionX);
  paramTTree->Branch("CellPositionY",&paramStruct.CellPositionY);
  paramTTree->Branch("CellPositionZ",&paramStruct.CellPositionZ);
  paramTTree->Branch("CellXelements",&paramStruct.CellXelements);
  paramTTree->Branch("CellYelements",&paramStruct.CellYelements);
  paramTTree->Branch("CellCrystalSizeX",&paramStruct.CellCrystalSizeX);
  paramTTree->Branch("CellCrystalSizeY",&paramStruct.CellCrystalSizeY);
  paramTTree->Branch("CellCrystalSizeZ",&paramStruct.CellCrystalSizeZ);
  paramTTree->Branch("CellCrystalPitchX",&paramStruct.CellCrystalPitchX);
  paramTTree->Branch("CellCrystalPitchY",&paramStruct.CellCrystalPitchY);
  paramTTree->Branch("CellMaterial",&paramStruct.CellMaterial);
  paramTTree->Branch("CellAirLayer",&paramStruct.CellAirLayer);
  paramTTree->Branch("CellHoleShape",&paramStruct.CellHoleShape);
  paramTTree->Branch("CellCrystalShape",&paramStruct.CellCrystalShape);
  paramTTree->Branch("CellCrystalCladding",&paramStruct.CellCrystalCladding);
  paramTTree->Branch("CellIntGapMaterial",&paramStruct.CellIntGapMaterial);
  paramTTree->Branch("CellStaggering",&paramStruct.CellStaggering);
  paramTTree->Branch("CellStaggeringAxis",&paramStruct.CellStaggeringAxis);
  paramTTree->Branch("CellStaggeringSize",&paramStruct.CellStaggeringSize);
  paramTTree->Branch("CellStaggeringParity",&paramStruct.CellStaggeringParity);
  paramTTree->Branch("CellStaggeringRemove",&paramStruct.CellStaggeringRemove);
  paramTTree->Branch("LAPPD_layers",&paramStruct.LAPPD_layers);
  paramTTree->Branch("LAPPD_materials",&paramStruct.LAPPD_materials);
  paramTTree->Branch("ecal_type",&paramStruct.ecal_type);
  paramTTree->Branch("container_volume",&paramStruct.container_volume);
  paramTTree->Branch("module_config_list",&paramStruct.module_config_list);
  paramTTree->Branch("crystalMaterialList",&paramStruct.crystalMaterialList);
}
