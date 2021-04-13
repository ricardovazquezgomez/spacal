
// a structure to hold all the config file keys
struct Data_t
{

  // ------------simple variables
  int ecal_position;
  double moduleZshift;
  int verbosity;
  // seed for random generator. if not given, default to -1,
  // which means the seed is computed randomly
  long int seed;
  // ask Geant4 to check if there are overlapping volumes - default = False
  bool checkOverlaps;
  // intensity of magnetic field, in Tesla - default = 0
  double B_field_intensity;
  // source can be defined by a single particle generator, with a gps file,
  // or by a ROOT file containing the a LHCb multi particle flux
  // choice between the two possibilities is stored in use_gps
  // (if false, multi particle flux used)
  bool use_gps;
  // start from event skipEvents in the tree - ignored if use_gps = true
  int skipEvents;
  // actual range: 12600 - 12629 - ignored if use_gps = true
  double Zshift;
  // production cut value in mm (if not given, default to 0.025 mm)
  double defaultCut;
  // key to set type of simulation.
  // 4 types are defined, and if chosen will overwrite the "output group" keys
  // and also "primaries" and "opticalMaterial" keys. If this key is not given instead, the
  // simulation type is considered "free" and all the keys above can be set freely by the user.
  // The simulation types are:
  // 0) only energy deposition
  // 1) full ray tracing
  // 2) hybrid simulation
  // 3) optical calibration
  int simulationType;
  // keys to control output. if none is specified, and sim type is free
  // so everything will be saved!
  int switchOnScintillation;
  int propagateScintillation;
  int switchOnCerenkov;
  int propagateCerenkov;
  bool saveAll;
  bool saveTree;
  bool saveShower;
  bool savePrimaries;
  bool savePhotons;
  bool savePhotonAbsPoint;
  bool saveStructure;
  bool saveSummary;
  bool saveSimple;
  bool savePhotonGen;
  bool saveEnergyPerModule;
  bool saveLAPPD;
  // number of primary particles shot for each /run/beamOn 1
  int primaries;
  // special key for optical calibration. ignored if primary in gps is not
  // opticalPhoton. Set to 0 if you want an instantaneous monocromatic
  // emission, with energy decided by gps. Leave commented you want
  // the optical emission to depend of the material where the point of emission is.
  int opticalMaterial;
  
  // world material: 1) Air (only choice)
  int world_material;
  // general dimension of air gaps, whenever two surfaces are touching in dry contact
  double gapSize;
  // extra gap to be added between the absorber and the interface
  double abs_interface_extraGap;
  // visualize or not the individual elements in viz mode. 1 = true, 0 = false
  int worldVisibility;
  int caloVisibility;
  int moduleVisibility;
  int holeVisibility;
  int crystalsVisibility;
  int interfaceVisibility;
  int gapsVisibility;
  int readoutVisibility;
  int absorberVisibility;
  int esrVisibility;
  int lgVisibility;
  int lappdVisibility;
  // special key for showing containers
  bool containerVisibility;
  // wire frame or solid Visualization 1 = wireFrame, 0 = solid
  int wireFrame;
  // 0)  Calorimeter center in z = 0, 1) Calorimeter center such
  // that calorimeter starts in z = 0
  int calorimeter_position;
  // calorimeter size, if the calo construction modality is 2
  double calorimeter_size_x;
  double calorimeter_size_y;
  double calorimeter_size_z;
  // In standard mode, the calorimeter is made of an array of identical modules
  // The array is in x-y, no space between modules
  int modules_nx;
  int modules_ny;
  // modules in the center of the calorimeter can be removed, leaving a hole
  // for the beam pipe the size of the hole is specified in terms of number
  // of missing modules, in x and y using the keys pipe_modules_nx and pipe_modules_ny
  int pipe_modules_nx;
  int pipe_modules_ny;
  // z length for both interfaces
  double InterfaceSizeZ;
  // 0 = air 1 = glue
  int gap_abs_interface_material;
  int gap_interface_readout_material;
  // size of both readout volumes
  double ReadoutSizeZ;
  // The absorber volume. these options are self explanatory
  double AbsSizeX;
  double AbsSizeY;
  double AbsSizeZ;
  double AbsPositionX;
  double AbsPositionY;
  double AbsPositionZ;
  // absorber material: 1) Brass 2) Tungsten alloy 3) Lead 4)
  // Iron 5) Aluminium 6) CopperTungstenAlloy 7) Pure Tungsten, 18.5 g/cm3
  // 8)  Pure Tungsten, 19.1 g/cm3 9) Air Killer 10) Air
  int AbsMaterial;
  // reflectivity of the absorber
  double absorber_reflectivity;
  double absorber_specular_lobe;
  double absorber_specular_spike;
  double absorber_backscatter   ;
  double absorber_sigma_alpha   ;
  // fraction of Tungsten in the alloy (only for option 6)
  double W_fraction;
  // material of optical light guides, if present
  // 0 = air 1 = plexiglass
  int cone_material;
  // 0 = no esr 1 = esr
  bool esr_on_cones;
  // 0 = no teflon 1 = teflon
  bool teflon_on_cones;
  // probability for a optical photon to cross ESR - default = 0
  double esrTransmittance;
  // module longitudinal separation
  // Type can be
  // 0 = nothing (air)
  // 1 = aluminization
  // 2 = reflector (esr)
  int cell_separation_type;
  // in mm
  double cell_separator_position;
  // ignored if cell_separation_type != 2
  double separation_thickness;
  // ignored if cell_separation_type != 2
  int separation_material;
  // user specific
  // readout type. 0) 1cm air + PMTs 1) test beam 2018 2) test beam 2019
  // 3)  1cm air + 1 big PMT 4) 0.1 mm air + 1 big squared PMT
  // 5) no readout (for hybrid simulation) - default = 5
  int readoutType;
  // multiply the abs length of the crystal material by a factor.
  double abs_length_scale_factor;
  // multiply the abs length of the cladding materials by a factor.
  double cladding_abs_length_scale_factor;
  // multiply the abs length of PLEX material by a factor.
  double plex_abs_length_scale_factor;
  // user light yield
  double user_lightyield;
  // depolishing of lateral and exit surfaces (comment out for no depolishing)
  double crystal_lateral_depolishing;
  double crystal_exit_depolishing;
  // special flag for optical calibration
  // int optical_calibration_module_type;
  // fraction of inner and outer cladding thickness, relative to diameter
  double crystal_inner_cladding_fraction;
  double crystal_outer_cladding_fraction;
  // user set time limit for events 
  double event_time_cut;

  double LAPPD_thickness;
  bool SingleSeparationLayer;

  bool esr_on_positive_exit;
  bool esr_on_negative_exit;

  // -------------strings 
  std::string input_filename;
  std::string input_treename;
  std::string fileName;
  std::string volume_file;
  std::string light_guide_file;
  std::string ecal_map_file;
  std::string ecal_map_histo;
  // gps macro file. this will be ignored if the program is run
  // in "exec with gps" mode (see above)
  std::string gps_instructions_file;
  // keys to control at which position the optical photons are killed
  // and recorded as "detected"
  // Photons are considered detected if they enter the logging_volume
  // and optionally the user can also impose that they have to be
  // entering it coming from pre_volume
  std::string logging_volume;
  std::string pre_volume;
  std::string AbsName;

  // -------------vectors
  // configure files for different modules
  std::vector<std::string> module_config_list;
  // ecal module type
  std::vector<int> ecal_type;
  // ecal region volumes, CAD file names
  std::vector<std::string> container_volume;
  // cell dimensions and stuff
  std::vector<std::string>   CellNames;
  std::vector<double>   CellPositionX;
  std::vector<double>   CellPositionY;
  std::vector<double>   CellPositionZ;
  std::vector<int>      CellXelements;
  std::vector<int>      CellYelements;
  std::vector<double>   CellCrystalSizeX;
  std::vector<double>   CellCrystalSizeY;
  std::vector<double>   CellCrystalSizeZ;
  std::vector<double>   CellCrystalPitchX;
  std::vector<double>   CellCrystalPitchY;
  std::vector<int>      CellMaterial;
  // thickness of layer of air between hole and crystal [mm]
  std::vector<double>   CellAirLayer  ;
  // shape of holes. 0 = squared, 1 = round. DEFAULT = 0. If = 1, x coordinate interpreted as diameter (y is ignored)
  std::vector<int>   CellHoleShape  ;
  // shape of crystals. 0 = squared, 1 = round. DEFAULT = 0. If = 1, x coordinate interpreted as diameter (y is ignored)
  std::vector<int>   CellCrystalShape  ;
  // double cladding. 0 = no, 1 = yes. DEFAULT = 0
  std::vector<bool>   CellCrystalCladding  ;
  // 1 = air 2 = optical grease
  std::vector<int>      CellIntGapMaterial;
  // STAGGERING
  // elements can be simply aligned, or there can be a row
  // staggering either in x or in y
  // if left commented out, there is no row staggering
  // 0 = no, 1 = yes
  std::vector<int>      CellStaggering;
  // 0 = x, 1 = y
  std::vector<int>      CellStaggeringAxis;
  // displacement in mm
  std::vector<double>   CellStaggeringSize  ;
  // 0 = even rows, 1 = odd rows
  std::vector<int>      CellStaggeringParity;
  // 1 = staggering means remove an element
  std::vector<int>      CellStaggeringRemove;
  //Thickness of LAPPD layers
  std::vector<double>   LAPPD_layers;
  //Materials of LAPPD layers
  std::vector<int>   LAPPD_materials;
  std::vector<int> crystalMaterialList;


  //---------------------------//
  // POINTERS FOR READING
  //---------------------------//

  // strings 
  std::string *v_input_filename = 0;
  std::string *v_input_treename = 0;
  std::string *v_fileName = 0;
  std::string *v_volume_file = 0;
  std::string *v_light_guide_file = 0;
  std::string *v_ecal_map_file = 0;
  std::string *v_ecal_map_histo = 0;
  std::string *v_gps_instructions_file = 0;
  std::string *v_logging_volume = 0;
  std::string *v_pre_volume = 0;
  std::string *v_AbsName = 0;

  // std::vector<std::string> module_config_list;
  // std::vector<std::string> container_volume;
  std::vector<int>           *v_ecal_type = 0;
  std::vector<double>        *v_CellPositionX = 0;
  std::vector<double>        *v_CellPositionY = 0;
  std::vector<double>        *v_CellPositionZ = 0;
  std::vector<int>           *v_CellXelements = 0;
  std::vector<int>           *v_CellYelements = 0;
  std::vector<double>        *v_CellCrystalSizeX = 0;
  std::vector<double>        *v_CellCrystalSizeY = 0;
  std::vector<double>        *v_CellCrystalSizeZ = 0;
  std::vector<double>        *v_CellCrystalPitchX = 0;
  std::vector<double>        *v_CellCrystalPitchY = 0;
  std::vector<int>           *v_CellMaterial = 0;
  std::vector<double>        *v_CellAirLayer   = 0;
  std::vector<int>           *v_CellHoleShape   = 0;
  std::vector<int>           *v_CellCrystalShape   = 0;
  std::vector<bool>          *v_CellCrystalCladding   = 0;
  std::vector<int>           *v_CellIntGapMaterial = 0;
  std::vector<int>           *v_CellStaggering = 0;
  std::vector<int>           *v_CellStaggeringAxis = 0;
  std::vector<double>        *v_CellStaggeringSize   = 0;
  std::vector<int>           *v_CellStaggeringParity = 0;
  std::vector<int>           *v_CellStaggeringRemove = 0;
  std::vector<double>        *v_LAPPD_layers = 0;
  std::vector<int>           *v_LAPPD_materials = 0;
  std::vector<int>           *v_crystalMaterialList = 0;
};



struct relevant_data_t
{
  char fileName[200]; // not relevant at all, but necessary for debugging
  double defaultCut;
  int ecal_position;
  char logging_volume[200];
  char pre_volume[200];
  double gapSize;
  double abs_interface_extraGap;
  int gap_abs_interface_material;
  int gap_interface_readout_material;
  double AbsSizeZ;
  double AbsPositionX;
  double AbsPositionY;
  double AbsPositionZ;
  int AbsMaterial;
  double absorber_reflectivity;
  double absorber_specular_lobe;
  double absorber_specular_spike;
  double absorber_backscatter;
  double absorber_sigma_alpha;
  double W_fraction;
  double esrTransmittance;
  int cell_separation_type;
  double cell_separator_position;
  double separation_thickness;
  int separation_material;
  double abs_length_scale_factor;
  double crystal_lateral_depolishing;
  double crystal_exit_depolishing;
  double crystal_inner_cladding_fraction;
  double crystal_outer_cladding_fraction;
  bool esr_on_positive_exit;
  bool esr_on_negative_exit;
  int sections; // only 1 or 2
  // these are vectors but could be just one value  
  double CellCrystalSizeZ[2];
  int CellMaterial[2];
  double CellAirLayer[2];
  int CellHoleShape[2];
  int CellCrystalShape[2];
  bool CellCrystalCladding[2];
  double CellIntGapMaterial[2];
};


void readParametersFillStruc(TTree* parameters, Data_t &footer);
void fillRelevant(relevant_data_t &relevant,Data_t &parameters);
int compareRelevant(relevant_data_t &deposition,relevant_data_t &calibration);

