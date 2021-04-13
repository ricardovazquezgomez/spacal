#include "../include/Data.hh"

void readParametersFillStruc(TTree* parameters, Data_t &footer)
{
  TBranch *b_moduleZshift;
  TBranch *b_verbosity;
  TBranch *b_fileName;
  TBranch *b_ecal_position;
  TBranch *b_volume_file;
  TBranch *b_light_guide_file;
  TBranch *b_ecal_map_file;
  TBranch *b_ecal_map_histo;
  TBranch *b_seed;
  TBranch *b_checkOverlaps;
  TBranch *b_B_field_intensity;
  TBranch *b_use_gps;
  TBranch *b_input_filename;
  TBranch *b_input_treename;
  TBranch *b_skipEvents;
  TBranch *b_Zshift;
  TBranch *b_gps_instructions_file;
  TBranch *b_defaultCut;
  TBranch *b_simulationType;
  TBranch *b_switchOnScintillation;
  TBranch *b_propagateScintillation;
  TBranch *b_switchOnCerenkov;
  TBranch *b_propagateCerenkov;
  TBranch *b_saveAll;
  TBranch *b_saveTree;
  TBranch *b_saveShower;
  TBranch *b_savePrimaries;
  TBranch *b_savePhotons;
  TBranch *b_savePhotonAbsPoint;
  TBranch *b_saveStructure;
  TBranch *b_saveSummary;
  TBranch *b_saveSimple;
  TBranch *b_savePhotonGen;
  TBranch *b_saveEnergyPerModule;
  TBranch *b_saveLAPPD;
  TBranch *b_primaries;
  TBranch *b_opticalMaterial;
  TBranch *b_logging_volume;
  TBranch *b_pre_volume;
  TBranch *b_world_material;
  TBranch *b_gapSize;
  TBranch *b_abs_interface_extraGap;
  TBranch *b_worldVisibility;
  TBranch *b_caloVisibility;
  TBranch *b_moduleVisibility;
  TBranch *b_holeVisibility;
  TBranch *b_crystalsVisibility;
  TBranch *b_interfaceVisibility;
  TBranch *b_gapsVisibility;
  TBranch *b_readoutVisibility;
  TBranch *b_absorberVisibility;
  TBranch *b_esrVisibility;
  TBranch *b_lgVisibility;
  TBranch *b_lappdVisibility;
  TBranch *b_containerVisibility;
  TBranch *b_wireFrame;
  TBranch *b_calorimeter_position;
  TBranch *b_calorimeter_size_x;
  TBranch *b_calorimeter_size_y;
  TBranch *b_calorimeter_size_z;
  TBranch *b_modules_nx;
  TBranch *b_modules_ny;
  TBranch *b_pipe_modules_nx;
  TBranch *b_pipe_modules_ny;
  TBranch *b_InterfaceSizeZ;
  TBranch *b_gap_abs_interface_material;
  TBranch *b_gap_interface_readout_material;
  TBranch *b_ReadoutSizeZ;
  TBranch *b_AbsName;
  TBranch *b_AbsSizeX;
  TBranch *b_AbsSizeY;
  TBranch *b_AbsSizeZ;
  TBranch *b_AbsPositionX;
  TBranch *b_AbsPositionY;
  TBranch *b_AbsPositionZ;
  TBranch *b_AbsMaterial;
  TBranch *b_absorber_reflectivity;
  TBranch *b_absorber_specular_lobe;
  TBranch *b_absorber_specular_spike;
  TBranch *b_absorber_backscatter   ;
  TBranch *b_absorber_sigma_alpha   ;
  TBranch *b_W_fraction;
  TBranch *b_cone_material;
  TBranch *b_esr_on_cones;
  TBranch *b_teflon_on_cones;
  TBranch *b_esrTransmittance;
  TBranch *b_cell_separation_type;
  TBranch *b_cell_separator_position;
  TBranch *b_separation_thickness;
  TBranch *b_separation_material;
  TBranch *b_readoutType;
  TBranch *b_abs_length_scale_factor;
  TBranch *b_cladding_abs_length_scale_factor;
  TBranch *b_plex_abs_length_scale_factor;
  TBranch *b_user_lightyield;
  TBranch *b_crystal_lateral_depolishing;
  TBranch *b_crystal_exit_depolishing;
  TBranch *b_crystal_inner_cladding_fraction;
  TBranch *b_crystal_outer_cladding_fraction;
  TBranch *b_event_time_cut;
  TBranch *b_LAPPD_thickness;
  TBranch *b_SingleSeparationLayer;
  TBranch *b_esr_on_positive_exit;
  TBranch *b_esr_on_negative_exit;
  TBranch *b_CellNames = 0;
  TBranch *b_CellPositionX = 0;
  TBranch *b_CellPositionY = 0;
  TBranch *b_CellPositionZ = 0;
  TBranch *b_CellXelements = 0;
  TBranch *b_CellYelements = 0;
  TBranch *b_CellCrystalSizeX = 0;
  TBranch *b_CellCrystalSizeY = 0;
  TBranch *b_CellCrystalSizeZ = 0;
  TBranch *b_CellCrystalPitchX = 0;
  TBranch *b_CellCrystalPitchY = 0;
  TBranch *b_CellMaterial = 0;
  TBranch *b_CellAirLayer = 0;
  TBranch *b_CellHoleShape = 0;
  TBranch *b_CellCrystalShape = 0;
  TBranch *b_CellCrystalCladding = 0;
  TBranch *b_CellIntGapMaterial = 0;
  TBranch *b_CellStaggering = 0;
  TBranch *b_CellStaggeringAxis = 0;
  TBranch *b_CellStaggeringSize = 0;
  TBranch *b_CellStaggeringParity = 0;
  TBranch *b_CellStaggeringRemove = 0;
  TBranch *b_LAPPD_layers = 0;
  TBranch *b_LAPPD_materials = 0;
  TBranch *b_ecal_type = 0;
  TBranch *b_container_volume = 0;
  TBranch *b_module_config_list = 0;
  TBranch *b_crystalMaterialList = 0;



  parameters->SetBranchAddress("moduleZshift",&footer.moduleZshift,&b_moduleZshift);
  parameters->SetBranchAddress("verbosity",&footer.verbosity,&b_verbosity);
  parameters->SetBranchAddress("ecal_position",&footer.ecal_position,&b_ecal_position);
  parameters->SetBranchAddress("checkOverlaps",&footer.checkOverlaps,&b_checkOverlaps);
  parameters->SetBranchAddress("B_field_intensity",&footer.B_field_intensity,&b_B_field_intensity);
  parameters->SetBranchAddress("use_gps",&footer.use_gps,&b_use_gps);
  parameters->SetBranchAddress("skipEvents",&footer.skipEvents,&b_skipEvents);
  parameters->SetBranchAddress("Zshift",&footer.Zshift,&b_Zshift);
  parameters->SetBranchAddress("defaultCut",&footer.defaultCut,&b_defaultCut);
  parameters->SetBranchAddress("simulationType",&footer.simulationType,&b_simulationType);
  parameters->SetBranchAddress("switchOnScintillation",&footer.switchOnScintillation,&b_switchOnScintillation);
  parameters->SetBranchAddress("propagateScintillation",&footer.propagateScintillation,&b_propagateScintillation);
  parameters->SetBranchAddress("switchOnCerenkov",&footer.switchOnCerenkov,&b_switchOnCerenkov);
  parameters->SetBranchAddress("propagateCerenkov",&footer.propagateCerenkov,&b_propagateCerenkov);
  parameters->SetBranchAddress("saveAll",&footer.saveAll,&b_saveAll);
  parameters->SetBranchAddress("saveTree",&footer.saveTree,&b_saveTree);
  parameters->SetBranchAddress("saveShower",&footer.saveShower,&b_saveShower);
  parameters->SetBranchAddress("savePrimaries",&footer.savePrimaries,&b_savePrimaries);
  parameters->SetBranchAddress("savePhotons",&footer.savePhotons,&b_savePhotons);
  parameters->SetBranchAddress("savePhotonAbsPoint",&footer.savePhotonAbsPoint,&b_savePhotonAbsPoint);
  parameters->SetBranchAddress("saveStructure",&footer.saveStructure,&b_saveStructure);
  parameters->SetBranchAddress("saveSummary",&footer.saveSummary,&b_saveSummary);
  parameters->SetBranchAddress("saveSimple",&footer.saveSimple,&b_saveSimple);
  parameters->SetBranchAddress("savePhotonGen",&footer.savePhotonGen,&b_savePhotonGen);
  parameters->SetBranchAddress("saveEnergyPerModule",&footer.saveEnergyPerModule,&b_saveEnergyPerModule);
  parameters->SetBranchAddress("saveLAPPD",&footer.saveLAPPD,&b_saveLAPPD);
  parameters->SetBranchAddress("primaries",&footer.primaries,&b_primaries);
  parameters->SetBranchAddress("opticalMaterial",&footer.opticalMaterial,&b_opticalMaterial);
  parameters->SetBranchAddress("world_material",&footer.world_material,&b_world_material);
  parameters->SetBranchAddress("gapSize",&footer.gapSize,&b_gapSize);
  parameters->SetBranchAddress("abs_interface_extraGap",&footer.abs_interface_extraGap,&b_abs_interface_extraGap);
  parameters->SetBranchAddress("worldVisibility",&footer.worldVisibility,&b_worldVisibility);
  parameters->SetBranchAddress("caloVisibility",&footer.caloVisibility,&b_caloVisibility);
  parameters->SetBranchAddress("moduleVisibility",&footer.moduleVisibility,&b_moduleVisibility);
  parameters->SetBranchAddress("holeVisibility",&footer.holeVisibility,&b_holeVisibility);
  parameters->SetBranchAddress("crystalsVisibility",&footer.crystalsVisibility,&b_crystalsVisibility);
  parameters->SetBranchAddress("interfaceVisibility",&footer.interfaceVisibility,&b_interfaceVisibility);
  parameters->SetBranchAddress("gapsVisibility",&footer.gapsVisibility,&b_gapsVisibility);
  parameters->SetBranchAddress("readoutVisibility",&footer.readoutVisibility,&b_readoutVisibility);
  parameters->SetBranchAddress("absorberVisibility",&footer.absorberVisibility,&b_absorberVisibility);
  parameters->SetBranchAddress("esrVisibility",&footer.esrVisibility,&b_esrVisibility);
  parameters->SetBranchAddress("lgVisibility",&footer.lgVisibility,&b_lgVisibility);
  parameters->SetBranchAddress("lappdVisibility",&footer.lappdVisibility,&b_lappdVisibility);
  parameters->SetBranchAddress("containerVisibility",&footer.containerVisibility,&b_containerVisibility);
  parameters->SetBranchAddress("wireFrame",&footer.wireFrame,&b_wireFrame);
  parameters->SetBranchAddress("calorimeter_position",&footer.calorimeter_position,&b_calorimeter_position);
  parameters->SetBranchAddress("calorimeter_size_x",&footer.calorimeter_size_x,&b_calorimeter_size_x);
  parameters->SetBranchAddress("calorimeter_size_y",&footer.calorimeter_size_y,&b_calorimeter_size_y);
  parameters->SetBranchAddress("calorimeter_size_z",&footer.calorimeter_size_z,&b_calorimeter_size_z);
  parameters->SetBranchAddress("modules_nx",&footer.modules_nx,&b_modules_nx);
  parameters->SetBranchAddress("modules_ny",&footer.modules_ny,&b_modules_ny);
  parameters->SetBranchAddress("pipe_modules_nx",&footer.pipe_modules_nx,&b_pipe_modules_nx);
  parameters->SetBranchAddress("pipe_modules_ny",&footer.pipe_modules_ny,&b_pipe_modules_ny);
  parameters->SetBranchAddress("InterfaceSizeZ",&footer.InterfaceSizeZ,&b_InterfaceSizeZ);
  parameters->SetBranchAddress("gap_abs_interface_material",&footer.gap_abs_interface_material,&b_gap_abs_interface_material);
  parameters->SetBranchAddress("gap_interface_readout_material",&footer.gap_interface_readout_material,&b_gap_interface_readout_material);
  parameters->SetBranchAddress("ReadoutSizeZ",&footer.ReadoutSizeZ,&b_ReadoutSizeZ);
  parameters->SetBranchAddress("AbsSizeX",&footer.AbsSizeX,&b_AbsSizeX);
  parameters->SetBranchAddress("AbsSizeY",&footer.AbsSizeY,&b_AbsSizeY);
  parameters->SetBranchAddress("AbsSizeZ",&footer.AbsSizeZ,&b_AbsSizeZ);
  parameters->SetBranchAddress("AbsPositionX",&footer.AbsPositionX,&b_AbsPositionX);
  parameters->SetBranchAddress("AbsPositionY",&footer.AbsPositionY,&b_AbsPositionY);
  parameters->SetBranchAddress("AbsPositionZ",&footer.AbsPositionZ,&b_AbsPositionZ);
  parameters->SetBranchAddress("AbsMaterial",&footer.AbsMaterial,&b_AbsMaterial);
  parameters->SetBranchAddress("absorber_reflectivity",&footer.absorber_reflectivity,&b_absorber_reflectivity);
  parameters->SetBranchAddress("absorber_specular_lobe",&footer.absorber_specular_lobe,&b_absorber_specular_lobe);
  parameters->SetBranchAddress("absorber_specular_spike",&footer.absorber_specular_spike,&b_absorber_specular_spike);
  parameters->SetBranchAddress("absorber_backscatter",&footer.absorber_backscatter   ,&b_absorber_backscatter   );
  parameters->SetBranchAddress("absorber_sigma_alpha",&footer.absorber_sigma_alpha   ,&b_absorber_sigma_alpha   );
  parameters->SetBranchAddress("W_fraction",&footer.W_fraction,&b_W_fraction);
  parameters->SetBranchAddress("cone_material",&footer.cone_material,&b_cone_material);
  parameters->SetBranchAddress("esr_on_cones",&footer.esr_on_cones,&b_esr_on_cones);
  parameters->SetBranchAddress("teflon_on_cones",&footer.teflon_on_cones,&b_teflon_on_cones);
  parameters->SetBranchAddress("esrTransmittance",&footer.esrTransmittance,&b_esrTransmittance);
  parameters->SetBranchAddress("cell_separation_type",&footer.cell_separation_type,&b_cell_separation_type);
  parameters->SetBranchAddress("cell_separator_position",&footer.cell_separator_position,&b_cell_separator_position);
  parameters->SetBranchAddress("separation_thickness",&footer.separation_thickness,&b_separation_thickness);
  parameters->SetBranchAddress("separation_material",&footer.separation_material,&b_separation_material);
  parameters->SetBranchAddress("readoutType",&footer.readoutType,&b_readoutType);
  parameters->SetBranchAddress("abs_length_scale_factor",&footer.abs_length_scale_factor,&b_abs_length_scale_factor);
  parameters->SetBranchAddress("plex_abs_length_scale_factor",&footer.plex_abs_length_scale_factor,&b_plex_abs_length_scale_factor);
  parameters->SetBranchAddress("user_lightyield",&footer.user_lightyield,&b_user_lightyield);
  parameters->SetBranchAddress("crystal_lateral_depolishing",&footer.crystal_lateral_depolishing,&b_crystal_lateral_depolishing);
  parameters->SetBranchAddress("crystal_exit_depolishing",&footer.crystal_exit_depolishing,&b_crystal_exit_depolishing);
  parameters->SetBranchAddress("crystal_inner_cladding_fraction",&footer.crystal_inner_cladding_fraction,&b_crystal_inner_cladding_fraction);
  parameters->SetBranchAddress("crystal_outer_cladding_fraction",&footer.crystal_outer_cladding_fraction,&b_crystal_outer_cladding_fraction);
  parameters->SetBranchAddress("event_time_cut",&footer.event_time_cut,&b_event_time_cut);
  parameters->SetBranchAddress("LAPPD_thickness",&footer.LAPPD_thickness,&b_LAPPD_thickness);
  parameters->SetBranchAddress("SingleSeparationLayer",&footer.SingleSeparationLayer,&b_SingleSeparationLayer);
  parameters->SetBranchAddress("esr_on_positive_exit",&footer.esr_on_positive_exit,&b_esr_on_positive_exit);
  parameters->SetBranchAddress("esr_on_negative_exit",&footer.esr_on_negative_exit,&b_esr_on_negative_exit);
  parameters->SetBranchAddress("cladding_abs_length_scale_factor",&footer.cladding_abs_length_scale_factor,&b_cladding_abs_length_scale_factor);
  parameters->SetBranchAddress("fileName",&footer.v_fileName);
  parameters->SetBranchAddress("volume_file",&footer.v_volume_file);
  parameters->SetBranchAddress("light_guide_file",&footer.v_light_guide_file);
  parameters->SetBranchAddress("ecal_map_file",&footer.v_ecal_map_file);
  parameters->SetBranchAddress("ecal_map_histo",&footer.v_ecal_map_histo);
  parameters->SetBranchAddress("input_filename",&footer.v_input_filename);
  parameters->SetBranchAddress("input_treename",&footer.v_input_treename);
  parameters->SetBranchAddress("gps_instructions_file",&footer.v_gps_instructions_file,&b_gps_instructions_file);
  parameters->SetBranchAddress("logging_volume",&footer.v_logging_volume);
  parameters->SetBranchAddress("pre_volume",&footer.v_pre_volume);
  parameters->SetBranchAddress("AbsName",&footer.v_AbsName);
  parameters->SetBranchAddress("CellPositionX",&footer.v_CellPositionX,&b_CellPositionX);
  parameters->SetBranchAddress("CellPositionY",&footer.v_CellPositionY,&b_CellPositionY);
  parameters->SetBranchAddress("CellPositionZ",&footer.v_CellPositionZ,&b_CellPositionZ);
  parameters->SetBranchAddress("CellXelements",&footer.v_CellXelements,&b_CellXelements);
  parameters->SetBranchAddress("CellYelements",&footer.v_CellYelements,&b_CellYelements);
  parameters->SetBranchAddress("CellCrystalSizeX",&footer.v_CellCrystalSizeX,&b_CellCrystalSizeX);
  parameters->SetBranchAddress("CellCrystalSizeY",&footer.v_CellCrystalSizeY,&b_CellCrystalSizeY);
  parameters->SetBranchAddress("CellCrystalSizeZ",&footer.v_CellCrystalSizeZ,&b_CellCrystalSizeZ);
  parameters->SetBranchAddress("CellCrystalPitchX",&footer.v_CellCrystalPitchX,&b_CellCrystalPitchX);
  parameters->SetBranchAddress("CellCrystalPitchY",&footer.v_CellCrystalPitchY,&b_CellCrystalPitchY);
  parameters->SetBranchAddress("CellMaterial",&footer.v_CellMaterial,&b_CellMaterial);
  parameters->SetBranchAddress("CellAirLayer",&footer.v_CellAirLayer,&b_CellAirLayer);
  parameters->SetBranchAddress("CellHoleShape",&footer.v_CellHoleShape,&b_CellHoleShape);
  parameters->SetBranchAddress("CellCrystalShape",&footer.v_CellCrystalShape,&b_CellCrystalShape);
  parameters->SetBranchAddress("CellCrystalCladding",&footer.v_CellCrystalCladding,&b_CellCrystalCladding);
  parameters->SetBranchAddress("CellIntGapMaterial",&footer.v_CellIntGapMaterial,&b_CellIntGapMaterial);
  parameters->SetBranchAddress("CellStaggering",&footer.v_CellStaggering,&b_CellStaggering);
  parameters->SetBranchAddress("CellStaggeringAxis",&footer.v_CellStaggeringAxis,&b_CellStaggeringAxis);
  parameters->SetBranchAddress("CellStaggeringSize",&footer.v_CellStaggeringSize,&b_CellStaggeringSize);
  parameters->SetBranchAddress("CellStaggeringParity",&footer.v_CellStaggeringParity,&b_CellStaggeringParity);
  parameters->SetBranchAddress("CellStaggeringRemove",&footer.v_CellStaggeringRemove,&b_CellStaggeringRemove);
  parameters->SetBranchAddress("LAPPD_layers",&footer.v_LAPPD_layers,&b_LAPPD_layers);
  parameters->SetBranchAddress("LAPPD_materials",&footer.v_LAPPD_materials,&b_LAPPD_materials);
  parameters->SetBranchAddress("ecal_type",&footer.v_ecal_type,&b_ecal_type);
  parameters->SetBranchAddress("crystalMaterialList",&footer.v_crystalMaterialList,&b_crystalMaterialList);
      
  

 

}

void fillRelevant(relevant_data_t &relevant,Data_t &parameters)
{ 
  // fill
  strncpy(relevant.fileName, parameters.v_fileName->c_str(), sizeof(relevant.fileName));
  // std::cout << "config   = " << std::string(relevant.fileName) << std::endl;
  // std::cout << "config2  = " << std::string(parameters.v_fileName->c_str()) << std::endl;
  relevant.defaultCut = parameters.defaultCut;
  relevant.ecal_position = parameters.ecal_position;
  strncpy(relevant.logging_volume, parameters.v_logging_volume->c_str(), sizeof(relevant.logging_volume));
  strncpy(relevant.pre_volume, parameters.v_pre_volume->c_str(), sizeof(relevant.pre_volume));
  relevant.gapSize = parameters.gapSize;
  relevant.abs_interface_extraGap = parameters.abs_interface_extraGap;
  relevant.gap_abs_interface_material = parameters.gap_abs_interface_material;
  relevant.gap_interface_readout_material = parameters.gap_interface_readout_material;
  relevant.AbsSizeZ = parameters.AbsSizeZ;
  relevant.AbsPositionX = parameters.AbsPositionX;
  relevant.AbsPositionY = parameters.AbsPositionY;
  relevant.AbsPositionZ = parameters.AbsPositionZ;
  relevant.AbsMaterial = parameters.AbsMaterial;
  relevant.absorber_reflectivity = parameters.absorber_reflectivity;
  relevant.absorber_specular_lobe = parameters.absorber_specular_lobe;
  relevant.absorber_specular_spike = parameters.absorber_specular_spike;
  relevant.absorber_backscatter = parameters.absorber_backscatter;
  relevant.absorber_sigma_alpha = parameters.absorber_sigma_alpha;
  relevant.W_fraction = parameters.W_fraction;
  relevant.esrTransmittance = parameters.esrTransmittance;
  relevant.cell_separation_type = parameters.cell_separation_type;
  relevant.cell_separator_position = parameters.cell_separator_position;
  relevant.separation_thickness = parameters.separation_thickness;
  relevant.separation_material = parameters.separation_material;
  relevant.abs_length_scale_factor = parameters.abs_length_scale_factor;
  relevant.crystal_lateral_depolishing = parameters.crystal_lateral_depolishing;
  relevant.crystal_exit_depolishing = parameters.crystal_exit_depolishing;
  relevant.crystal_inner_cladding_fraction = parameters.crystal_inner_cladding_fraction;
  relevant.crystal_outer_cladding_fraction = parameters.crystal_outer_cladding_fraction;
  relevant.esr_on_positive_exit = parameters.esr_on_positive_exit;
  relevant.esr_on_negative_exit = parameters.esr_on_negative_exit;

  // find number of sections 
  std::vector<double> zPos;
  
  for(int i = 0 ; i < parameters.v_CellPositionZ->size(); i++)
  { 
    bool zCellIsAlreadyThere = false;
    for(unsigned int iz = 0; iz < zPos.size(); iz++)
    {
      if(parameters.v_CellPositionZ->at(i) == zPos[iz])
      {
        zCellIsAlreadyThere = true;
      }
    }
    if(!zCellIsAlreadyThere)
    {
      zPos.push_back(parameters.v_CellPositionZ->at(i));
    }
  }
  // std::cout << "sections = " << zPos.size() << std::endl;
  relevant.sections = zPos.size(); 

  relevant.CellCrystalSizeZ[0] = parameters.v_CellCrystalSizeZ->at(0);
  relevant.CellMaterial[0] = parameters.v_CellMaterial->at(0);
  relevant.CellAirLayer[0] = parameters.v_CellAirLayer->at(0);
  relevant.CellHoleShape[0] = parameters.v_CellHoleShape->at(0);
  relevant.CellCrystalShape[0] = parameters.v_CellCrystalShape->at(0);
  relevant.CellCrystalCladding[0] = parameters.v_CellCrystalCladding->at(0);
  relevant.CellIntGapMaterial[0] = parameters.v_CellIntGapMaterial->at(0);
  
  if(zPos.size() == 2)
  {
    relevant.CellCrystalSizeZ[1]    = parameters.v_CellCrystalSizeZ   ->at(parameters.v_CellCrystalSizeZ->size() - 1)   ;
    relevant.CellMaterial[1]        = parameters.v_CellMaterial       ->at(parameters.v_CellMaterial->size() - 1)         ;
    relevant.CellAirLayer[1]        = parameters.v_CellAirLayer       ->at(parameters.v_CellAirLayer->size() - 1)        ;
    relevant.CellHoleShape[1]       = parameters.v_CellHoleShape      ->at(parameters.v_CellHoleShape->size() - 1)       ;
    relevant.CellCrystalShape[1]    = parameters.v_CellCrystalShape   ->at(parameters.v_CellCrystalShape->size() - 1)    ;
    relevant.CellCrystalCladding[1] = parameters.v_CellCrystalCladding->at(parameters.v_CellCrystalCladding->size() - 1)  ;
    relevant.CellIntGapMaterial[1]  = parameters.v_CellIntGapMaterial ->at(parameters.v_CellIntGapMaterial->size() - 1)  ;
  }
}

bool AreSame(double &a, double &b)
{
  double EPSILON = 0.0001; // so 0.1 micron since in our sim distances are in mm
  return fabs(a - b) < EPSILON;
}

// compare relevant entries in two relevant_data_t structures, one taken from the 
// energy deposition file, one from the calibration file coupled to that module 
int compareRelevant(relevant_data_t &deposition,relevant_data_t &calibration)
{
  std::cout << "Deposition config   = " << std::string(deposition.fileName) << std::endl;
  std::cout << "Calibration config  = " << std::string(calibration.fileName) << std::endl;

  int ret = 0;
  
  if(!AreSame(deposition.defaultCut,calibration.defaultCut)) {ret++; std::cout << "WARNING: incompatible defaultCut found!" << std::endl;}
  if(strcmp(deposition.logging_volume,calibration.logging_volume) != 0) {ret++; std::cout << "WARNING: incompatible logging_volume found!" << std::endl;}
  if(strcmp(deposition.pre_volume,calibration.pre_volume) != 0) {ret++; std::cout << "WARNING: incompatible pre_volume found!" << std::endl;}
  if(!AreSame(deposition.gapSize,calibration.gapSize)) {ret++; std::cout << "WARNING: incompatible gapSize found!" << std::endl;}
  if(!AreSame(deposition.abs_interface_extraGap,calibration.abs_interface_extraGap)) {ret++; std::cout << "WARNING: incompatible abs_interface_extraGap found!" << std::endl;}
  if(deposition.gap_abs_interface_material != calibration.gap_abs_interface_material) {ret++; std::cout << "WARNING: incompatible gap_abs_interface_material found!" << std::endl;}
  if(deposition.gap_interface_readout_material != calibration.gap_interface_readout_material) {ret++; std::cout << "WARNING: incompatible gap_interface_readout_material found!" << std::endl;}
  if(!AreSame(deposition.AbsSizeZ,calibration.AbsSizeZ)) {ret++; std::cout << "WARNING: incompatible AbsSizeZ found!" << std::endl;}
  if(!AreSame(deposition.AbsPositionX,calibration.AbsPositionX)) {ret++; std::cout << "WARNING: incompatible AbsPositionX found!" << std::endl;}
  if(!AreSame(deposition.AbsPositionY,calibration.AbsPositionY)) {ret++; std::cout << "WARNING: incompatible AbsPositionY found!" << std::endl;}
  if(!AreSame(deposition.AbsPositionZ,calibration.AbsPositionZ)) {ret++; std::cout << "WARNING: incompatible AbsPositionZ found!" << std::endl;}
  if(deposition.AbsMaterial != calibration.AbsMaterial) {ret++; std::cout << "WARNING: incompatible AbsMaterial found!" << std::endl;}
  if(!AreSame(deposition.absorber_reflectivity,calibration.absorber_reflectivity)) {ret++; std::cout << "WARNING: incompatible absorber_reflectivity found!" << std::endl;}
  if(!AreSame(deposition.absorber_specular_lobe,calibration.absorber_specular_lobe)) {ret++; std::cout << "WARNING: incompatible absorber_specular_lobe found!" << std::endl;}
  if(!AreSame(deposition.absorber_specular_spike,calibration.absorber_specular_spike)) {ret++; std::cout << "WARNING: incompatible absorber_specular_spike found!" << std::endl;}
  if(!AreSame(deposition.absorber_backscatter,calibration.absorber_backscatter)) {ret++; std::cout << "WARNING: incompatible absorber_backscatter found!" << std::endl;}
  if(!AreSame(deposition.absorber_sigma_alpha,calibration.absorber_sigma_alpha)) {ret++; std::cout << "WARNING: incompatible absorber_sigma_alpha found!" << std::endl;}
  
  // check W_fraction only if AbsMaterial == 6
  if(deposition.AbsMaterial == 6)
  {
    if(!AreSame(deposition.W_fraction,calibration.W_fraction)) {ret++; std::cout << "WARNING: incompatible W_fraction found!" << std::endl;}
  }

  if(deposition.cell_separation_type != calibration.cell_separation_type) {ret++; std::cout << "WARNING: incompatible cell_separation_type found!" << std::endl;}
  // these need to be checked only if separation type == 2
  if(deposition.cell_separation_type == 2)
  {
    if(!AreSame(deposition.esrTransmittance,calibration.esrTransmittance)) {ret++; std::cout << "WARNING: incompatible esrTransmittance found!" << std::endl;}
    if(!AreSame(deposition.cell_separator_position,calibration.cell_separator_position)) {ret++; std::cout << "WARNING: incompatible cell_separator_position found!" << std::endl;}
    if(!AreSame(deposition.separation_thickness,calibration.separation_thickness)) {ret++; std::cout << "WARNING: incompatible separation_thickness found!" << std::endl;}
    if(deposition.separation_material != calibration.separation_material) {ret++; std::cout << "WARNING: incompatible separation_material found!" << std::endl;}
  }
  
  if(!AreSame(deposition.abs_length_scale_factor,calibration.abs_length_scale_factor)) {ret++; std::cout << "WARNING: incompatible abs_length_scale_factor found!" << std::endl;}
  if(!AreSame(deposition.crystal_lateral_depolishing,calibration.crystal_lateral_depolishing)) {ret++; std::cout << "WARNING: incompatible crystal_lateral_depolishing found!" << std::endl;}
  if(!AreSame(deposition.crystal_exit_depolishing,calibration.crystal_exit_depolishing)) {ret++; std::cout << "WARNING: incompatible crystal_exit_depolishing found!" << std::endl;}
  if(deposition.esr_on_positive_exit != calibration.esr_on_positive_exit) {ret++; std::cout << "WARNING: incompatible esr_on_positive_exit found!" << std::endl;}
  if(deposition.esr_on_negative_exit != calibration.esr_on_negative_exit) {ret++; std::cout << "WARNING: incompatible esr_on_negative_exit found!" << std::endl;}
  
  if(deposition.sections != calibration.sections) {ret++; std::cout << "WARNING: incompatible number of sections found!" << std::endl;}
  // compare only relevant part of 2 entries arrays
  if(!AreSame(deposition.CellCrystalSizeZ[0],calibration.CellCrystalSizeZ[0])) {ret++; std::cout << "WARNING: incompatible CellCrystalSizeZ[0] found!" << std::endl;}
  if(deposition.CellMaterial[0] != calibration.CellMaterial[0]) {ret++; std::cout << "WARNING: incompatible CellMaterial[0] found!" << std::endl;}
  if(!AreSame(deposition.CellAirLayer[0],calibration.CellAirLayer[0])) {ret++; std::cout << "WARNING: incompatible CellAirLayer[0] found!" << std::endl;}
  if(deposition.CellHoleShape[0] != calibration.CellHoleShape[0]) {ret++; std::cout << "WARNING: incompatible CellHoleShape[0] found!" << std::endl;}
  if(deposition.CellCrystalShape[0] != calibration.CellCrystalShape[0]) {ret++; std::cout << "WARNING: incompatible CellCrystalShape[0] found!" << std::endl;}
  if(deposition.CellCrystalCladding[0] != calibration.CellCrystalCladding[0]) {ret++; std::cout << "WARNING: incompatible CellCrystalCladding[0] found!" << std::endl;}
  if(!AreSame(deposition.CellIntGapMaterial[0],calibration.CellIntGapMaterial[0])) {ret++; std::cout << "WARNING: incompatible CellIntGapMaterial[0] found!" << std::endl;}
  if(deposition.sections == calibration.sections)
  {
    if(deposition.sections == 2)
    {
      if(!AreSame(deposition.CellCrystalSizeZ[1],calibration.CellCrystalSizeZ[1])) {ret++; std::cout << "WARNING: incompatible CellCrystalSizeZ[1] found!" << std::endl;}
      if(deposition.CellMaterial[1] != calibration.CellMaterial[1]) {ret++; std::cout << "WARNING: incompatible CellMaterial[1] found!" << std::endl;}
      if(!AreSame(deposition.CellAirLayer[1],calibration.CellAirLayer[1])) {ret++; std::cout << "WARNING: incompatible CellAirLayer[1] found!" << std::endl;}
      if(deposition.CellHoleShape[1] != calibration.CellHoleShape[1]) {ret++; std::cout << "WARNING: incompatible CellHoleShape[1] found!" << std::endl;}
      if(deposition.CellCrystalShape[1] != calibration.CellCrystalShape[1]) {ret++; std::cout << "WARNING: incompatible CellCrystalShape[1] found!" << std::endl;}
      if(deposition.CellCrystalCladding[1] != calibration.CellCrystalCladding[1]) {ret++; std::cout << "WARNING: incompatible CellCrystalCladding[1] found!" << std::endl;}
      if(!AreSame(deposition.CellIntGapMaterial[1],calibration.CellIntGapMaterial[1])) {ret++; std::cout << "WARNING: incompatible CellIntGapMaterial[1] found!" << std::endl;}
    }
  }

  // check cladding fraction only if there is cladding 
  if(deposition.CellCrystalCladding[0])
  {
    if(!AreSame(deposition.crystal_inner_cladding_fraction,calibration.crystal_inner_cladding_fraction)) {ret++; std::cout << "WARNING: incompatible crystal_inner_cladding_fraction found!" << std::endl;}
    if(!AreSame(deposition.crystal_outer_cladding_fraction,calibration.crystal_outer_cladding_fraction)) {ret++; std::cout << "WARNING: incompatible crystal_outer_cladding_fraction found!" << std::endl;}
  }
  
  
  return ret;

}