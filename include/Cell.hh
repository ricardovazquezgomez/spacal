// Class for cell elements
// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#ifndef Cell_h
#define Cell_h 1

#include <iostream>
#include <string>
#include <fstream>
#include <utility>

#include "ConfigFile.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"
// #include "DetectorParameterisation.hh"

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"


struct row_t
{
  int elements;
  std::vector<G4double> hole_dimension_x;
  std::vector<G4double> hole_dimension_y;
  std::vector<G4double> hole_dimension_z;
  std::vector<G4double> hole_pos_x;
  std::vector<G4double> hole_pos_y;
  std::vector<G4double> hole_pos_z;

  std::vector<G4double> crystal_dimension_x;
  std::vector<G4double> crystal_dimension_y;
  std::vector<G4double> crystal_dimension_z;
  std::vector<G4double> crystal_pos_x;
  std::vector<G4double> crystal_pos_y;
  std::vector<G4double> crystal_pos_z;

};

// class for Module Cells
class Cell
{
public:

  //! ctor
  Cell ();

  //! dtor
  ~Cell ();

  //! construct method
  G4VPhysicalVolume* Construct () ;

  //! other methods
  // generic identification of cell
  void           SetID(int aNum){id = aNum;};
  G4int          GetID(){return id;}
  void           SetName(G4String aString){name = aString;};
  G4String       GetName(){return name;};

  // position in absorber volume
  void           SetPosition(G4double x, G4double y ,G4double z){pos_x = x;pos_y = y; pos_z = z;};
  G4double       GetPositionX(){return pos_x;};
  G4double       GetPositionY(){return pos_y;};
  G4double       GetPositionZ(){return pos_z;};
  G4double       GetSizeX(){return size_x;};
  G4double       GetSizeY(){return size_y;};
  G4double       GetSizeZ(){return size_z;};

  // number of elements in x and y (i.e. number of holes)
  void           SetXelements(G4int aNum){xElements = aNum;};
  G4int          GetXelements(){return xElements;};
  void           SetYelements(G4int aNum){yElements = aNum;};
  G4int          GetYelements(){return yElements;};

  // methods for row-staggering
  void           SetStaggering(G4int aNum){staggering = aNum;};
  G4int          GetStaggering(){return staggering;};
  void           SetStaggeringAxis(G4int aNum){staggeringAxis = aNum;};
  G4int          GetStaggeringAxis(){return staggeringAxis;};
  void           SetStaggeringSize(G4double aNum){staggeringSize = aNum;};
  G4double       GetStaggeringSize(){return staggeringSize;};
  void           SetStaggeringParity(G4int aNum){staggeringParity = aNum;};
  G4int          GetStaggeringParity(){return staggeringParity;};
  void           SetStaggeringRemove(G4int aNum){staggeringRemove = aNum;};
  G4int          GetStaggeringRemove(){return staggeringRemove;};

  // methods for crystals
  // air layer between crystal and hole walls
  void           SetAirLayer(G4double aNum){airLayer = aNum;};
  G4double       GetAirLayer(){return airLayer;};
  // crystal dimensions
  void           SetNominalCrystalDimensions (G4double x, G4double y, G4double z)
                                             {crystal_nominal_size_x = x;
                                              crystal_nominal_size_y = y;
                                              crystal_nominal_size_z = z;};
  G4double       GetRealCrystalDimensionX(){return crystal_real_size_x;};
  G4double       GetRealCrystalDimensionY(){return crystal_real_size_y;};
  G4double       GetRealCrystalDimensionZ(){return crystal_real_size_z;};
  // pitch
  void           SetCrystalPitch(G4double x,G4double y)
                                {crystal_pitch_x = x;
                                 crystal_pitch_y = y;};
  G4double       GetCrystalPitchX(){return crystal_pitch_x;};
  G4double       GetCrystalPitchY(){return crystal_pitch_y;};
  // materials
  void           SetCrystalMaterial (G4int mat,double user_lightyield,double scaleFactor,double cladding_abs_length_scale_factor);
  G4Material*    GetCrystalMaterial(){return CrystalMaterial;};
  void           SetCrystalMaterialIdentifier(G4int mat);
  G4int          GetCrystalMaterialIdentifier(){return CrystalMaterialIdentifier;};
  // material of the gap between front-back faces of the fibers and either outside or other cell
  void           SetExtGapMaterial (G4int mat);
  G4Material*    GetExtGapMaterial(){return ExternalGapMaterial;};
  void           SetIntGapMaterial (G4int mat);
  G4Material*    GetIntGapMaterial(){return InternalGapMaterial;};
  // colors
  G4Colour       GetCrystalColor(){return crystalColor;};
  G4Colour       GetOuterCladdingColor(){return OuterCladdingColor;};
  G4Colour       GetInnerCladdingColor(){return InnerCladdingColor;};
  G4Colour       GetExtGapMaterialColor(){return ExternalGapMaterialColor;};
  G4Colour       GetIntGapMaterialColor(){return InternalGapMaterialColor;};

  // hole dimensions
  G4double       GetHoleDimensionX(){return hole_size_x;};
  G4double       GetHoleDimensionY(){return hole_size_y;};
  G4double       GetHoleDimensionZ(){return hole_size_z;};

  G4double       GetGapSize(){return gapSize;};
  void           SetGapSize(G4double aNum){gapSize = aNum;};
  
  // hole and crystal shapes
  void           SetHoleShape(G4int aNum){holeShape = aNum;};
  void           SetCrystalShape(G4int aNum){crystalShape = aNum;};
  void           SetCrystalCladding(G4bool aNum){crystalCladding = aNum;};
  G4int          GetHoleShape(){return holeShape;};
  G4int          GetCrystalShape(){return crystalShape;};
  G4bool         GetCrystalCladding(){return crystalCladding;};

  G4Material*    GetInnerCladdingMaterial(){return InnerCladdingMaterial;};
  G4Material*    GetOuterCladdingMaterial(){return OuterCladdingMaterial;};
  G4double       GetInnerCladdingDiameter(){return cladding_inner_diameter;};
  G4double       GetOuterCladdingDiameter(){return cladding_outer_diameter;};
  G4double       GetCoreDiameter(){return cladding_core_diameter;};

  void           SetCrystalInnerCladdingFraction(G4double aNum){crystal_inner_cladding_fraction = aNum;};
  void           SetCrystalOuterCladdingFraction(G4double aNum){crystal_outer_cladding_fraction = aNum;};
  G4double       GetInnerCladdingFraction(){return crystal_inner_cladding_fraction;};
  G4double       GetOuterCladdingFraction(){return crystal_outer_cladding_fraction;};



  //layers of the cell
  // vector of layers
  std::vector<row_t> layer;
  // method to produce the cell structure
  void MakeCellStruture();
  void CalculateCellDimensions();
  G4int GetNumberOfLayers(){return layer.size();};
  row_t GetLayer(G4int id){return layer[id];};

private:

  // identifiers
  G4int id;
  std::string name;

  // dimensions, material etc
  // position in absorber volume
  G4double pos_x;             // position in x [mm]
  G4double pos_y;             // position in y [mm]
  G4double pos_z;             // position in z [mm]
  G4double size_x;             // size in x [mm]
  G4double size_y;             // size in y [mm]
  G4double size_z;             // size in z [mm]
  // elements in x and y
  G4int xElements;
  G4int yElements;

  G4int holeShape;
  G4int crystalShape;
  G4bool crystalCladding;


  G4double gapSize;

  G4int    staggering;
  G4int    staggeringAxis;
  G4double staggeringSize;
  G4int    staggeringParity;
  G4int    staggeringRemove;

  //crystals
  G4double airLayer;
  G4double crystal_nominal_size_x;             // cry dimensions in x [mm]
  G4double crystal_nominal_size_y;             // cry dimensions in y [mm]
  G4double crystal_nominal_size_z;             // cry dimensions in z [mm]
  G4double crystal_real_size_x;             // cry dimensions in x [mm]
  G4double crystal_real_size_y;             // cry dimensions in y [mm]
  G4double crystal_real_size_z;             // cry dimensions in z [mm]
  G4double crystal_pos_x;
  G4double crystal_pos_y;
  G4double crystal_pos_z;
  G4double crystal_pitch_x;            // cry pitch in x [mm]
  G4double crystal_pitch_y;            // cry pitch in y [mm]
  // material
  // G4Material* AbMaterial;
  G4Material* CrystalMaterial;
  G4Material* OuterCladdingMaterial;
  G4Material* InnerCladdingMaterial;
  G4int       CrystalMaterialIdentifier;
  G4Material* ExternalGapMaterial;
  G4Material* InternalGapMaterial;
  G4Colour crystalColor;
  G4Colour InnerCladdingColor;
  G4Colour OuterCladdingColor;
  G4Colour ExternalGapMaterialColor;
  G4Colour InternalGapMaterialColor;

  G4double hole_size_x;
  G4double hole_size_y;
  G4double hole_size_z;

  G4double cladding_core_diameter;
  G4double cladding_inner_diameter;
  G4double cladding_outer_diameter;

  G4double crystal_inner_cladding_fraction;
  G4double crystal_outer_cladding_fraction;

} ;

#endif /*Cell_h*/
