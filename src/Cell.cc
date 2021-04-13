// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#include "Cell.hh"

#include "CreateTree.hh"

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

Cell::Cell()
{
  // constructor
  // doesn't actually do anything, all will be set by the methods
}

Cell::~Cell(){}


void Cell::SetCrystalMaterialIdentifier(G4int mat)
{
  CrystalMaterialIdentifier = mat;
}

void Cell::SetCrystalMaterial(G4int mat,double user_lightyield,double scaleFactor,double cladding_abs_length_scale_factor)
{

  SetCrystalMaterialIdentifier(mat);

  // set cladding materials, no matter what 
  OuterCladdingMaterial = MyMaterials::FlurPoly(cladding_abs_length_scale_factor);
  InnerCladdingMaterial = MyMaterials::PLEX(cladding_abs_length_scale_factor);
  

  CrystalMaterial = NULL ;
  if      ( mat == 1 ) CrystalMaterial = MyMaterials::Quartz () ;
  else if ( mat == 2 ) CrystalMaterial = MyMaterials::SiO2_Ce () ;
  else if ( mat == 3 ) CrystalMaterial = MyMaterials::DSB_Ce () ;
  else if ( mat == 4 ) CrystalMaterial = MyMaterials::LuAG_Ce () ;
  else if ( mat == 5 ) CrystalMaterial = MyMaterials::YAG_Ce (user_lightyield,scaleFactor) ;
  else if ( mat == 6 ) CrystalMaterial = MyMaterials::GAGG_Ce_Mg(user_lightyield,scaleFactor) ;
  else if ( mat == 7 ) CrystalMaterial = MyMaterials::Water() ;
  else if ( mat == 8 ) CrystalMaterial = MyMaterials::GFAG(user_lightyield,scaleFactor) ;
  else if ( mat == 9 ) CrystalMaterial = MyMaterials::GAGG_very_fast(user_lightyield,scaleFactor) ;
  else if ( mat == 10 ) CrystalMaterial = MyMaterials::GYAGG(user_lightyield,scaleFactor) ;
  else if ( mat == 11 ) CrystalMaterial = MyMaterials::GAGG_slow(user_lightyield,scaleFactor) ;
  else if ( mat == 12 ) CrystalMaterial = MyMaterials::Polystyrene(user_lightyield,scaleFactor) ;
  else if ( mat == 13 ) CrystalMaterial = MyMaterials::Lead_crystal() ;
  else if ( mat == 14 ) CrystalMaterial = MyMaterials::PureTungsten2_crystal() ;
  else if ( mat == 15 ) CrystalMaterial = MyMaterials::Al_crystal() ;
  else if ( mat == 16 ) CrystalMaterial = MyMaterials::Cu_crystal() ;
  else if ( mat == 17 ) CrystalMaterial = MyMaterials::GAGG_ILM(user_lightyield,scaleFactor) ;

  else
  {
    G4cerr << "<Cell>: Invalid fiber material specifier " << mat << G4endl ;
    exit (-1) ;
  }
  // G4cout << "Cell "<< name << " - Crystal material : "<< CrystalMaterial << G4endl ;

  G4Colour  white   (1.00, 1.00, 1.00) ;  // white
  G4Colour  grey    (0.50, 0.50, 0.50) ;  // grey
  G4Colour  black   (0.00, 0.00, 0.00) ;  // black
  G4Colour  red     (1.00, 0.00, 0.00) ;  // red
  G4Colour  green   (0.00, 1.00, 0.00) ;  // green
  G4Colour  blue    (0.00, 0.00, 1.00) ;  // blue
  G4Colour  cyan    (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  air     (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  magenta (1.00, 0.00, 1.00) ;  // magenta
  G4Colour  yellow  (1.00, 1.00, 0.00) ;  // yellow
  G4Colour  brass   (0.80, 0.60, 0.40) ;  // brass
  G4Colour  brown   (0.70, 0.40, 0.10) ;  // brown
  
  // and colors 
  OuterCladdingColor = brown;
  InnerCladdingColor = green;

  // Fibres Materials: 1) Quartz 2) SiO2:Ce 3) DSB:Ce 4) LuAG 5) YAG 6) GAGG 7) Water (insteed Polystyrene)
  if      ( mat == 1 ) crystalColor = grey;
  else if ( mat == 2 ) crystalColor = magenta;
  else if ( mat == 3 ) crystalColor = blue;
  else if ( mat == 4 ) crystalColor = green;
  else if ( mat == 5 ) crystalColor = yellow;
  else if ( mat == 6 ) crystalColor = red;
  else if ( mat == 7 ) crystalColor = black;
  else if ( mat == 8 ) crystalColor = blue;
  else if ( mat == 17 ) crystalColor = black;
  else // just set a default
  {
    crystalColor = brass;
  }

  // SetCrystalColor(mat);
}

void Cell::SetExtGapMaterial(G4int mat)
{
  ExternalGapMaterial = NULL ;
  if      ( mat == 1 ) ExternalGapMaterial = MyMaterials::Air () ;
  else if ( mat == 2 ) ExternalGapMaterial = MyMaterials::OpticalGrease () ;

  else
  {
    G4cerr << "<Cell>: Invalid front gap material specifier " << mat << G4endl ;
    exit (-1) ;
  }
  // G4cout << "Cell "<< name << " - front gap material : "<< ExternalGapMaterial << G4endl ;

  G4Colour  air     (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  yellow  (1.00, 1.00, 0.00) ;  // yellow
  G4Colour  red     (1.00, 0.00, 0.00) ;  // red
  G4Colour  green   (0.00, 1.00, 0.00) ;  // green
  G4Colour  blue    (0.00, 0.00, 1.00) ;  // blue

  if      ( mat == 1 ) ExternalGapMaterialColor = red;
  else if ( mat == 2 ) ExternalGapMaterialColor = blue;
  else // just set a default
  {
    ExternalGapMaterialColor = air;
  }
}

void Cell::SetIntGapMaterial(G4int mat)
{
  InternalGapMaterial = NULL ;
  if      ( mat == 1 ) InternalGapMaterial = MyMaterials::Air () ;
  else if ( mat == 2 ) InternalGapMaterial = MyMaterials::OpticalGrease () ;

  else
  {
    G4cerr << "<Cell>: Invalid back gap material specifier " << mat << G4endl ;
    exit (-1) ;
  }
  // G4cout << "Cell "<< name << " - back gap material : "<< InternalGapMaterial << G4endl ;

  G4Colour  air     (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  yellow  (1.00, 1.00, 0.00) ;  // yellow
  G4Colour  red     (1.00, 0.00, 0.00) ;  // red
  G4Colour  green   (0.00, 1.00, 0.00) ;  // green
  G4Colour  blue    (0.00, 0.00, 1.00) ;  // blue

  if      ( mat == 1 ) InternalGapMaterialColor = red;
  else if ( mat == 2 ) InternalGapMaterialColor = blue;
  else // just set a default
  {
    InternalGapMaterialColor = air;
  }
}


void Cell::MakeCellStruture()
{
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << xElements << " " <<yElements<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // set hole and crystal sizes
  // holes have the section of crystals, plus twice the air layer
  hole_size_x = crystal_nominal_size_x + 2.0 * airLayer;
  hole_size_y = crystal_nominal_size_y + 2.0 * airLayer;
  hole_size_z = crystal_nominal_size_z;

  // one special case, on the basis of crystal and hole shape: 
  // if crystal shape is round and hole shape is rectangular (crazy!)
  // then hole diameter needs to be crystal diagonal + 2*airLayer

  if(holeShape == 1 && crystalShape == 0)
  {
    G4double diagonal = sqrt(pow(crystal_nominal_size_x,2) + pow(crystal_nominal_size_y,2));
    hole_size_x = diagonal + 2.0 * airLayer;
    hole_size_y = diagonal + 2.0 * airLayer; // not needed for round hole, but assigned here for consistency
    hole_size_z = crystal_nominal_size_z;
  }

  // crystals have dimensions define by the user, in case the calorimeter is not split in two longitudinally
  // if it is split in two, then it can be a single absorber with each hole filled with two crystals, or
  // actually two abs separated by either air or a reflector.
  // since there cannot be more than 2 abs, at least for now, the discriminant is given by the z position of the
  // cell in the absorber space. If it is = 0, then it's not split in two.

  if(pos_z == 0)
  {
    // if not split in two, the actual size of the crystal in z is the nominal one
    crystal_real_size_x = crystal_nominal_size_x ;
    crystal_real_size_y = crystal_nominal_size_y ;
    crystal_real_size_z = crystal_nominal_size_z ;
    // and the position is just in the center of the hole
    crystal_pos_x = 0.0 ;
    crystal_pos_y = 0.0 ;
    crystal_pos_z = 0.0 ;
  }
  else
  {
    // otherwise, we need a air gap from the end of the crystal (towards the center of the module)
    // to the separation. so the crystal needs to be shorter
    crystal_real_size_x = crystal_nominal_size_x ;
    crystal_real_size_y = crystal_nominal_size_y ;
    crystal_real_size_z = crystal_nominal_size_z - gapSize ;
    // and the xy position is just in the center of the hole
    crystal_pos_x = 0.0 ;
    crystal_pos_y = 0.0 ;
    // but the z position has to be shifted, so that the face where there is no gap is at the end of the hole
    if(pos_z > 0)
    {
      crystal_pos_z = + 0.5 *gapSize ;
    }
    else
    {
      crystal_pos_z = -  0.5 *gapSize ;
    }
  }

  // cladding
  // define the cladding elements diameters, only if cladding is used 
  if(crystalCladding)
  {
    cladding_outer_diameter = crystal_real_size_x; // as usual, only x matters in round fibers, and it's the diameter
    cladding_inner_diameter = crystal_real_size_x - 2.0*(crystal_real_size_x*crystal_outer_cladding_fraction); // as per kuraray definition
    cladding_core_diameter  = cladding_inner_diameter - 2.0*(crystal_real_size_x*crystal_inner_cladding_fraction); // as per kuraray definition
  }

  // find out if there is staggering, and what it means
  if(staggering) // staggering has to be done
  {
    //group holes in layers. layers are then shifted to make the row-staggering
    if(staggeringAxis == 0) // staggering along x axis
    {

      for(int iY = 0; iY < yElements; iY++)
      {
        row_t temp_layer;
        double stagShift = 0.0;
        int layerElements = xElements;
        if(staggeringRemove == 1) // if staggering means "remove 1 hole, and shift"
        {
          if((iY % 2) == staggeringParity) // if this row is staggered
          {
            layerElements = layerElements - 1; // remove one element from layer, because of staggering
          }
        }
        else // staggering means shift without removing holes
        {
          if((iY % 2) == staggeringParity) // if this row is staggered
          {
            stagShift = staggeringSize;
          }
        }
        temp_layer.elements = layerElements;
        // set hole positions
        for(int iX = 0; iX < xElements; iX++)
        {
          // set the position of the centers of each hole in the row (crystals will be just
          // in the center of the holes)
          temp_layer.hole_pos_x.push_back( stagShift +(iX*crystal_pitch_x) - crystal_pitch_x*((layerElements - 1)/2.0));
          temp_layer.hole_pos_y.push_back((iY*crystal_pitch_y) - crystal_pitch_y*((yElements - 1)/2.0));
          temp_layer.hole_pos_z.push_back(0.0);

          // and set the other parameters
          temp_layer.hole_dimension_x.push_back(hole_size_x);
          temp_layer.hole_dimension_y.push_back(hole_size_y);
          temp_layer.hole_dimension_z.push_back(hole_size_z);

          temp_layer.crystal_pos_x.push_back(crystal_pos_x);
          temp_layer.crystal_pos_y.push_back(crystal_pos_y);
          temp_layer.crystal_pos_z.push_back(crystal_pos_z);
          temp_layer.crystal_dimension_x.push_back(crystal_real_size_x);
          temp_layer.crystal_dimension_y.push_back(crystal_real_size_y);
          temp_layer.crystal_dimension_z.push_back(crystal_real_size_z);
        }

        layer.push_back(temp_layer);
      }

    }

    else if(staggeringAxis == 1)// staggering along y axis
    {

      for(int iX = 0; iX < xElements; iX++)
      {
        row_t temp_layer;
        double stagShift = 0.0;
        int layerElements = yElements;
        if(staggeringRemove == 1) // if staggering means "remove 1 hole, and shift"
        {
          if((iX % 2) == staggeringParity) // if this row is staggered
          {
            layerElements = layerElements - 1; // remove one element from layer, because of staggering
          }
        }
        else // staggering means shift without removing holes
        {
          if((iX % 2) == staggeringParity) // if this row is staggered
          {
            stagShift = staggeringSize;
          }
        }
        temp_layer.elements = layerElements;

        // set hole positions
        for(int iY = 0; iY < yElements; iY++)
        {
          // set the position of the centers of each hole in the row (crystals will be just
          // in the center of the holes)
          temp_layer.hole_pos_x.push_back((iX*crystal_pitch_x) - crystal_pitch_x*((xElements - 1)/2.0));
          temp_layer.hole_pos_y.push_back( stagShift +(iY*crystal_pitch_y) - crystal_pitch_y*((layerElements - 1)/2.0));
          temp_layer.hole_pos_z.push_back(0.0);

          temp_layer.hole_dimension_x.push_back(hole_size_x);
          temp_layer.hole_dimension_y.push_back(hole_size_y);
          temp_layer.hole_dimension_z.push_back(hole_size_z);

          temp_layer.crystal_pos_x.push_back(crystal_pos_x);
          temp_layer.crystal_pos_y.push_back(crystal_pos_y);
          temp_layer.crystal_pos_z.push_back(crystal_pos_z);
          temp_layer.crystal_dimension_x.push_back(crystal_real_size_x);
          temp_layer.crystal_dimension_y.push_back(crystal_real_size_y);
          temp_layer.crystal_dimension_z.push_back(crystal_real_size_z);
        }

        layer.push_back(temp_layer);
      }


    }
    else
    {
      G4cerr << "<Cell>: Invalid staggering axis " << staggeringAxis << G4endl ;
      exit (-1) ;
    }
    G4cout << "Cell "<< name << " - staggering over axis : "<< staggeringAxis << G4endl ;
  }
  else // no row-staggering
  {
    for(int iY = 0; iY < yElements; iY++)
    {
      row_t temp_layer;
      temp_layer.elements = xElements;
      for(int iX = 0; iX < xElements; iX++)
      {
        temp_layer.hole_pos_x.push_back((iX*crystal_pitch_x) - crystal_pitch_x*((xElements - 1)/2.0));
        temp_layer.hole_pos_y.push_back((iY*crystal_pitch_y) - crystal_pitch_y*((yElements - 1)/2.0) );
        temp_layer.hole_pos_z.push_back(0.0);

        temp_layer.hole_dimension_x.push_back(hole_size_x);
        temp_layer.hole_dimension_y.push_back(hole_size_y);
        temp_layer.hole_dimension_z.push_back(hole_size_z);

        temp_layer.crystal_pos_x.push_back(crystal_pos_x);
        temp_layer.crystal_pos_y.push_back(crystal_pos_y);
        temp_layer.crystal_pos_z.push_back(crystal_pos_z);
        temp_layer.crystal_dimension_x.push_back(crystal_real_size_x);
        temp_layer.crystal_dimension_y.push_back(crystal_real_size_y);
        temp_layer.crystal_dimension_z.push_back(crystal_real_size_z);
      }

      layer.push_back(temp_layer);
    }
  }

  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // for(int iLay = 0; iLay < layer.size(); iLay++)
  // {
  //   for(int i = 0 ; i < layer[iLay].elements; i++)
  //   {
  //     G4cout << "(" << layer[iLay].hole_pos_x[i] << ","
  //                   << layer[iLay].hole_pos_y[i] << ","
  //                   << layer[iLay].hole_pos_z[i] << ")  ";
  //   }
  //   G4cout << G4endl;
  //
  // }
  //
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
  // G4cout << "-----------------------------------------------"<< G4endl;
}


void Cell::CalculateCellDimensions()
{
  // take first elements of holes position x and y
  G4double current_min_x = layer[0].hole_pos_x[0];
  G4double current_min_y = layer[0].hole_pos_y[0];
  G4double current_max_x = layer[0].hole_pos_x[0];
  G4double current_max_y = layer[0].hole_pos_y[0];
  // run on layers vector
  for(int i = 0 ; i < layer.size(); i++)
  {
    // run on all elements
    for(int j = 0 ; j < layer[i].elements ; j++)
    {
      if(layer[i].hole_pos_x[j] < current_min_x)
      {
        current_min_x = layer[i].hole_pos_x[j];
      }
      if(layer[i].hole_pos_y[j] < current_min_y)
      {
        current_min_y = layer[i].hole_pos_y[j];
      }

      if(layer[i].hole_pos_x[j] > current_max_x)
      {
        current_max_x = layer[i].hole_pos_x[j];
      }
      if(layer[i].hole_pos_y[j] > current_max_y)
      {
        current_max_y = layer[i].hole_pos_y[j];
      }
    }
  }

  size_x = (current_max_x - current_min_x) + crystal_pitch_x ;
  size_y = (current_max_y - current_min_y) + crystal_pitch_y ;
  size_z = hole_size_z;
}
