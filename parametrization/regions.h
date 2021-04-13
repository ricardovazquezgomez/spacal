// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

struct region_t
{
  // definition of region structure
  Float_t min[3];
  Float_t max[3];
  int type;
  int sections;
  float separation_z;
  float position_x;
  float position_y;
  float position_z;
};

// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

int FindMaterial(TTree* elementTree)
{
  int         elementID;
  int         elementType;
  int         elementMaterial;
  Float_t     elementX;
  Float_t     elementY;
  Float_t     elementZ;
  Float_t     elementDX;
  Float_t     elementDY;
  Float_t     elementDZ;
  Float_t     elementSeparationZ;
  int         elementSections;
  TBranch    *b_elementID;
  TBranch    *b_elementType;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  TBranch    *b_elementSeparationZ;
  TBranch    *b_elementSections;
  elementTree->SetBranchAddress("ID"           ,&elementID          ,&b_elementID);
  elementTree->SetBranchAddress("type"         ,&elementType        ,&b_elementType);
  elementTree->SetBranchAddress("material"     ,&elementMaterial    ,&b_elementMaterial);
  elementTree->SetBranchAddress("x"            ,&elementX           ,&b_elementX);
  elementTree->SetBranchAddress("y"            ,&elementY           ,&b_elementY);
  elementTree->SetBranchAddress("z"            ,&elementZ           ,&b_elementZ);
  elementTree->SetBranchAddress("dx"           ,&elementDX          ,&b_elementDX);
  elementTree->SetBranchAddress("dy"           ,&elementDY          ,&b_elementDY);
  elementTree->SetBranchAddress("dz"           ,&elementDZ          ,&b_elementDZ);
  elementTree->SetBranchAddress("separation_z" ,&elementSeparationZ ,&b_elementSeparationZ);
  elementTree->SetBranchAddress("sections"     ,&elementSections    ,&b_elementSections);
  elementTree->GetEvent(0);
  return elementMaterial;
}

int FindMaterialOld(TTree* elementTree)
{
  int         elementID;
  int         elementMaterial;
  Double_t       elementX;
  Double_t       elementY;
  Double_t       elementZ;
  Double_t       elementDX;
  Double_t       elementDY;
  Double_t       elementDZ;
  TBranch    *b_elementID;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  elementTree->SetBranchAddress("ID"       ,&elementID      ,&b_elementID);
  elementTree->SetBranchAddress("material" ,&elementMaterial,&b_elementMaterial);
  elementTree->SetBranchAddress("x"        ,&elementX       ,&b_elementX);
  elementTree->SetBranchAddress("y"        ,&elementY       ,&b_elementY);
  elementTree->SetBranchAddress("z"        ,&elementZ       ,&b_elementZ);
  elementTree->SetBranchAddress("dx"       ,&elementDX      ,&b_elementDX);
  elementTree->SetBranchAddress("dy"       ,&elementDY      ,&b_elementDY);
  elementTree->SetBranchAddress("dz"       ,&elementDZ      ,&b_elementDZ);
  elementTree->GetEvent(0);
  return elementMaterial;
}

std::vector<region_t> CalculateRegions(TTree* elementTree)
{
  // routine to calculate the boundaries of regions loaded from simulation output file
  int         elementID;
  int         elementType;
  int         elementMaterial;
  Float_t     elementX;
  Float_t     elementY;
  Float_t     elementZ;
  Float_t     elementDX;
  Float_t     elementDY;
  Float_t     elementDZ;
  Float_t     elementSeparationZ;
  int         elementSections;
  TBranch    *b_elementID;
  TBranch    *b_elementType;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  TBranch    *b_elementSeparationZ;
  TBranch    *b_elementSections;
  elementTree->SetBranchAddress("ID"           ,&elementID          ,&b_elementID);
  elementTree->SetBranchAddress("type"         ,&elementType        ,&b_elementType);
  elementTree->SetBranchAddress("material"     ,&elementMaterial    ,&b_elementMaterial);
  elementTree->SetBranchAddress("x"            ,&elementX           ,&b_elementX);
  elementTree->SetBranchAddress("y"            ,&elementY           ,&b_elementY);
  elementTree->SetBranchAddress("z"            ,&elementZ           ,&b_elementZ);
  elementTree->SetBranchAddress("dx"           ,&elementDX          ,&b_elementDX);
  elementTree->SetBranchAddress("dy"           ,&elementDY          ,&b_elementDY);
  elementTree->SetBranchAddress("dz"           ,&elementDZ          ,&b_elementDZ);
  elementTree->SetBranchAddress("separation_z" ,&elementSeparationZ ,&b_elementSeparationZ);
  elementTree->SetBranchAddress("sections"     ,&elementSections    ,&b_elementSections);
  std::vector<region_t> regions;
  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);

    region_t temp_region;
    temp_region.min[0] =  elementX - 0.5*elementDX;    // min x
    temp_region.min[1] =  elementY - 0.5*elementDY;    // min y
    temp_region.min[2] =  elementZ - 0.5*elementDZ;    // min z
    temp_region.max[0] =  elementX + 0.5*elementDX;    // max x
    temp_region.max[1] =  elementY + 0.5*elementDY;    // max y
    temp_region.max[2] =  elementZ + 0.5*elementDZ;    // max z
    temp_region.type       = elementType;
    temp_region.sections   = elementSections;
    temp_region.separation_z   = elementSeparationZ;
    temp_region.position_x     = elementX;
    temp_region.position_y     = elementY;
    temp_region.position_z     = elementZ;
    regions.push_back(temp_region);
  }
  return regions;
}


void CopyStructure(TTree* elementTree,TTree* outTree)
{
  // routine to calculate the boundaries of regions loaded from simulation output file
  int         elementID;
  int         elementType;
  int         elementMaterial;
  Float_t     elementX;
  Float_t     elementY;
  Float_t     elementZ;
  Float_t     elementDX;
  Float_t     elementDY;
  Float_t     elementDZ;
  Float_t     elementSeparationZ;
  int         elementSections;
  TBranch    *b_elementID;
  TBranch    *b_elementType;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  TBranch    *b_elementSeparationZ;
  TBranch    *b_elementSections;
  elementTree->SetBranchAddress("ID"           ,&elementID          ,&b_elementID);
  elementTree->SetBranchAddress("type"         ,&elementType        ,&b_elementType);
  elementTree->SetBranchAddress("material"     ,&elementMaterial    ,&b_elementMaterial);
  elementTree->SetBranchAddress("x"            ,&elementX           ,&b_elementX);
  elementTree->SetBranchAddress("y"            ,&elementY           ,&b_elementY);
  elementTree->SetBranchAddress("z"            ,&elementZ           ,&b_elementZ);
  elementTree->SetBranchAddress("dx"           ,&elementDX          ,&b_elementDX);
  elementTree->SetBranchAddress("dy"           ,&elementDY          ,&b_elementDY);
  elementTree->SetBranchAddress("dz"           ,&elementDZ          ,&b_elementDZ);
  elementTree->SetBranchAddress("separation_z" ,&elementSeparationZ ,&b_elementSeparationZ);
  elementTree->SetBranchAddress("sections"     ,&elementSections    ,&b_elementSections);

  Int_t    out_ID       ;
  Int_t    out_Type     ;
  Int_t    out_Material ;
  Float_t  out_X        ;
  Float_t  out_Y        ;
  Float_t  out_Z        ;
  Float_t  out_DX       ;
  Float_t  out_DY       ;
  Float_t  out_DZ       ;
  Float_t  out_SeparationZ       ;
  Int_t  out_Sections       ;
  outTree   ->Branch("ID"       , &out_ID          , "ID/I");
  outTree   ->Branch("type"     , &out_Type        , "type/I");
  outTree   ->Branch("material" , &out_Material    , "material/I");
  outTree   ->Branch("x"        , &out_X           , "x/F");
  outTree   ->Branch("y"        , &out_Y           , "y/F");
  outTree   ->Branch("z"        , &out_Z           , "z/F");
  outTree   ->Branch("dx"       , &out_DX          , "dx/F");
  outTree   ->Branch("dy"       , &out_DY          , "dy/F");
  outTree   ->Branch("dz"       , &out_DZ          , "dz/F");
  outTree   ->Branch("separation_z" , &out_SeparationZ , "separation_z/F");
  outTree   ->Branch("sections" , &out_Sections , "sections/I");
  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);
    out_ID      = elementID;
    out_Type    = elementType;
    out_Material= elementMaterial;
    out_X       = elementX;
    out_Y       = elementY;
    out_Z       = elementZ;
    out_DX      = elementDX;
    out_DY      = elementDY;
    out_DZ      = elementDZ;
    out_SeparationZ = elementSeparationZ;
    out_Sections = elementSections;
    outTree->Fill();
  }
}

TH2I* ComputeDetectorMap(TTree* moduleTree,TTree* cellTree, int front_back) // front = 0, z negative, back = 1, z positive
{
  TH2I *detector_ID_map = NULL;
  Float_t cell_bin_size = 5.; // in mm
  // set var and adresses
  // MODULE TTREE
  int         moduleID;
  int         moduleType;
  int         moduleMaterial;
  Float_t     moduleX;
  Float_t     moduleY;
  Float_t     moduleZ;
  Float_t     moduleDX;
  Float_t     moduleDY;
  Float_t     moduleDZ;
  Float_t     moduleSeparationZ;
  int         moduleSections;
  TBranch    *b_moduleID;
  TBranch    *b_moduleType;
  TBranch    *b_moduleMaterial;
  TBranch    *b_moduleX;
  TBranch    *b_moduleY;
  TBranch    *b_moduleZ;
  TBranch    *b_moduleDX;
  TBranch    *b_moduleDY;
  TBranch    *b_moduleDZ;
  TBranch    *b_moduleSeparationZ;
  TBranch    *b_moduleSections;
  moduleTree->SetBranchAddress("ID"           ,&moduleID          ,&b_moduleID);
  moduleTree->SetBranchAddress("type"         ,&moduleType        ,&b_moduleType);
  moduleTree->SetBranchAddress("material"     ,&moduleMaterial    ,&b_moduleMaterial);
  moduleTree->SetBranchAddress("x"            ,&moduleX           ,&b_moduleX);
  moduleTree->SetBranchAddress("y"            ,&moduleY           ,&b_moduleY);
  moduleTree->SetBranchAddress("z"            ,&moduleZ           ,&b_moduleZ);
  moduleTree->SetBranchAddress("dx"           ,&moduleDX          ,&b_moduleDX);
  moduleTree->SetBranchAddress("dy"           ,&moduleDY          ,&b_moduleDY);
  moduleTree->SetBranchAddress("dz"           ,&moduleDZ          ,&b_moduleDZ);
  moduleTree->SetBranchAddress("separation_z" ,&moduleSeparationZ ,&b_moduleSeparationZ);
  moduleTree->SetBranchAddress("sections"     ,&moduleSections    ,&b_moduleSections);
  // CELL TTREE
  int         cellID;
  int         cellType;
  int         cellMaterial;
  Float_t     cellX;
  Float_t     cellY;
  Float_t     cellZ;
  Float_t     cellDX;
  Float_t     cellDY;
  Float_t     cellDZ;
  Float_t     cellSeparationZ;
  int         cellSections;
  TBranch    *b_cellID;
  TBranch    *b_cellType;
  TBranch    *b_cellMaterial;
  TBranch    *b_cellX;
  TBranch    *b_cellY;
  TBranch    *b_cellZ;
  TBranch    *b_cellDX;
  TBranch    *b_cellDY;
  TBranch    *b_cellDZ;
  TBranch    *b_cellSeparationZ;
  TBranch    *b_cellSections;
  cellTree->SetBranchAddress("ID"           ,&cellID          ,&b_cellID);
  cellTree->SetBranchAddress("type"         ,&cellType        ,&b_cellType);
  cellTree->SetBranchAddress("material"     ,&cellMaterial    ,&b_cellMaterial);
  cellTree->SetBranchAddress("x"            ,&cellX           ,&b_cellX);
  cellTree->SetBranchAddress("y"            ,&cellY           ,&b_cellY);
  cellTree->SetBranchAddress("z"            ,&cellZ           ,&b_cellZ);
  cellTree->SetBranchAddress("dx"           ,&cellDX          ,&b_cellDX);
  cellTree->SetBranchAddress("dy"           ,&cellDY          ,&b_cellDY);
  cellTree->SetBranchAddress("dz"           ,&cellDZ          ,&b_cellDZ);
  cellTree->SetBranchAddress("separation_z" ,&cellSeparationZ ,&b_cellSeparationZ);
  cellTree->SetBranchAddress("sections"     ,&cellSections    ,&b_cellSections);


  // determine map size
  Float_t xMapMin = +INFINITY;
  Float_t xMapMax = -INFINITY;
  Float_t yMapMin = +INFINITY;
  Float_t yMapMax = -INFINITY;
  Float_t cell_size_x = cell_bin_size ; // fixed to 5 mm for the histogram
  Float_t cell_size_y = cell_bin_size ; // fixed to 5 mm for the histogram
  // find min and max in x an y
  for(int i = 0 ; i < moduleTree->GetEntries() ; i++)
  {
    moduleTree->GetEvent(i);
    Float_t elMinX =  moduleX - 0.5*moduleDX;    // min x
    Float_t elMinY =  moduleY - 0.5*moduleDY;    // min y
    Float_t elMinZ =  moduleZ - 0.5*moduleDZ;    // min z
    Float_t elMaxX =  moduleX + 0.5*moduleDX;    // max x
    Float_t elMaxY =  moduleY + 0.5*moduleDY;    // max y
    Float_t elMaxZ =  moduleZ + 0.5*moduleDZ;    // max z
    if(elMinX  < xMapMin)
    {
      xMapMin = elMinX;
    }
    if(elMinY  < yMapMin)
    {
      yMapMin = elMinY;
    }
    if(elMaxX  > xMapMax)
    {
      xMapMax = elMaxX;
    }
    if(elMaxY  > yMapMax)
    {
      yMapMax = elMaxY;
    }
  }

  // calc number of bins
  int nBinsX = (int) ((xMapMax - xMapMin)/cell_size_x);
  int nBinsY = (int) ((yMapMax - yMapMin)/cell_size_y);

  // std::cout << "xMapMin = " << xMapMin << std::endl;
  // std::cout << "xMapMax = " << xMapMax << std::endl;
  // std::cout << "yMapMin = " << yMapMin << std::endl;
  // std::cout << "yMapMax = " << yMapMax << std::endl;
  // std::cout << "nBinsX  = " << nBinsX  << std::endl;
  // std::cout << "nBinsY  = " << nBinsY  << std::endl;

  // create the histo map

  std::stringstream sname;
  sname << "detector_ID_map_";
  if(front_back == 0)
  {
    sname << "front";
  }
  else
  {
    sname << "back";
  }
  detector_ID_map = new TH2I(sname.str().c_str(),sname.str().c_str(),nBinsX,xMapMin,xMapMax,nBinsY,yMapMin,yMapMax);
  detector_ID_map->GetXaxis()->SetTitle("x[mm]");
  detector_ID_map->GetYaxis()->SetTitle("y[mm]");

  long int counter = 0;
  for(int iMod = 0 ; iMod < moduleTree->GetEntries() ; iMod++) // run on modules
  {
    moduleTree->GetEvent(iMod);
    // for each modules, run on cells
    for(int iCell = 0 ; iCell < cellTree->GetEntries() ; iCell++) // run on modules
    {
      // find cell center, in absolute x and y
      cellTree->GetEvent(iCell);

      // take only cells in positive or negative side wrt module separation, according to front_back
      bool keepEvent = false;
      if(front_back == 0) // front, so z negative, so cellZ must be more negative than module separation
      {
        if(cellZ < moduleSeparationZ)
        {
          keepEvent = true;
        }
      }
      else // back, so z positive, so cellZ must be more positive than module separation
      {
        if(cellZ > moduleSeparationZ)
        {
          keepEvent = true;
        }
      }
      // take only cell events of this module type
      if(cellType != moduleType)
      {
        keepEvent = false;
      }

      if(keepEvent)
      {

        // fill all histogram points from center - 0.5 size to center + 0.5 size with same number
        // calc how many cell bits of 5x5 mm2 are in this cell
        int nCellBinsX = (int) ((cellDX)/cell_size_x);
        int nCellBinsY = (int) ((cellDY)/cell_size_y);
        for(int iBin = 0; iBin < nCellBinsX; iBin++)
        {
          for(int jBin = 0; jBin < nCellBinsY; jBin++)
          {
            Float_t bin_center_x = moduleX + cellX - ((nCellBinsX-1)/2.0)*cell_size_x + iBin*cell_size_x;
            Float_t bin_center_y = moduleY + cellY - ((nCellBinsY-1)/2.0)*cell_size_y + jBin*cell_size_y;
            detector_ID_map->Fill(bin_center_x,bin_center_y,counter);
          }
        }

      }
      counter++; // counter has to increase for each cell, regardless of keep or not
    }

  }





  return detector_ID_map;
}


TH2I* ComputeModuleMap(TTree* elementTree)
{
  TH2I *histo_modules = NULL;
  int         elementID;
  int         elementType;
  int         elementMaterial;
  Float_t     elementX;
  Float_t     elementY;
  Float_t     elementZ;
  Float_t     elementDX;
  Float_t     elementDY;
  Float_t     elementDZ;
  Float_t     elementSeparationZ;
  int         elementSections;
  TBranch    *b_elementID;
  TBranch    *b_elementType;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  TBranch    *b_elementSeparationZ;
  TBranch    *b_elementSections;
  elementTree->SetBranchAddress("ID"           ,&elementID          ,&b_elementID);
  elementTree->SetBranchAddress("type"         ,&elementType        ,&b_elementType);
  elementTree->SetBranchAddress("material"     ,&elementMaterial    ,&b_elementMaterial);
  elementTree->SetBranchAddress("x"            ,&elementX           ,&b_elementX);
  elementTree->SetBranchAddress("y"            ,&elementY           ,&b_elementY);
  elementTree->SetBranchAddress("z"            ,&elementZ           ,&b_elementZ);
  elementTree->SetBranchAddress("dx"           ,&elementDX          ,&b_elementDX);
  elementTree->SetBranchAddress("dy"           ,&elementDY          ,&b_elementDY);
  elementTree->SetBranchAddress("dz"           ,&elementDZ          ,&b_elementDZ);
  elementTree->SetBranchAddress("separation_z" ,&elementSeparationZ ,&b_elementSeparationZ);
  elementTree->SetBranchAddress("sections"     ,&elementSections    ,&b_elementSections);

  Float_t xMapMin = +INFINITY;
  Float_t xMapMax = -INFINITY;
  Float_t yMapMin = +INFINITY;
  Float_t yMapMax = -INFINITY;
  Float_t module_size_x = 0. ;
  Float_t module_size_y = 0. ;
  // find min and max in x an y
  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);

    // all modules with the same dimensions! so just take the first
    if(i==0)
    {
      module_size_x = elementDX;
      module_size_y = elementDY;
    }

    Float_t elMinX =  elementX - 0.5*elementDX;    // min x
    Float_t elMinY =  elementY - 0.5*elementDY;    // min y
    Float_t elMinZ =  elementZ - 0.5*elementDZ;    // min z
    Float_t elMaxX =  elementX + 0.5*elementDX;    // max x
    Float_t elMaxY =  elementY + 0.5*elementDY;    // max y
    Float_t elMaxZ =  elementZ + 0.5*elementDZ;    // max z


    if(elMinX  < xMapMin)
    {
      xMapMin = elMinX;
    }
    if(elMinY  < yMapMin)
    {
      yMapMin = elMinY;
    }
    if(elMaxX  > xMapMax)
    {
      xMapMax = elMaxX;
    }
    if(elMaxY  > yMapMax)
    {
      yMapMax = elMaxY;
    }

  }


  // calc number of bins
  int nBinsX = (int) ((xMapMax - xMapMin)/module_size_x);
  int nBinsY = (int) ((yMapMax - yMapMin)/module_size_y);

  // std::cout << "xMapMin = " << xMapMin << std::endl;
  // std::cout << "xMapMax = " << xMapMax << std::endl;
  // std::cout << "yMapMin = " << yMapMin << std::endl;
  // std::cout << "yMapMax = " << yMapMax << std::endl;
  // std::cout << "nBinsX  = " << nBinsX  << std::endl;
  // std::cout << "nBinsY  = " << nBinsY  << std::endl;

  // create the histo map

  histo_modules = new TH2I("module_ID_map","module_ID_map",nBinsX,xMapMin,xMapMax,nBinsY,yMapMin,yMapMax);
  histo_modules->GetXaxis()->SetTitle("x[mm]");
  histo_modules->GetYaxis()->SetTitle("y[mm]");
  // run on modules again, and fill the histo_modules
  long int counter = 0;
  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);
    histo_modules->Fill(elementX,elementY,counter);
    counter++;
  }

  return histo_modules;


}


TH2I* ComputeElementMap(TTree* elementTree,TString name,int front_back, int type, int &counter) // front = 0, z negative, back = 1, z positive. -1 means ignore (for modules)
{
  TH2I *histo = NULL;
  int         elementID;
  int         elementType;
  int         elementMaterial;
  Float_t     elementX;
  Float_t     elementY;
  Float_t     elementZ;
  Float_t     elementDX;
  Float_t     elementDY;
  Float_t     elementDZ;
  Float_t     elementSeparationZ;
  int         elementSections;
  TBranch    *b_elementID;
  TBranch    *b_elementType;
  TBranch    *b_elementMaterial;
  TBranch    *b_elementX;
  TBranch    *b_elementY;
  TBranch    *b_elementZ;
  TBranch    *b_elementDX;
  TBranch    *b_elementDY;
  TBranch    *b_elementDZ;
  TBranch    *b_elementSeparationZ;
  TBranch    *b_elementSections;
  elementTree->SetBranchAddress("ID"           ,&elementID          ,&b_elementID);
  elementTree->SetBranchAddress("type"         ,&elementType        ,&b_elementType);
  elementTree->SetBranchAddress("material"     ,&elementMaterial    ,&b_elementMaterial);
  elementTree->SetBranchAddress("x"            ,&elementX           ,&b_elementX);
  elementTree->SetBranchAddress("y"            ,&elementY           ,&b_elementY);
  elementTree->SetBranchAddress("z"            ,&elementZ           ,&b_elementZ);
  elementTree->SetBranchAddress("dx"           ,&elementDX          ,&b_elementDX);
  elementTree->SetBranchAddress("dy"           ,&elementDY          ,&b_elementDY);
  elementTree->SetBranchAddress("dz"           ,&elementDZ          ,&b_elementDZ);
  elementTree->SetBranchAddress("separation_z" ,&elementSeparationZ ,&b_elementSeparationZ);
  elementTree->SetBranchAddress("sections"     ,&elementSections    ,&b_elementSections);

  Float_t xMapMin = +INFINITY;
  Float_t xMapMax = -INFINITY;
  Float_t yMapMin = +INFINITY;
  Float_t yMapMax = -INFINITY;
  Float_t element_size_x = 0. ;
  Float_t element_size_y = 0. ;
  // find min and max in x an y
  bool firstFound = false;
  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);


    if(type != -1)
    {
      if(elementType != type)
      {
        continue;
      }
    }

    // all elements with the same dimensions! so just take the first
    if(!firstFound)
    {
      element_size_x = elementDX;
      element_size_y = elementDY;
      firstFound = true;
    }

    Float_t elMinX =  elementX - 0.5*elementDX;    // min x
    Float_t elMinY =  elementY - 0.5*elementDY;    // min y
    Float_t elMinZ =  elementZ - 0.5*elementDZ;    // min z
    Float_t elMaxX =  elementX + 0.5*elementDX;    // max x
    Float_t elMaxY =  elementY + 0.5*elementDY;    // max y
    Float_t elMaxZ =  elementZ + 0.5*elementDZ;    // max z


    if(elMinX  < xMapMin)
    {
      xMapMin = elMinX;
    }
    if(elMinY  < yMapMin)
    {
      yMapMin = elMinY;
    }
    if(elMaxX  > xMapMax)
    {
      xMapMax = elMaxX;
    }
    if(elMaxY  > yMapMax)
    {
      yMapMax = elMaxY;
    }

  }


  // calc number of bins
  int nBinsX = (int) round((xMapMax - xMapMin)/element_size_x);
  int nBinsY = (int) round((yMapMax - yMapMin)/element_size_y);

  // std::cout << "xMapMin = " << xMapMin << std::endl;
  // std::cout << "xMapMax = " << xMapMax << std::endl;
  // std::cout << "yMapMin = " << yMapMin << std::endl;
  // std::cout << "yMapMax = " << yMapMax << std::endl;
  // std::cout << "nBinsX  = " << nBinsX  << std::endl;
  // std::cout << "nBinsY  = " << nBinsY  << std::endl;

  // create the histo map

  histo = new TH2I(name,name,nBinsX,xMapMin,xMapMax,nBinsY,yMapMin,yMapMax);
  histo->GetXaxis()->SetTitle("x[mm]");
  histo->GetYaxis()->SetTitle("y[mm]");
  // run on modules again, and fill the histo_modules

  for(int i = 0 ; i < elementTree->GetEntries() ; i++)
  {
    elementTree->GetEvent(i);
    bool keepFB = false;
    bool keepType = false;
    if(front_back == -1)
    {
      keepFB = true;
    }
    else
    {
      if(front_back == 0) // front, so z negative, so cellZ must be more negative than module separation
      {
        if(elementZ < elementSeparationZ)
        {
          keepFB = true;
        }
      }
      else // back, so z positive, so cellZ must be more positive than module separation
      {
        if(elementZ > elementSeparationZ)
        {
          keepFB = true;
        }
      }
    }
    if(type == -1)
    {
      keepType = true;
    }
    else
    {
      if(elementType == type)
      {
        keepType = true;
      }
    }
    if(keepType && keepFB)
    {
      histo->Fill(elementX,elementY,counter);
      counter++;
    }
  }

  return histo;


}

struct maps_t
{
  std::vector<int> moduleTypes;
  // lots of maps
  std::map<int, int> moduleNumberMap;         // map from module type to number of modules
  std::map<int, int> moduleNumberOfBeforeMap; // map from module type to number of other modules
  std::map<int, int> moduleSectionsMap;       // map from module type to number of segments
  std::map<int, float> moduleMinMap;          // map from module type to min z
  std::map<int, float> moduleMaxMap;          // map from module type to max z
  std::map<int, float> moduleSeparationMap;   // map from module type to z of separation
  std::map<int, int> mapOfCellMaps; // map from showerModuleType to i-index of cell_map_per_module (then j is 0 or 1)

  TH2I *module_ID_map;
  TH2I ***cell_map_per_module;
};


// for the moment, this is working only for the calculateChiFactors program
int readStructure(TFile *_file0, maps_t &maps)
{
  // calc regions
  std::stringstream feedbackString;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| READING SPACAL STRUCTURE                                      |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  // TFile *_file0 = TFile::Open(listInputFiles[0].c_str());
  // TFile *inputFile = TFile::Open(inputFileName.c_str());
  TTree *modules = (TTree *)_file0->Get("modules");
  TTree *absorbers = (TTree *)_file0->Get("absorbers");
  TTree *cells = (TTree *)_file0->Get("cells");
  TTree *holes = (TTree *)_file0->Get("holes");
  TTree *fibres = (TTree *)_file0->Get("fibres");
  // TTree *shower = (TTree *)_file0->Get("shower");
  _file0->cd("Configuration");
  TNamed SeedNameD("Seed", ((TNamed *)gDirectory->Get("Seed"))->GetTitle());
  TNamed HostNameD("Hostname", ((TNamed *)gDirectory->Get("Hostname"))->GetTitle());
  TNamed PWDNameD("PWD", ((TNamed *)gDirectory->Get("PWD"))->GetTitle());
  TNamed ConfigNameD("ConfigFile", ((TNamed *)gDirectory->Get("ConfigFile"))->GetTitle());
  TNamed GpsNameD("GpsFile", ((TNamed *)gDirectory->Get("GpsFile"))->GetTitle());
  TNamed PrimariesNameD("Primaries", ((TNamed *)gDirectory->Get("Primaries"))->GetTitle());
  _file0->cd();

  // produce the module ID TH2I map
  int det_counter = 0;
  maps.module_ID_map = ComputeElementMap(modules, "module_ID_map", -1, -1, det_counter);
  // and the detector id map
  // TH2I *detector_ID_map_front = ComputeDetectorMap(modules,cells,0);
  // TH2I *detector_ID_map_back  = ComputeDetectorMap(modules,cells,1);

  std::vector<region_t> modules_regions = CalculateRegions(modules);
  std::vector<int> moduleTypes;
  std::vector<int> moduleSections;
  std::vector<float> moduleSeparationZ;
  std::vector<float> modulePositionZ;
  // find how many module types are involved
  for (int iMod = 0; iMod < modules_regions.size(); iMod++)
  {
    bool moduleTypeIsAlreadyThere = false;
    for (unsigned int im = 0; im < moduleTypes.size(); im++)
    {
      if (modules_regions[iMod].type == moduleTypes[im])
      {
        moduleTypeIsAlreadyThere = true;
      }
    }
    if (!moduleTypeIsAlreadyThere)
    {
      moduleTypes.push_back(modules_regions[iMod].type);
      // also get sections
      moduleSections.push_back(modules_regions[iMod].sections);
      moduleSeparationZ.push_back(modules_regions[iMod].separation_z);
      modulePositionZ.push_back(modules_regions[iMod].position_z);
    }
  }

  maps.moduleTypes = moduleTypes;

  

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    int modCounter = 0;
    for (int iMod = 0; iMod < modules_regions.size(); iMod++)
    {
      if (modules_regions[iMod].type == moduleTypes[i])
      {
        modCounter++;
      }
    }
    maps.moduleNumberMap.insert(std::make_pair(moduleTypes[i], modCounter));
  }

  // sum all modules
  // std::map<int,int>::iterator it_moduleNumberMap = moduleNumberMap.begin();

  int totModules = 0;
  for (std::map<int, int>::iterator it = maps.moduleNumberMap.begin(); it != maps.moduleNumberMap.end(); ++it)
  {
    maps.moduleNumberOfBeforeMap.insert(std::make_pair(it->first, totModules));
    totModules += it->second;
  }
  feedbackString << "Total number of modules                 = " << totModules << std::endl;

  int NofModuleTypes = moduleTypes.size();
  // int NofSpacalSegments = moduleTypes.size();
  feedbackString << "Number of module types                 = " << NofModuleTypes << std::endl;

  
  // produce 2 cell maps for each module type
  
  maps.cell_map_per_module = new TH2I **[NofModuleTypes];
  for (int i = 0; i < NofModuleTypes; i++)
  {
    maps.cell_map_per_module[i] = new TH2I *[2];
    int det_counter = 0;
    for (int j = 0; j < 2; j++)
    {
      std::stringstream cname;
      cname << "cell_map_";
      if (j == 0) // front
      {
        cname << "front_";
      }
      else // back
      {
        cname << "back_";
      }
      cname << moduleTypes[i];
      // if sections == 1, then keep always front_back
      int fb = j;
      if (moduleSections[i] == 1)
      {
        fb = -1;
      }
      maps.cell_map_per_module[i][j] = ComputeElementMap(cells, cname.str(), fb, moduleTypes[i], det_counter);
      maps.mapOfCellMaps.insert(std::make_pair(moduleTypes[i], i));
    }
  }

  // vector of found
  std::vector<std::vector<bool>> foundLimit;

  std::vector<region_t> cells_regions = CalculateRegions(cells);

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    maps.moduleSectionsMap.insert(std::make_pair(moduleTypes[i], moduleSections[i]));
    maps.moduleSeparationMap.insert(std::make_pair(moduleTypes[i], moduleSeparationZ[i] + modulePositionZ[i]));
    std::vector<bool> foundLimitForModule;
    foundLimitForModule.push_back(false);
    foundLimitForModule.push_back(false);
    foundLimit.push_back(foundLimitForModule);
  }

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    for (int iCell = 0; iCell < cells_regions.size(); iCell++)
    {
      if (foundLimit[i][0] && foundLimit[i][1])
      {
        break;
      }
      // look for a cell of this type
      if (cells_regions[iCell].type == moduleTypes[i])
      {
        if (maps.moduleSectionsMap.at(moduleTypes[i]) == 1)
        {
          // only 1 sections, just take max and min
          foundLimit[i][1] = true;
          maps.moduleMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
          foundLimit[i][0] = true;
          maps.moduleMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
        }
        else
        {
          // find mean
          float meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2]) / 2.0;
          if (meanCellZ > 0)
          {
            if (foundLimit[i][1] == false)
            {
              foundLimit[i][1] = true;
              maps.moduleMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
            }
          }
          else
          {
            if (foundLimit[i][0] == false)
            {
              foundLimit[i][0] = true;
              maps.moduleMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
            }
          }
        }
      }
    }
  }


  int last_cell_front;

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    feedbackString << "Type " << moduleTypes[i] << " has " << maps.moduleSectionsMap.at(moduleTypes[i]) << " sections " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << maps.moduleNumberMap.at(moduleTypes[i]) << " modules " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << maps.moduleNumberOfBeforeMap.at(moduleTypes[i]) << " modules before" << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " min " << maps.moduleMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " max " << maps.moduleMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " separation at " << maps.moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
    for(int j = 0; j < maps.moduleSectionsMap.at(moduleTypes[i]); j++)
    {
      feedbackString << "Type " << moduleTypes[i] << " section " << j << " has " << maps.cell_map_per_module[i][j]->GetEntries() << " cells " << std::endl;
      feedbackString << "Type " << moduleTypes[i] << " section " << j << " max cell id is " << maps.cell_map_per_module[i][j]->GetBinContent(maps.cell_map_per_module[i][j]->GetMaximumBin()) << std::endl;

      if(j == 0){
        last_cell_front = maps.cell_map_per_module[i][j]->GetBinContent(maps.cell_map_per_module[i][j]->GetMaximumBin());
      }
      

    }
    
  }
  
  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;
  return last_cell_front;

}
