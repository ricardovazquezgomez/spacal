
# **How to setup and run an ECAL simulation**

*M. Pizzichemi (marco.pizzichemi@cern.ch)*

This document will guide you through setting up and running a full hybrid simulation of an ECAL prototype.


## **Download and compile the simulation package**

Create a new folder where you will download and compile the Simulation package. In this document, we will call this folder **$SIMULATION**. You can actually set it as a variable in bash, if you prefer. This will make following this tutorial much easier. Suppose that the full path of the folder you chose is **/path/to/my/folder**, then

```
SIMULATION=/path/to/my/folder
```

will set a temporary variable in your bash terminal. Since it's temporary, so you will lose it if you close the terminal or logout, so if you prefer you can also save the value somewhere.

Initialize the folder and get the latest source code from this GitLab repository: here the example uses the SSH authentication cloning method, but you may choose another method depending on your GIT setup. For LXPLUS, the KRB5 method is recommended.

```
git init
git remote add origin ssh://git@gitlab.cern.ch:7999/spacal-rd/spacal-simulation.git
git fetch
git pull origin master
```

Compile the entire simulation package by running the **deploy.sh** script. If you are compiling on your local machine, you need to have a working version of ROOT and GEANT4 installed, and the proper environment variables set up. More specifically, these software packages have been used for development and testing:

```
g++/gcc 9.20
Geant4 10.06.p01
ROOT v6.20.02
Python 3.5.2
```

Then, from the same folder where you downloaded the source code, simply do

```
./deploy.sh
```

If you are compiling on LXPLUS, simply do:

```
./deploy.sh lxplus
```

In this latter case, the script will source the appropriate environment variables for you in the terminal before compiling (in particular, it will setup the LCG_97python3 environment).

The executables are now created in the new folder **./build**, and in particular the main Geant4 simulation executable will be

```
$SIMULATION/build/FibresCalo
```

If you run the executable without arguments, it will provide you with a short explanation on the different running modes:

```
Syntax for exec with gps, seed, and flux: FibresCalo <configuration file> <output file without extension> <gps file> <random seed> <flux file>
Syntax for exec with gps and seed: FibresCalo <configuration file> <output file without extension> <gps file> <random seed>
Syntax for exec with gps: FibresCalo <configuration file> <output file without extension> <gps file>
Syntax for exec: FibresCalo <configuration file> <output file>
Syntax for viz: FibresCalo <configuration file>
```

In this tutorial we will use the first running mode on LXPLUS.


## **Prepare simulation files**

You can now prepare the files needed to perform the simulation.

### **Setting of main config**
You need a base config file that will specify the ECAL configuration. You can take as example the

```
$SIMULATION/configuration/spacal_full_ecal.cfg
```

file. The configuration is explained in the file itself, however in brief this file describes a full ECAL based on 12x12 cm^2 section modules, divided in 6 regions. The regions differ for the module technology or structure.

It is important to underline here that whenever a simulation is based on a multiple module type structure, the .cfg config files needed are more than one. One config file, which we will call **main config**, will define the general properties of the calorimeter, and provide information on the **module config** files, which instead will describe the individual module types. This modular approach to config files is based on the idea that the module config files can be used to drive a full simulation for the characterizion of individual prototypes, and then can be reused without any modification as part of a full ECAL simulation.

The example is already set to hybrid simulation type (`simulationType = 2`). You can modify the ECAL configuration by altering the following part:

```
### ROOT file with x-y position of modules and types
# and name of TH2F histogram in ecal_map_file
ecal_map_file  = ECAL_LS4_option1.root   
ecal_map_histo = histo_ls4_option1       

### Configuration files for each module region
# list of module types corresponding to regions
# list of CAD drawing files corresponding to regions
modules_config   = |spacal_Pb_gagg_12x12cm.cfg|spacal_Pb_gagg_12x12cm.cfg|spacal_Pb_gagg_12x12cm.cfg|spacal_Pb_gagg_12x12cm.cfg|spacal_Pb_gagg_12x12cm.cfg|spacal_W_gagg_12x12cm.cfg|
ecal_type        = |1|2|3|4|5|6|
container_volume = |vol_1.obj|vol_2.obj|vol_3.obj|vol_4.obj|vol_5.obj|vol_6.obj|

### calorimeter dimensions in mm
calorimeter_size_x = 8000.     
calorimeter_size_y = 8000.     
calorimeter_size_z = 250.      
```

Pay attention to this **VERY IMPORTANT** suggestion:

```
NB/1: In general, it's better to specify ALL the file names in the config files with
the full path, when running on LXPLUS (for example /eos/experiment/spacal/Simulations/file.cfg)

For the sake of simplicity, the full paths are not specified in this guide.
Please pay attention to add them when you setup your real simulation.
```

You then need to provide a **map of the module positions and types**, under the form of a TH2F histogram. You can take as example the file

```
$SIMULATION/ecal_configurations/ECAL_LS4_option1.root  
```

You also need to provide the name of the TH2F to be used. The x-y positions of the bin centers in this histogram will be used as the coordinates of the centers of the modules in the ECAL simulation, while the bin content will define the module type.

You then need to provide a **configuration file**, a **type number**, and a **CAD drawing file** for each of the regions you want to include in the simulation. It is not necessary to include all the regions defined in ecal_map_file (this file could define for example 6 regions, but you could prefer to simulate only 2, as we will do in this example). It is however important that you assign the correct type number to the region config files, i.e. the same type number defined by the TH2F bin contents. Furthermore, you need to provide a proper CAD drawing for each region: **no checks** will be done by the simulation code to ensure that the modules are created within the correct CAD volume. You need to be sure of what you are doing!

You can use the files provided by the simulation package:

```
$SIMULATION/configurations/spacal_Pb_gagg_12x12cm.cfg
$SIMULATION/configurations/spacal_W_gagg_12x12cm.cfg
$SIMULATION/cad/vol_5.obj
$SIMULATION/cad/vol_6.obj
```

and for the sake of this example, set:

```
calorimeter_position = 0
ecal_map_file        = ECAL_LS4_option1.root   
ecal_map_histo       = histo_ls4_option1
modules_config       = |spacal_Pb_gagg_12x12cm.cfg|spacal_W_gagg_12x12cm.cfg|
ecal_type            = |5|6|
container_volume     = |vol_5.obj|vol_6.obj|
```

but remember, use the full paths for the files specified by the ecal_map_file, modules_config and container_volume keys.

**THE REST OF THE MAIN CONFIG FILE WILL JUST BE IGNORED BY THE SIMULATION.** Now you need to set the module config files

### **Setting of module config files**

The two example config files are already set for you (files spacal_W_gagg_12x12cm.cfg and spacal_Pb_gagg_12x12cm.cfg). Region 6 (the innermost) is structured as a SPACAL module based on W absorber and GAGG fibers, split in two regions 4 + 10 cm, with a Moliere radius of approximately 1.5 cm. Each module is split in 8x8 cells (in both front and back part). Region 5 is made of SPACAL modules longitudinally split in a 4 cm and a 10 cm section, GAGG fibers, but Pb absorber. As a consequence, the Moliere radius is around 3 cm and the cell size has been adjusted accordingly: each module is split in 4x4 cells in both front and back side. You can leave the original files as they are.

### **Setting up optical Calibrations files**

You will need run an optical calibration of each different ECAL module type defined in your simulation. In this example, we have defined two modules types, so we will need two optical calibrations. If the calibrations have already been produced, you can of course skip this step, but you need to be sure they have been produced with exactly the same crystal and absorber configurations. Again, be careful, the **calorimeter_position** key in the optical calibration files needs to match the **calorimeter_position** key in the **main config** file!

As explained in the READ.me file, the optical calibrations require a minimal configuration file, derived in this case from the corresponding module config file. You can take as templates the files

```
$SIMULATION/configurations/single_Pb_gagg_12x12cm.cfg
$SIMULATION/configurations/single_W_gagg_12x12cm.cfg
```  

If you have not made any modifications to the example module config files, these files are already ok.



### **Particle source**

You will need to provide a definition of particle source for your simulation. As examplained in the config files, you can either use a single particle source (`useGPS = true`) or take a flux of particle as source. Here we will assume that you will choose this latter case (for the former, you can check the README.md file). You need to set

```
useGPS = false
deltaZ      = -12690.5  ### actual range: 12600 - 12629                         - ignored if useGPS = true
```

in reality, deltaZ depends on the condition you used to generate the flux file we will specify later. In general, this value corresponds to -12620 minus the half length of the spacal absorber in z (in our case 70.5). You also need a gps file, in this case you can take the very simple file

```
$SIMULATION/gps/flux_gps.mac
```

that you will use later to start the simulation.



## **Visualize geometry**

To ensure everything is ok, you can visualize the geometry by running the simulation in visualization mode. You need to have a file like

```
$SIMULATION/vis.mac
```

in the same folder where you start the simulation. Be careful to the setting of visualization options in the config files. In particular, these keys:

```
### Visualization
# visualize or not the individual elements in viz mode. 1 = true, 0 = false
worldVisibility       = 0         # default = 0
caloVisibility        = 1         # default = 0
moduleVisibility      = 1         # default = 0
holeVisibility        = 0         # default = 1
crystalsVisibility    = 0         # default = 1
interfaceVisibility   = 0         # default = 1
gapsVisibility        = 0         # default = 1
readoutVisibility     = 0         # default = 1
absorberVisibility    = 0         # default = 0
containerVisibility   = 0         # default = 1
esrVisibility         = 0         # default = 1
lgVisibility          = 0         # default = 1
# wire frame or solid Visualization 1 = wireFrame, 0 = solid
wireFrame = 1
```

control which volumes are shown by the Geant4 visualization toolkit. If the simulation aims to build a very large calorimeter, the number of elements can easily diverge (in a single 12x12 cm SPACAL module of the innermost region, for example, the crystals volumes are more than 10k). The Geant4 visualization toolkit might therefore easily crash if you try to visualize all the volumes. With the options above it's safe, you will only visualize the calorimeter and the individual module boxes. Still, even in this minimal configuration, you will need to wait some minute if you are attempting to visualize the entire ECAL, with all the 6 regions filled with modules (this is not the case in this example, were we build only region 5 and tegion 6).

Now you can run

```
$SIMULATION/build/FibresCalo spacal_full_ecal.cfg
```

Check that the geometry makes sense, then move forward.



## **Prepare the simulation**

We will now prepare the simulation jobs. If you want to see how a simulation can be run locally, check out the README.md file. Here, for the rest of this document we will assume that you are setting up and running the simulation on the LXPLUS cluster. For what follows (in particular, to run the python scripts that create the jobs), you need to source the LCG_97python3 environment:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
```

Remember that the simulation consists of several steps. To summarize:

```
1) Optical calibration
2) Energy deposition by primary particles
3) Ray tracing propagation of Cerenkov photons
4) Parametrized propagation of Scintillation photons
5) Pulse formation
```

Step 5 needs to run on the combined results of step 4 and 3. Step 4, in turns, needs to run on the combined results of step 1 and 2. Finally, step 3 is usually performed together with step 2. As a result, in this HowTo we propose the following work flow: we group step 2 and 3, and execute them in parallel to step 1. Then, when these first tasks are all completed, we perform together steps 4 and 5, and save only the output of step 5 (i.e. the detector pulses). The reason for this can be found in the size of the output of step 4: these files contain the information on all the optical photons that reach the detectors, and are therefore by far the biggest output files produced by the simulation. Saving them to disk might quickly become unfeasible.


### **Optical calibration (step 1)**

As stated before, you need an optical calibration for each ECAL region. We assume you are now on LXPLUS, and you have created a **howto** folder in your workspace (meaning the AFS workspace, since the Condor jobs cannot be launched from the EOS folders). We will call this folder **$BASESIM**. You have also copied the configuration files in a subfolder, **$BASESIM/files**. Since the LHCb flux file can be quite big, we also assume that you left it where it is on EOS, so in

```
/eos/experiment/spacal/users/zhangy/FluxGammaB2KstGamma_MB_Merged.root
```

We can now create 3 other subfolders in **hotwo** (the first one with two subfolders)

```
mkdir -p OpticalCalibration/region6
mkdir -p OpticalCalibration/region5
mkdir HybridSimulation
mkdir PulseFormation

```  

We also need to prepare a folder for the output files. This needs to be in the EOS space: you will not have space on your personal AFS workspace for this. Let's say that you create a folder in **EOS** that we will refer to as **$OUTSIM** . You also created under this one all the subfolders as in the step above.

Assuming you have properly set the configuration files, you can now move to each of the 2 OpticalCalibration subfolders in **$BASESIM** and use the python scripts to prepare the Condor jobs.

```
cd $BASESIM/OpticalCalibration/region6
python3 $SIMULATION/parametrization/prepareOpticalCalibration.py --config $BASESIM/files/single_W_gagg_12x12cm.cfg --baseFolderJobs jobs --baseFolderOut $OUTSIM/OpticalCalibration/region6 --build $SIMULATION/build --queue workday --primaries 2000000 --jobs 1000 --xmin -0.5 --xmax 0.5 --ymin -0.5 --ymax 0.5 --zmin -70.5 --zmax 70.5 --xn 1 --yn 1 --zn 564 --emin 1.5 --emax 2.7 --en 24
cd $BASESIM/OpticalCalibration/region5
python3 $SIMULATION/parametrization/prepareOpticalCalibration.py --config $BASESIM/files/single_Pb_gagg_12x12cm.cfg --baseFolderJobs jobs --baseFolderOut $OUTSIM/OpticalCalibration/region5 --build $SIMULATION/build --queue workday --primaries 2000000 --jobs 1000 --xmin -0.5 --xmax 0.5 --ymin -0.5 --ymax 0.5 --zmin -70.5 --zmax 70.5 --xn 1 --yn 1 --zn 564 --emin 1.5 --emax 2.7 --en 24

```

See the README.md file for an explanation of the script parameters.


### **Geant4 simulation (step 2 and 3)**

We can now setup the Condor jobs for the energy deposition and Cerenkov propagation.

```
cd $BASESIM/HybridSimulation
python3 $SIMULATION/parametrization/prepareChain.py --baseFolderOut $OUTSIM/HybridSimulation --config $BASESIM/files/spacal_full_ecal.cfg --baseGPS /$BASESIM/files/flux_gps.mac --build $SIMULATION/build --listEnergy 1 --baseFolderJobs jobs --listCalibrations useless.data --listEvents 1 --events 1000 --listQueue workday --onlyEnDepo --useFlux /eos/experiment/spacal/users/zhangy/FluxGammaB2KstGamma_MB_Merged.root
```

## **Run the simulation**

We can now start the optical calibration and the Geant4 simulations.

For the **optical calibration**

```
cd $BASESIM/OpticalCalibration/region6/jobs
condor_submit jobs.sub

cd $BASESIM/OpticalCalibration/region5/jobs
condor_submit jobs.sub
```

For the **Geant4 simulation**
```
cd $BASESIM/HybridSimulation
./runAll_jobs.sh
```

Once all the **optical calibration** jobs are finished, move to the optical calibration output folders and produce the calibration files

```
cd $OUTSIM/OpticalCalibration/region6
condor_submit jobs.sub
$SIMULATION/build/mergeCalibration -f ./ -p out -o calibration_region6.data

cd $OUTSIM/OpticalCalibration/region5
condor_submit jobs.sub
$SIMULATION/build/mergeCalibration -f ./ -p out -o calibration_region5.data
```

this will produce two calibration binary files, that we will need for parametrized propagation.


## **Generate output pulses (steps 4 and 5)**

We can now perform the last two steps, i.e. the parametrized propagation of Scintillation photons and the pulse formation. You need a configuration file for pulse formation: take the file

```
$SIMULATION/signalFormation/SignalConfigFile.cfg
```

place it in **files** subfolder, rename it **SigConfig_FL4.cfg**. Then, modify just one line, i.e. set

```
UnifLossFac     = 4
```

You are now ready to prepare the last step:

```
cd $BASESIM/PulseFormation
python3 $SIMULATION/scripts/PrepareHybridAndSignalFormation_custom.py --baseFolderJobs jobs --baseFolderOut $OUTSIM/PulseFormation/ --baseFolderIn $OUTSIM/HybridSimulation/ --buildFolder $SIMULATION/build/ --config $BASESIM/files/SigConfig_FL4.cfg  --inputFolderKey GeV --inputFilePrefix out --listPFrac 0.3  --delete_hybrid  --listCalibrations $OUTSIM/OpticalCalibration/region5/calibration_region5.data $OUTSIM/OpticalCalibration/region6/calibration_region6.data --listTypes 5 6 --queue tomorrow --requestDisk 100GB
```
This command will prepare jobs doing hybrid propagation, readout simulation and time estimation of each cell.
The parameter ```--delete_hybrid``` deletes the output of the hybrid propagation at the end of the job, in order to save space. Additionally, the parameter ```--delete_signal``` can be used acting on the output of the readout simulation.

Then run the simulation with the script returned by the python code. For example

```
#########
Done. Please run /afs/cern.ch/work/m/mpizzich/simulations/jobs/howto/PulseFormation/jobs/submitALL_script.sh to submit all jobs. Bye!
#########
```








<!-- ### **Main config file description**


The first part of the config file deals with general simulation properties. The meaning of these keys is already explained in the config file comments (in general all lines starting with **#** are considered comments and are ignored by the simulation).

```
### remember that in general you can start a single, local, simulation from command line in 3 ways
# Syntax for exec with gps: FibresCalo <configuration file> <output file without extention> <gps file>
# Syntax for exec: FibresCalo <configuration file> <output file>
# Syntax for viz: FibresCalo <configuration file>

#####################
# GENERAL
#####################
### general configuration keys
# seed                  seed for random generator. if not given, default to -1, which means the seed is computed randomly
# checkOverlaps         ask Geant4 to check if there are overlapping volumes   - default = False
# B_field_intensity     intensity of magnetic field, in Tesla - default = 0
# defaultCut            production cut value in mm (if not given, default to 0.025 mm)
seed = -1
checkOverlaps = False
B_field_intensity = 0
defaultCut = 0.025
```

The primary particles used by the simulation are defined in the second and third section. The user need to choose whether to generate single primary particles (`useGPS = true`) or get the input particle flux from a root file.

```
#####################
# SOURCE
#####################
# source can be defined by a single particle generator, with a gps file,
# or by a ROOT file containing the a LHCb multi particle flux
# choice between the two possibilities is made by the key useGPS (if false, multi particle flux used)
# IMPORTANT - when using the multi particle flux, it is better to set calorimeter_position = 1 (see below)
useGPS      = true
input_file  = ./input_fakeGamma_fixedE1GeV.root     # ROOT file with flux                - ignored if useGPS = true
input_tree  = tree                                  # mane of ROOT TTree in input_file   - ignored if useGPS = true
skipEvents  = 0        # 0 #start from event skipEvents in the tree            - ignored if useGPS = true
deltaZ      = -12820.  ### actual range: 12600 - 12629                         - ignored if useGPS = true


#####################
# GPS files
#####################
### in general, GPS files are passed directly by command line, using the "exec with gps" mode, so this key can be omitted
### some example are in the gps folder
### gps_instructions_file             gps macro file. this will be ignored if the program is run in "exec with gps" mode (see above)
gps_instructions_file = gps_electron.mac
``` -->
