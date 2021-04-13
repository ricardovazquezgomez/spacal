---------------------------------------------------------------------------------------------------------

# **SPACAL SIMULATION FOR LHCb**

This code aims at studying the design of the SPACAL for LHCb U2 era.

It consists in a complete re-implementation of geometry and readout, to allow modular scale up of dimensions. It also implements a parametrization of optical transport to speed up simulation time. This modality is called **Hybrid-MC**.

For any questions, contact:
*Marco Pizzichemi (marco.pizzichemi@cern.ch)*

## **Instructions to compile the code**


To compile the simulation package on your local machine, execute

`./deploy.sh`

To compile the code on LXPLUS, execute

`./deploy.sh lxplus`




## **Geometry definition**

The geometry structure is modular. A module is defined in a configuration file, and can be placed in several positions in the simulation world coordinates. For the moment, only SPACAL modules are implemented, but work is undergoing for including also SHASHLIK geometry as well.

The simulation geometry is explained in the details in the configuration file examples that can be found in folder

```
./configuration
```

It has to be noted that there are two different ways to define the simulation geometry, as described in what follows.

### 1. Standalone modules

This modality is intended for testing and prototyping of different configurations. The user defines a module with an individual configuration file. See for example files **tb_spacal_W_gagg.cfg** or **spacal_W_gagg_12x12cm.cfg**. Optionally, the user can replicate the module on a grid (and also leave a hole in the center of the grid, to mimic a beampipe).

### 2. Full ECAL

In this modality, a full LHCb ECAL is simulated, or at least 1 or more regions of the ECAL. In this case, the user needs to provide more configuration files, as well as some auxiliary files to define the geometry. In particular:

```
1) A main config file: a config file that specifies the global characteristics of the ECAL (see for example the spacal_full_ecal.cfg file).
2) One or more module config files: these files are standard configuration files of SPACAL modules
3) A TH2F map of module types and positions (see an example in file ./ecal_configurations/ECAL_LS4_option1.root)
4) One or more CAD drawing of the region(s) the user wants to simulate (see the files .obj in ./cad folder)
```
For more details on the necessity for CAD drawings, see presentation in

`
./doc/UseOfCADregions.pdf
`

See file `spacal_full_ecal.cfg` for more explanations on how to setup a full ECAL simulation.


## **Running the simulation**

### Running locally

The simulation executable is located in the build folder

`./build/FibresCalo`

You can run locally the simulation in 3 modalities (they are also explained automatically if you run FibresCalo without arguments).

1. **Visualization only**:

`/path/to/build/FibresCalo <configuration file>`

2. **Standard exec**:

`/path/to/build/FibresCalo <configuration file> <output file>`

This modality will use the GPS file specified in the configuration file, at the key "gps_instructions_file".

3. **Exec with GPS**:

`/path/to/build/FibresCalo <configuration file> <output file> <gps file>`

This modality will use the GPS file specified in the command line argument <gps file>, ignoring the key "gps_instructions_file" written in the configuration file.

4. **Exec with GPS, and manually setting random seed**:

`/path/to/build/FibresCalo <configuration file> <output file without extension> <gps file> <random seed>`

This modality will use the GPS file specified in the command line argument <gps file>, ignoring the key "gps_instructions_file" written in the configuration file, and will also set the seed for random generator to the value specified by the user with the last command line argument.


```
ATTENTION: There is now the possibility to run the simulation using a multi particle source, that takes as input some ROOT files containing the LHCb particle flux. This can be activated setting

useGPS = false

in the config file. Other parameters will have to be properly set as well (see the config file).
In any case, the Exec procedure remains the same as explained above, but you will need to use a proper gps file (see for example gps/flux_gps.mac).
```


### Running full Hybrid-MC simulation on LXPLUS

A full Hybrid-MC simulation can be divided, to simplify the user activity, in 3 consecutive steps.

```
For all steps ahead, it's convenient to source the proper LXPLUS environment via LCG view:

source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh

Furthermore, we will define the main folder of this code as $CODEDIR

export CODEDIR=/path/to/git/folder

a working directory (aka the directory where lxplux jobs are created and submitted) as $WORKDIR

export WORKDIR=/path/to/work/folder

and an output directory (aka an EOS folder) as

export OUTDIR=/path/to/out/folder

```

1. **Prepare simulation geometry**:

The user should write a configuration file adapted to his/her needs, adapting the example found in $CODEDIR/configuration/ folder It is highly suggested to run some short test simulations locally. As the user will see in the configuration files, it can be useful to run such tests in simulationType = 0 or simulationType = 2 mode. Let's call the output of this step

```
mySpacal.cfg
```

If the user wants to run a full ECAL simulation, then in this step more than one config file will be needed, as long as some other auxiliary files, as explained in point 2 of "Geometry definition" section.


2. **Optical calibration**:

**NEW** The optical calibration file for individual modules does not need to the created manually anymore. The python script will generate automatically on the basis of the module configuration file, changing all the keys needed.

**NEW** The x,y,z min and max values do not need to be provided to `prepareOpticalCalibration.py` anymore, they are automatically calculated by the python script. Only the number of points in x,y,z needs to be set (if different from default)

**NEW** The `moduleZshift` key is now irrelevant to optical calibrations. The Hybrid Propagation step will automatically compensate for it, which means that an optical calibration produced with a given moduleZshift value (e.g. with the module centered in z = 0) can be used to propagate optical photons for a module whose energy deposition step has been performed with a different moduleZshift value (e.g. when the module is shifted in z, in a full ECAL configuration, to align the front faces of all sections).

The optical calibration needs to be executed just once per configuration, i.e. you need to perform it only if you change the mySpacal.cfg file. Otherwise, several hybrid simulations (described at STEP 3) can be performed with the output of the same optical calibration.

```
IMPORTANT: For the moment, there are some limitations:
1) Only one crystal material per SPACAL can be used in an optical calibration
2) No staggering should be used
```

The optical configuration run can be then prepared on LXPLUS as follows. Execute the script

`$CODEDIR/parametrization/prepareOpticalCalibration.py`

with the following command line arguments

```
--config            <value>    - module config file                (required)
--baseFolderJobs    <value>    - Calibration jobs folder           (required)
--baseFolderOut     <value>    - Calibration output folder         (required)
--build             <value>    - path to spacal build folder       (required)
--queue             <value>    - JobFlavour for condor             (required)
--primaries         <value>    - Number of opt. photons for each (x,y,z,energy) point (a good guess is 2000000) - default 2000000
--jobs              <value>    - Number of parallel jobs - default 1000
--xn                <value>    - Number of points in x scan - default 1
--yn                <value>    - Number of points in y scan - default 1
--zn                <value>    - Number of points in z scan - default 1
--emin              <value>    - Min optical photon energy value [eV] - default 1
--emax              <value>    - Max optical photon energy value [eV] - default 5
--en                <value>    - Number of points in energy scan - default 1
```

To give an example, the following command produces a calibration dataset for a module described by ./configurations/spacal_W_gagg_pitch1.67_12x12cm.cfg

```
python3 $CODEDIR/parametrization/prepareOpticalCalibration.py --config $CODEDIR/configurations/spacal_W_gagg_pitch1.67_12x12cm.cfg --baseFolderJobs jobs --baseFolderOut $OUTDIR/OpticalCalibrations/WGAGG  --build $CODEDIR/build --queue workday --primaries 2000000 --jobs 1000 --zn 604 --emin 1.5 --emax 2.7 --en 24
```

Typically one calibration point (i.e an individual x,y,z,energy choice) takes 30 sec to be completed, if the primaries shot are 100k. In this example, 2M optical photons are shot per calibration point, but divided in 1000 jobs, so 2000 optical photon per calibration point per job. Each job has 1 x 1 x 604 x 24 = 14496 calibration points (1 point in x, 1 point in y, 564 points in z, 24 energy points), so each job has 14496 x 2000 = 28992000 optical photons. As a reference, each job of the calibration example above lasted on average about 1h.

The prepareOpticalCalibration script will create a subfolder with the name given to --baseFolderJobs. It will also create the configuration file used by the optical calibration procedure, with a name created on the bases of the `--config`, adding a prefix `optiCali_`. At the same time, a `runAll_jobs.sh` will be created in the current folder: use it to submit the entire optical calibration with

`./runAll_jobs.sh`

Wait for the jobs to finish. The output files needs to be analyzed and **merged**, to create the calibration histograms. This can be performed using the program mergeCalibration, with the following command line parameters:

```
--folder            <value>    - path to the folder where the output files of optical calibration jobs are saved (--baseFolderOut in previous script)
--prefix            <value>    - prefix of optical calibration output files                 - default = out
--output            <value>    - name of calibration file output                            - default =  calibration.data
--pdfBins           <value>    - number of divisions for sampling of propagation time PDFs  - default = 5000
--pdfEnd            <value>    - end in ns of propagation time PDFs (begin is always 0ns)   - default = 10

N.B. the combination of pdfBins and pdfEnd determins the LSB of PDFs (default = 10ns/5000 = 2ps)
```

Pay attention to properly set the environment variables in your terminal in lxplus, before executing mergeCalibration:

```
ATTENTION: you need to properly set the environment variables to be able to run mergeCalibration on lxplus.
In particular, the same source command that is executed by ./deploy to compile the code need to be executed in the current terminal.

source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh

```

To give an example, the following command will merge the output of a calibration datased, if this output is located in current working directory and the file names all begin with "out":

```
$CODEDIR/build/mergeCalibration -f ./ -p out -o calibration.data
```

The calibration file will be named, in this example,

`calibration.data`

Store it somewhere. BTW, it's simple a binary file, so you cannot look inside with ROOT, but the histograms generated will be dumped in a `calibration_histograms.root` file, in order to allow some checks if needed.


3. **Hybrid-MC**:

In this step you will execute the simulation for generating the energy deposition maps and propagating the Cerenkov photons. Then you will perform the hybrid propagation of scintillation photons and finally you will analyze the output photons to produce "signals". Everything can be done with one script, namely

`./parametrization/prepareChain.py`

Go to a folder on lxplus, copy there the original mySpacal.cfg file here. At the same time, prepare a GPS file (let's call it gpsBase.mac) with the specifications of the beam you want to use. This file should include these commands (of course with the values you prefer) to allow you to choose the particle, position ,the angle, and the beam shape you want to simulate:

```
/gps/particle e-
/gps/pos/type Plane
/gps/pos/shape Square
/gps/pos/centre 1.5 1.5 -100  mm
/gps/pos/halfx 1  mm
/gps/pos/halfy 1  mm
/gps/direction -0.05236 -0.05236 0.9986
```

The /gps/energy and /run/beamOn commands will be written by the python script (you can omit them, or leave them in gpsBase.mac, they will be ignored). You can now run the script with the following commad line parameters:

```
--config            <value>    - SPACAL config file [mySpacal.cfg]         (required)
--baseFolderJobs    <value>    - Hybrid jobs folder                (required)
--baseFolderOut     <value>    - Hybrid output folder              (required)        
--baseGPS           <value>    - Base GPS file [gpsBase.mac]      (required)    
--build             <value>    - path to spacal build folder       (required)
--calibration       <value>    - Optical calibration binary file [calibration.data]
--events            <value>    - Total number of primary particles to shoot for a given energy   
--listEnergy        <value>    - List of energies to simulate      
--listEvents        <value>    - For each energy, how may primaries to shoot per job  
--listQueue         <value>    - For each energy, on which queue to submit the jobs
--listCalibrations  <value>    - list of optical calibration binary files [calibration.data]
--listTypes         <value>    - list of module type associated to a calibration file
--onlyEnDepo                   - Whether or not the program will stop at the outN.root files level, hence skipping the hybrid propagation.
```

Some parameters require further explanation. In particular the lists of energy, queues and events:

```
These lists and the parameter "events" determine how many jobs will be prepared.
For example, suppose you want to perform a simulation with energies (in GeV)

1 2 3 4 5 10 25 50 75 100

and shooting always 1000 electrons per energy. By your tests at STEP 1, you also know how long one event takes per each energy. So you want to run for example 25 primaries per job for energies from 1 GeV to 5 GeV, then 10 primaries for 10 GeV, 5 primaries for 25 GeV, and finally 2 primaries for 50, 75 and 100 GeV. You calculated that for 1 and 2 GeV the appropriate queue would be the "longlunch" one, while "workday" will be better suited for the others. You can therefore use the arguments like this:

--events 1000 --listEnergy 1 2 3 4 5 10 20 35 50 100 --listEvents 25 25 25 25 25 10 5 2 2 2 --listQueue longlunch longlunch workday workday workday workday workday workday workday workday

This will automatically create (1000/25) = 40 jobs (so 40 lines in the corresponding args.txt file) for the first 5 energies, (1000/10) = 100 jobs for 10 GeV and so on.
```

and the list of calibration files and types

```
These lists determine which optical calibration files to use, and what module type is associated with a specific calibration file.

As you can see in section "Geometry definition 2. Full ECAL" it is now possible to simulate more than one type of modules per simulation. The module types and their position in space is set via a TH2F histogram in a ROOT file, as explained in that section and in the config files.  
This of course forces the need to use more than one optical calibration in order to perform the hybrid propagation. The user will need to perform an optical calibration for each module type, then pass the optical calibration files as list to the flag --listCalibrations. Furthermore, the user will need to tell the module type of each calibration to the simulation program, by a list with the flag --listTypes. IMPORTANT: the lists need to be in the same order!

If the simulation is performed without using multiple modules, then the user can specify only one calibration file in --listCalibration, and omit the --listTypes flag.
```

The output of the simulation can be quite big. To explain it, we can say that the real hybrid simulation output is divided into 2 files per job:

```
outN.root           - output of energy deposition and Cerenkov propagation, for job N
hybridN.root        - output of parametrized scintillation propagation, energy deposition and Cerenkov propagation, for job N
```

Especially the second file can take a lot of disk space. Tipically, a scan of energies like the one above can take up to 6TB of disk space. The user can choose not to perform the hybrid propagation, by using the flag --onlyEnDepo.


1. From outN files, one can perform several (other) parametrized propagations of scintillation photons, on the basis of different optical calibrations (more refined ones, for example) or different scintillator properties. See for this section `Parametric propagation of optical photons` and `SPECIAL CASE: propagate with different scintillating material`
2. From hybridN files, the signals seen by detectors can be produced. See the section **Signals production - Readout simulation**

For reference, an example of full ./parametrization/prepareChain.py command is reported here.

```
python3 /afs/cern.ch/work/m/mpizzich/simulations/spacal/gitlab/parametrization/prepareChain.py --baseFolderOut /eos/experiment/spacal/Simulations/SpacalAngleScan2/rescale_real/Angle_3_3 --config /afs/cern.ch/work/m/mpizzich/simulations/jobs/SpacalAngleScan2/realisticSpacal_rescale2.43.cfg --baseGPS /afs/cern.ch/work/m/mpizzich/simulations/jobs/SpacalAngleScan2/gps_electron_3_3.mac --build /afs/cern.ch/work/m/mpizzich/simulations/spacal/gitlab/build --listEnergy 1 2 3 4 5 10 20 35 50 100 --baseFolderJobs Angle_3_3 --listCalibrations /eos/experiment/spacal/Simulations/SpacalAngleScan3/realistic_rescale/OpticalCalibrationOut/calibration.data --listEvents 25 25 25 25 25 10 5 2 2 2 --events 1000 --listQueue longlunch longlunch workday workday workday workday workday workday workday workday
```

At this point, a bash script named

`runAll_[baseFolderJobs].sh`

has been created in the folder. Of course [baseFolderJobs] is the value of the parameter --baseFolderJobs. Run it to submit all jobs to cluster, executing

```
./runAll_[baseFolderJobs].sh       
```

The output files (after some time...) will be in the folder selected with the --baseFolderOut parameter.


## **Hybrid-MC output**

As explained, the basic output of the simulation consists in the

```
outN.root           - output of energy deposition and Cerenkov propagation, for job N
hybridN.root        - output of parametrized scintillation propagation, energy deposition and Cerenkov propagation, for job N
```

files. On top of these, the user can produce the `pulses` seen by the detectors coupled to the simulated modules. Two approaches have been developed.




### Signals Production - Readout simulation

*L. Martinazzoli (loris.martinazzoli@cern.ch)*

The code  ```simReadout.cpp``` performs a GroupBy on the photons of the hybrid file grouping them by event number and simulates a photodetector producing a sampled pulse.

Use the following command to run it:

```
./simReadout -c SignalConfigFile.cfg -i InputHybridFile.root -o OutputFileName
```

The number of modules, cells and their geometry is automatically retrieved from the maps in the InputHybridFile.
The configuration file allows to set all the input and output name branches, as well as the digitization settings (gate, binning, delay, electronic noise, etc...).   

*Program Workflow:*
The program takes in input the hybrid file, groups the photons by event, module and cell and produce pulses.
1. Some photons are discarded to simulate losses in the coupling. This can be tuned in the configuration file through the parameter ```UnifLossFac```, defined as Number of photons exiting the fibres / Number of photons reaching the photodetector. For instance if ```UnifLossFac = 3```, 1 photon every 3 reaches the photodetector. 
2. The photodetector response is simulated. QE is taken into account and, for each photon detected, the timestamp of arrival is smeared according to the detector single photon time response (SPTR), an overall delay is added simulating cables and the pulse is produced. All the relevant parameters can be set in the configuration file.

The output of this program will contain 2 + 2 x #modules in the simulation branches. That is, for a single module simulation, it will contain 2 + 2 = 4 branches.
They are:
- ```Total_Light``` (int): for each event, total number or photons detected (or photoelectrons produced) in the whole calorimeter. Useful to understand the total amount of energy deposited.
- ```modulesHit``` (std::vector<int>): for each event, vector containing the ID of the modules detecting at least one photon during that event.
- ```modN_ph``` (std::vector<int>) where N is the ID of the module: for each event, vector with as many entries as the number of cells in the module N. Each entry is the number of photons detected in each cell.
- ```modN_pulse``` (std::vector<std::vector<float>>) where N is the ID of the module: for each event, vector with as many entries as the number of cells in the module N. Each entry is a std::vector<float> containing the simulated and digitised pulse.

To visualise the output file, it is possible to use the ROOT macro ```scripts/pulseVisualization.C```. It will show the pulses of an event of a module. These, along the input file, can be selected in the macro through 3 variables.

*Getting the timestamps: ApplyCFD.cpp:*
The program ```ApplyCFD.cpp``` allows to apply a constant fraction discriminator (CFD) technique to get a timestamp for each pulse.
It takes in input the output file of ```simReadout.cpp``` and creates an output .root file which contains all the information of the input file plus a branch of timestamps per branch of pulses. Additionally, it draws the pulse of an event for each branch of pulses.
It can be run with:

```
./ApplyCFD -c SignalConfigFile.cfg -i InputFileName.root -o OutputFileName -t 0.5
```
The flags stand for:
- *c* Configuration file
- *i* Input data file
- *o* Output data file
- *t* Peak fraction whereat the threshold is set. If the value is valid, it overrides that set in the configuration file. [optional]

Some configuration parameters required by the ```ApplyCFD``` are in common with ```simReadout```, therefore they can share the same configuration file. In fact, the digitization parameters are in common for the two programs.
Key parameters are:
- ```baseLineSamples``` : Number of samples used to calculate mean and sigma of the baseline
- ```peakFraction```    : Fraction of the pulse peak at which the timestamp will be calculated.
- ```convert_to_ns```   : All the timestamps and the plots will be converted in time units. The conversion factor between clock units and time is calculated from the digitizer parameters in the configuration file. OFF by default.


The output of this program will be a copy of the input tree where the ```modN_pulse``` branches are delete and replaced by:
- ```modN_t``` (std::vector<float>) where N is the ID of the module: for each event, vector with as many entries as the number of cells in the module N. Each entry is the timestamp obtained with CFD on the pulse. If the algorithm failed, for instance due to a pulse indistinguishable from the electronic noise, a negative timestamp is returned. Hence, *always filter out negative values when doing timing analysis*

## ***STEP BY STEP***

In case you are interested, the chain of G4 simulation, hybrid propagation, etc can be performed step by step. This allows the user to perform some special actions. Also, if the user is not deleting the output files generated by the various steps of the chain described above (in particular the outN.root and hybridN.root files), it is possible to precess them again with different parameters (see below).

### **From config file to outN.root files**

The first step of the chain is to run the Geant4 simulation that will produce a map of energy depositions and will propagate the Cerenkov photons. This can be done, starting from proper mySpacal.cfg and gpsBase.mac files, with the following command:

/path/to/your/spacal/build/folder/FibresCalo mySpacal.cfg out gpsBase.mac

This will produce the output file

```
out.root
```


### **Optical calibration**

If you want to run an optical calibration on local machine. First of all you need to produce a opticalConfiguration.cfg file, in the same way as described in point 2 of "Running full Hybrid-MC simulation on LXPLUS". Then, you have two possible strategies. The first one, strongly suggested, is that you start by using the prepareOpticalCalibration.py script, but pointing to local PC folders with the parameters

```
--baseFolderJobs
--baseFolderOut
```

Put whatever value you want for the

```
--queue
```

parameter, it will be ignored (but do put something, or the python script will fail). Finally, split the calibration in as many jobs as you like, but most likely you will want not to exceed the number of cores of you PC. So for example use

```
--jobs 4
```

Run the script. Now go the the folder selected with --baseFolderJobs, and manually run the jobs. You can do this by running run_script.sh followed by an entire line from the file args.txt. Otherwise, you can paste all together with some smart script. We paste here a possibility in bash, just for reference:

```
#!/bin/bash

for i in 0 1 2 3
do
  time ./run_script.sh /path/to/your/spacal/build/folder/FibresCalo /path/to/baseFolderJobs/optConfig.cfg ~/path/to/baseFolderOut/out$i /path/to/baseFolderJobs/gps.mac &> log_${i}.txt &
done

jobs
wait
```

When the jobs are done, merge them together with the mergeCalibration program, following the same procedure described at the end of point 2 of "Running full Hybrid-MC simulation on LXPLUS" section. You will obatin your

```
calibration.data
```

file.
A second strategy is to write your own gps file(s), choosing the points of photon emission and the number of primaries per point.

```
Be careful, you NEED to use a regular grid of x,y,z and energy points (i.e., within one coordinate, the distance between points must always be the same). Otherwise the fast algorithm on which propagateHybrid is based will not work (but it will not crash, it will just produce wrong results.. so beware)
```

The you can run the optical calibration as any Geant4 simulation (see the beginning of this readme). This procedure is quite dangerous, so please try to use the first strategy.


## **Parametric generation and propagation of scintillation photons**

You can now perform manually the parametric generation and transportation of scintillation photons, using this command

```
/path/to/your/spacal/build/folder/propagateHybrid -i out.root -c calibration.data -o hybrid.root
```

where obviously out.root is the output of step 1, and calibration.data the output of step 2.


## **SPECIAL CASE: propagate with different scintillating material**

On the basis of out.root data, the user can re-run propagateHybrid several time, using crystal materials that are different from the one used in the generation of out.root and calibration.data.

```
ATTENTION:
This approach makes sense if you want to test for example the same type of crystal but with different timing performance. For example, GaGG for different producers. In this case, the energy deposition, and the Cerenkov transportation simulated in first step will make sense. Also, the optical calibration will be still valid. Only the generation (time and energy) of scintillation photons will be changed. In this case, the procedure makes sense. Otherwise, it doesn't.
```

In order to do this, the user needs to provide a file that we call

```
material.root
```

that contains 5 entries that are usually found in all outN.root files. More specifically

```
TH1F                  enHisto_N                      - emission spectrum of the material
TH1F                  tHisto_N                       - scintillation time profile of the material
std::vector<int>      crystalMaterialList            - crystal material number
std::vector<float>    crystalLightYield              - crystal light yield
std::vector<float>    crystalResolutionScale         - crystal resolution scale (Geant4 parameter, check Geant5 documentation for an explanation)
```

Look into one of the outN.root files to understand how these entries look like.
Then the propagation can be performed with

 ```
/path/to/your/spacal/build/folder/propagateHybrid -i out.root -c calibration.data -o hybrid.root --external material.root
 ```


## **Production of signals - simple approach**

```
ATTENTION: this method does not work anymore!
```

From hybrid.root files, realistic "signals" seen by detectors can be produced. The

```
resolution_[DETECTOR]_N.root
```

file described at the end of the chain are produced with a simple data analysis, that can be run manually using the "readHybrid" program. The program assumes that all photons propagated out of the module by propagateHybrid hit a detector (so there is no coupling efficiency). It then sorts photons, event by event, on different detectors, assuming a detector for each "cell" defined in the original simulation configuration file (mySpacal.cfg). For each detectors, it then processes photons 1 by 1, applying the quantum efficiency filter (which here is flat, no energy dependence). The surviving photons for each detector are summed and saved as the "charge" seen by the detector. The times of arrival of these surviving photons is smeared by the detector timing resolution, using a simple Gaussian function. These smeared timestamps are sorted in chronological order, and finally the timestamp of the detector is computed as the average of the first N timestamps.

All this can be performed running readHybrid with the following parameters:

```
--input       <value>   - input file name (i.e. one of the hybridN.root files)
--output      <value>   - output file name
--seed        <value>   - random seed                                     - default = 0 , i.e. seed automatically computed via a TUUID object
--qe          <value>   - detector quantum efficiency (as a fraction)     - default = 0.07
--sptr        <value>   - detector time resolution (either TTS or SPTR) in ns  - default = 0.087
--photons     <value>   - compute the timsestamp as the average of first <value> photons detected   - default = 5
--events      <value>   - primary events in input file                    - default = -1 (if you leave default, program will calculate it. slower)
--detectors   <value>   - number of detectors                             - default = -1 (if you leave default, program will calculate it. slower)
--verbose     <value>   - verbose output. better keep this 0, i.e false, for condor jobs            - default = 0
```

Be careful, the values of "events" and "detectors" HAVE to be the ones actually saved in the input file. You cannot use this keys to simulate different detector pixellizations.
So for example

```
/path/to/your/spacal/build/folder/readHybrid -i hybrid.root -o resolution_DET.root -s 1 --qe 0.07 --photons 5 --sptr 0.232 --events 15 --detectors 162
```

will do the trick for a detector called "DET", with qe = 0.07 and SPTR = 0.232 ns, computing the timestamps for each detector as the average of the time of arrival of the first 5 photons detected, assuming a configuration with 162 detectors and on an input file with 15 events. The output of this command is already described in at the end of "Hybrid-MC output" section.

### **Test beam configuration**

*L. Capriotti (lorenzo.capriotti@cern.ch)*

The test beam configuration (configurations/testBeamConfig.cfg) has been updated with the possibility to introduce several layers of different materials acting as separation volume between the two absorber sections in case a separation is requested. This is ensured by setting

```
cell_separation_type = 2
```

in the cfg file. The new keys to consider are the following:

```
LAPPD_layers_thickness = |t1|t2|t3|t4|...|tN|
LAPPD_layers_materials = |m1|m2|m3|m4|...|mN|
```

They will create N layers, each of them with their respective thickness and material, between the two absorber sections. The possible choice for materials is listed under Absorber::SetLAPPDMaterials, which creates a vector of materials to be paired with the vector of thicknesses imported from the cfg file. Then, a boolean flag ensures that, if the size of the vector of thicknesses is zero, which means

```
LAPPD_layers_thickness = NULL
```

then both this key and the LAPPD_layers_materials key are ignored, and the positioning of a separation volume is regulated by the keys separation_thickness and separation_material. Note that the vector of materials does not have to be set to NULL for this to happen. The position of the center of the separation volume is still set by the cell_separator_position key, as in the case of single-layer separator.
Note that another flag has been introduced (saveLAPPD) which creates a TTree with observables of interest for performance studies of a LAPPD as separator.
