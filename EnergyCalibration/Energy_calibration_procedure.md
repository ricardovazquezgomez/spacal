# Energy Calibration for full ECAL

*M. Pizzichemi (marco.pizzichemi@cern.ch)*

## Introduction

Calculation of factors needed to pass from Photoelectrons seen by PMTs to real energy (MeV), with the HybridMC toolkit simulations. The example reported here is based on the energy calibration of 2 regions in a full ECAL simulation, but it can be easily adapted to derive a generic calibration procedure for any SPACAL module.

### Reference presentations for HybridMC

- [HybridMC toolkit presentation](https://cernbox.cern.ch/index.php/s/JG5WfSMzzZPK5CX) (10 March 2020)
- [HybridMC toolkit update](https://cernbox.cern.ch/index.php/s/Xhm1YgedlHqAzKk) (14 May 2020)
- [HybridMC radiation damage study](https://cernbox.cern.ch/index.php/s/oOXawQFbfU0sjlL) (17 Sept 2020)
- [CPU time slide](https://cernbox.cern.ch/index.php/s/5GncgYMEbRBokHU)

Documentation (in [GitLab](https://gitlab.cern.ch/spacal-rd/spacal-simulation)) is referred to in the presentations.


### Working directory on EOS

```
/eos/experiment/spacal/Simulations/EnergyCalibrationECAL
```

### Configuration files

The **main configuration files** are derived from this one in EOS:
```
/eos/experiment/spacal/Simulations/FluxSimulations/FluxDiGammaFromPi0_B2KPiPi0/spacal_full_ecal.cfg
```

In particular, the module steering files pointed by that are copied


```
cp /eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/Configurations/spacal_Pb_Polystirene_pitch1.67_12x12cm.cfg /eos/experiment/spacal/Simulations/EnergyCalibrationECAL/Configuration/
cp /eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/Configurations/spacal_W_gagg_pitch1.67_12x12cm.cfg /eos/experiment/spacal/Simulations/EnergyCalibrationECAL/Configuration/
```


The **optical calibration files** used for these ECAL simulations are here:

```
/eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/OpticalCalibration/WGagg/calibration_WGagg.data  
/eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/OpticalCalibration/PbPolystirene/calibration_PbPolystyrene.data
```

For this study, the modules will be analyzed individually, using therefore the individual configuration files. The optical calibrations used will be the same. The GPS files are produced ad hoc (using single electrons and photons, no flux), and stored in:
```
/eos/experiment/spacal/Simulations/EnergyCalibrationECAL/GPS/
```

folder. Finally, the signal formation file used is:

```
/eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/Configurations/SignalConfigFile_HPKR7899-20.cfg
```



## Strategy

A full hybrid simulation of both module types is performed. The information on real energy deposited is saved, as well as the photoelectrons generated on the detectors at the end of the hybrid chain. These two quantities are then compared to derive the calibration factors needed to reconstruct the energy deposited from photoelectrons.

## Simulation dataset creation

Individual modules are simulated. Electrons and photons of different energies (1-2-5-10-20-35-50-100 GeV) are shot in the center of the modules (squared beam, 16.7x16.7 mm$^2$ surface), the full hybrid simulation chain is performed.


Simulations are prepared in some AFS folder (for self reference, in /afs/cern.ch/work/m/mpizzich/simulations/jobs/EnergyCalibrationECAL), with commands like this one:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
python3 /afs/cern.ch/work/m/mpizzich/simulations/spacal/gitlab/parametrization/prepareChain.py --baseFolderOut /eos/experiment/spacal/Simulations/EnergyCalibrationECAL/WGAGG/Electrons --config /eos/experiment/spacal/Simulations/EnergyCalibrationECAL/Configuration/spacal_W_gagg_pitch1.67_12x12cm.cfg --baseGPS /eos/experiment/spacal/Simulations/EnergyCalibrationECAL/GPS/gps_electron_wgagg.mac --build /afs/cern.ch/work/m/mpizzich/simulations/spacal/gitlab/build --listEnergy 1 2 5 10 20 35 50 100 --baseFolderJobs jobs --listCalibrations /eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/OpticalCalibration/WGagg/calibration_WGagg.data --listEvents 200 100 40 20 10 5 4 2 --events 1000 --listQueue workday workday workday workday workday workday workday workday --pulse --pulseConfig /eos/experiment/spacal/Simulations/Note_WGagg_PbPolystirene/Configurations/SignalConfigFile_HPKR7899-20.cfg
```

For simplicity these commands are usually saved in a bash script. The simulations are then submitted to the cluster with

```
./runAll_jobs.sh
```
## Module dataset production

The output of the simulations consists in the outN files and the OutTrigdN files. The former contain the information on actual energy deposited in the module, the latter the photo-electrons generated on the detectors coupled to the module.

The **N** number is conserved, meaning for example that the file out123.root corresponds to OutTrigd_123.root. Within these files, the event order is conserved as well. Thanks to this, it is easy to analyze pair of files and extract the information need. This is done with the program calculateChiFactors (in the spacal-simulation gitlab), with following parameters:

```
calculateChiFactors
		[-e | --efile]   <input energy deposition file name>   
		[-p | --pfile]   <input pulse file name>   
		[-o | --output]  <output file name>    
		[--verbose]      <verbosity level - default = 0>
```

The output is a TTree, with following entries:

```
primaryEnergy_mev       // energy of primary particle
primaryPDGID            // PDG ID of primary particle
endep_total_mev         // total energy deposited in the module (including absorber)
endep_front_mev         // energy deposited in front section of the module (including absorber)
endep_back_mev          // energy deposited in back section of the module (including absorber)
endep_cry_total_mev     // energy deposited in crystals of whole module
endep_cry_front_mev     // energy deposited in crystals of front section of the module
endep_cry_back_mev      // energy deposited in crystals of back section of the module
photoElectrons_front    // photo-electrons generated on detectors coupled to front section of the module
photoElectrons_back     // photo-electrons generated on detectors coupled to back section of the module
```

This has been run on the entire output of the simulations, producing a **chiN.root** file for each {outN,OutTrigdN} pair. To give an example of how to run the calculateChiFactors program on all the dataset just produced, take a look into the `ComputeChi` subfolder. Here you will find a simple structure to submit jobs to HTCondor, assuming that the output of the previous step is in a folder

```
/eos/experiment/spacal/Simulations/EnergyCalibrationECAL/
```

with the structure

```
.
├── Configuration
├── GPS
├── PbPoly
│   ├── Electrons
│   │   ├── 100.0GeV
│   │   ├── 10.0GeV
│   │   ├── 1.0GeV
│   │   ├── 20.0GeV
│   │   ├── 2.0GeV
│   │   ├── 35.0GeV
│   │   ├── 50.0GeV
│   │   └── 5.0GeV
│   └── Gammas
│       ├── 100.0GeV
│       ├── 10.0GeV
│       ├── 1.0GeV
│       ├── 20.0GeV
│       ├── 2.0GeV
│       ├── 35.0GeV
│       ├── 50.0GeV
│       └── 5.0GeV
├── WGAGG
    ├── Electrons
    │   ├── 100.0GeV
    │   ├── 10.0GeV
    │   ├── 1.0GeV
    │   ├── 20.0GeV
    │   ├── 2.0GeV
    │   ├── 35.0GeV
    │   ├── 50.0GeV
    │   └── 5.0GeV
    └── Gammas
        ├── 100.0GeV
        ├── 10.0GeV
        ├── 1.0GeV
        ├── 20.0GeV
        ├── 2.0GeV
        ├── 35.0GeV
        ├── 50.0GeV
        └── 5.0GeV
```

This produces all the **chiN.root** files. Now you need to merge them into a dataset that summarizes all the information for each module. In order to do this, you can use hadd and adapt the `merge.sh` script that you find in this folder. For example:

```
cd /eos/experiment/spacal/Simulations/EnergyCalibrationECAL
./merge.sh
```

This will produce the final ROOT files, in this case namely:

```
/eos/experiment/spacal/Simulations/EnergyCalibrationECAL/WGAGG/calibrationWGAGG.root
/eos/experiment/spacal/Simulations/EnergyCalibrationECAL/PbPoly/calibrationPbPoly.root
```

## Calibration factors computation

Both summary files can be analized with the same procedure. Producing the plots below amounts to just some Draw and Fit commands run on the branches explained above, so those commands are not reported here.

The results are shown for WGAGG ones.

- photo-electrons in **front section** plotted against real energy deposited in the entire front part of the module (so in crystals AND absorber). Points from simulations with primary <span style="color:red">gammas are in red</span>, while points from primary <span style="color:blue">electrons in blue</span>.

![](https://codimd.web.cern.ch/uploads/upload_a8e4437df34137bd4451437b64dba1f2.png)

- photo-electrons in **back section** plotted against real energy deposited in the entire back part of the module (so in crystals AND absorber). Points from simulations with primary <span style="color:red">gammas are in red</span>, while points from primary <span style="color:blue">electrons in blue</span>.

![](https://codimd.web.cern.ch/uploads/upload_c967840081cff3aeaec8cddbe9e705be.png)


- since the plots follow the same behavior, the points are analyzed together. A **simple line regression** is made to derive a relation from photo-electrons measured to energy deposited in a section. More sophisticated analysis are possible of course

![](https://codimd.web.cern.ch/uploads/upload_c6d13350694ec15eae4c8e82b0d282a9.png)

![](https://codimd.web.cern.ch/uploads/upload_ce48ebee8b0dec3b4251b1c53e568a4f.png)

:::info
Notice that the lines do not seem compatible with crossing in the origin. This point has to be discussed (I would expect to reconstruct 0 MeV for 0 photoelectrons). Possibly the relation is not completely linear (not totally surprising, since the shower center of gravity moves as a function of primary energy, resulting in an average response from a section that changes because of the effective attenuation curve).
:::
- Assuming the simple linear relation is ok, the reconstructed energy deposited in a module can be reconstructed simply by

$$
E_{reco} = (P_f \cdot m_f + q_f) + (P_b \cdot m_b + q_b)  
$$

where $P_x$ are the photo-electrons measured, $m_x$ and $q_x$ are coefficients of the linear regressions above, and of course $x$ is $f$ for the front, $b$ for the back.

- We can test this by plotting the reconstructed energy versus the MC truth


![](https://codimd.web.cern.ch/uploads/upload_3e05d3dc6ab650d3f435b58d60c17c10.png)


A very rough linear fit shows angular coefficient compatible with 1, within 0.2%.

- Same test is reported for Pp+Poly modules

![](https://codimd.web.cern.ch/uploads/upload_b53e49d3e00cd645d20979963b972316.png)




## Summary

In summary, when in the final pulse files a certain number of photoelectrons are found, the total energy in a cell/module can be reconstructed with

$$
E_{reco} = (P_f \cdot m_f + q_f) + (P_b \cdot m_b + q_b)  
$$

using the coefficients are reported here:

| Module  | Section | m         | q        |
| ------- | ------- | --------- | -------- |
| W+GAGG  | Front   | 0.02566   | -69.35   |
| W+GAGG  | Back    | 0.0341725 | 47.6619  |
| Pb+Poly | Front   | 0.144837  | -51.3804 |
| Pb+Poly | Back    | 0.229883  | 199.982  |
| W+Poly  | Front   | 0.19795   | -3.092   |
| W+Poly  | Back    | 0.26426   | 71.5739  |

Obviously one can use only half of the equation above, if interested in reconstructing the energy in just front or just back section.
