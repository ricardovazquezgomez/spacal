# These are instructions how to generate B0 -> K+pi-pi0 or B+ -> K+pi0 events (without pileup), possibly changing the number of events
# Check gauss_upgrade_pgun.py, ReadSim.py, flux_from_LHCbXgenPi02GammaGamma.py for comments on B0 -> K+pi-pi0 or B+ -> K+pi0 decays,
# and make sure that necessary lines are uncommented

# modify gauss_upgrade_pgun.py: change LHCbApp().EvtMax and output_name, and run:
lb-run Gauss/v52r2 gaudirun.py gauss_upgrade_pgun.py

# modify ReadSim.py: modify DaVinci().Input and DaVinci().TupleFile
# to fill the timing information, we need to modify the DecayTreeTupleMC as done in DaVinciDev_v45r4
/eos/experiment/spacal/users/anlp/pgun/DaVinciDev_v45r4/run gaudirun.py ReadSim.py

# convert the generator output to our flux file format 
# (original script: /eos/experiment/spacal/users/zhangy/flux_from_LHCbXgenPi02GammaGamma.py is modified to work with both B+ and B0)
# modify flux_from_LHCbXgenPi02GammaGamma.py: change fi, fs
python2 flux_from_LHCbXgenPi02GammaGamma.py

# finally, one might want to filter events to remove photons outside the interesting regions
