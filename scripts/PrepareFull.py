###### USAGE:
### 		python3 prepareChain.py [args]
### see README.md for details on [args]
# Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

import math
import os
import stat
import sys
import argparse
import subprocess
from subprocess import Popen, PIPE, STDOUT
import shutil
from shutil import copyfile
def main(args):
  # bla

  #parsing args
  parser = argparse.ArgumentParser(description='Python script to prepare optical calibration run')
  parser.add_argument('--baseFolderJobs'   , default=''    , help='Output folder',required=True)
  parser.add_argument('--baseFolderOut'   , default=''    , help='Output folder',required=True)
  parser.add_argument('--config'                         , help='Base config file',required=True)
  parser.add_argument('--baseGPS'  ,   default=''  ,      help='Number of optical photons generated for each point',required=True)
  parser.add_argument('--configSF'                         , help='Base config file for signal formation',required=True)
  parser.add_argument('--build'                          , help='path to spacal build folder',required=True)
  # parser.add_argument('--calibration'      , default=""        , help='optical calibration file',required=True)
  parser.add_argument('--events'      , default=""        , help='optical calibration file',required=True)
  parser.add_argument('--delete_out', action='store_true', help='do not save out* files at the end of job')
  parser.add_argument('--delete_hybrid', action='store_true', help='do not save hybrid* files at the end of job')
  parser.add_argument('--delete_signal', action='store_true', help='do not save simReadout output files at the end of job')
  parser.add_argument('--delete_CFD', action='store_true', help='do not save ApplyCFD output files at the end of job')
#   parser.add_argument('--onlyEnDepo', action='store_true', help='Stop at the out* files stage')
  # parser.add_argument('--simple', action='store_true', help='produce the simple signal outputs')
  parser.add_argument("--listEnergy",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=int,
                       default=[1],  # default if nothing is provided
  )
  parser.add_argument("--listEvents",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=int,
                       default=[100],  # default if nothing is provided
  )
  parser.add_argument("--listQueue",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=str,
                       default=["workday"],  # default if nothing is provided
  )
  parser.add_argument("--listCalibrations",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=str,
                       default=[""],  # default if nothing is provided
  )
  parser.add_argument("--listTypes",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=int,
                       default=[0],  # default if nothing is provided
  )


  args = parser.parse_args()

  # assign variables
  baseFolderOut   = args.baseFolderOut
  baseFolderJobs  = args.baseFolderJobs
  configFile      = args.config
  configFileSF    = args.configSF
  baseGPS         = args.baseGPS
  build           = args.build
  events          = int(args.events)
  listEnergy      = args.listEnergy
  listEvents      = args.listEvents
  listQueue       = args.listQueue
  listCalibrations= args.listCalibrations
  listTypes       = args.listTypes
  deleteOut       = args.delete_out
  deleteHybrid    = args.delete_hybrid
  deleteSignal    = args.delete_signal
  deleteCFD       = args.delete_CFD


#   onlyEnDepo      = args.onlyEnDepo

  # simpleAnaysis   = args.simple
  #check lists size
  if len(listEnergy) != len(listEvents):
      print("ERROR! Lists of energies and events must have the same length!!! Aborting...")
      sys.exit(1)
  if len(listEnergy) != len(listQueue):
      print("ERROR! Lists of energies and queues must have the same length!!! Aborting...")
      sys.exit(1)
  if len(listCalibrations) != len(listTypes):
      print("ERROR! Lists of calibration files and types must have the same length!!! Aborting...")
      sys.exit(1)

  ## calc args list lengths
  listArgLines = []
  for i in listEvents:
      if (events % i) != 0:
           print("ERROR! Number of events is not divisible by events per job!!! Aborting...")
           sys.exit(1)
      ArgLines = int(events / i)
      listArgLines.append(ArgLines)
  listGevFolders = []
  for i in listEnergy:
      listGevFolders.append(str(i)+"GeV")

  # if onlyEnDepo and deleteFiles:
      # print ("ERROR! Both deleteFiles and onlyEnDepo!!! Just a waste of CPU time! Aborting...")
      # sys.exit(1)

  ### abs paths

  script_identifier = baseFolderJobs
  baseFolderJobs   = os.path.abspath(baseFolderJobs)
  baseFolderOut    = os.path.abspath(baseFolderOut)
  build            = os.path.abspath(build)
  configFile       = os.path.abspath(configFile)
  configFileSF     = os.path.abspath(configFileSF)
  baseGPS          = os.path.abspath(baseGPS)
  listCalibrations = list(map(os.path.abspath, listCalibrations))

  print("")
  print("########################################")
  print("# FEEDBACK                             #")
  print("########################################")
  print("")
  print("baseFolderOut     = %s" % baseFolderOut)
  print("baseFolderJobs    = %s" % baseFolderJobs)
  print("configFile        = %s" % configFile)
  print("configFileSignal  = %s" % configFileSF)
  print("baseGPS           = %s" % baseGPS)
  print("build             = %s" % build)
  print("listEnergy        = %s" % listEnergy)
  print("listEvents        = %s" % listEvents)
  print("listCalibrations  = %s" % listCalibrations)
  print("listTypes         = %s" % listTypes)
  print("events            = %s" % events)
  print("listArgLines      = %s" % listArgLines)
  print("deleteOut         = %s" % deleteOut   )
  print("deleteHybrid      = %s" % deleteHybrid)
  print("deleteSignal      = %s" % deleteSignal)
  print("deleteCFD         = %s" % deleteCFD   )

  ## make job folders in a subfolder here
  if not os.path.exists(baseFolderJobs):
      print("Creating job folder at %s" %(baseFolderJobs))
      os.mkdir(baseFolderJobs)
  else:
      print("Job folder exists at %s" %(baseFolderJobs))

  baseFolderJobs = os.path.abspath(baseFolderJobs)
  #create output folder
  if not os.path.exists(baseFolderOut):
      print("Creating output folder at %s" %(baseFolderOut))
      os.mkdir(baseFolderOut)
  else:
      print("Output folder exists at %s" %(baseFolderOut))

  ## create subfolders for energies
  index = 0
  for e in listEnergy:

      newFolder = baseFolderJobs + "/" + str(listEnergy[index]) + "GeV"
      if not os.path.exists(newFolder):
          print("Creating output folder at %s" %(newFolder))
          os.mkdir(newFolder)
      else:
          print("Output folder exists at %s" %(newFolder))
      outDir = newFolder + "/out"
      errDir = newFolder + "/err"
      logDir = newFolder + "/log"
      if not os.path.exists(outDir):
          os.mkdir(outDir)
      else:
          print("Out folder exists at %s" %(outDir))
      if not os.path.exists(errDir):
          os.mkdir(errDir)
      else:
          print("Err folder exists at %s" %(errDir))
      if not os.path.exists(logDir):
          os.mkdir(logDir)
      else:
          print("Log folder exists at %s" %(errDir))
      newOutFolder = baseFolderOut + "/" + str(listEnergy[index]) + "GeV"
      if not os.path.exists(newOutFolder):
          os.mkdir(newOutFolder)
      else:
          print("Output subfolder exists at %s" %(newOutFolder))
      #create condor files...

      #########################
      # RUN SCRIPT            #
      #########################
      run_script = newFolder + "/run_script.sh"
      with open(run_script, "w") as output:

          output.write("#!/bin/bash\n\n")
          output.write("#parse input, they are written in args.txt\n\n")
          output.write("COMMAND=$1\n")
          output.write("TEMPLATE=$2\n")
          output.write("OUTPUT=$3\n")
          output.write("GPS=$4\n")
          output.write("PROPAGATE=$5\n")
          output.write("CALIBRATION=$6\n")
          output.write("HYBRID=$7\n")
          output.write("SEED=$8\n")
          output.write("SIGNAL=$9\n")
          output.write("TRIG=${10}\n")
          output.write("CONFIG_SF=${11}\n")
          output.write("OUTPUT_SF1=${12}\n")
          output.write("OUTPUT_CFD=${13}\n")

          output.write("set --\n\n")
          output.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh\n")
          # old sources, left here for reference
          # output.write("source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
          # output.write("source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh\n")
          # output.write("source /cvmfs/geant4.cern.ch/geant4/10.3/x86_64-slc6-gcc49-opt/bin/geant4.sh\n\n")
          output.write("# do the MC\n")
          output.write("$COMMAND $TEMPLATE output_temp $GPS $SEED\n\n")
          output.write("# do the hybrid propagation\n")
          output.write("$PROPAGATE -i output_temp.root -c \"$CALIBRATION\" -o hybrid_temp.root --photonSeed $SEED --type ")
          for tIndex in range(len(listTypes)):
              output.write("%d"  %listTypes[tIndex])
              if tIndex != (len(listTypes)-1):
                  output.write(",")
              else:
                  output.write(" ")
          output.write("\n\n")
          output.write("### Do the GroupBy\n")
          output.write("$SIGNAL -c $CONFIG_SF -i hybrid_temp.root -o output_sf1\n\n")
          output.write("### Do the pulse triggering\n")
          output.write("$TRIG   -c $CONFIG_SF -i output_sf1.root -o output_sf2\n" )
          output.write("\n\n")
          # copy output files
          if deleteOut    != True:
            output.write("echo \"Copying output file...\"\n\n")
            output.write("xrdcp --nopbar output_temp.root $OUTPUT.root\n\n")
          if deleteHybrid != True:
            output.write("echo \"Copying hybrid file...\"\n\n")
            output.write("xrdcp --nopbar hybrid_temp.root $HYBRID\n\n")
          if deleteSignal != True:
            output.write("echo \"Copying pulse file...\"\n\n")
            output.write("xrdcp --nopbar output_sf1.root $OUTPUT_SF1\n\n")
          if deleteCFD    != True:
            output.write("echo \"Copying CFD file...\"\n\n")
            output.write("xrdcp --nopbar output_sf2.root $OUTPUT_CFD\n\n")
          output.close()
      #and make executable
      st = os.stat(run_script)
      os.chmod(run_script, st.st_mode | stat.S_IEXEC)

      #########################
      # JOBS SCRIPT           #
      #########################
      jobsub = newFolder + "/jobs.sub"
      with open(jobsub, "w") as output:
          output.write("executable              = run_script.sh \n")
          output.write("output                  = out/out.$(ClusterId).$(ProcId).out \n")
          output.write("error                   = err/err.$(ClusterId).$(ProcId).err \n")
          output.write("log                     = log/log.$(ClusterId).log  \n")
          output.write("transfer_output_files   = \"\" \n")
          output.write("+JobFlavour             = \"%s\" \n" %(listQueue[index]))
          output.write("queue arguments from args.txt \n")
          output.close()

      #########################
      # GPS                  #
      #########################
      gpsFile = newFolder + "/gps.mac"
      gpsFile = os.path.abspath(gpsFile)

      to_be_removed = ['/gps/energy' ,'/run/beamOn']
      with open(baseGPS, "r") as input:
        with open(gpsFile, "w") as output:
            for line in input:
                if not any(bad_word in line for bad_word in to_be_removed):
                    output.write(line)
            output.write('\n')
            output.write("/gps/energy %d GeV\n" % listEnergy[index] )
            output.write("/run/beamOn %d\n" % listEvents[index] )
            output.close()
      #########################
      # ARGS                  #
      #########################
      argsFile = newFolder + "/args.txt"
      with open(argsFile, "w") as output:
          for i in range(listArgLines[index]):
              output.write("%s/FibresCalo " % build )                # output.write("COMMAND=$1\n")
              output.write("%s " % configFile )                      # output.write("TEMPLATE=$2\n")
              output.write("%s/out%d "%(newOutFolder,i))             # output.write("OUTPUT=$3\n")
              output.write("%s "% gpsFile)                           # output.write("GPS=$4\n")
              output.write("%s/propagateHybrid " % build)            # output.write("PROPAGATE=$5\n")
              for cIndex in range(len(listCalibrations)):
                  output.write("%s" % listCalibrations[cIndex])
                  if cIndex != (len(listCalibrations)-1):
                    output.write(",")
                  else:
                    output.write(" ")
              output.write("%s/hybrid%d.root " %(newOutFolder,i) )   # output.write("HYBRID=$7\n")
              output.write("%d " %(i+1) )   # SEED, cannot be 0!
              output.write("%s/simReadout " % build )                # output.write("SIGNAL=$9\n")
              output.write("%s/ApplyCFD   " % build )                # output.write("TRIG=$10\n")
              output.write("%s " % configFileSF)                     # output.write("CONFIG_SF=$11\n")                                                    
              output.write("%s/signal%d.root " %(newOutFolder,i) )   # output.write("OUTPUT_SF1=$12\n")
              output.write("%s/CFD%d.root "    %(newOutFolder,i) )   # output.write("OUTPUT_CFD=$13\n")
              output.write("\n")

          output.close()
      index = index +1
  #
  # now do the submit script
  subScript = "runAll_" + script_identifier + ".sh"
  with open(subScript, "w") as output:
      output.write("#!/bin/bash\n\n")
      output.write("for i in ")
      for i in listEnergy:
          output.write("%dGeV " % i)
      output.write("\n")
      output.write("do\n")
      output.write("  cd %s/$i\n" % baseFolderJobs)
      output.write("  condor_submit jobs.sub\n")
      output.write("  cd - \n")
      output.write("done\n")
      output.close()
  #and make executable
  st = os.stat(subScript)
  os.chmod(subScript, st.st_mode | stat.S_IEXEC)


  ### and it's all...
  print("Done. Please run %s to submit all jobs. Bye!" % subScript)
  sys.exit()




if __name__ == "__main__":
   main(sys.argv[1:])
