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
import random
from shutil import copyfile
def main(args):
  # bla

  #parsing args
  parser = argparse.ArgumentParser(description='Python script to prepare optical calibration run')
  parser.add_argument('--baseFolderJobs'   , default=''    , help='Output folder',required=True)
  parser.add_argument('--baseFolderOut'   , default=''    , help='Output folder',required=True)
  parser.add_argument('--config'                         , help='Base config file',required=True)
  parser.add_argument('--baseGPS'  ,   default=''  ,      help='Number of optical photons generated for each point',required=True)
  parser.add_argument('--build'                          , help='path to spacal build folder',required=True)
  # parser.add_argument('--calibration'      , default=""        , help='optical calibration file',required=True)
  parser.add_argument('--events'      , default=""        , help='optical calibration file',required=True)
  # parser.add_argument('--delete', action='store_true', help='delete out* and hybrid* files at the end of job')
  parser.add_argument('--onlyEnDepo', action='store_true', help='Stop at the out* files stage')
  parser.add_argument('--useFlux', default='', help='Use flux file as primary')
  # parser.add_argument('--simple', action='store_true', help='produce the simple signal outputs')
  parser.add_argument("--listEnergy",    # name on the CLI - drop the `--` for positional/required parameters
                       nargs="*",  # 0 or more values expected => creates a list
                       type=float,
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
  parser.add_argument('--local', action='store_true', help='Jobs will be run locally')
  parser.add_argument('--pulse', action='store_true', help='Perform pulse formation')
  parser.add_argument('--pulseConfig', help='Pulse formation config file', default='')
  ## decide what to copy.
  # if the simulation is just energy deposition(+cerenkov) then the user will never want to get rid of those files
  # being this the baseline, energy depo are kept by default, and you need to specify that you don't want them  
  parser.add_argument('--discardEnergyDepositions', action='store_true', help='Copy EnergyDepositions files from LXPLUS node to EOS storage',default=False)
  parser.add_argument('--keepHybrid', action='store_true', help='Copy Hybrid files from LXPLUS node to EOS storage', default=False)
  parser.add_argument('--keepGroupBy', action='store_true', help='Copy GroupBy files from LXPLUS node to EOS storage',default=False)
  parser.add_argument('--requestDisk' )

  args = parser.parse_args()

  # assign variables
  baseFolderOut   = args.baseFolderOut
  baseFolderJobs  = args.baseFolderJobs
  configFile      = args.config
  baseGPS         = args.baseGPS
  build           = args.build
  events          = int(args.events)
  listEnergy      = args.listEnergy
  listEvents      = args.listEvents
  listQueue       = args.listQueue
  listCalibrations= args.listCalibrations
  listTypes       = args.listTypes
  local           = args.local
  pulse           = args.pulse
  pulseConfigFile = args.pulseConfig
  discardEnergyDepositions = args.discardEnergyDepositions
  keepHybrid            = args.keepHybrid           
  keepGroupBy           = args.keepGroupBy     
  requestDisk     = args.requestDisk     

  # deleteFiles     = args.delete
  onlyEnDepo      = args.onlyEnDepo
   
  # check for useless choices
  if(onlyEnDepo == True):
    if(discardEnergyDepositions == True):
      print("ERROR! You chose to perform only energy deposition (and eventually Cerenkov) but you set energy deposition output to be discared. You will end up with no output! Aborting...")
      sys.exit(1)
  
  # real variable to save of not energy dep
  keepEnergyDepositions = True
  if(discardEnergyDepositions == True):
    keepEnergyDepositions = False


  fluxBaseFile    = ''

  if(args.useFlux == ''):
      useFlux = False
  else:
      useFlux = True
      fluxBaseFile = args.useFlux
  # useFlux         = args.useFlux
  if(args.pulseConfig == ''):
    pulseConfigFile = 'dummyConfigFile.cfg'
    
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
  baseGPS          = os.path.abspath(baseGPS)
  
  request_disk = False
  if requestDisk is not None:
      request_disk = True
      request_disk_space = requestDisk



  if(useFlux == True):
      fluxBaseFile     = os.path.abspath(fluxBaseFile)
  listCalibrations = list(map(os.path.abspath, listCalibrations))

  # create folder for flux files if needed
  baseFolderFlux = baseFolderJobs # default
  if(useFlux == True):
      baseFolderFlux = baseFolderJobs + "/fluxFiles"
      baseFolderFlux = os.path.abspath(baseFolderFlux)

  print("")
  print("########################################")
  print("# FEEDBACK                             #")
  print("########################################")
  print("")
  print("baseFolderOut     = %s" % baseFolderOut)
  print("baseFolderJobs    = %s" % baseFolderJobs)
  if(useFlux):
      print("baseFolderFlux    = %s" % baseFolderFlux)
      print("fluxBaseFile      = %s" % fluxBaseFile)
  print("configFile        = %s" % configFile)
  print("baseGPS           = %s" % baseGPS)
  print("build             = %s" % build)
  print("listEnergy        = %s" % listEnergy)
  print("listEvents        = %s" % listEvents)
  print("listCalibrations  = %s" % listCalibrations)
  print("listTypes         = %s" % listTypes)
  print("events            = %s" % events)
  print("listArgLines      = %s" % listArgLines)
  # print("deleteFiles       = %s" % deleteFiles)
  print("onlyEnDepo        = %s" % onlyEnDepo)
  print("useFlux           = %s" % useFlux)
  print("local             = %s" % local)
  print("pulse             = %s" % pulse)
  print("pulseConfigFile   = %s" % pulseConfigFile)
  print("keepEnergyDepositions = %s" % keepEnergyDepositions)
  print("keepHybrid            = %s" % keepHybrid)
  print("keepGroupBy           = %s" % keepGroupBy)

  ## make job folders in a subfolder here
  if not os.path.exists(baseFolderJobs):
      print("Creating job folder at %s" %(baseFolderJobs))
      os.mkdir(baseFolderJobs)
  else:
      print("Job folder exists at %s" %(baseFolderJobs))

  #create output folder
  if not os.path.exists(baseFolderOut):
      print("Creating output folder at %s" %(baseFolderOut))
      os.mkdir(baseFolderOut)
  else:
      print("Output folder exists at %s" %(baseFolderOut))


  if(useFlux):
      ## create flux folder if needed
      if not os.path.exists(baseFolderFlux):
          print("Creating flux files folder at %s" %(baseFolderFlux))
          os.mkdir(baseFolderFlux)
      else:
          print("Flux files folder exists at %s" %(baseFolderFlux))

      print("### SPLITTING INPUT FLUX FILE" )
      ## split input flux file
      fluxExec = build + "/splitFluxFile"

      cmd = [fluxExec,fluxBaseFile,str(events),baseFolderFlux]
      cmdString = ' '.join(cmd)
      print (cmdString)
      subprocess.Popen(cmd).wait()
      # count flux files extracted
      countdir = 0
      for path in os.listdir(baseFolderFlux):
          if os.path.isfile(os.path.join(baseFolderFlux, path)):
              countdir += 1
      # force --events and --listEvents to this value
      # and force listArgLines to be recalculated
      listArgLines = []
      listEvents   = [1]
      events       = countdir
      print("#### FORCING PARAMETERS BECAUSE OF USING FLUX ####")
      print("listEvents        = %s" % listEvents)
      print("events            = %s" % events)
      
      for i in listEvents:
          if (events % i) != 0:
               print("ERROR! Number of events is not divisible by events per job!!! Aborting...")
               sys.exit(1)
          ArgLines = int(events / i)
          listArgLines.append(ArgLines)
      print("listArgLines      = %s" % listArgLines)

      

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
      
      #temp file names 
      output_temp="output_temp"
      hybrid_temp="hybrid_temp"
      groupby_temp="output_f1"
      tr_temp="output_f2"
      if local == True:
          output_temp="$OUTPUT"
          hybrid_temp="$HYBRID"
          groupby_temp="$OUTGROUPBY"
          tr_temp="$OUTTR"
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
          # if useFlux == True:
          output.write("FLUX=$9\n")
          # if pulse == True:
          output.write("PULSEGROUPBY=${10}\n")
          output.write("PULSETR=${11}\n")
          output.write("PULSECONFIG=${12}\n")
          output.write("OUTGROUPBY=${13}\n")
          output.write("OUTTR=${14}\n")
          # if simpleAnaysis == True:
              # output.write("RESOLUTION=$8\n\n")
              # output.write("COUNTER=$9\n\n")
          output.write("set --\n\n")
          if local == False:
              output.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh\n")
          # old sources, left here for reference
          # output.write("source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
          # output.write("source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh\n")
          # output.write("source /cvmfs/geant4.cern.ch/geant4/10.3/x86_64-slc6-gcc49-opt/bin/geant4.sh\n\n")
          output.write("# do the MC\n")
          output.write("$COMMAND $TEMPLATE %s $GPS $SEED " %(output_temp))
          if useFlux == True:
              output.write("$FLUX")
          output.write("\n\n")
          if onlyEnDepo != True:
            output.write("# do the hybrid propagation\n")
            output.write("$PROPAGATE -i %s.root -c \"$CALIBRATION\" -o %s.root --photonSeed $SEED --type " %(output_temp,hybrid_temp))
            for tIndex in range(len(listTypes)):
                output.write("%d"  %listTypes[tIndex])
                if tIndex != (len(listTypes)-1):
                    output.write(",")
                else:
                    output.write(" ")
            output.write("\n\n")
          if pulse == True:
              output.write("### Do the GroupBy\n")  
              output.write("$PULSEGROUPBY -c $PULSECONFIG -i %s.root -o %s\n\n" %(hybrid_temp,groupby_temp ))
              # print("$PULSEGROUPBY -c $PULSECONFIG -i %s.root -o %s\n\n" %(hybrid_temp,groupby_temp ))
              output.write("### Do the pulse triggering\n")
              output.write("$PULSETR -c $PULSECONFIG -i %s.root -o %s \n" %(groupby_temp,tr_temp ) )
              # print("$PULSETR -c $PULSECONFIG -i %s.root -o %s \n" %(groupby_temp,tr_temp ) )
          # if the run is not just in local pc
          if local == False:
              # copy output files, if user wants it
              if keepEnergyDepositions == True:
                output.write("echo \"Copying output file...\"\n\n")
                output.write("xrdcp --nopbar %s.root $OUTPUT.root\n\n" %(output_temp) )
              # then the rest, if it's performed
              if onlyEnDepo != True:
                  if keepHybrid == True:
                    output.write("echo \"Copying hybrid file...\"\n\n")
                    output.write("xrdcp --nopbar %s.root $HYBRID.root\n\n" %(hybrid_temp) )
                  if pulse == True:
                    if keepGroupBy == True:
                      output.write("echo \"Copying GroupBy file...\"\n\n")
                      output.write("xrdcp --nopbar %s.root $OUTGROUPBY.root\n\n" %(groupby_temp))
                    # pulses always saved if pulses are produced (otherwise, why doing them?)
                    output.write("echo \"Copying Trigger file...\"\n\n")
                    output.write("xrdcp --nopbar %s.root $OUTTR.root\n\n" %(tr_temp))
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
          if(request_disk == True):
              output.write("request_disk            = " + request_disk_space + " \n")
          output.write("+JobFlavour             = \"%s\" \n" %(listQueue[index]))
          output.write("queue arguments from args.txt \n")
          output.close()

      #########################
      # GPS                  #
      #########################
      gpsFile = newFolder + "/gps.mac"
      gpsFile = os.path.abspath(gpsFile)


      if useFlux == True:
          # just write a new gps file, one line
          with open(gpsFile, "w") as output:
              output.write("/run/beamOn %d\n" % listEvents[index] )
              output.close()
      else:
          # modify base gps
          to_be_removed = ['/gps/energy' ,'/run/beamOn']
          with open(baseGPS, "r") as input:
            with open(gpsFile, "w") as output:
                for line in input:
                    if not any(bad_word in line for bad_word in to_be_removed):
                        output.write(line)
                output.write('\n')
                output.write("/gps/energy %s GeV\n" % str(listEnergy[index]) )
                output.write("/run/beamOn %d\n" % listEvents[index] )
                output.close()
      #########################
      # ARGS                  #
      #########################
      argsFile = newFolder + "/args.txt"
      startOfRandom = random.randint(1, 100000)
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
              output.write("%s/hybrid%d " %(newOutFolder,i) )   # output.write("HYBRID=$7\n")
              output.write("%d " %(i+startOfRandom) )   # SEED, cannot be 0!
              # if useFlux == True:
                  ## compose flux file name
              fluxFileName = baseFolderFlux + "/flux_" + str(i) + ".root";
              output.write("%s " %(fluxFileName) )      # FLUX FILE
              # if pulse == True:
              output.write("%s/simReadout " % build )                    #output.write("PULSEGROUPBY=$10\n")
              output.write("%s/ApplyCFD " % build )                      #output.write("PULSETR=$11\n")
              output.write("%s " % pulseConfigFile )                     #output.write("PULSECONFIG=$12\n")
              output.write("%s/OutGroupd_%d " %(newOutFolder,i) )   #output.write("OUTGROUPBY=$13\n")
              output.write("%s/OutTrigd_%d " %(newOutFolder,i) )    #output.write("OUTTR=$14\n")
              output.write("\n")
          output.close()
      index = index +1
  #
  # now do the submit script

  if local == True:
    subScript = "runLocal_" + script_identifier + ".py"
    with open(subScript, "w") as output:
      output.write("#!/usr/bin/python3\n")
      output.write("# -*- coding: utf-8 -*-\n")
      output.write("\n")
      output.write("import math\n")
      output.write("import os\n")
      output.write("import stat\n")
      output.write("import sys\n")
      output.write("import argparse\n")
      output.write("import subprocess\n")
      output.write("from subprocess import Popen, PIPE, STDOUT\n")
      output.write("import threading\n")
      output.write("import time\n")
      output.write("import multiprocessing\n")
      output.write("\n")
      output.write("def worker(path,line,count,en):\n")
      output.write("  \n")
      output.write("  comm = path + '/jobs/' + en + 'GeV/run_script.sh'\n")
      output.write("  cmd = [comm]\n")
      output.write("  for i in line.split():\n")
      output.write("    cmd.append(i)\n")
      output.write("  # print(cmd)\n")
      output.write("  logName = path + '/jobs/' + en + 'GeV/log/log_' + str(count) +  '.log'\n")
      output.write("  log = open(logName, 'w')\n")
      output.write("  subprocess.Popen(cmd,stdout = log,stderr=None).wait()\n")
      output.write("  log.close()\n")
      output.write("  return\n")
      output.write("\n")
      output.write("def main(argv):\n")
      output.write("  #parsing args\n")
      output.write("  parser = argparse.ArgumentParser(description='Python script to start jobs in parallel')\n")
      output.write("  parser.add_argument('--num' , default=-1, help='Number of PCs in the farm')\n")
      output.write("  parser.add_argument('--here', default=-1, help='Part of jobs to run on this pc')\n")
      output.write("  parser.add_argument('--cores', default=8, help='Number of cores to use')\n")
      output.write("  args = parser.parse_args()\n")
      output.write("  numberOfPCs = int(args.num)\n")
      output.write("  thisPC      = int(args.here)\n")
      output.write("  cores       = int(args.cores)\n")
      output.write("  splitJobs = False\n")
      output.write("  if numberOfPCs == -1:\n")
      output.write("    if thisPC != -1:\n")
      output.write("      print ('ERROR: if you specify --num you should specify also --here  ')\n")
      output.write("      sys.exit()\n")
      output.write("  if thisPC == -1:\n")
      output.write("    if numberOfPCs != -1:\n")
      output.write("      print ('ERROR: if you specify --num you should specify also --here  ')\n")
      output.write("      sys.exit()\n")
      output.write("  # now both -1 or both a number \n")
      output.write("  if thisPC != -1:\n")
      output.write("    if numberOfPCs != -1:\n")
      output.write("      splitJobs = True\n")
      output.write("  pool = multiprocessing.Pool(cores) \n")
      output.write("  # preparing \n")
      output.write("  energies = ['")
      countEnergies = 0
      for e in listEnergy:
        countEnergies = countEnergies + 1
        if countEnergies == len(listEnergy):
          output.write("%s']\n" % str(e) )
        else:
          output.write("%s','" % str(e) )
      output.write("  path = os.path.abspath('./')\n")
      output.write("  # take only 1/3 of the jobs if lineNum is specified \n")
      output.write("  processList = []\n")
      output.write("  for en in energies:\n")
      output.write("    filepath = path + '/jobs/' + en + 'GeV/args.txt'\n")
      output.write("    countLine = 0\n")
      output.write("    lines = len(open(filepath).readlines())\n")
      output.write("    limit = lines / numberOfPCs\n")
      output.write("    uplimit = limit * thisPC \n")
      output.write("    downlimit = limit * (thisPC -1) \n")
      output.write("    with open(filepath) as fp:\n")
      output.write("      for line in fp:\n")
      output.write("        if splitJobs == True:\n")
      output.write("          if countLine >= downlimit:\n")
      output.write("            if countLine < uplimit:\n")
      output.write("              pool.apply_async(worker, args=(path,line,countLine,en))\n")
      output.write("              #print(countLine)\n")
      output.write("        else:\n")
      output.write("          pool.apply_async(worker, args=(path,line,countLine,en))\n")
      output.write("        countLine = countLine + 1\n")
      output.write("  pool.close()\n")
      output.write("  pool.join()\n")
      output.write("\n")
      output.write("if __name__ == '__main__':\n")
      output.write("  main(sys.argv[1:])\n")
  else:
    subScript = "runAll_" + script_identifier + ".sh"
    with open(subScript, "w") as output:
        output.write("#!/bin/bash\n\n")
        output.write("for i in ")
        for i in listEnergy:
            output.write("%sGeV " % str(i))
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
