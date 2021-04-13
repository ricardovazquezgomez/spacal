###### USAGE:
### 		python3 prepareOpticalCalibration.py [args]
### see README.md for details on [args]
# Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

import math
import os
import stat
import sys
import argparse
import subprocess
import random
from subprocess import Popen, PIPE, STDOUT
import shutil
from shutil import copyfile
from decimal import *
import csv
import ntpath


def main(args):

    #parsing args
    parser = argparse.ArgumentParser(description='Python script to prepare optical calibration run')
    parser.add_argument('--config'                         , help='Base config file',required=True)
    parser.add_argument('--baseFolderJobs'   , default=''    , help='Output folder',required=True)
    parser.add_argument('--baseFolderOut'     , default=''    , help='Output folder',required=True)
    parser.add_argument('--build'                          , help='path to spacal build folder',required=True)
    parser.add_argument('--queue'            , default=""        , help='JobFlavour for condor',required=True)
    parser.add_argument('--primaries'        , default='2000000'  , help='Number of optical photons generated for each (x,y,z,energy) point')
    parser.add_argument('--jobs'             , default='1000'    , help='Number of parallel jobs')
    # parser.add_argument('--xmin'             , default='0'       , help='Min x of crystals [mm]')
    # parser.add_argument('--xmax'             , default='0'       , help='Max x of crystals [mm]')
    parser.add_argument('--xn'               , default='1'       , help='Number of points in x scan')
    # parser.add_argument('--ymin'             , default='0'       , help='Min y of crystals [mm]')
    # parser.add_argument('--ymax'             , default='0'       , help='Max y of crystals [mm]')
    parser.add_argument('--yn'               , default='1'       , help='Number of points in y scan')
    # parser.add_argument('--zmin'             , default='0'       , help='Min z of crystals [mm]')
    # parser.add_argument('--zmax'             , default='0'       , help='Max z of crystals [mm]')
    parser.add_argument('--zn'               , default='1'       , help='Number of points in z scan')
    parser.add_argument('--emin'             , default='1'       , help='Min optical photon energy value [eV]')
    parser.add_argument('--emax'             , default='5'       , help='Max optical photon energy value [eV]')
    parser.add_argument('--en'               , default='1'       , help='Number of points in energy scan')
    parser.add_argument('--local', action='store_true', help='Jobs will be run locally')
    # parser.add_argument('--strategy'   , default='0'       , help='Simulation strategy',required=True)
    # parser.add_argument('--optConfig'                      , help='Optical config file',required=True)
    # parser.add_argument('--crystalMaterial' , default='6'       , help='1) Quartz 2) SiO2:Ce 3) DSB:Ce 4) LuAG 5) YAG 6) GAGG 7) Water')
    args = parser.parse_args()

    # assign variables
    # jobFlavour = ""
    configFile      = args.config
    baseFolderOut   = args.baseFolderOut
    baseFolderJobs  = args.baseFolderJobs
    buildFolder     = args.build

    # xmin        = float(args.xmin)
    # xmax        = float(args.xmax)
    xn          = int  (args.xn)
    # ymin        = float(args.ymin)
    # ymax        = float(args.ymax)
    yn          = int  (args.yn)
    # zmin        = float(args.zmin)
    # zmax        = float(args.zmax)
    zn          = int  (args.zn)
    emin        = float(args.emin)
    emax        = float(args.emax)
    en          = int  (args.en)

    queue       = args.queue
    primaries   = int(args.primaries)
    jobs        = int(args.jobs)

    local       = args.local

    #############################################
    ### PARSE AND MODIFY ORIGINAL CONFIG FILE ###
    #############################################
    
    ## DICTIONARIES
    # dictionary of strings to read 
    stringDict = {
      "pos_z"            : ['cell_pos_z']  ,
      "size_z"           : ['cell_crystal_size_z']  ,
      "size_x"           : ['cell_crystal_size_x']  ,
      "size_y"           : ['cell_crystal_size_y']  ,
      "crystal_material" : ['cell_crystal_material'],
      "air_layer"        : ['cell_air_layer'],
      "int_gap_material" : ['cell_int_gap_material'],
      "cry_shape"        : ['cell_crystal_shape'],
      "hole_shape"       : ['cell_hole_shape'],
      "cladding"         : ['cell_crystal_cladding'],
      "staggering"       : ['cell_staggering'],
      "staggering_axis"  : ['cell_staggering_axis'],
      "staggering_size"  : ['cell_staggering_size'],
      "staggering_parity": ['cell_staggering_parity'],
      "staggering_remove": ['cell_staggering_remove']
    }
    # dictionary with stored values  
    valueDict = {
      "name"             : [],
      "pos_x"            : [],
      "pos_y"            : [],
      "pos_z"            : [],
      "x_elements"       : [],
      "y_elements"       : [],
      "size_z"           : [],
      "size_x"           : [],
      "size_y"           : [],
      "pitch_x"          : [],
      "pitch_y"          : [],
      "crystal_material" : [],
      "air_layer"        : [],
      "int_gap_material" : [],
      "cry_shape"        : [],
      "hole_shape"       : [],
      "cladding"         : [],
      "staggering"       : [],
      "staggering_axis"  : [],
      "staggering_size"  : [],
      "staggering_parity": [],
      "staggering_remove": []
    }
    # dictionary with all key names (also the ones that we don't need to read) 
    allDict = {
      "name"             : ["cell_name"],
      "pos_x"            : ["cell_pos_x"],
      "pos_y"            : ["cell_pos_y"],
      "pos_z"            : ["cell_pos_z"],
      "x_elements"       : ["cell_x_elements"],
      "y_elements"       : ["cell_y_elements"],
      "size_z"           : ["cell_crystal_size_z"],
      "size_x"           : ["cell_crystal_size_x"],
      "size_y"           : ["cell_crystal_size_y"],
      "pitch_x"          : ["cell_crystal_pitch_x"],
      "pitch_y"          : ["cell_crystal_pitch_y"],
      "crystal_material" : ["cell_crystal_material"],
      "air_layer"        : ["cell_air_layer"],
      "int_gap_material" : ["cell_int_gap_material"],
      "cry_shape"        : ["cell_crystal_shape"],
      "hole_shape"       : ["cell_hole_shape"],
      "cladding"         : ["cell_crystal_cladding"],
      "staggering"       : ["cell_staggering"],
      "staggering_axis"  : ["cell_staggering_axis"],
      "staggering_size"  : ["cell_staggering_size"],
      "staggering_parity": ["cell_staggering_parity"],
      "staggering_remove": ["cell_staggering_remove"]
    }

    ## RUN ON ORIGINAL CONFIG, EXTRACT THE INFO
    with open(configFile, "r") as input:    
      for line in input:
        # inDict = False
        for key in stringDict:
          if any(key_word in line for key_word in stringDict[key]):
            # check if not a comment
            # print(line)
            if not (line.strip().startswith('#')):
              # print(line)
              # delete everything before = included (hence +1)
              line = line[line.index('=')+1:].strip()
              # remove trailing and leading |
              line = line.rstrip('|').lstrip('|')
              # print(line)
              arr = line.split('|')
              arr = list(map(float, arr)) 
              unique_words = set(arr)
              sections = len(unique_words)
              # print(sections)
              valueDict[key].append(arr[0])
              if sections == 2:
                valueDict[key].append(arr[len(arr)-1])
        # now if it was in dict 

      # print(valueDict)
      input.close()
    
    ## NOW PROPERLY FILL THE DICTIONARY
    # first assumption: the number of sections (1 or 2)
    # is derived from len(valueDict["pos_z"])
    moduleSections = len(valueDict["pos_z"])
    if(moduleSections < 1 or moduleSections > 2):
      print("ERROR invalid number of sections found = %i\n" %(moduleSections))
      sys.exit()
    # else:
    #   print("Number of sections found = %i\n" %(moduleSections))
    ### one by one, on values dictionary (no easy way to solve this with a smart loop..)
    
    ## "name"             
    valueDict["name"].append(0)
    if moduleSections == 2:
      valueDict["name"].append(1)
    
    ## "pos_x"   
    valueDict["pos_x"].append(0)
    if moduleSections == 2:
      valueDict["pos_x"].append(0)      
    
    ## "pos_y"   
    valueDict["pos_y"].append(0)
    if moduleSections == 2:
      valueDict["pos_y"].append(0)            
    
    ## "pos_z"   
    # stay as it is 
    
    ## "x_elements"
    valueDict["x_elements"].append(1)
    if moduleSections == 2:
      valueDict["x_elements"].append(1)

    ## "y_elements"    
    valueDict["y_elements"].append(1)
    if moduleSections == 2:
      valueDict["y_elements"].append(1)   

    ## "size_z" 
    # stay as it is

    ## "size_x" 
    # needs to be same size of moduleSections
    # if it is, nothing to do
    # otherwise, enlarge by copy
    if not(len(valueDict["size_x"]) == moduleSections):    
      valueDict["size_x"].append(valueDict["size_x"][0])

    ## "size_y" 
    # needs to be same size of moduleSections
    # if it is, nothing to do
    # otherwise, enlarge by copy
    if not(len(valueDict["size_y"]) == moduleSections):    
      valueDict["size_y"].append(valueDict["size_y"][0])  
    
    ## "pitch_x" 
    valueDict["pitch_x"].append(0)
    if moduleSections == 2:
      valueDict["pitch_x"].append(0)   

    ## "pitch_y" 
    valueDict["pitch_y"].append(0)
    if moduleSections == 2:
      valueDict["pitch_y"].append(0)   

    ## "crystal_material" 
    # needs to be same size of moduleSections
    # if it is, nothing to do
    # otherwise, enlarge by copy
    if not(len(valueDict["crystal_material"]) == moduleSections):    
      valueDict["crystal_material"].append(valueDict["crystal_material"][0])  

    ## "air_layer"        
    if not(len(valueDict["air_layer"]) == moduleSections):    
      valueDict["air_layer"].append(valueDict["air_layer"][0])  
    
    ## "int_gap_material" 
    if not(len(valueDict["int_gap_material"]) == moduleSections):    
      valueDict["int_gap_material"].append(valueDict["int_gap_material"][0])  
    
    ### LAST 3, crystal/hole shape, cladding, could be not specified 
    ### (they would be 0 0 0 by default), so fill them only if already not empty
    
    ## "cry_shape"
    if(len(valueDict["cry_shape"]) > 0):
      if not(len(valueDict["cry_shape"]) == moduleSections):    
        valueDict["cry_shape"].append(valueDict["cry_shape"][0])      
    
    ## "hole_shape"  
    if(len(valueDict["hole_shape"]) > 0):
      if not(len(valueDict["hole_shape"]) == moduleSections):    
        valueDict["hole_shape"].append(valueDict["hole_shape"][0])    

    ## "cladding"
    if(len(valueDict["cladding"]) > 0):
      if not(len(valueDict["cladding"]) == moduleSections):    
        valueDict["cladding"].append(valueDict["cladding"][0])   

    ## "staggering"
    if(len(valueDict["staggering"]) > 0):
      if not(len(valueDict["staggering"]) == moduleSections):    
        valueDict["staggering"].append(valueDict["staggering"][0])   

    ## "staggering_axis"  
    if(len(valueDict["staggering_axis"]) > 0):
      if not(len(valueDict["staggering_axis"]) == moduleSections):    
        valueDict["staggering_axis"].append(valueDict["staggering_axis"][0])  

    ## "staggering_size"  
    if(len(valueDict["staggering_size"]) > 0):
      if not(len(valueDict["staggering_size"]) == moduleSections):    
        valueDict["staggering_size"].append(valueDict["staggering_size"][0])   

    ## "staggering_parity"
    if(len(valueDict["staggering_parity"]) > 0):
      if not(len(valueDict["staggering_parity"]) == moduleSections):    
        valueDict["staggering_parity"].append(valueDict["staggering_parity"][0])   

    ## "staggering_remove" 
    if(len(valueDict["staggering_remove"]) > 0):
      if not(len(valueDict["staggering_remove"]) == moduleSections):    
        valueDict["staggering_remove"].append(valueDict["staggering_remove"][0])    
    
    ## END OF DICT FILL
    print(valueDict)

    ## RUN AGAIN AND REPLACE THE LINES 
    ## ON NEW OPTICALI FILE 
    optiBasename = os.path.basename(configFile)
    optiCaliFile = "optiCali_" + optiBasename
    to_be_changed = ['simulationType']
    with open(configFile, "r") as input:   
      with open(optiCaliFile, "w") as output: 
        for line in input:
          inDict = False
          for key in allDict:
            if any(key_word in line for key_word in allDict[key]):
              # check if not a comment
              if not (line.strip().startswith('#')):
                inDict = True
                # don't write this line and write the new version
                # if the corresponding value in valueDict is not empty (which should not happen)
                if(len(valueDict[key]) > 0):
                  output.write(allDict[key][0])
                  output.write(" = |")
                  for value in valueDict[key]:
                    output.write(str(value))
                    output.write("|")
                  output.write("\n")
          # now if it was not in dict 
          if inDict == False:
            if not any(key_word in line for key_word in to_be_changed):
              output.write(line)
            else:
              if not (line.strip().startswith('#')):
                output.write('simulationType = 3\n')
              else:
                output.write(line)
        #
        output.close()
      input.close()

    #############################################
    ### END OF PARSE AND MODIFY               ###
    #############################################

    #############################################
    ### FIND MODULE LIMITS IN XYZ             ###
    #############################################
    ## rather than relying on the user to input them
    ## read them from the config file 
    ## this can be calculated from 
    # absorber_size_z
    # absorber_pos_z
    # cell_separation_type
    # separation_thickness
    # moduleZshift
    # LAPPD_layers_thickness
    ## so read all of them...
    limitsNamesDict = {
      "absorber_size_z"        : ["absorber_size_z"       ],
      "absorber_pos_z"         : ["absorber_pos_z"        ],
      "cell_separation_type"   : ["cell_separation_type"  ],
      "separation_thickness"   : ["separation_thickness"  ],
      "moduleZshift"           : ["moduleZshift"          ],
      "LAPPD_layers_thickness" : ["LAPPD_layers_thickness"]
    }
    limitsValuesDict = {
      "absorber_size_z"        : [],
      "absorber_pos_z"         : [],
      "cell_separation_type"   : [],
      "separation_thickness"   : [],
      "moduleZshift"           : [],
      "LAPPD_layers_thickness" : []
    }
    with open(optiCaliFile, "r") as input:
      for line in input:
        for key in limitsNamesDict:
          if any(key_word in line for key_word in limitsNamesDict[key]):
            if not (line.strip().startswith('#')):
              line = line[line.index('=')+1:].strip()
              # remove trailing and leading |
              line = line.rstrip('|').lstrip('|')
              # print(line)
              arr = line.split('|')
              arr = list(map(float, arr)) 
              unique_words = set(arr)
              limitsValuesDict[key] = arr
            
      input.close()
    # print(limitsValuesDict)  
    ## calculate limits
    # total length of active part is
    totalLengthActive = float(limitsValuesDict["absorber_size_z"][0])
    if(int(limitsValuesDict["cell_separation_type"][0]) == 2):
      if(len(limitsValuesDict["LAPPD_layers_thickness"]) > 0):
        for parts in limitsValuesDict["LAPPD_layers_thickness"]:
          totalLengthActive += float(parts)
      else:
        totalLengthActive += float(limitsValuesDict["separation_thickness"][0])
    # add a default for missing keys (at least the ones that could be left not specified by the user)
    if(len(limitsValuesDict["moduleZshift"]) == 0):
      limitsValuesDict["moduleZshift"] = [0.]

    # now compute center in world coordinates
    moduleCenter = float(limitsValuesDict["absorber_pos_z"][0]) + float(limitsValuesDict["moduleZshift"][0])
    # and module min and max z
    moduleMinZ = moduleCenter - totalLengthActive/2.0
    moduleMaxZ = moduleCenter + totalLengthActive/2.0
    # for x and y, data is already acquired, because it's only the single crystal size 
    # only problem is to choose the greatest among the 2 values, in case they are different (they shouldn't, unless the user is [CENSORED])
    crySizeX = max(valueDict["size_x"])
    crySizeY = max(valueDict["size_y"])
    # and center is always 0,0, so
    cryMinX = - crySizeX/2.0
    cryMaxX = + crySizeX/2.0
    cryMinY = - crySizeY/2.0
    cryMaxY = + crySizeY/2.0

    

    ### calc x y z e steps
    ### NOW OVERWRITING xmin xmax etc 
    xmin = cryMinX
    xmax = cryMaxX
    ymin = cryMinY
    ymax = cryMaxY
    zmin = moduleMinZ
    zmax = moduleMaxZ

    xstep = round((xmax - xmin) / xn,6)
    ystep = round((ymax - ymin) / yn,6)
    zstep = round((zmax - zmin) / zn,6)
    estep = round((emax - emin) / en,6)

    ### compute lists
    xList = []
    yList = []
    zList = []
    eList = []
    getcontext().prec = 2
    for i in range(xn):
        xList.append( round(xmin + (xstep/2.0)  + i*xstep,6) )
    for i in range(yn):
        yList.append( round(ymin + (ystep/2.0)  + i*ystep,6) )
    for i in range(zn):
        zList.append( round(zmin + (zstep/2.0)  + i*zstep,6) )
    for i in range(en):
        eList.append( round(emin + (estep/2.0)  + i*estep,6) )

    primaries_per_point_per_job = primaries / jobs
    points_per_job      = xn*yn*zn*en
    primaries_per_job   = points_per_job*primaries_per_point_per_job

    print("")
    print("########################################")
    print("# FEEDBACK                             #")
    print("########################################")
    print("")
    print("Total active length [mm]      = %f" %totalLengthActive)
    print("Module center in z [mm]       = %f" %moduleCenter)
    print("Crystal min in x [mm]         = %f" %cryMinX)
    print("Crystal max in x [mm]         = %f" %cryMaxX)
    print("Crystal min in y [mm]         = %f" %cryMinY)
    print("Crystal max in y [mm]         = %f" %cryMaxY)
    print("Module min in z [mm]          = %f" %moduleMinZ)
    print("Module max in z [mm]          = %f" %moduleMaxZ)
    print("Base configFile               = %s" % configFile)
    print("optiCaliFile                  = %s" % optiCaliFile)
    print("baseFolderOut                 = %s" % baseFolderOut)
    print("baseFolderJobs                = %s" % baseFolderJobs)
    print("build                         = %s" % buildFolder)
    print("xmin                          = %f" % xmin )
    print("xmax                          = %f" % xmax )
    print("xn                            = %d" % xn   )
    print("ymin                          = %f" % ymin )
    print("ymax                          = %f" % ymax )
    print("yn                            = %d" % yn   )
    print("zmin                          = %f" % zmin )
    print("zmax                          = %f" % zmax )
    print("zn                            = %d" % zn   )
    print("emin                          = %f" % emin )
    print("emax                          = %f" % emax )
    print("en                            = %d" % en   )
    print("xstep                         = %f" % xstep   )
    print("ystep                         = %f" % ystep   )
    print("zstep                         = %f" % zstep   )
    print("estep                         = %f" % estep   )
    print("xList                         = %s" % xList   )
    print("yList                         = %s" % yList   )
    print("zList                         = %s" % zList   )
    print("eList                         = %s" % eList   )
    print("queue                         = %s" % queue   )
    print("primaries                     = %d" % primaries   )
    print("jobs                          = %d" % jobs   )
    print("Primaries per point per job   = %d" % primaries_per_point_per_job  )
    print("Primaries per job             = %d" % primaries_per_job  )
    print("local                         = %s" % local)


    if not os.path.exists(baseFolderJobs):
        print("Creating job folder at %s" %(baseFolderJobs))
        os.mkdir(baseFolderJobs)
    else:
        print("Job folder exists at %s" %(baseFolderJobs))
    script_identifier = baseFolderJobs
    baseFolderJobs = os.path.abspath(baseFolderJobs)
    #create output folder
    if not os.path.exists(baseFolderOut):
        print("Creating output folder at %s" %(baseFolderOut))
        os.mkdir(baseFolderOut)
    else:
        print("Output folder exists at %s" %(baseFolderOut))
    baseFolderOut = os.path.abspath(baseFolderOut)

    ### do condor folders
    outDir = baseFolderJobs + "/" + "out"
    errDir = baseFolderJobs + "/" + "err"
    logDir = baseFolderJobs + "/" + "log"
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
        print("Log folder exists at %s" %(logDir))

    ### copy config file
    configFileJob = baseFolderJobs + "/optConfig.cfg"
    # no more copy configFile, but optiCaliFile
    shutil.copy(optiCaliFile,configFileJob)

    #temp file names 
    output_temp="output_temp"
    if local == True:
        output_temp="$OUTPUT"

    #########################
    # RUN SCRIPT            #
    #########################
    run_script = baseFolderJobs + "/run_script.sh"
    with open(run_script, "w") as output:
        output.write("#!/bin/bash\n\n")
        output.write("#parse input, they are written in args.txt\n\n")
        output.write("COMMAND=$1\n")
        output.write("TEMPLATE=$2\n")
        output.write("OUTPUT=$3\n")
        output.write("GPS=$4\n")
        output.write("SEED=$5\n")
        output.write("set --\n\n")
        if local == False:
            output.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh\n")
        # old sources, left here for reference
        # output.write("source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
        # output.write("source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh\n")
        # output.write("source /cvmfs/geant4.cern.ch/geant4/10.3/x86_64-slc6-gcc49-opt/bin/geant4.sh\n\n")
        output.write("# do the mc\n")
        output.write("$COMMAND $TEMPLATE %s $GPS $SEED\n\n" %(output_temp))
        if local == False:
            output.write("echo -e \"Copying result to: $OUTPUT.root\"\n\n")
            output.write("xrdcp --nopbar output_temp.root $OUTPUT.root\n\n")
        output.close()
    #and make executable
    st = os.stat(run_script)
    os.chmod(run_script, st.st_mode | stat.S_IEXEC)

    #########################
    # JOBS.SUB            #
    #########################
    jobsub = baseFolderJobs + "/jobs.sub"
    with open(jobsub, "w") as output:
        output.write("executable              = run_script.sh \n")
        output.write("output                  = out/out.$(ClusterId).$(ProcId).out \n")
        output.write("error                   = err/err.$(ClusterId).$(ProcId).err \n")
        output.write("log                     = log/log.$(ClusterId).log  \n")
        output.write("transfer_output_files   = \"\" \n")
        output.write("+JobFlavour             = \"%s\" \n" %(queue))
        output.write("queue arguments from args.txt \n")
        output.close()

    #########################
    # GPS                  #
    #########################
    ###
    gpsFile = baseFolderJobs + "/gps.mac"
    gpsFile = os.path.abspath(gpsFile)
    with open(gpsFile, "w") as output:
        output.write("/gps/particle               opticalphoton\n")
        output.write("/gps/pos/type               Point\n")
        output.write("/gps/ang/type               iso\n")
        output.write('\n')
        ### loop on energies
        for e in eList:
            ### loop on x
            for x in xList:
                ### loop on y
                for y in yList:
                    ### loop on z
                    for z in zList:
                        output.write("/gps/energy  %f eV\n" % e)
                        output.write("/gps/pos/centre  %f %f %f mm\n" %(x,y,z))
                        output.write("/run/beamOn  %d\n" %(primaries_per_point_per_job))
                        output.write('\n')
        # output.write("/gps/energy %d GeV\n" % listEnergy[index] )
        # output.write("/run/beamOn %d\n" % listEvents[index] )
        output.close()


    #########################
    # ARGS                  #
    #########################
    # random number to start seeds 
    startOfRandom = random.randint(1, 100000)
    argsFile = baseFolderJobs + "/args.txt"
    with open(argsFile, "w") as output:
        for i in range(jobs):
            output.write("%s/FibresCalo " % buildFolder )                # output.write("COMMAND=$1\n")
            output.write("%s " % configFileJob )                      # output.write("TEMPLATE=$2\n")
            output.write("%s/out%d "%(baseFolderOut,i))             # output.write("OUTPUT=$3\n")
            output.write("%s "% gpsFile)                           # output.write("GPS=$4\n")
            output.write("%d "% (i+startOfRandom))                           # output.write("GPS=$4\n")
            output.write('\n')
        output.close()
    

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
        output.write("def worker(path,line,count):\n")
        output.write("  \n")
        output.write("  comm = '%s/run_script.sh'\n" %(baseFolderJobs))
        output.write("  cmd = [comm]\n")
        output.write("  for i in line.split():\n")
        output.write("    cmd.append(i)\n")
        output.write("  # print(cmd)\n")
        output.write("  logName = '%s/log/log_' + str(count) +  '.log'\n" %(baseFolderJobs))
        output.write("  log = open(logName, 'w')\n")
        output.write("  subprocess.Popen(cmd,stdout = log,stderr=log).wait()\n")
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
        # output.write("  energies = ['")
        # countEnergies = 0
        # for e in listEnergy:
        #   countEnergies = countEnergies + 1
        #   if countEnergies == len(listEnergy):
        #     output.write("%s']\n" % e )
        #   else:
        #     output.write("%s','" % e )
        output.write("  path = os.path.abspath('./')\n")
        output.write("  # take only 1/3 of the jobs if lineNum is specified \n")
        output.write("  processList = []\n")
        # output.write("  for en in energies:\n")
        output.write("  filepath = '%s/args.txt'\n" %(baseFolderJobs))
        output.write("  countLine = 0\n")
        output.write("  lines = len(open(filepath).readlines())\n")
        output.write("  limit = lines / numberOfPCs\n")
        output.write("  uplimit = limit * thisPC \n")
        output.write("  downlimit = limit * (thisPC -1) \n")
        output.write("  with open(filepath) as fp:\n")
        output.write("    for line in fp:\n")
        output.write("      if splitJobs == True:\n")
        output.write("        if countLine >= downlimit:\n")
        output.write("          if countLine < uplimit:\n")
        output.write("            pool.apply_async(worker, args=(path,line,countLine))\n")
        output.write("            #print(countLine)\n")
        output.write("      else:\n")
        output.write("        pool.apply_async(worker, args=(path,line,countLine))\n")
        output.write("      countLine = countLine + 1\n")
        output.write("  pool.close()\n")
        output.write("  pool.join()\n")
        output.write("\n")
        output.write("if __name__ == '__main__':\n")
        output.write("  main(sys.argv[1:])\n")
    else:
      subScript = "runAll_" + script_identifier + ".sh"
      with open(subScript, "w") as output:
        output.write("#!/bin/bash\n\n")
        # output.write("for i in ")
        # for i in listEnergy:
        #     output.write("%dGeV " % i)
        # output.write("\n")
        # output.write("do\n")
        output.write("cd %s\n\n" % baseFolderJobs)
        output.write("condor_submit jobs.sub\n\n")
        output.write("cd - \n")
        # output.write("done\n")
        output.close()
      #and make executable
      st = os.stat(subScript)
      os.chmod(subScript, st.st_mode | stat.S_IEXEC)

    sys.exit(0)




if __name__ == "__main__":
   main(sys.argv[1:])
