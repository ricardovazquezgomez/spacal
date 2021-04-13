###### USAGE:
### 		python3 runCalc.py [args]
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

import threading
import time
import multiprocessing


def worker(baseFolder,buildFolder,listEnergies,detector,angle_x,angle_y):
    ##
    print(listEnergies)
    ## delet previous summary and plots file if they exist
    fileSummary = 'summary_' + detector + ".txt"
    removePath = os.path.join(baseFolder, fileSummary)
    if os.path.exists(removePath):
      os.remove(removePath)
    else:
      print("No need to remove summary file")



    ### prepare command
    executable = os.path.join(buildFolder, "calculateTimeAndEnergy")

    for energy in listEnergies:
        ##log file
        logName = 'log_' + detector + '.log'
        logFile = os.path.join(baseFolder, logName)
        log = open(logFile, 'w')
        ## cmd
        cmd = [executable]

        ## energy folder
        energyFolder = baseFolder + "/" + str(energy) + "GeV/"

        ## remove plots files
        filePlots = 'plots_' + detector + ".root"
        removePlotsPath = os.path.join(energyFolder, filePlots)
        if os.path.exists(removePlotsPath):
          os.remove(removePlotsPath)
        else:
          print("No need to remove plots file")

        ## remove small files
        # get all files startin with resolution_
        prefixResolution = 'resolution_' + detector
        allFiles = [i for i in os.listdir(energyFolder) if os.path.isfile(os.path.join(energyFolder,i)) and prefixResolution in i]
        # delete files smaller than 10k
        # print(allFiles)
        for file in allFiles:
            fileNameRes = os.path.join(energyFolder,file)
            if os.path.getsize(fileNameRes) < 10 * 1024:
                os.remove(fileNameRes)

        ## command for analysis
        cmd.append('-f')
        cmd.append(energyFolder)
        cmd.append('-i')
        filePrefix = "resolution_" + detector
        cmd.append(filePrefix)
        cmd.append('-o')
        outputFilePrefix = "plots_" + detector
        outputFilePath = os.path.join(energyFolder, outputFilePrefix)
        cmd.append(outputFilePath)
        cmd.append('--energy')
        cmd.append(str(energy))
        cmd.append('--angle_x')
        cmd.append(str(angle_x))
        cmd.append('--angle_y')
        cmd.append(str(angle_y))
        # print(cmd)
        ### exectute
        subprocess.Popen(cmd,stdout = log,stderr=None).wait()
        log.close()
        print ("Detector %s energy %d analysis done" %(detector,energy) )
        cmd = []

        plotsTxtFile = outputFilePath + ".txt"
        last_line = ""
        with open(plotsTxtFile, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        summary_file_out = "summary_" + detector + ".txt"
        targetFile = os.path.join(baseFolder, summary_file_out)
        with open(targetFile, "a") as myfile:
            myfile.write(last_line)
            myfile.write("\n")

    return


def main(args):

    ###parsing command line args
    parser = argparse.ArgumentParser(description='Python script to automatically analyze resolution data files')
    parser.add_argument('--baseFolder'   , help='Output folder',required=True)
    parser.add_argument('--build'        , help='Build folder' ,required=True)
    parser.add_argument('--angle_x'      , help='Angle X'      ,required=True)
    parser.add_argument('--angle_y'      , help='Angle Y'      ,required=True)
    parser.add_argument("--energies",    # 1 2 3 4 5 10 25 50 75 100
                         nargs="*",      # 0 or more values expected => creates a list
                         type=int,
                         required=True
    )
    parser.add_argument("--detectors",    # lhcbPMT H13700 H13543  HPK
                         nargs="*",       # 0 or more values expected => creates a list
                         type=str,
                         required=True
    )
    args = parser.parse_args()

    # assign variables
    baseFolder      = args.baseFolder
    buildFolder     = args.build
    listEnergies    = args.energies
    listDetectors   = args.detectors
    angle_x         = args.angle_x
    angle_y         = args.angle_y

    ### FEEDBACK

    print("")
    print("########################################")
    print("# FEEDBACK                             #")
    print("########################################")
    print("")
    print("baseFolder        = %s" % baseFolder)
    print("buildFolder       = %s" % buildFolder)
    print("Energies          = %s" % listEnergies)
    print("Detectors         = %s" % listDetectors)
    print("")

    ### move to target folder
    if os.path.exists(baseFolder) == False:
        print("ERROR! baseFolder does not exist! Aborting...")
        sys.exit(1)
    ### start analysis. detectors run in parallel
    baseFolder = os.path.abspath(baseFolder) ## get abs paths
    buildFolder = os.path.abspath(buildFolder) ## get abs paths
    processList = []
    for detector in listDetectors:
        proc = multiprocessing.Process(target=worker, args=(baseFolder,buildFolder,listEnergies,detector,angle_x,angle_y) )
        processList.append(proc)

    #start processes
    for i in range(len(processList)):
        processList[i].start()
        # processList[i].join()





    ### and it's all...
    # print("Done. Bye!")
    # sys.exit()




if __name__ == "__main__":
     main(sys.argv[1:])
