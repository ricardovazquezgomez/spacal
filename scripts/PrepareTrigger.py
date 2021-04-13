#!/usr/bin/python3

##############################################################################
### Python script to prepare the jobs to analyse the Hybrid SPACAL simulation
### Based on the work of M. Pizzichemi, CERN,  marco.pizzichemi@cern.ch
#####################
### USAGE EXAMPLE
### python3 ../PrepareSignalFormation.py  --baseFolderJobs . --baseFolderOut . --baseFolderIn data/ --config ../ConfigFile.cfg --build build/ --inputPrefix hybrid --listFilters ../Filters.root
#########

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
  

    #parsing args
    parser = argparse.ArgumentParser(description='Python script to prepare the hybrid analysis')
    parser.add_argument('--baseFolderJobs'   , default=''    , help='Jobs folder',required=True)
    parser.add_argument('--baseFolderOut'   , default=''    , help='Base folders for the output files', required=True)
    parser.add_argument('--baseFolderIn'   , default=''    , help='Base folders for the input files', required=True)
    parser.add_argument('--buildFolder'                          , help='path to op build',required=True)
    parser.add_argument('--config'                         , help='Base config file',required=True)
    parser.add_argument('--queue'   , default='microcentury')
    parser.add_argument('--output1'   , default='OutTrigd')
    parser.add_argument('--inputFolderKey'   ,  default='GeV', help='Keyword of the input folders, to be looked for in the BaseFolderIn', required=True)
    parser.add_argument('--inputFilePrefix'   ,  default='hybrid', help='Prefix of the input data files', required=True)
    parser.add_argument("--listFilters",    # name on the CLI - drop the `--` for positional/required parameters
                        nargs="*",  # 0 or more values expected => creates a list
                        type=str,
                        default=[],
                        help='Filters File!'
    )
    parser.add_argument("--listPFrac",    # name on the CLI - drop the `--` for positional/required parameters
                        nargs="*",  # 0 or more values expected => creates a list
                        type=str,
                        default=[],
                        help='List of thresholds for the CFD technique'
    )
    parser.add_argument('--AccountingGroup'   ,  default='', help='Accounting Group for these jobs')
    parser.add_argument('--delete', default='False', action='store_true')
    parser.add_argument('--saveInJobF', default='False', action='store_true')
    parser.add_argument('--verbose', default='False', action='store_true')





    args = parser.parse_args()









    baseFolderJobs  = args.baseFolderJobs
    baseFolderOut   = args.baseFolderOut
    baseFolderIn    = args.baseFolderIn
    buildFolder     = args.buildFolder
    config          = args.config
    queue           = args.queue
    inputFoldKey    = args.inputFolderKey
    inputPrefix     = args.inputFilePrefix
    listFilters     = args.listFilters
    output1         = args.output1
    accountingGroup = args.AccountingGroup
    peakFracs       = args.listPFrac
    deleteFiles     = args.delete
    saveinJobF      = args.saveInJobF
    verbose         = args.verbose





    ### Get the absolute path
    baseFolderJobs  = os.path.abspath(baseFolderJobs)
    baseFolderOut   = os.path.abspath(baseFolderOut)
    baseFolderIn    = os.path.abspath(baseFolderIn)
    # baseFolderOut   = list(map(os.path.abspath, baseFolderOut))
    # baseFolderIn    = os.path.abspath(baseFolderIn)
    buildFolder     = os.path.abspath(buildFolder)
    config          = os.path.abspath(config)
    # for i in range(len(listFilters)):
    #     listFilters[i] = os.path.abspath(listFilters[i])
    listFilters    = list(map(os.path.abspath, listFilters))









    # ### Check that the number of output folders corresponds to the input ones!
    # if not (len(baseFolderIn) == len(baseFolderOut)):
    #     print ("ERROR! number of input (%i) and output (%i) folders is different! Aborting..." % (len(baseFolderIn), len(baseFolderOut)))
    #     sys.exit(1)

    ### add the ending '/' if the paths do not end with it
    if not (baseFolderJobs[-1] == '/'):
        baseFolderJobs += '/'
    if not (baseFolderIn[-1] == '/'):
        baseFolderIn += '/'
    if not (baseFolderOut[-1] == '/'):
        baseFolderOut += '/'
    if not (buildFolder[-1] == '/'):
        buildFolder += '/'




    ### Look in the input basefolder to find folders matching with the expected suffix
    inputFolderList = []
    # inputSuffix     = []
    for file in os.listdir(baseFolderIn):
        if inputFoldKey in file:
            # inputSuffix.append(file)
            # inputFolderList.append(baseFolderIn + file)
            inputFolderList.append(file)



    # for i in range(len(inputFolderList)):
    #     if not (inputFolderList[i][-1] == '/'):
    #         inputFolderList[i] += '/'



    if (saveinJobF == True):
        print("Saving the output in the Job folder.\n")
        baseFolderOut = baseFolderJobs


    print("")
    print("########################################")
    print("# F E E D B A C K                      #")
    print("########################################")
    print("")
    print("baseFolderJobs          = %s" % baseFolderJobs)
    print("baseFolderOut           = %s" % baseFolderOut)
    print("baseFolderIn            = %s" % baseFolderIn)
    print("InputFolderList         = %s" % inputFolderList)
    print("buildFolder             = %s" % buildFolder)
    print("configFile              = %s" % config)
    print("queue                   = %s" % queue)
    print("inputPrefix             = %s" % inputPrefix)
    print("InputFolderKeyword      = %s" % inputFoldKey)
    print("listFilters             = %s" % listFilters)
    print("peakFracs               = %s" % peakFracs)
    print("deleteFiles             = %s" % deleteFiles)
    print("saveInJobF              = %s" % saveinJobF)
    if not (peakFracs):
        print("")
        print("CFD Threshold will be taken from the configuration file")
    print("")
    print("For more information use --verbose")
    print("")









    jobslist = []



    print ("#####################################")
    ### Job base folder
    if not os.path.exists(baseFolderJobs):
        print("Creating base job folder at %s" %(baseFolderJobs))
        os.mkdir(baseFolderJobs)
    else:
        print("Job folder exists at %s" %(baseFolderJobs))

    ### Output base folder
    if not os.path.exists(baseFolderOut):
        print("Creating base Output folder at %s" %(baseFolderOut))
        os.mkdir(baseFolderOut)
    else:
        print("Output folder exists at %s" %(baseFolderOut))
    print ("#####################################\n")


    #####################################
    ### FOR loop on the input folders ###
    for inputFolder in inputFolderList:


        ###########################################################
        ### Scan the Input folder and take all the files to analyse
        print ("###################")
        currentFolderIn = baseFolderIn + inputFolder + "/"
        print ("Opening folder %s" %(currentFolderIn))
        suffixes = []

        inputfiles  = []
        inputsuffix = []
        for file in os.listdir(currentFolderIn):
            if file.startswith(inputPrefix):
                inputfiles.append(file)
                inputsuffix.append(file[int(len(inputPrefix)): int(file.rfind('.'))])
        suffixes.append(inputsuffix)

        if not (inputfiles):   ### Check that you actually found some files.
            print ("ERROR! No input files found. Skipping the folder...")
            # sys.exit(1)
            continue

        print ("Found %i file(s) to analyse." % len(inputfiles))
        if (verbose == True): 
            print ("Input files:", inputfiles)
        if (verbose == True):
            print ("Input suffixes:", inputsuffix)
        




        ###################################
        ### Create Jobs related folders ###
        # print(baseFolderJobs)
        newFolder = baseFolderJobs + inputFolder + "/"

        if not os.path.exists(newFolder):
            print("Creating job folder at %s" %(newFolder))
            os.mkdir(newFolder)
        else:
            print("Job folder exists at %s" %(newFolder))

        jobsFolderJobs    = newFolder + "jobs/"
        if not os.path.exists(jobsFolderJobs):
            os.mkdir(jobsFolderJobs)
        else:
            print("Jobs folder exists at %s" %(jobsFolderJobs))

        ### Create other Jobs-related folders needed
        outDir = jobsFolderJobs + "out"
        errDir = jobsFolderJobs + "err"
        logDir = jobsFolderJobs + "log"
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
 
        ### Create Jobs related folders ###
        ###################################

        ##############################
        ### Create Output folders ###
        # print(baseFolderJobs)
        currentFolderOut = baseFolderOut + inputFolder + "/"
        if not os.path.exists(currentFolderOut):
            print("Creating output folder at %s" %(currentFolderOut))
            os.mkdir(currentFolderOut)
        else:
            print("Output folder exists at %s" %(currentFolderOut))
        ### Create Output folders ###
        ##############################



        #########################
        # R U N   S C R I P T   #
        #########################
        run_script = jobsFolderJobs + "run_script.sh"
        with open(run_script, "w") as output:

            output.write("#!/bin/bash\n\n")
            output.write("#parse input, they are written in args.txt\n\n")
            output.write("CONFIG_F=$1\n")
            output.write("INPUT_F=$2\n")
            output.write("OUTPUT_F1=$3\n")
            # output.write("FILTERS=$5\n")

            output.write("\nset --\n\n")
            output.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh\n")
            output.write("### Do the pulse triggering\n")
            if (peakFracs):
                for peakFrac in peakFracs:
                    output.write( buildFolder + "ApplyCFD   -c $CONFIG_F -i $INPUT_F -o ${OUTPUT_F1}_CFD" + str(float(peakFrac)) + " -t " + str(float(peakFrac)) + "\n")
            else:
                output.write( buildFolder + "ApplyCFD   -c $CONFIG_F -i $INPUT_F -o $OUTPUT_F1\n")



            output.close()
        #and make executable
        st = os.stat(run_script)
        os.chmod(run_script, st.st_mode | stat.S_IEXEC)













        #########################
        # J O B   S U B         #
        #########################
        jobsub = jobsFolderJobs + "jobs.sub"
        with open(jobsub, "w") as output:
            output.write("executable              = " + jobsFolderJobs + "run_script.sh \n")
            output.write("output                  = " + jobsFolderJobs + "out/out.$(ClusterId).$(ProcId).out \n")
            output.write("error                   = " + jobsFolderJobs + "err/err.$(ClusterId).$(ProcId).err \n")
            output.write("log                     = " + jobsFolderJobs + "log/log.$(ClusterId).log  \n")
            output.write("+JobFlavour             = \"%s\" \n" %queue)
            if (accountingGroup):
                output.write("+AccountingGroup        = \"" + accountingGroup + "\"\n")
            output.write("queue arguments from " + jobsFolderJobs + "args.txt \n")

            output.close()
        jobslist.append(jobsub)






        #########################
        # A R G S               #
        #########################
        argsFile = jobsFolderJobs + "args.txt"
        with open(argsFile, "w") as output:
            for i, inputFile in enumerate(inputfiles):

                output.write("%s " % config)                                                    # CONFIG_F=$1
                output.write("%s " % (currentFolderIn  + inputFile) )                           # INPUT_F=$3
                output.write("%s " % (currentFolderOut + output1 + "_" + inputsuffix[i]) )      # OUTPUT_F1=$4

                output.write("\n")
            output.close()
        print ("###################\n")
    ### FOR loop on the input folders ###
    #####################################









    #########################
    # S U B M I T           #
    #########################
    subScript = baseFolderJobs + "submitALL_script.sh"
    with open(subScript, "w") as output:
        output.write("#!/bin/bash\n\n")
        for job in jobslist:
            output.write("echo " + job + "\n")
            output.write("condor_submit " + job + "\n")
        output.close()
    #and make executable
    st = os.stat(subScript)
    os.chmod(subScript, st.st_mode | stat.S_IEXEC)  

    ### and it's all...
    print("")
    print("#########")
    print("Done. Please run %s to submit all jobs. Bye!" % subScript)
    print("#########")
    sys.exit()
















if __name__ == "__main__":
   main(sys.argv[1:])

