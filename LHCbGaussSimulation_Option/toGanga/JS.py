AppDt = GaudiExec()
AppDt.directory = '/afs/cern.ch/user/z/zhangy/workdir/EcalUpgrade/GaussDev_v52r2'
#AppDt.platform  = 'x86_64-centos7-gcc7-opt'
AppDt.platform  = 'x86_64-centos7-gcc62-opt'

def create_job(Name,Application,OptsFiles,Outputs,evtsPJ=100,nJob=100): 

   Application.options =OptsFiles
   job=Job(
	   name = Name,
	   application  = Application,
	   splitter     = GaussSplitter( eventsPerJob=evtsPJ, numberOfJobs = nJob, firstEventNumber=0),
	   inputsandbox = [],
	   outputfiles  = Outputs,#[DiracFile(OutputTuple)],
	   backend      = Dirac(),
	   #inputdata    = BKQuery(bkk_directory,dqflag=['OK','UNCHECKED']).getDataset()
	   )
   job.submit()

#
#create_job(Name="Gen_11102432",
#        Application=AppDt,
#        OptsFiles=["gauss_upgrade_signal.py","11102432.py"],
#        Outputs=[LocalFile("*.root"),LocalFile("*.sim"),LocalFile("*.xgen")]
#        )
#create_job(Name="MC_MB",
#        Application=AppDt,
#        OptsFiles=["gauss_upgrade_MB.py"],
#        Outputs=[LocalFile("*.root")],
#        evtsPJ=1000,
#        nJob=100
#        )


create_job(Name="Gen_15102430",
        Application=AppDt,
        OptsFiles=["gauss_upgrade_signal.py","15102430.py"],
        Outputs=[LocalFile("*.root"),LocalFile("*.sim"),LocalFile("*.xgen")]
        )
