#!/usr/bin/env python
import sys
import os, subprocess
from glob import glob
import optparse
import ast
import re
import multiprocessing
from time import gmtime, strftime

# crabControl.py -- submit CRAB3 jobs and retrieve job information
#
# usage: crabControl.py --config datasamples.cfg
#
# comments @ viktor.kutzner@cern.ch

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def ShellExec(command):

    print color.BOLD + color.BLUE + command + color.END
    if not options.dryrun:
        os.system(command)

def readConfig(configFileName):
    
    samples = {}
    cfg_data = ""
    
    with open(configFileName, 'r') as f:
        cfg_data = f.read()

    currentSection = ""
    CMSConfigParameters = None
    CMSConfigParameterAdd = None
    addFiles = None
    mergeOutputFiles = None
    inputDBS = "global"
    
    for line in cfg_data.split("\n"):
       
        if len(line)>0 and line[0] == "#":
            continue
        
        if "[data]" in line:                  currentSection = "data"
        if "[MC]" in line:                    currentSection = "MC"
        if "workArea" in line:                workArea = line.split("=")[-1].strip()
        if "outputDatasetTagPrefix" in line:  outputDatasetTagPrefix = line.split("=")[-1].strip()
        if "storageSite" in line:             storageSite = line.split("=")[-1].strip()
        if "lumiMask" in line:                lumiMask = line.split("=")[-1].strip()
        if "runRange" in line:                runRange = line.split("=")[-1].strip()
        if "CMSConfigFile" in line:           CMSConfigFile = line.split("=")[-1].strip()
        if "CMSConfigPlugin" in line:         CMSConfigPlugin = line.split("=")[-1].strip()
        if "CMSConfigParameters" in line:     CMSConfigParameters = line.split("CMSConfigParameters = ")[-1].strip()
        if "CMSConfigParameterAdd" in line:   CMSConfigParameterAdd = line.split("CMSConfigParameterAdd = ")[-1].strip()
        if "addFiles" in line:                addFiles = line.split("addFiles = ")[-1].strip()
        if "mergeOutputFiles" in line:        mergeOutputFiles = line.split("mergeOutputFiles = ")[-1].strip()
                
        try:
            if "=" not in line and "[" not in line and len(line.split())>0:
                shortname = line.split(":")[-1].strip()
                samples[shortname] = {}
                samples[shortname]["datasetpath"] = line.split(":")[0]
                samples[shortname]["type"] = currentSection
                samples[shortname]["workArea"] = workArea
                samples[shortname]["outputDatasetTagPrefix"] = outputDatasetTagPrefix
                samples[shortname]["storageSite"] = storageSite
                samples[shortname]["lumiMask"] = lumiMask
                samples[shortname]["runRange"] = runRange
                samples[shortname]["CMSConfigFile"] = os.path.abspath(CMSConfigFile)
                samples[shortname]["CMSConfigPlugin"] = CMSConfigPlugin
                samples[shortname]["CMSConfigParameters"] = CMSConfigParameters
                samples[shortname]["CMSConfigParameterAdd"] = CMSConfigParameterAdd
                samples[shortname]["addFiles"] = addFiles
                samples[shortname]["mergeOutputFiles"] = mergeOutputFiles
                               
        except:
            print "Missing / malformed parameters in configuration file"
            quit()

    return samples

logo = "          ,__,\n     (/__/\oo/\__(/\n       _/\/__\/\_\n        _/    \_   "

def crabSubmit(samples):

    for sample in samples:
        tagname = samples[sample]["outputDatasetTagPrefix"]+"_"+sample
        cmd = "source /cvmfs/cms.cern.ch/crab3/crab.sh; crab submit -c crabsubmit_%s.py" % tagname
        ShellExec(cmd)

    print "\n" + logo + "done! \n"
    print " Check crab jobs with --status option."


def buildCrabSubmitScripts(samples):
    
    for sample in samples:

        item = samples[sample]

        tagname = item["outputDatasetTagPrefix"]+"_"+sample
        
        isData = False
        if item["type"] == "data":
            isData = True
            
        f = open('crabsubmit_%s.py' % tagname, 'w')
        
        f.write("from CRABClient.UserUtilities import config, getUsernameFromSiteDB\n")
        f.write("config = config()\n\n")
        f.write("config.General.requestName = '%s'\n" % tagname)
        f.write("config.General.workArea = '%s'\n" % item["workArea"])
        f.write("config.General.transferOutputs = True\n")
        f.write("config.General.transferLogs = True\n\n")
        
        f.write("config.JobType.pluginName = '%s'\n" % item["CMSConfigPlugin"])
        f.write("config.JobType.psetName = '%s'\n" % item["CMSConfigFile"])
        if item["CMSConfigParameters"]:
            if item["CMSConfigParameterAdd"]:
                item["CMSConfigParameters"] = item["CMSConfigParameters"] + ", " + item["CMSConfigParameterAdd"]
            f.write("config.JobType.pyCfgParams = [%s]\n\n" % item["CMSConfigParameters"])
        if item["addFiles"]:
            f.write("config.JobType.inputFiles = [%s]\n\n" % item["addFiles"])

        f.write("config.Data.inputDataset = '%s'\n" % item["datasetpath"])
        f.write("config.Data.inputDBS = 'global'\n")
        #~ f.write("config.Data.inputDBS = 'phys03'\n")
        if isData:
            f.write("config.Data.splitting = 'LumiBased'\n")
            f.write("config.Data.unitsPerJob = 15\n")
            f.write("config.Data.lumiMask = '%s'\n" % item["lumiMask"])
            f.write("config.Data.runRange = '%s'\n" % item["runRange"])
        else:
            f.write("config.Data.splitting = 'FileBased'\n")
            f.write("config.Data.unitsPerJob = 2\n")
        f.write("config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())\n")
        f.write("config.Data.publication = True\n")
        f.write("config.Data.outputDatasetTag = '%s'\n\n" % tagname)
        f.write("config.Site.storageSite = 'T2_PT_NCG_Lisbon'\n")
        
        f.close()

    print "all crab submit scripts written"


def mergeCrabJobOutput(samples):

    if options.cores == -1:
        options.cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(options.cores)

    print "\n" + logo + "merging crab job output... \n"

    # merge all datasets individually:
    commands = []
    for sample in samples:

        print "Merging", sample

        logfile = samples[sample]["workArea"] + "/crab_" + samples[sample]["outputDatasetTagPrefix"] + "_" + sample + "/crab.log"
        process = subprocess.Popen(['grep', '-m 1', 'Submitting', logfile], stdout=subprocess.PIPE)
        out, err = process.communicate()
                                    
        info = ast.literal_eval(out.split("Submitting")[-1].strip())

        if "lfn" in info:
            outputfolder = options.folderprefix + info["lfn"] + samples[sample]["datasetpath"].split("/")[1] + "/" + info["workflow"].split("crab_")[-1]

            outputname = info["edmoutfiles"][0].split(".root")[0]

            # merge all output root files:
                       
            files = ast.literal_eval(samples[sample]["mergeOutputFiles"])
            for ifile in files:
                cmd = "hadd -f9 " + options.outputfolder + "/" + sample + ifile + " " + outputfolder + "/*/*/" + files[ifile]
                commands.append(cmd)
            
    # run actual merge commands in parallel mode:
    results = pool.map(ShellExec, commands)
    
    print "\n" + logo + "done! \n"


def mergeDataOutput(samples):

    # merge all data into allData.root:
    datalist = []
    files = ""
    for sample in samples:
        if samples[sample]["type"] == "data":
            datalist.append(sample)
            files = ast.literal_eval(samples[sample]["mergeOutputFiles"])
                
    commands = []
    for ifile in files:
        cmd = "hadd -k -f9 " + options.outputfolder + "/" + "allData%s" % ifile + " "
        for sample in datalist:
            cmd += options.outputfolder + "/" + sample + ifile + " "
        commands.append(cmd)

    if options.cores == -1:
        options.cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(options.cores)
    results = pool.map(ShellExec, commands)
    
    print "\n" + logo + "done! \n"


def getCrabStatus(samples):
    
    status = {}

    print "Time: ", strftime("%Y-%m-%d %H:%M:%S", gmtime()), ", skim configuration:", options.config, "\n"
    print logo + "searching crab jobs... \n"
    print "Sample                                | Status       | Finished           | Failed             | Running            "
    print "--------------------------------------+--------------+--------------------+--------------------+--------------------"

    jobStatusDict = {}

    for sample in samples:

        cmd = "source /cvmfs/cms.cern.ch/crab3/crab.sh; crab status " + samples[sample]["workArea"] + "/crab_" + samples[sample]["outputDatasetTagPrefix"] + "_" + sample
        
        process = subprocess.Popen(['sh', '-c', cmd], stdout=subprocess.PIPE)
        out, err = process.communicate()
        if options.debug: print out, err
        status[sample] = [out, err]

        try:
            taskstatus = status[sample][0].split("Task status:\t\t\t")[1].split("\n")[0]
        except: taskstatus = ""
            
        try:
            jobsstatus = status[sample][0].split("Jobs status:\t\t\t")[1].split("\n\nPublication status:")[0]
        except: jobsstatus = ""

        if "finished" in jobsstatus:
            jobsfinished = jobsstatus.split("finished      ")[1].split("\n")[0]
        else: jobsfinished = ""

        if "failed" in jobsstatus:
            jobsfailed = jobsstatus.split("failed        ")[1].split("\n")[0]
        else: jobsfailed = ""

        if "running" in jobsstatus:
            jobsrunning = jobsstatus.split("running       ")[1].split("\n")[0]
        else: jobsrunning = ""

        while len(sample)<38:       sample += " "
        while len(taskstatus)<13:   taskstatus += " "
        while len(jobsfinished)<19: jobsfinished += " "
        while len(jobsfailed)<19:   jobsfailed += " "
        while len(jobsrunning)<19:  jobsrunning += " "

        print sample[:38] + "| " + taskstatus[:13] + "| " + jobsfinished[:19] + "| " + jobsfailed[:19] + "| " + jobsrunning[:19]
        
        jobStatusDict[sample] = jobsstatus

    print "\nResubmit failed tasks with --resubmit option. When finished, merge all job output files with the --merge and optionally with the --outputfolder option.\n"
    print "For detailed job statuses use " + color.UNDERLINE + "http://dashb-cms-job.cern.ch/dashboard/templates/task-analysis/" + color.END + "\n"
    return jobStatusDict


def crabJobCommand(samples, crabcommand):

    for sample in samples:
        cmd = "source /cvmfs/cms.cern.ch/crab3/crab.sh; crab %s " % crabcommand + samples[sample]["workArea"] + "/crab_" + samples[sample]["outputDatasetTagPrefix"] + "_" + sample
        ShellExec(cmd)


def collectLumiJSON(samples):

    lumis = ""
    for sample in samples:    
        if samples[sample]["type"] == "data":
            jsonfile = (glob(samples[sample]["workArea"] + "/crab_" + samples[sample]["outputDatasetTagPrefix"] + "_" + sample + "/results/processedLumis.json") )
            if len(jsonfile)>0:
                with open(jsonfile[0], 'r') as f:
                    lumis += f.read().split("{")[-1].split("}")[0] + ","
       
    f = open('%s/allProcessedLumis.json' % options.outputfolder, 'w')
    f.write('{' + lumis + '}')
    f.close()


if __name__ == "__main__":

    parser = optparse.OptionParser()
    
    parser.add_option("--config", action="store", type="string", default="", dest="config")
    parser.add_option("--folderprefix", action="store", type="string", default="/gstore/t3cms/store/user/giles/", dest="folderprefix")
    parser.add_option("--outputfolder", action="store", type="string", default=".", dest="outputfolder")
    parser.add_option("--cores", action="store", type="int", default=-1, dest="cores")
    parser.add_option("--dryrun", action="store_true", dest="dryrun")
    parser.add_option("--submit", action="store_true", dest="submit")
    parser.add_option("--merge", action="store_true", dest="merge")
    parser.add_option("--mergeData", action="store_true", dest="mergeData")
    parser.add_option("--status", action="store_true", dest="status")
    parser.add_option("--resubmit", action="store_true", dest="resubmit")
    parser.add_option("--report", action="store_true", dest="report")
    parser.add_option("--createSubmitScripts", action="store_true", dest="createSubmitScripts")
    parser.add_option("--debug", action="store_true", dest="debug")
    (options, pfade) = parser.parse_args()

    sys.dont_write_bytecode = True

    if options.config == "":
        print "Submit jobs with: ./crabControl.py --config datasets.cfg --submit"
        quit()

    samples = readConfig(options.config)
    if options.debug:
        print "Using samples:", samples
    
    if options.submit:
        buildCrabSubmitScripts(samples)
        crabSubmit(samples)

    if options.merge:
        mergeCrabJobOutput(samples)
        mergeDataOutput(samples)

    if options.mergeData:
        mergeDataOutput(samples)

    if options.status:
        getCrabStatus(samples)

    if options.resubmit:
        crabJobCommand(samples, "resubmit")

    if options.report:
        crabJobCommand(samples, "report")
        collectLumiJSON(samples)

    if options.createSubmitScripts:
        buildCrabSubmitScripts(samples)
