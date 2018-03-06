#!/usr/bin/env python
import sys
import os, subprocess
from glob import glob
import optparse
import multiprocessing

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


if __name__ == "__main__":

    parser = optparse.OptionParser()
    
    parser.add_option("--cores", action="store", type="int", default=-1, dest="cores")
    parser.add_option("--edmfolder", action="store", type="string", default="/eos/user/k/kandroso/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_miniGen/", dest="edmfolder")
    parser.add_option("--fwoutputfolder", action="store", type="string", default="/eos/user/g/gstrong/cms/n-tuples/hh_v3/", dest="fwoutputfolder")
    parser.add_option("--debug", action="store_true", dest="debug")
    parser.add_option("--dryrun", action="store_true", dest="dryrun")
    (options, pfade) = parser.parse_args()

    sys.dont_write_bytecode = True
    
    if options.edmfolder == "":
        print "run analysis with: ./runItalianSkimAnalyser.py --edmfolder /path/to/edm-files"
        print "uses all available cpu cores, control with --cores option."
        quit()
    if options.fwoutputfolder == "":
        options.fwoutputfolder = options.edmfolder + "/fwoutput"
        
    print "Output files will be in", options.fwoutputfolder
    
    if not os.path.exists(options.fwoutputfolder):
        os.makedirs(options.fwoutputfolder)
    
    if options.cores == -1:
        options.cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(options.cores)
        
    if not os.path.exists(options.fwoutputfolder):
        os.makedirs(options.fwoutputfolder)
    
    commands = []
    
    edmfiles = glob (options.edmfolder + "/*.root")
    for edmfile in edmfiles:
        
        # ignore inverted selection
        #if "Inv" in edmfile: continue
        #if "selection" in edmfile: continue
       
        # ignore non-edm files
        if "monitor" in edmfile: continue

        #if 'Tau' not in edmfile: continue
        
        outputfile = options.fwoutputfolder + "/" + edmfile.split("/")[-1]
        
        more = ""
        if options.debug:
            more = "debug=True"
        
        commands.append( "$(echo $CMSSW_BASE)/bin/$(echo $SCRAM_ARCH)/italian_skim_analyser inputFiles=%s outputFile=%s monitorFile=none %s" % (edmfile, outputfile, more ) )
       
    results = pool.map(ShellExec, commands)
