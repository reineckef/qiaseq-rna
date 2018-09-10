'''
main script to run Qiagen speRNA pipeline
creates detailed report on read and tag counts, primer performance, and fragment lengths
M Nezami Ranjbar, Qiagen, Aug 2018
'''

__version__ = 2.0
__author__  = 'mohammad.nezami@qiagen.com'


#-----------------------------------------------------------------------------#
import sys
import argparse
import joblib
from datetime import datetime

# local
import trim
import align
import qc
import vc
from param import *


#-----------------------------------------------------------------------------#
def main(args):
    # logo
    print("\n" + sectionSeparator)
    print("\t\t" + "QIAGEN speRna")

    # initialize
    numCores = min(args.numCores, joblib.cpu_count())

    #---------- parse run list
    cfgs = {}

    for line in open(args.runList, "r"):
        # skip header
        if line.startswith("#"):
            continue
        
        # get values    
        vals = line.strip().split("\t")
        (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = vals[:10]
        
        # read set
        readset = runname + "." + samplename
        
        # config parameter
        cfgs[readset] = (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform)

    #---------- run modules one at a time
    # trimming
    print(sectionSeparator)
    trim.run(cfgs, numCores)

    # alignment and post processing       
    print(sectionSeparator)
    align.run(cfgs, numCores)

    # primer qc and filtering 
    print(sectionSeparator)
    qc.run(cfgs, numCores) 

    # varinat calling
    print(sectionSeparator)
    vc.run(cfgs, numCores, args.vc, args.fc)

    # clean-up
    print(sectionSeparator + "\n")
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " completed.\n")


#-----------------------------------------------------------------------------#
def parseArg(flag):     
    parser = argparse.ArgumentParser(description = "RNAscan 4-in-1 pipeline")
    parser.add_argument("--runList", default=None, help="list of readsets with annotations, file", required=True)
    parser.add_argument("--numCores", type=int, default=1, help="number of CPU threads to use in parallel")
    parser.add_argument("--ge", action = "store_true", help="run gene expression")
    parser.add_argument("--vc", action = "store_true", help="run snp/indel calling")
    parser.add_argument("--fc", action = "store_true", help="run fusion calling")

    if flag:
        return parser.parse_args()
    else:
        parser.print_help()        
    
   
#-----------------------------------------------------------------------------------
if __name__ == "__main__":
    # check if any arguments provided:
    if len(sys.argv) == 1:
        parseArg(False)
        
    # parse arguments and run main function
    args = parseArg(True)
    
    # run pipeline
    main(args)