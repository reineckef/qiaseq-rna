import os
import subprocess
import joblib
from datetime import datetime

from param import *
import paramStar

'''
STAR MAPQ (5th field in SAM):

255 = uniquely mapped reads
3 = read maps to 2 locations
2 = read maps to 3 locations
1 = reads maps to 4-9 locations
0 = reads maps to 10 or more locations
'''


#------------------------------------------------------------------------------#
# load/unload STAR index to/from memory (if required)
def loadStarIndex(loadIndex):
    # load
    if loadIndex:
        cmd = "mkdir -p tmp && STAR --genomeDir " + starIndex + " --outFileNamePrefix tmp/ --genomeLoad LoadAndExit > /dev/null"
        subprocess.check_call(cmd, shell = True)
    # unload
    else:
        try:
            cmd = "mkdir -p tmp && STAR --genomeDir " + starIndex + " --outFileNamePrefix tmp/ --genomeLoad Remove > /dev/null && rm -r tmp"
            subprocess.check_call(cmd, shell = True)
        except:
            print("Warning: unload STAR Index not successful, most likely in use by another run."    )
 

#------------------------------------------------------------------------------#
# RNA alignment to genome
def align(cfg, numCores):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename
    
    # get STAR parameters    
    cmdParam = ''
    for key, value in paramStar.__dict__.items():
        if not key.startswith("__"):
            cmdParam += "--" + key + " " + str(value) + " "
    
    # STAR alignment
    cmd = "STAR " + cmdParam + "--runThreadN " + str(numCores) + " --genomeDir " + starIndex + " --readFilesIn " + readset + "/trim.R1.fq " + readset + "/trim.R2.fq " + "--outFileNamePrefix " + readset + "/align. > " + readset + "/align.star.log 2>&1"
    subprocess.check_call(cmd, shell = True)
    
    subprocess.check_call("mv " + readset + "/align.Aligned.out.bam " + readset + "/align.star.bam", shell = True)
    


#------------------------------------------------------------------------------#
# sort BAM files
def sortBam(cfg, numCores):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools sort -o " + readset + "/" + readset + ".bam -T " + readset + "/tmp.bam " + "-@ " + str(numCores) + " " + readset + "/align.star.bam > /dev/null 2>&1"
    subprocess.check_call(cmd, shell = True)


#------------------------------------------------------------------------------#
# index BAM files
def indexBam(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools index " + readset + "/" + readset + ".bam"
    subprocess.check_call(cmd, shell = True)
    
    
#------------------------------------------------------------------------------#
# convert BAM to SAM
def bamToSam(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools view -h -F 256 " + readset + "/align.star.bam > " + readset + "/align.star.sam"
    subprocess.check_call(cmd, shell = True)    
    
    subprocess.check_call("rm " + readset + "/align.star.bam", shell = True)
    
    
#------------------------------------------------------------------------------#
# run all functions in this module
def run(cfgs, numCores):
    numJobs = min(len(cfgs), numCores)

    # load STAR index
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " load STAR index")
    loadStarIndex(True)
    
    # align RNA reads
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " alignment")
    for (readset, cfg) in cfgs.items():
        align(cfg, numCores)
        
    # unload STAR index
    loadStarIndex(False) 
    
    # sort BAM
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " BAM/SAM processing")
    for (readset, cfg) in cfgs.items():
        sortBam(cfg, numCores)
    
    # index BAM
    joblib.Parallel(n_jobs = numJobs)(joblib.delayed(indexBam)(cfgs[readset]) for readset in cfgs)
    
    # convert BAM to SAM
    joblib.Parallel(n_jobs = numJobs)(joblib.delayed(bamToSam)(cfgs[readset]) for readset in cfgs)