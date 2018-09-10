import os
import subprocess
from datetime import datetime
from param import *


#-----------------------------------------------------------------------------#
dnaComplementTranslation = str.maketrans("ATGC", "TACG")
def revComp(seq):
    ''' provided by John Dicarlo '''
    seq = seq[::-1]
    return seq.translate(dnaComplementTranslation)
   

#-----------------------------------------------------------------------------#
def primer(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg
   
    # check if panel is already sorted
    nameout = panelId + ".primers.txt"
    if os.path.isfile(nameout):
        return
        
    # sort primers
    namein = panelId + ".txt"
    fileout = open(nameout, "w")
    header = ["gene.id", "gene.symbol", "gene.strand", "chrom", "loc5", "loc3", "strand", "primer", "target.loc", "rna.distance", "pd.id", "offTarget.genes", "notMappale"]
    fileout.write("\t".join(header) + "\n")
    
    # regular primers
    firstPrimer = True
    for line in open(namein, "r"):
        if line.startswith("#"):
            continue
            
        row = line.split("\n").split("\t")
        
        controlType = row[11]
        if controlType > 0:
            continue
        
        vec = list(row[1:8])
        loc = vec[5]
        primer = row[9]
        key = (primer, loc)
        targetLoc = str(row[12])
        targetDist = str(row[13])
        pdId = str(row[15])
        offTargetGenes = "" if row[10] == None else row[10]
        notMappable = "" if row[16] == None else row[16]
        
        if targetLoc in ["0", "None"]:
            targetLoc = "-1"
            targetDist = ""
    
        if firstPrimer:
            firstPrimer = False
            vecLast = vec
            primerLast = primer
            locLast = loc
            keyLast = (primerLast, locLast)
            pdIdLast = pdId
            offTargetGenesLast = offTargetGenes
            notMappableLast = notMappable
            targetLocLsat = targetLoc
            targetLocs = targetLoc
            targetDists = targetDist
        elif key == keyLast:
            if targetLoc == targetLocLsat:
                targetDists += distanceSeparator + targetDist
            else:
                targetLocs += targetSeparator + targetLoc
                targetDists += targetSeparator + targetDist
                targetLocLsat = targetLoc
        else:
            vecLast.extend([primerLast, targetLocs, targetDists, pdIdLast, offTargetGenesLast, notMappableLast])
            fileout.write("\t".join((str(x) for x in vecLast)) + "\n")
            
            targetLocs = targetLoc
            targetDists = targetDist
            
            vecLast = vec
            primerLast = primer
            locLast = loc
            keyLast = (primerLast, locLast)
            pdIdLast = pdId
            offTargetGenesLast = offTargetGenes
            notMappableLast = notMappable
            targetLocLsat = targetLoc
         
    vec.extend([primer, targetLocs, targetDists, pdId, offTargetGenes, notMappable])           # last one
    fileout.write("\t".join((str(x) for x in vec)) + "\n")      
    
    # control primers at bottom
    firstPrimer = True
    for line in open(namein, "r"):
        row = line.split("\n").split("\t")
        
        controlType = row[11]
        if controlType == 0:
            continue
        
        vec = list(row[1:8])
        loc = vec[5]
        primer = row[9]
        key = (primer, loc)
        targetLoc = str(row[12])
        targetDist = str(row[13])
        pdId = str(row[15])
        offTargetGenes = "" if row[10] == None else row[10]
        notMappable = "" if row[16] == None else row[16]
        
        if targetLoc in ["0", "None"]:
            targetLoc = "-1"
            targetDist = ""
    
        if firstPrimer:
            firstPrimer = False
            vecLast = vec
            primerLast = primer
            locLast = loc
            keyLast = (primerLast, locLast)
            pdIdLast = pdId
            offTargetGenesLast = offTargetGenes
            notMappableLast = notMappable
            targetLocLsat = targetLoc
            targetLocs = targetLoc
            targetDists = targetDist
        elif key == keyLast:
            if targetLoc == targetLocLsat:
                targetDists += distanceSeparator + targetDist
            else:
                targetLocs += targetSeparator + targetLoc
                targetDists += targetSeparator + targetDist
                targetLocLsat = targetLoc
        else:
            vecLast.extend([primerLast, targetLocs, targetDists, pdIdLast, offTargetGenesLast, notMappableLast])
            fileout.write("\t".join((str(x) for x in vecLast)) + "\n")
            
            targetLocs = targetLoc
            targetDists = targetDist
            
            vecLast = vec
            primerLast = primer
            locLast = loc
            keyLast = (primerLast, locLast)
            pdIdLast = pdId
            offTargetGenesLast = offTargetGenes
            notMappableLast = notMappable
            targetLocLsat = targetLoc
         
    vec.extend([primer, targetLocs, targetDists, pdId, offTargetGenes, notMappable])           # last one
    fileout.write("\t".join((str(x) for x in vec)) + "\n")      
    
    fileout.close()    

   
#-----------------------------------------------------------------------------#
def makePair(fastq, speSide):
    # get fastq name
    fastqName = fastq.split(".fastq")[0]
    
    # output files
    mtSide = "2" if speSide == "1" else "1"
    fastqSpe = open(fastqName + "_R" + speSide + ".fastq", "w")
    fastqUmi = open(fastqName + "_R" + mtSide + ".fastq", "w")
    
    # loop over input fastq
    count = 0
    for line in open(fastq, "r"):
        if count == 1:
            seq = line.strip("\n")
            fastqSpe.write(revComp(seq) + "\n")
            
            # STAR compatibility: max Lread of 650 (paired)
            seq_ = seq
            seqLen = len(seq)
            if seqLen > 325:
                delta = 650-seqLen
                seq_ = seq[:delta]
            fastqUmi.write(seq_ + "\n")
        elif count == 3:
            bq = line.strip("\n")
            fastqSpe.write(bq[::-1] + "\n") 
            
            # STAR compatibility: max Lread of 650 (paired)
            bq_ = bq
            if seqLen > 325:
                bq_ = bq[:delta]
            fastqUmi.write(bq_ + "\n")
        else:
            # just save to disk
            fastqSpe.write(line)
            fastqUmi.write(line)
    
        # increment counter
        count += 1
        if count == 4:
            count = 0
        
    # clean up
    fastqSpe.close()
    fastqUmi.close()
    
    # clean up
    os.remove(fastq)


#-----------------------------------------------------------------------------#
def prep(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg
    
    # read set related preprocessing
    readSet = runname + "." + samplename
    samplepath = runpath + "/" + runname
    if not os.path.exists(samplepath):
        samplepath = runpath    
    subprocess.check_call("mkdir -p " + readSet, shell = True)   
    
    if platform == "2":
        # cmd = "mv " + samplepath + "/" + samplename + "*.fastq " + readset + "/" + samplename + ".fastq"
        cmd = "cp " + samplepath + "/" + samplename + "*.fastq " + readset + "/" + samplename + ".fastq"
        subprocess(cmd, shell = True)
        makePair(readset + "/" + samplename + ".fastq", speSide)
    else:
        fastq1 = samplepath + "/" + samplename + "*_R1*.fastq"
        fastq2 = samplepath + "/" + samplename + "*_R2*.fastq"        
        # cmd1 = "if [ -e " + fastq1 + " ]; then mv " + fastq1 + " " + readSet + "/" + samplename + "_R1.fastq; else mv " + fastq1 + ".gz " + readSet + "/" + samplename + "_R1.fastq.gz; fi"
        # cmd2 = "if [ -e " + fastq2 + " ]; then mv " + fastq2 + " " + readSet + "/" + samplename + "_R2.fastq; else mv " + fastq2 + ".gz " + readSet + "/" + samplename + "_R2.fastq.gz; fi"            
        cmd1 = "if [ -e " + fastq1 + " ]; then cp " + fastq1 + " " + readSet + "/" + samplename + "_R1.fastq; else cp " + fastq1 + ".gz " + readSet + "/" + samplename + "_R1.fastq.gz; fi"
        cmd2 = "if [ -e " + fastq2 + " ]; then cp " + fastq2 + " " + readSet + "/" + samplename + "_R2.fastq; else cp " + fastq2 + ".gz " + readSet + "/" + samplename + "_R2.fastq.gz; fi"    
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)

    
#-----------------------------------------------------------------------------#
def trim(cfg, numCores):
    '''
    trim and filter fastq
    '''

    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg

    # read set
    readset = runname + "." + samplename
    
    # SPE side related parameters
    if speSide == "1":
        speSideAdapter = r1universal
        umiSideAdapter = r2universal
        primer3BasesR1 = "-1"
        primer3BasesR2 = "8"
    else:
        speSideAdapter = r2universal
        umiSideAdapter = r1universal
        primer3BasesR2 = "-1"
        primer3BasesR1 = "8"
    
    fastq1 = readset + "/" + samplename + "_R1.fastq"
    if not os.path.isfile(fastq1):
        fastq1 += ".gz"
    fastq2 = readset + "/" + samplename + "_R2.fastq"
    if not os.path.isfile(fastq2):
        fastq2 += ".gz"
    
    # run trimming command
    cmd = "python3 " + readTrimmer + " --seqtype rna --check_primer_side --no_tagnames --primer_col " + str(primerCol)
    cmd += " --ncpu " + str(numCores) +  " --tag_seperator \"" + readIdSeparator + "\" --min_primer_side_len 30 --min_umi_side_len 30 --poly_tail_primer_side polyA "
    cmd += " --r1 " + fastq1 + " --r2 " + fastq2 + " --out_r1 " + readset + "/trim.R1.fq --out_r2 " + readset + "/trim.R2.fq " + " --out_metrics " + readset + "/trim.read.txt"
    cmd += " --primer_file " + panelId + ".primers.txt" + " --primer3_bases_R1 " + primer3BasesR1 + " --primer3_bases_R2 " + primer3BasesR2    
    cmd += " --common_seq_len " + str(len(umiSideAdapter)) + " --umi_len " + str(umiLen) + " --umi_filter_max_lowQ_bases " + str(umiMaxLowqBases) + " --umi_filter_max_Ns " + str(umiMaxNs)
    if platform == "1":
        cmd += " --is_nextseq "        
    elif platform == "2":
        cmd += " --trim_custom_seq_adapter "
    if speSide == "1":
        if len(r1universal):
            cmd += " --custom_seq_adapter " + speSideAdapter
    else:
        if len(r2universal):
            cmd += " --custom_seq_adapter " + speSideAdapter
        cmd += " --is_r2_primer_side"
    cmd += " > " + readset + "/trim.read.log"
    subprocess.check_call(cmd, shell = True)
    
    # clean up
    subprocess.check_call("rm -f " + readset + "/*.fastq " + readset + "/*.fastq.gz", shell = True)


#-----------------------------------------------------------------------------#
def run(cfgs, numCores):
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " preparing reads")

    # file transfer and decompression
    for readset in sorted(cfgs.keys(), reverse = True):
        prep(cfgs[readset])
    #joblib.Parallel(n_jobs = numJobs)(joblib.delayed(trim.fastqPrep)(cfgs[readSet], isDebug) for readSet in cfgs)
    
    # trim and tag reads
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " trim and tag"    )
    for (readSet, cfg) in cfgs.items():
        # prepare primer file
        primer(cfg)
        
        # trim
        trim(cfg, numCores)
    
    # clean up
    subprocess.check_call("rm -f *.index.cache", shell = True)