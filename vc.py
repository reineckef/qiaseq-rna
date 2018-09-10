import subprocess
import joblib
from datetime import datetime

import bed
from param import *


#-----------------------------------------------------------------------------#
dnaComplementTranslation = str.maketrans("ATGC", "TACG")
def revComp(seq):
    ''' provided by John Dicarlo '''
    seq = seq[::-1]
    return seq.translate(dnaComplementTranslation)
    
    
#-----------------------------------------------------------------------------#
def ionfq(cfgs, speSide):
    '''
    function to convert trimmed fastq to a correct format, only for Ion readset
    '''

    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg
        
    if platform != "2":
        return            
    
    readset = runname + "." + samplename
    
    fileout = open(readset + "/" + readset + "trim.fq", "w")
    if speSide == "1":
        filein = open(readset + "/trim.R1.fq", "r")        
    else:
        filein = open(readset + "/trim.R2.fq", "r")
    
    count = 0
    for line in filein:
        if count == 0 or count == 2:
            fileout.write(line)
        elif count == 1:
            fileout.write(revComp(line.strip("\n")) + "\n")
        else:
            fileout.write(line.strip("\n")[::-1] + "\n")
            
        count += 1
        if count == 4:
            count = 0
    
    fileout.close()


#------------------------------------------------------------------------------#
def qimera(cfgs, numCores):
    pass


#------------------------------------------------------------------------------#
# convert BAM to SAM
def samToBam(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools view -1 -hb  " + readset + "/qc.clip.sam > " + readset + "/qc.clip.bam 2> " + readset + "/vc.samtobam.log"
    subprocess.check_call(cmd, shell = True)    
    
    subprocess.check_call("rm " + readset + "/qc.clip.sam", shell = True)
    
    
#------------------------------------------------------------------------------#
# sort BAM files
def sortBam(cfg, numCores):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools sort -o " + readset + "/vc.clip.bam -T " + readset + "/tmp.bam " + "-@ " + str(numCores) + " " + readset + "/qc.clip.bam > /dev/null 2>&1"
    subprocess.check_call(cmd, shell = True)

    subprocess.check_call("rm " + readset + "/qc.clip.bam", shell = True)


#------------------------------------------------------------------------------#
# index BAM files
def indexBam(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename

    cmd = "samtools index " + readset + "/vc.clip.bam"
    subprocess.check_call(cmd, shell = True)
    

#------------------------------------------------------------------------------#
def geneCov(gene, genePrimers, fragLen):
    '''
    get gene coverage using exon models and a max fragment length
    based on step01.py by John Dicarlo
    '''
    
    # intit bed coverage
    bedCovOneGene = []
    bedTrackSet = set()
    bedWarnings = []
    
    # loop over RNAs
    for rnaId in gene:
        rnaLen = 0
        bedExons = []
        
        # get exons
        firstExon = True
        for (geneName, strand, chrom, exonStart, exonEnd) in gene[rnaId]:
            exonStart = int(exonStart)
            exonEnd = int(exonEnd)
            
            if firstExon:
                geneLocL = exonStart
                geneLocR = exonEnd
                firstExon = False
            else:
                geneLocL = min(exonStart, geneLocL)
                geneLocR = max(exonEnd, geneLocR)
            
            rnaLen += exonEnd - exonStart
            bedExons.append((chrom, exonStart, exonEnd))
        
        bedExons.sort()
        
        # init coverage for this RNA      
        bedCovOneRna = []
        
        # loop over primers, make RNA coverage BED tracks
        for (chrom, locDna5, locDna3, strand, primer) in genePrimers:
            locDna5 = int(locDna5)
            locDna3 = int(locDna3)
            strand = int(strand)
            
            # get primer RNA loc3
            exonsLen =  0
            locRna3 = None
            for (chrom, locL, locR) in bedExons:
                if locL <= locDna3 < locR:
                    locRna3 = exonsLen + locDna3 - locL
                    break
                exonsLen += (locR - locL)
            # check if primer match RNA
            if locRna3 == None:
                continue
            
            # get loc5 on RNA
            primerLen = len(primer)
            if strand == 0:
                locRna5 = locRna3 - primerLen + 1
            else:
                locRna5 = locRna3 + primerLen - 1
                
            # get DNA position of end of fragment
            if strand == 0:
                locRnaEnd = min(locRna5 + fragLen - 1, rnaLen - 1)
            else:
                locRnaEnd = max(locRna5 - fragLen + 1, 0)
            locL_ =  0
            locR_ =  0
            locDnaEnd = None
            for (chrom, locL, locR) in bedExons:
                locR_ += (locR - locL)
                if locL_ <= locRnaEnd < locR_:
                    locDnaEnd = locL + locRnaEnd - locL_
                    break
                locL_ = locR_
            if locDnaEnd == None:
                raise Exception()
                
            # frag coverage region
            if strand == 0:
                bedDelete = [(chrom, geneLocL, locDna3+1),(chrom, locDnaEnd+1, geneLocR)]
            else:
                bedDelete = [(chrom, geneLocL, locDnaEnd),(chrom, locDna3    , geneLocR)]
            bedCov = bed.subtract(bedExons, bedDelete)
            
            # save coverage across whole RNA
            bedCovOneRna.extend(bedCov)
        
            # make subtraction bed for full frag, including primer
            if strand == 0:
                bedDelete = [(chrom, geneLocL, locDna5  ),(chrom, locDnaEnd+1, geneLocR)]
            else:
                bedDelete = [(chrom, geneLocL, locDnaEnd),(chrom, locDna5 + 1, geneLocR)]
                
            # do bed subtraction to get enrichment frag
            bedFrag = bed.subtract(bedExons, bedDelete)
            bedFrag = bed.merge(bedFrag) # should not do anything
            
            # get size of enrichment frag (might be less than fragLen at ends of RNA)
            bpFrag = sum((x[2]-x[1] for x in bedFrag))
    
            # convert bedFrag to a one-row bed
            bedLocL = bedFrag[0][1]
            bedLocR = bedFrag[-1][2]
            if strand == 0:
                bedStrand = "+"
                bedThickStart = locDna3 + 1
                bedThickStop  = bedLocR
            else:
                bedStrand = "-"
                bedThickStart = bedLocL
                bedThickStop  = locDna3
            if bedThickStart >= bedThickStop:
                bedWarnings.append((chrom, locDna5, locDna3, strand, primer, geneName, rnaId, bedThickStart, bedThickStop))
            numBlocks = len(bedFrag)
            blockSizes  = ",".join([str(x[2]-x[1]   ) for x in bedFrag])
            blockStarts = ",".join([str(x[1]-bedLocL) for x in bedFrag])
            bedScore = 0
            bedOne = (chrom, bedLocL, bedLocR, geneName, bedScore, bedStrand, bedThickStart, bedThickStop, 0, numBlocks, blockSizes, blockStarts)
            # bedOne = (chrom, bedLocL, bedLocR, bpFrag, bedScore, bedStrand, bedThickStart, bedThickStop, 0, numBlocks, blockSizes, blockStarts)
            bedTrackSet.add(bedOne)
        
        # update BED for all RNAs coverage
        bedCovOneGene.extend(bedCovOneRna) 
            
    # post processing
    bedCovOneGene.sort()
    bedTrackSet = list(bedTrackSet)
    bedTrackSet.sort()

    return (bedCovOneGene, bedTrackSet, bedWarnings)


#------------------------------------------------------------------------------#
# run variant calling, smCounter V2
def roiBed(cfg):
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename
    
    inPrimers = panelId + ".primers.txt"
    inBam = readset + "/vc.clip.bam"

    #---------- E: get gencode basic exons (GRCh38)
    genes = {}
    firstLine = True
    for line in open(exonModels, "r"):
        # skip header
        if firstLine:
            firstLine = False
            continue
        
        (geneId, rnaId, geneName, strand, chrom, codingStart, codingEnd, exonStart, exonEnd) = line.strip().split("\t")
        
        if chrom == "MT":
            chrom = "M"
        if not "CHR" in chrom.upper():
            chrom = "chr" + chrom
        
        exonStart = str(int(exonStart) - 1)
        val = (geneName, strand, chrom, exonStart, exonEnd)
        
        if geneId in genes:
            if rnaId in genes[geneId]:
                (genes[geneId][rnaId]).append(val)
            else:
                genes[geneId][rnaId] = [val]
        else:
            genes[geneId] = {}
            genes[geneId][rnaId] = [val]    
        
    #---------- P: get primers (GRCh38)
    primers = {}
    firstLine = True
    for line in open(inPrimers, "r"):
        # skip header
        if firstLine:
            firstLine = False
            continue
        
        vals = line.strip("\n").split("\t")
        (geneId, geneName, geneStrand, chrom, loc5, loc3, strand, primer) = vals[:8]
        
        if geneId == "":
            continue    
        val = (chrom, loc5, loc3, strand, primer)
        
        if geneId in primers:
            (primers[geneId]).append(val)
        else:
            primers[geneId] = [val]        
    
    #---------- EP: overlap of exons and primers (GRCh38)
    fileout1 = open(readset + "/vc." + panelId + ".exons.bed", "w")
    fileout2 = open(readset + "/vc." + panelId + ".exons.track.bed", "w")
    fileoutwarn = open(readset + "/vc.gencov.log", "w")
    
    # loop over genes
    for geneId in genes:
        # check if any primer covers this gene
        if not geneId in primers:
            continue
            
        # get info of one gene
        gene = genes[geneId]
        for rnaId in gene:
            geneName = gene[rnaId][0][0]
            break
        
        # get primers info on the gene
        (bedCovOneGene, bedTrackSet, bedWarnings) = geneCov(gene, primers[geneId], vcFragLenMax)
    
        # write primer/RNA tracks for this gene         
        for vec in bedCovOneGene:
            outVec = [str(x) for x in vec]
            outVec.append(geneName)
            fileout1.write("\t".join(outVec) + "\n") 
    
        for vec in bedTrackSet:
            fileout2.write("\t".join((str(x) for x in vec)) + "\n")
            
        # write warnings
        for item in bedWarnings:
            fileoutwarn.write("Warning: " + "\t".join([str(x) for x in item]) + "\n")
    fileout1.close()
    fileout2.close()
    fileoutwarn.close()
    
    cmd = "bedtools sort -i " + readset + "/vc." + panelId + ".exons.bed | bedtools merge -c 4 -o distinct -i - > " + readset + "/vc." + panelId + ".exons.tmp"
    cmd += " && mv " + readset + "/vc." + panelId + ".exons.tmp " + readset + "/vc." + panelId + ".exons.bed"
    subprocess.check_call(cmd, shell = True)    
    
    #---------- EPB: overlap of BAM and EP (GRCh38)
    # get BAM coverage
    bamFileName = inBam.split("/")[-1].split(".bam")[:-1][0]
    cmd = "bedtools genomecov -bg -split -ibam " + inBam + " -g " + refGenome + " > " + readset + "/" + bamFileName + ".gencov.bed"
    subprocess.check_call(cmd, shell = True)
    
    # filter low-coverage regions
    fileout = open(readset + "/" + bamFileName + ".gencov.filtered.bed", "w")
    for line in open(readset + "/" + bamFileName + ".gencov.bed", "r"):
        (chrom, locL, locR, cov) = line.strip().split("\t")
        if float(cov) >= vcDepthMin:
            fileout.write(line)
    fileout.close()
    
    cmd = "bedtools sort -i " + readset + "/" + bamFileName + ".gencov.filtered.bed" + " | bedtools merge -i - > " + readset + "/" + bamFileName + ".gencov.filtered.tmp"
    cmd += " && mv " + readset + "/" + bamFileName + ".gencov.filtered.tmp " + readset + "/" + bamFileName + ".gencov.filtered.bed"
    subprocess.check_call(cmd, shell = True)
    
    # overlap BAM coverage with primer-exon BED
    cmd = "bedtools intersect -a " + readset + "/" + bamFileName + ".gencov.filtered.bed -b " + readset + "/vc." + panelId + ".exons.bed" + " > " + readset + "/" + bamFileName + ".gencov.filtered." + panelId + ".exons.bed"
    subprocess.check_call(cmd, shell = True)
    
    cmd = "bedtools sort -i " + readset + "/" + bamFileName + ".gencov.filtered." + panelId + ".exons.bed" + " | bedtools merge -i - > " + readset + "/" + bamFileName + ".gencov.filtered." + panelId + ".exons.tmp"
    cmd += " && mv " + readset + "/" + bamFileName + ".gencov.filtered." + panelId + ".exons.tmp " + readset + "/" + bamFileName + ".gencov.filtered." + panelId + ".exons.bed"
    subprocess.check_call(cmd, shell = True)
    

#------------------------------------------------------------------------------#
def vc(cfg, numCores):
    ''' 
    run variant calling, smCounter V2
    '''
    
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg
    
    # read set
    readset = runname + "." + samplename

    cmd = "python2 " + smCounter + " --isRna --bkgErrorDistSimulation " + bkgErrorDistSimulation + " --runPath " + readset + " --bamFile vc.clip.bam --bedTarget vc.clip.gencov.filtered." + panelId + ".exons.bed --outPrefix " + readset + " --nCPU " + str(numCores) + " --minBQ 25 --minMQ 50 --hpLen 8 --mismatchThr 6.0 --primerDist 2 --mtThreshold 0.8 --primerSide " + speSide + " --refGenome " + refGenome + " --repBed " + repBed + " --srBed " + srBed + " > " + readset + "/vc.snpindel.log 2>&1"
    subprocess.check_call(cmd, shell = True)


#------------------------------------------------------------------------------#
def run(cfgs, numCores, runVc, runFc):
    # fusion calling
    if runFc:
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " fusion calling")
        joblib.Parallel(n_jobs = numCores)(joblib.delayed(ionfq)(cfgs[readset]) for readset in cfgs)
        qimera(cfgs, numCores)
    else:
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " fusion calling skipped")

    # variant calling
    if runVc:
        # BAM/SAM processing
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " SAM/BAM processing")
        joblib.Parallel(n_jobs = numCores)(joblib.delayed(samToBam)(cfgs[readset]) for readset in cfgs)    
        for cfg in cfgs.values():
            sortBam(cfg, numCores)
        joblib.Parallel(n_jobs = numCores)(joblib.delayed(indexBam)(cfgs[readset]) for readset in cfgs)

        # make BED of exons (for variant calling)
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " ROI BED")
        joblib.Parallel(n_jobs = numCores)(joblib.delayed(roiBed)(cfgs[readset]) for readset in cfgs)

        # run snp/indel calling
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " snp/indel calling")
        for cfg in cfgs.values():
            vc(cfg, numCores)
    else:
        print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " snp/indel calling skipped")
