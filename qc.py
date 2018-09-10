import sys
import subprocess
import joblib
import multiprocessing
from collections import defaultdict
from datetime import datetime

from param import *


# SAM flags
sam_flag_read_is_paired = 1
sam_flag_mapped_in_proper_pair = 2
sam_flag_read_is_unmapped = 4
sam_flag_is_reverse = 16
sam_flag_is_read1 = 64
sam_flag_is_read2 = 128
sam_flag_not_primary = 256
sam_flag_supplementary = 2048

# SAM columns
sam_qname = 0
sam_flag  = 1
sam_rname = 2
sam_pos   = 3
sam_mapq  = 4
sam_cigar = 5
sam_rname_mate = 6
sam_pos_mate = 7
sam_tlen  = 8
sam_seq   = 9
sam_bq    = 10

# extra columns (would be translated to SAM tags at the end) considering first 11 standard SAM columns, index from end
sam_tag_umi = -7
sam_tag_primer_id = -6
sam_tag_primer_info = -5
sam_tag_revalign = -4
sam_tag_cigar = -3
sam_tag_aend = -2
sam_tag_ontarget = -1

# counter columns
index_read_supplementary = 0
index_not_proper_pair = 1
index_spe_not_mapped = 2
index_not_same_chrom = 3
index_not_same_loc = 4
index_reverse_align = 5
index_low_mapq = 6
index_pass_filter = 7
index_on_target = 8
index_multi_align = 9
numCountFields = 10

# cigar op to char
cigarOpChar = ["M", "I", "D", "N", "S", "H"]

# segments for length binning
segmentLength, fragmentLengthMin, fragmentLengthMax = 25, 25, 250    
numSegments = (fragmentLengthMax-fragmentLengthMin) // segmentLength + 1
segmentsIdx = range(numSegments)

# number of columns for each off-target site (loc, reads, percentage)
offSiteFields = 3
    
    
#------------------------------------------------------------------------------#
def getPrimers(panelId):
    #initialize
    primers = []
    
    for line in open(panelId + ".primers.txt", "r"):
        # skip header
        if line.startswith("#"):
            continue
        
        # get values 
        (geneId, geneSymbol, geneStrand, chrom, loc5, loc3, strand, primer, targets, distances) = line.strip("\n").split("\t")[0:10]
        
        # no target site and distance provided
        if len(targets) == 0 or targets == "None":
            targets = "-1"
            distances = -1
        elif len(distances) == 0:
            if not targets == "-1":
                distances = -1
        
        # update primer dictionary
        primers.append((geneId, geneSymbol, geneStrand, chrom, loc5, loc3, strand, primer, targets, distances))        
    
    return primers


#------------------------------------------------------------------------------#
def getSuppTagAlignment(read):
    suppAlignments = []
    for tag in read[11:]:
        if "SA:Z:" in tag:
            # get all supplementary alignments info
            suppAligns = tag.split(":")[-1]
            
            for suppAlign in suppAligns.strip(";").split(";"):
                chrom, loc, strand, cigar, mapq, nm = suppAlign.split(",")
                suppAlignments.append((chrom, int(loc), cigar))
    
    return suppAlignments


#------------------------------------------------------------------------------#
def decode(cigar):
    '''
    functions to calculate coverage and align lengths based on CIGAR
    '''
    # initialize
    idx = -1
    out = []
    
    # validate
    if not cigar.isalnum():
        return out
    
    # calculate fragment and gap lengths
    for i in range(len(cigar)):
        x = cigar[i]
        
        # if a letter
        if x.isalpha():
            # if a match
            if x == "M":
                out.append((0, int(cigar[idx+1:i])))
            # if an insertion
            elif x == "I":
                out.append((1, int(cigar[idx+1:i])))
            # if a deletion
            elif x == "D":
                out.append((2, int(cigar[idx+1:i])))
            # if a gap    
            elif cigar[i] == "N":
                out.append((3, int(cigar[idx+1:i])))
            # if a soft clip
            elif x == "S":
                out.append((4, int(cigar[idx+1:i])))
            # if a hard clip
            elif x == "H":
                out.append((5, int(cigar[idx+1:i])))
    
            # update index
            idx = i            
    
    return out
    

#------------------------------------------------------------------------------#
def coverage(roiLen, matchBp, jumpBp):
    '''
    calculate coverage based on coordinates from alignment
    '''
    # initialize
    coverage = 0
    idx = 0
    numMatch = len(matchBp)
    numJump = len(jumpBp)
    alnLen = sum(matchBp)
    gapLen = sum(jumpBp)

    for n in range(numMatch):
        coverage += matchBp[n]
                
        if coverage > roiLen:
            idx = -(n+1)
            break
        elif n < numJump:
            coverage += jumpBp[n]
            if coverage > roiLen:
                idx = n+1
                break
            
    # target lies in middle of a gap, consider sum of following matches
    if idx > 0:
        coverage = 0
        for k in range(idx, numMatch):
            coverage += matchBp[k]
    # probably not covered, check max possible coverage to ensure
    elif idx == 0:
            coverage = alnLen + gapLen - roiLen
    # target lies in middle of a match, sum the match and following ones
    else:
        coverage -= roiLen
        for k in range(-idx, numMatch):
            coverage += matchBp[k]
    
    return coverage 


#------------------------------------------------------------------------------#
def check(read):
    '''
    check if read is on and calculate coverage of target based on RNA and DNA coordinates
    '''
    # initialize for read length and other parameter calculations
    readId = read[sam_qname]
    alnChrom = read[sam_rname]
    alnLoc = read[sam_pos]
    primerIds = read[sam_tag_primer_id]
    
    # initialize for length calculations (when primary alignment)
    rnaBp = 0    
    dnaBp = []
    alnPrimLoc = ""
    readType = 0    
    
    # check if CIGAR is OK
    cigar = read[sam_cigar]
    if not cigar.isalnum():
        return primerIds[0], readType, (alnChrom, alnLoc, alnPrimLoc, rnaBp, dnaBp)    

    # get supplementary alignment info
    suppAlignments = getSuppTagAlignment(read[:-7])
    
    # get primer info
    primerInfo = read[sam_tag_primer_info]
    
    # loop over possible primer indexes matching the same primer sequence
    alnLocOnTarg = None
    alnPrimLocOnTarg = None
    alnChromOnTarg = None
    
    for primerIdx in range(len(primerIds)):
        primerId = primerIds[primerIdx]
        
        # get primer info for each primer index
        (chrom, loc5, loc3, strand, primerLen, target, distance) = primerInfo[primerIdx]
            
        # check if any supplementary alignments
        alnLocs = [alnLoc]
        alnChroms = [alnChrom]
        cigars = [cigar]

        for (alnChrom, alnLoc, cigar) in suppAlignments:
            alnChroms.append(alnChrom)
            alnLocs.append(alnLoc)
            cigars.append(cigar)
                            
        # if different locus than primer, mark as off-target, skip to next primer
        if chrom not in alnChroms:
            continue
        
        for idx in range(len(alnLocs)):
            # calculate fragment and gap lengths
            cigarVec = decode(cigars[idx])
            matchBp = [bp for (op, bp) in cigarVec if op == 0]
            jumpBp = [bp for (op, bp) in cigarVec if op == 3]
            numMatch = len(matchBp)
            numJump = len(jumpBp)
            alnLen = sum(matchBp)
            gapLen = sum(jumpBp)

            # adjust alignment positions and alignment lengths based on direction
            alnLoc = alnLocs[idx]
            if strand == 0:
                alnPrimLoc = alnLoc
                alnLoc3Bases = loc3 - alnPrimLoc
            else:
                alnPrimLoc = alnLoc + alnLen + gapLen
                alnLoc3Bases = alnPrimLoc - loc3
            
            # check if primer loc3 on this segment
            if numMatch > 1 and not (minDistPos < alnLoc3Bases < primerLen + minDistPos):
                # check next segment, possible over-junction primer
                # first see if there is gap
                jumpBpL = 0
                jumpBpR = 0
                if numJump > 0:
                    jumpBpL = jumpBp[0]
                    jumpBpR = jumpBp[-1]
                
                # if primer aligned to next segment    
                if strand == 0:
                    alnPrimLoc += matchBp[0] + jumpBpL
                    alnLoc3Bases = loc3 - alnPrimLoc + matchBp[0]
                else:
                    alnPrimLoc -= (matchBp[-1] + jumpBpR)
                    alnLoc3Bases = alnPrimLoc - loc3 + matchBp[-1]
                
            # update totals
            rnaBp += alnLen
            
            # check if on-target
            if minDistPos < alnLoc3Bases < primerLen + minDistPos and alnLen - alnLoc3Bases > followBasesMin:
                # mark as on-target
                readType = 1

                # primer, chrom and location
                primerIdOn = primerId
                alnLocOnTarg = alnLoc
                alnPrimLocOnTarg = alnPrimLoc
                alnChromOnTarg = alnChroms[idx]
                loc3on = loc3
                
                # adjust direction
                if strand == 1:
                    matchBp = matchBp[::-1]
                    jumpBp = jumpBp[::-1] 
                
                # length of last segment
                alnLenLast = matchBp[-1]                
        
                # check if target is provided
                if target == "-1":
                    continue
                    
                targets = target.split(targetSeparator)
                numTargets = len(targets)
            
                # calculate coverage for each target site
                for targetIdx in range(numTargets):
                    # check if target is valid
                    if not targets[targetIdx].isdigit():
                        continue

                    roiLen = int(targets[targetIdx]) - alnPrimLoc     # alnPrimLoc instead of loc5 because of starting softclips usually caused by over-junction primers
                    if strand == 1:
                        roiLen = -roiLen
                        
                    dnaBp.append(coverage(roiLen, matchBp, jumpBp))           

    if readType:
        primerId = primerIdOn
        loc3 = loc3on
    else:
        primerId = primerIds[0]    
    
    return primerId, readType, loc3, (alnChrom, alnLoc, alnPrimLoc, rnaBp, dnaBp)


#------------------------------------------------------------------------------#
def clip(qname, pos, cigar, alignmentReoriented, basesGenomeNeeded):
    '''
    clip oligo, i.e. primer, UMI, adapter from a single read
    
    based on a work from John Dicarlo for DNA reads
    '''
    
    # init return flag for 100% softclip
    keepFlag = True
    
    # use CIGAR to get number of read bases in the primer binding region that are now aligned to genome
    cigarNew = cigar
    if alignmentReoriented:
        cigar.reverse()
    readLen = sum((x[1] for x in cigar if x[0] in [0, 1, 4]))
    basesGenome = 0
    basesRead = 0
    basesGap = 0
    for idx in range(0, len(cigar)):
        (op, bases) = cigar[idx]
        if op in [4, 5]:    # soft/hard clip - only the read - ignore to get bases of the read aligned
            continue
        elif op == 3:       # gap containing N's in CIGAR, intronic regions in RNA alignments
            basesGenome += bases
            basesGap += bases
        elif op == 1:       # insertion to reference - only the read
            basesRead += bases
        elif op == 0:       # genome and read in tandem
            bases = min(bases, basesGenomeNeeded - basesGenome) # don't go past the primer 3' binding design site
            basesRead   += bases
            basesGenome += bases
        elif op == 2:       # deletion from reference - only the genome
            basesGenome += bases
        else:
            raise Exception("unexpected CIGAR code")
        if basesGenome >= basesGenomeNeeded:
            break
            
    # debug check
    if basesRead < 0:
        raise Exception("bad CIGAR accounting!")
        
    # if last op was a deletion, move the read clipping one base forward - this could be a misaligned primer base - also forces wipeout of deletion in S walk below
    if op == 2 and (idx + 1) < len(cigar):
        (op, bases) = cigar[idx+1]
        if op == 0:
            basesGenome += 1
            basesRead   += 1
        elif op == 1:
            basesRead   += 1
    
    # if first is hard clip:
    (op, bases) = cigar[0]
    if op == 5:
        cigar = cigar[1:]
    
    # if read bases to clip, modify the cigar
    if basesRead > 0:
        # increase the soft clip amount, if already some soft clipping
        (op, bases) = cigar[0]
        if op == 4:
            cigar[0] = (op, basesRead+bases)
        # otherwise push a new code onto the cigar
        else:
            cigar.insert(0,(4, basesRead))
        
        # walk the cigar, subtract the M,D,I bases overwritten by the new S bases
        basesToSubtract = basesRead
        for idx in range(1, len(cigar)):
            (op, bases) = cigar[idx]
            if op == 2:  # deletion from reference does not help with cigar accounting
                cigar[idx] = (op, 0)
            elif op == 3:
                continue
            elif bases <= basesToSubtract:
                cigar[idx] = (op, 0)
                basesToSubtract -= bases
            else:
                basesNew = bases - basesToSubtract
                cigar[idx] = (op, basesNew)
                basesToSubtract = 0
            if basesToSubtract == 0:
                break
    
        # debug check
        if basesToSubtract > 0:
            raise Exception("ERROR: entire read soft clipped!", qname)
            
        # remove any zero-base ops
        cigarNew = []
        for val in cigar:
            (op, bases) = val
            if bases > 0:
                cigarNew.append(val)
                
        # collapse same adjacent codes (to handle example such as 24S105S and 20S1000N)
        basesGapRemoved = 0
        cigar = cigarNew
        cigarNew = []
        cigarNew.append(cigar[0])
        for idx in range(1,len(cigar)):
            val = cigar[idx]
            (op,  bases)  = val
            (op_, bases_) = cigarNew[-1]
            if op == 3 and op_ == 4:
                basesGapRemoved += bases
                continue
            elif op != op_:
                cigarNew.append(val)
            else:
                bases = bases + bases_
                cigarNew[-1] = (op, bases)
        
        # debug check on valid cigar string
        numCigarBases = sum((bases for (op, bases) in cigarNew if op in [0, 1, 4]))
        if numCigarBases != readLen:
            raise Exception("new cigar does not match read length!", read.qname, read.cigarstring, alignmentReoriented, cigarNew, basesGenomeNeeded)
    
        # identify 100% clipped reads - need to drop these (and PE mate read)
        ops = [x0 for (x0, x1) in cigarNew]
        if not 0 in ops:
            keepFlag = False
        else:
            # reverse the new cigar if negative strand
            if alignmentReoriented:
                cigarNew.reverse()
            
            # update the read alignment start position, if it has changed from the cigar edit (aend must be auto recalculated by pysam)
            if not alignmentReoriented:
                pos += basesGenome + (basesGapRemoved - basesGap)
                
    return cigarNew, pos, keepFlag


#------------------------------------------------------------------------------#
def filtr(chunk, speSide):
    # initialize
    readCounts = [0] * numCountFields
    speSide1 = (speSide == "1")
    speSide2 = not speSide1
    primerqc = []
    readClips = []
    
    #---------- loop over reads in a chunk
    for reads in chunk:
        # initialize
        speNotMap = True
        offTarget = True
        notProperPair = True        
        rdPrim1 = []
        rdPrim2 = []
        rdSupp1 = []
        rdSupp2 = []
        alignGroup = []        
        
        #----- loop over a collection of alignments with same read Id (alignment group)
        for idx in range(len(reads)):
            rd = reads[idx]

            #----- get read info, calculate some necessary values
            flag = int(rd[sam_flag])
            rd[sam_flag] = flag
            pos = int(rd[sam_pos])
            rd[sam_pos] = pos
            rdrev = True if flag & sam_flag_is_reverse else False
            cigarVec = decode(rd[sam_cigar])
            aend = pos + sum([bases for (op, bases) in cigarVec if op in [0, 3]])
            isRead1 = flag & sam_flag_is_read1
            
            # add SAM tags, mark alignment as off-target, change status downstream if otherwise
            rd.extend([rdrev, cigarVec, aend, 0])
                        
            if flag & sam_flag_supplementary:
                readCounts[index_read_supplementary] += 1
                if flag & sam_flag_read_is_unmapped:
                    continue
                if isRead1:
                    rdSupp1.append(rd)
                else:
                    rdSupp2.append(rd)
                continue
        
            if flag & sam_flag_mapped_in_proper_pair:
                notProperPair = False
            
            #----- check if mapped and on-target
            if flag & sam_flag_read_is_unmapped:
                continue
            else:
                if (speSide1 and isRead1) or (speSide2 and not isRead1):
                    speNotMap = False
                    
                    #----- check if on-target
                    primerId, readType, loc3, fragmentInfo = check(rd)
                    alignGroup.append((readType, fragmentInfo))
                    
                    # mark on-target and adjust primer Id based on mapping 
                    if readType:
                        rd[-1] = 1
                        rd[sam_tag_primer_id] = [primerId]
                        primerIdOnTarget = primerId
                        offTarget = False                        
                        chromSpe, posSpe, posMateSpe = rd[sam_rname], pos, rd[sam_pos_mate]
                        
                # update primeray alignment list
                if isRead1:
                    rdPrim1.append(rd)
                else:
                    rdPrim2.append(rd)
        
        #---------- general filtering 
        if notProperPair:
            readCounts[index_not_proper_pair] += 1
        
        if speNotMap:
            readCounts[index_spe_not_mapped] += 1
            primerId = rd[sam_tag_primer_id][0]
        else:
            readCounts[index_pass_filter] += 1         
        
        sameChrom = False
        notSameLoc = False
        reverseReads = True
        lowMapq = True        
        rdToClipSpe = []
        rdToClipUmi = []

        if speSide1:
            rdPrimSpe = rdPrim1
            rdPrimUmi = rdPrim2
            rdSuppSpe = rdSupp1
            rdSuppUmi = rdSupp2
        else:
            rdPrimSpe = rdPrim2
            rdPrimUmi = rdPrim1
            rdSuppSpe = rdSupp2
            rdSuppUmi = rdSupp1
        
        numAlignSpe = len(rdPrimSpe)
        numAlignUmi = len(rdPrimUmi)            
        numAlignSuppSpe = len(rdSuppSpe)
        numAlignSuppUmi = len(rdSuppUmi)

        if offTarget:
            readCounts[index_multi_align] += numAlignSpe            
            primerqc.append((primerId, rd[sam_tag_umi], alignGroup))
        else:            
            readCounts[index_on_target] += 1
            readCounts[index_multi_align] += (numAlignSpe - 1)
            primerqc.append((primerIdOnTarget, rd[sam_tag_umi], alignGroup))

            for idx in range(numAlignSpe):
                rd1 = rdPrimSpe[idx]
                chrom1 = rd1[sam_rname]
                pos1 = rd1[sam_pos]
                if chrom1 == chromSpe and pos1 == posSpe:
                    aend1 = rd1[sam_tag_aend]
                    rd1rev = rd1[sam_tag_revalign]
                    loc1 = aend1 - 1 if rd1rev else pos1
                    mapq1 = rd1[sam_mapq]
                    rdToClipSpe.append(rd1)
                    break
            
            chrom2 = None
            rd2rev = None
            loc2 = None
            mapq2 = None
            
            for idx in range(numAlignUmi):
                rd2 = rdPrimUmi[idx]
                chrom2 = rd2[sam_rname]
                pos2 = rd2[sam_pos]
                if pos2 == int(posMateSpe):
                    aend2 = rd2[sam_tag_aend]
                    rd2rev = rd2[sam_tag_revalign]
                    loc2 = aend2 - 1 if rd2rev else pos2
                    mapq2 = rd2[sam_mapq]
                    rdToClipUmi.append(rd2)
                    rd2[sam_tag_ontarget] = 1
                    rd2[sam_tag_primer_id] = [primerIdOnTarget]
                    break
                
            # not mapped to same chrom
            if chrom1 == chrom2:
                sameChrom = True
            else:                    
                # reads not mapped to same locus
                if loc2 == None or abs(loc1 - loc2) > maxSameLocDist:
                    notSameLoc = True

            # odd alignment orientation
            if rd1rev != rd2rev:
                reverseReads = False

            # if R1 or R2 not uniquely mapped 
            if mapq1 == "255" and mapq2 == "255":
                lowMapq = False

            # collect on-target supplementary alignments
            if numAlignSuppSpe:
                suppAlignments = getSuppTagAlignment(rd1[:sam_tag_umi])
                for idxSupp in range(numAlignSuppSpe):
                    rdsupp = rdSupp1[idxSupp]
                    if (rdsupp[sam_rname], rdsupp[sam_pos], rdsupp[sam_cigar]) in suppAlignments:
                        rdsupp[-1] = 1
                        rdsupp[sam_tag_primer_id] = [primerIdOnTarget]
                        rdToClipSpe.append(rdsupp)
            
            if len(rdToClipUmi) and numAlignSuppUmi:
                suppAlignments = getSuppTagAlignment(rd2[:sam_tag_umi])
                for idxSupp in range(numAlignSuppUmi):
                    rdsupp = rdSupp2[idxSupp]
                    if (rdsupp[sam_rname], rdsupp[sam_pos], rdsupp[sam_cigar]) in suppAlignments:
                        rdsupp[-1] = 1
                        rdsupp[sam_tag_primer_id] = [primerIdOnTarget]
                        rdToClipUmi.append(rdsupp)
        
            #---------- clipping reads (for variant calling)
            readClipSpe, readClipUmi = [], []
                        
            #----- clip primer at 5' of SPE side
            rdToClipSpe3 = []
            for idx in range(len(rdToClipSpe)):
                rd = rdToClipSpe[idx]
                alignmentReoriented = rd[sam_tag_revalign]                    
                pos = rd[sam_pos]
                
                if rd[sam_rname] == chromSpe:
                    direction = "1" if alignmentReoriented else "0"
                    
                    basesGenomeNeeded = rd[sam_tag_aend] - loc3 if alignmentReoriented else loc3 + 1 - pos
                    if basesGenomeNeeded > 0:                     
                        cigarNew, pos, keepRead = clip(rd[sam_qname], pos, rd[sam_tag_cigar], alignmentReoriented, basesGenomeNeeded)
                        
                        if keepRead:
                            # update cigar tag, cigar field, and position
                            rd[sam_tag_cigar] = cigarNew                    
                            rd[sam_cigar] = "".join([str(cigarNew[idx][1]) + cigarOpChar[cigarNew[idx][0]] for idx in range(len(cigarNew))])
                            rd[sam_pos] = pos
                            # rd[sam_tlen] = int(rd[sam_tlen])
                            # rd[sam_tlen] = rd[sam_tlen] - 
                            # update set of reads for 3' clipping
                            rdToClipSpe3.append(rd)
                            
                        continue
                            
                # update set of reads for 3' clipping
                rdToClipSpe3.append(rd)
                        
            #----- clip 3' ends
            edgeUmi = {}
            if len(rdToClipSpe3):
                if len(rdToClipUmi):
                    # clip 3' end of UMI side reads
                    for idx in range(len(rdToClipUmi)):                        
                        rd = rdToClipUmi[idx]
                        alignmentReoriented = not rd[sam_tag_revalign]
                        rdChrom = rd[sam_rname]
                        pos = rd[sam_pos]
                        aend = rd[sam_tag_aend]
                        
                        if rdChrom in edgeUmi:
                            edgeUmi[rdChrom].append((pos, aend))
                        else:
                            edgeUmi[rdChrom] = [(pos, aend)]
                        
                        if rdChrom == chromSpe:
                            basesGenomeNeeded = aend - loc3 if alignmentReoriented else loc3 + 1 - pos
                            
                            if basesGenomeNeeded > 0:                     
                                cigarNew, pos, keepRead = clip(rd[sam_qname], pos, rd[sam_tag_cigar], alignmentReoriented, basesGenomeNeeded)                            
                                if keepRead:
                                    # update cigar tag, cigar field, and position
                                    rd[sam_tag_cigar] = cigarNew                    
                                    rd[sam_cigar] = "".join([str(cigarNew[idx][1]) + cigarOpChar[cigarNew[idx][0]] for idx in range(len(cigarNew))])
                                    rd[sam_pos] = pos
                                    readClipUmi.append(rd)                            
                                continue                                    
                            
                        # update set of reads for 3' clipping
                        readClipUmi.append(rd)
                
                    # clip 3' end of SPE side reads
                    for idx in range(len(rdToClipSpe3)):
                        rd = rdToClipSpe3[idx]
                        alignmentReoriented = not rd[sam_tag_revalign]
                        rdChrom = rd[sam_rname]
                        pos = rd[sam_pos]
                        
                        keepRead = True
                        for (chromUmi, edges) in edgeUmi.items():
                            if rdChrom == chromUmi:
                                # find corresponding edge (right or left)
                                if alignmentReoriented:
                                    loc = max([x[1] for x in edges])
                                    basesGenomeNeeded = rd[sam_tag_aend] - loc
                                else:
                                    loc = min([x[0] for x in edges])
                                    basesGenomeNeeded = loc + 1 - pos
                                
                                if basesGenomeNeeded > 0:                     
                                    cigarNew, pos, keepRead = clip(rd[sam_qname], pos, rd[sam_tag_cigar], alignmentReoriented, basesGenomeNeeded)                                
                                    if keepRead:
                                        # update cigar tag, cigar field, and position
                                        rd[sam_tag_cigar] = cigarNew                    
                                        rd[sam_cigar] = "".join([str(cigarNew[idx][1]) + cigarOpChar[cigarNew[idx][0]] for idx in range(len(cigarNew))])
                                        rd[sam_pos] = pos                            
                                        # update set of reads for 3' clipping
                                        readClipSpe.append(rd)
                                        # already saved
                                        keepRead = False                            
                                    continue
                                
                        # update set of reads for 3' clipping, if needed
                        if keepRead:
                            readClipSpe.append(rd)                        
                else:
                    readClipSpe = rdToClipSpe3
                    
            if len(readClipSpe):
                for rdClip in (readClipSpe + readClipUmi):
                    rdClip[sam_flag] = str(rdClip[sam_flag])
                    rdClip[sam_pos] = str(rdClip[sam_pos])
                    rdClip = rdClip[:sam_tag_umi] + [umiTag + ":Z:" + rd[sam_tag_umi] + "-" + direction, primerIdTag + ":i:" + str(rd[sam_tag_primer_id][0]), onTargetTag + ":i:" + str(rd[sam_tag_ontarget])]
                    readClips.append("\t".join(rdClip))

        # count different categories
        if not sameChrom:           
            readCounts[index_not_same_chrom] += 1
        
        if notSameLoc:
            readCounts[index_not_same_loc] += 1
        
        if reverseReads:
            readCounts[index_reverse_align] += 1
        
        if lowMapq:
            readCounts[index_low_mapq] += 1     
        
    return readClips, primerqc, readCounts


#------------------------------------------------------------------------------#
def analyze(primerId, primerInfo, fragments):    
    '''
    process one primer
    '''
    
    # initialization
    (geneId, geneSymbol, geneStrand, chrom, loc5, loc3, strand, primer, target, distance) = primerInfo
    
    timeStart = datetime.now()
    
    if target == '-1':
        noTarget = True
        numTargets = 0
    else:
        noTarget = False        
        targets = target.split(targetSeparator)
        numTargets = len(targets)
        
        dnaCoverageRead = [0] * numTargets
        dnaCoverageMinBp = [0] * numTargets
        
        distances = distance.split(targetSeparator)
        if numTargets != len(distances):
            raise Exception("Error: target " + target + " does not match distance info " + distances)
        
        rnaDistances = []
        numRnaDistances = []
        for idx in range(numTargets):
            z = [int(x) for x in distances[idx].split(distanceSeparator)]
            rnaDistances.append(z)
            numRnaDistances.append(len(z))

        rnaCoverageRead = [[0] * numRnaDistances[x] for x in range(numTargets)]
        rnaCoverageMinBp = [[0] * numRnaDistances[x] for x in range(numTargets)]
    
    readCount = 0
    readOnTargetCount = 0    
    multiAlignCount = 0
    
    offSiteRead = defaultdict(lambda: 0)
    umis = []
    umisOnTarget = []
    fragmentLengths = [0] * numSegments
    
    #---------- main loop
    for umi, alignGroup in fragments:
        # update total read count
        readCount += 1        

        # initialize
        isOnTarget = False
        numAlign = len(alignGroup)
        readOnTargetOne = 0
        
        # loop over all alignments for same read Id
        for (readType, fragmentInfo) in alignGroup:
            
            (alnChrom, alnLoc, alnPrimLoc, rnaBp, dnaBp) = fragmentInfo            
        
            multiAlignCount += 1
            
            # if off-target
            if readType == 0:
                offtargetSite = alnChrom + ":" + str(alnLoc) if alnPrimLoc == "" else alnChrom + ":" + str(alnPrimLoc)
                offSiteRead[offtargetSite] += 1
                continue
            
            #----- if on-target
            umisOnTarget.append(umi)                
            readOnTargetOne = 1
            
            for idx in range(rnaBp//segmentLength):
                fragmentLengths[idx] += 1 
           
            # if no target from DB for this primer
            if noTarget:
                continue                
            
            for targetIdx in range(numTargets):
                numRnaDistance = numRnaDistances[targetIdx]
                
                #----- RNA-based target coverage      
                for rnaDistanceIdx in range(numRnaDistance):
                    rnaDistance = rnaDistances[targetIdx][rnaDistanceIdx]
                    
                    # count reads that cover target
                    if rnaDistance != -1 and rnaBp >= rnaDistance:
                        rnaCoverageRead[targetIdx][rnaDistanceIdx] += 1
                        
                        # count reads that cover target (RNA) with a min bp over
                        if rnaBp >= (rnaDistance + followBasesMin):
                            rnaCoverageMinBp[targetIdx][rnaDistanceIdx] += 1
                    
                #----- DNA-based target coverage
                if dnaBp[targetIdx] >= followBasesMin:
                    dnaCoverageRead[targetIdx] += 1
                    dnaCoverageMinBp[targetIdx] += dnaBp[targetIdx]

        # update on-target read count (adjust for more than on-target SPE-side alignment because of multialignments)
        readOnTargetCount += readOnTargetOne
        
        # collect UMI for this alignment group
        umis.append(umi)
    
    #---------- summarize            
    readOnTargetPct = round(100.0 * readOnTargetCount/readCount, 2) if readCount > 0 else 0.00        
    
    # get MT counts
    umiCount = len(set(umis))
    umiOnTargetCount = len(set(umisOnTarget))
    umiOnTargetPct = round(100.0 * umiOnTargetCount/umiCount, 2) if umiCount > 0 else 0.00        
    
    # make histogram based on length bins
    fragmentLengthHistogram = [round(100.0 * fragmentLengths[x]/readOnTargetCount, 1) for x in range(numSegments)] if readOnTargetCount > 0 else [0.0] * numSegments  
    
    # top off-target sites
    offSiteTop = [""] * offSiteTopMax * offSiteFields
    offSiteCount = len(offSiteRead)
    offSiteReadSort = sorted(offSiteRead.items(), key = lambda x: x[1], reverse = True)
    offSiteReadTotal = sum([x[1] for x in offSiteReadSort])
    for offSiteIdx in range(min(offSiteCount, offSiteTopMax)):
        (offSite, offSiteReadCount) = offSiteReadSort[offSiteIdx]
        offSiteTop[offSiteFields*offSiteIdx] = offSite
        offSiteTop[offSiteFields*offSiteIdx+1] = offSiteReadCount
        offSiteTop[offSiteFields*offSiteIdx+2] = round(100.0*offSiteReadCount/offSiteReadTotal, 1) if offSiteReadTotal > 0 else ""
        
    # prepare output
    outVec = [int(primerId)]
    outVec.extend(primerInfo)
    outVec.extend([readCount, readOnTargetCount, readOnTargetPct, umiCount, umiOnTargetCount, umiOnTargetPct])
    
    dnaCoveragesRead, dnaCoverageRatios, dnaCoveragesAve, dnaCoverageAveMax, rnaDistanceMin, dnaCoverageRatioMinDistance, rnaDistanceMax, dnaCoverageRatioMaxDistance = [""] * 8
    rnaCoveragesRatiosAll, rnaCoverageRatioMinDistance, rnaCoverageRatioMaxDistance, rnaCoveragesMinBpRatiosAll = [""] * 4
    
    #----- if target provided
    firstTarget = True
    dnaCoverageAveMax = 0 if numTargets else ""
    
    # loop over targets 
    for targetIdx in range(numTargets):
        #----- DNA-based coverage
        dnaCoverageRead_ = dnaCoverageRead[targetIdx]
        dnaCoverageAve = 0
        dnaCoverageRatio = 0.0
        dnaCoverageAve = round(dnaCoverageMinBp[targetIdx]/dnaCoverageRead_, 1) if dnaCoverageRead_ > 0 else 0
        if dnaCoverageAveMax < dnaCoverageAve:
            dnaCoverageAveMax = dnaCoverageAve
        dnaCoverageRatio = round(100.0 * dnaCoverageRead_/readOnTargetCount, 1) if readOnTargetCount > 0 else 0.0
        
        #----- RNA-based coverage        
        # loop over distances for each primer-target pair
        rnaCoverageRatios = []
        rnaCoverageMinBpRatios = []
        firstDistance = True        
        numRnaDistance = numRnaDistances[targetIdx]
        
        for rnaDistanceIdx in range(numRnaDistance):
            rnaDistance = rnaDistances[targetIdx][rnaDistanceIdx]        
            rnaCoverageRatio = round(100.0 * rnaCoverageRead[targetIdx][rnaDistanceIdx]/readOnTargetCount, 1) if readOnTargetCount > 0 else 0.0
            rnaCoverageRatios.append(str(rnaCoverageRatio))
            
            rnaCoverageMinBpRatio = round(100.0 * rnaCoverageMinBp[targetIdx][rnaDistanceIdx]/readOnTargetCount, 1) if readOnTargetCount > 0 else 0.0
            rnaCoverageMinBpRatios.append(str(rnaCoverageMinBpRatio))

            # find maximum and minimum distance(s)
            if firstDistance:
                firstDistance = False
                rnaDistanceMin = rnaDistance
                rnaCoverageRatioMinDistance = rnaCoverageRatio
                dnaCoverageRatioMinDistance = dnaCoverageRatio
                rnaDistanceMax = rnaDistance
                rnaCoverageRatioMaxDistance = rnaCoverageRatio
                dnaCoverageRatioMaxDistance = dnaCoverageRatio
            else:
                if rnaDistanceMin > rnaDistance:
                    rnaDistanceMin = rnaDistance
                    rnaCoverageRatioMinDistance = rnaCoverageRatio
                    dnaCoverageRatioMinDistance = dnaCoverageRatio
                if rnaDistanceMax < rnaDistance:
                    rnaDistanceMax = rnaDistance
                    rnaCoverageRatioMaxDistance = rnaCoverageRatio
                    dnaCoverageRatioMaxDistance = dnaCoverageRatio
                    
        rnaCoveragesRatios = distanceSeparator.join(rnaCoverageRatios)
        rnaCoveragesMinBpRatios = distanceSeparator.join(rnaCoverageMinBpRatios)
        
        # assemble output variables for multiple targets/distances
        if firstTarget:
            dnaCoveragesRead = str(dnaCoverageRead_)
            dnaCoverageRatios = str(dnaCoverageRatio)
            dnaCoveragesAve = str(dnaCoverageAve)
            rnaCoveragesRatiosAll = str(rnaCoveragesRatios)
            rnaCoveragesMinBpRatiosAll = str(rnaCoveragesMinBpRatios)
            firstTarget = False
        else:
            dnaCoveragesRead += targetSeparator + str(dnaCoverageRead_)
            dnaCoverageRatios += targetSeparator + str(dnaCoverageRatio)
            dnaCoveragesAve += targetSeparator + str(dnaCoverageAve)
            rnaCoveragesRatiosAll += targetSeparator + str(rnaCoveragesRatios)
            rnaCoveragesMinBpRatiosAll += targetSeparator + str(rnaCoveragesMinBpRatios)
            
    outVec.extend([dnaCoveragesRead, dnaCoverageRatios, dnaCoveragesAve, dnaCoverageAveMax, rnaDistanceMin, dnaCoverageRatioMinDistance, rnaDistanceMax, dnaCoverageRatioMaxDistance])
    outVec.extend(fragmentLengthHistogram)
    outVec.extend([rnaCoveragesRatiosAll, rnaCoverageRatioMinDistance, rnaCoverageRatioMaxDistance, rnaCoveragesMinBpRatiosAll, multiAlignCount])
    outVec.append(offSiteCount)
    outVec.extend(offSiteTop)
    outVec = [str(x) for x in outVec]

    # save time
    logTime = datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "\t" + str(datetime.now() - timeStart) + "\tprimer Id " + str(primerId)
    
    return outVec, logTime


#------------------------------------------------------------------------------#
def qc(cfg, numCores):
    ''' 
    filtering steps includ:
        - general read quality metrics
        - on/off-target
        - clipping oligos based on alignment positions
        
        it collects alignments with same read Id in one list
        then it makes a chunk which is a list of these alignment groups
        finally, submits a list of chunks to threads, default is one chunk per thread
    '''
    # get parameter settings
    (runpath, runname, samplename, speSide, r1universal, r2universal, umiLen, umiOffset, panelId, platform) = cfg    
    
    # read set
    readset = runname + "." + samplename
    
    # I/O 
    filein = open(readset + "/align.star.sam", "r")
    fileoutsam = open(readset + "/qc.clip.sam", "w")
    fileoutqc = open(readset + "/qc.primer.txt", "w")
    fileoutrc = open(readset + "/qc.read.txt", "w")
    fileoutfilterlog = open(readset + "/qc.filter.log", "w")
    fileoutqclog = open(readset + "/qc.primer.log", "w")
    
    # get primers info
    primers = getPrimers(panelId)    
    numPrimers = len(primers)
    
    # initialize
    numJobs = numCores * numJobsScale
    numJobsMinus = numJobs - 1
    readPairCounts = [0] * numCountFields
    readOfftargetCounts = [0] * numPrimers
    readMultialignCounts = [0] * numPrimers
    numReadsInChunk = 0
    numReadPairsTotal = 0
    numAlignmentsTotal = 0
    numChunksInBucket = 0
    chunk = []
    bucket = []
    fragments = [[] for idx in range(numPrimers)]
    
    #---------- loop over SAM
    # print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " filter and clip")
    
    firstRead = True
    for line in filein:
        # header
        if line[0] == "@":
            fileoutsam.write(line)
            continue
        
        # get values
        read = line.strip().split("\t")
        
        # initialize
        dropReadPair = False
        flag = int(read[sam_flag])        
        (readId, umi, primerId) = read[sam_qname].split(readIdSeparator)[:3]
        if primerId == "-1":
            continue
        read[0] = readId
        primerIds = []
        primerInfo = []        
        for idx in primerId.split(","):
            primerIdInt = int(idx)
            primerIds.append(primerIdInt)
            (geneId, geneSymbol, geneStrand, chrom, loc5, loc3, strand, primer, target, distance) = primers[primerIdInt]
            primerInfo.append((chrom, int(loc5), int(loc3), int(strand), len(primer), target, distance))
        read.extend([umi, primerIds, primerInfo])
        
        # crash if read is not paired
        if not flag & sam_flag_read_is_paired:
            raise Exception("read not paired: " + readId)
        
        # crash if secondary (neither primary nor supplementary)
        if flag & sam_flag_not_primary and not flag & sam_flag_supplementary:
            raise Exception("secondary alignment found: " + readId)
        
        # update total number of alignments
        # numAlignmentsTotal += 1
        
        if firstRead:
            readIdLast = readId
            reads = [read]            
            firstRead = False
            continue
        
        if readId == readIdLast:
            reads.append(read)
            continue

        # update total number of alignment groups
        numReadPairsTotal += 1
        
        # update chunk
        if numReadsInChunk < chunkSize:            
            chunk.append(reads)
            numReadsInChunk += 1
        # update bucket
        elif numChunksInBucket < numJobsMinus:
            bucket.append(chunk)
            numChunksInBucket += 1
            chunk = [reads]
            numReadsInChunk = 1
        # run last bucket
        else:
            # add last chunk
            bucket.append(chunk)
            
            timeStart = datetime.now()

            pool = multiprocessing.Pool(processes = numCores)
            output = [pool.apply_async(filtr, args = (bucket[idx], speSide)) for idx in range(numJobs)]
            pool.close()
            pool.join()                
            results = [p.get() for p in output]
            
            timeEnd = datetime.now()
            
            # process results
            for (readClips, primerqc, readCounts) in results:
                # update read/alignment counts (after map reduce step)
                for idx in range(numCountFields):
                    readPairCounts[idx] += readCounts[idx]
            
                # save clipped reads to disk
                for readClip in readClips:
                    fileoutsam.write(readClip + "\n")
                
                # collect fragment info to perform primer qc next
                for (primerId, umi, alignGroup) in primerqc:
                    fragments[primerId].append((umi, alignGroup))
                    
            # reset
            fileoutfilterlog.write("\t".join([datetime.now().strftime("%Y/%m/%d %H:%M:%S"), str(timeEnd-timeStart), str(datetime.now()-timeEnd), str(datetime.now()-timeStart), str(numReadPairsTotal-1)]) + "\n")
            results = []            
            bucket = []
            numChunksInBucket = 0
            chunk = [reads]
            numReadsInChunk = 1
            
        #----- initialize for next group (with same read Id)
        readIdLast = readId
        reads = [read]

    #---------- process last chunk, if any chunk(s) yet to be processe
    if len(reads):
        timeStart = datetime.now()
        numReadPairsTotal += 1
        if len(bucket):
            chunk.append(reads)
            bucket.append(chunk)
            pool = multiprocessing.Pool(processes = numCores)
            output = [pool.apply_async(filtr, args = (bucket[idx], speSide)) for idx in range(len(bucket))]
            pool.close()
            pool.join()    
            results = [p.get() for p in output]            
        elif len(chunk):
            chunk.append(reads)
            results = [filtr(chunk, speSide)]
        else:        
            results = [filtr([reads], speSide)]
        
        # process last results
        for (readClips, primerqc, readCounts) in results:
            # update read/alignment counts (after map reduce step)
            for idx in range(numCountFields):
                readPairCounts[idx] += readCounts[idx]
        
            # save clipped reads to disk
            for readClip in readClips:
                fileoutsam.write(readClip + "\n")
            
            # collect fragment info to perform primer qc next
            for (primerId, umi, alignGroup) in primerqc:
                fragments[primerId].append((umi, alignGroup))
        fileoutfilterlog.write("\t".join([datetime.now().strftime("%Y/%m/%d %H:%M:%S"), str(datetime.now()-timeStart), str(numReadPairsTotal)]) + "\n")
        results = []
    fileoutsam.close()
    fileoutfilterlog.close()
    
    #---------- primer qc
    # print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " primer qc")
    
    brk = []
    headSeg = []
    for i in segmentsIdx:
        brk.append((i+1)*segmentLength)
        headSeg.append(str(brk[i]) + "bp_pct")
    
    # make and print(header)
    header = ["primer_index", "gene_Id", "gene_symbol", "gene_strand", "chrom", "loc5", "loc3", "strand", "primer", "target_sites", "dist_to_targets_RNA", "reads", "ontarget", "ontarget_pct", "umi", "umi_ontarget", "umi_ontarget_pct", "covering_target_DNA", "covering_target_DNA_pct", "avg_following_bp", "max_avg_following_bp","min_RNA_distance", "covering_target_DNA_with_min_RNA_dist_pct", "max_RNA_distance", "covering_target_DNA_with_max_RNA_dist_pct"]
    header.extend(headSeg)
    header.extend(["target_dist_RNA_pct", "min_target_dist_RNA_pct", "max_target_dist_RNA_pct", "target_dist_RNA_with_min_bp_coverage_pct", "multialignments", "off_target_sites"])
    headerOffSiteTitle = ["offtarget_site." + str(idx + 1) for idx in range(offSiteTopMax)]
    headerOffSiteCount = ["offtarget_count." + str(idx + 1) for idx in range(offSiteTopMax)]
    headerOffSitePct = ["offtarget_pct." + str(idx + 1) for idx in range(offSiteTopMax)]
    headerOffSite = [""] * offSiteTopMax * offSiteFields
    headerOffSite[::3] = headerOffSiteTitle
    headerOffSite[1::3] = headerOffSiteCount
    headerOffSite[2::3] = headerOffSitePct
    header.extend(headerOffSite)
    
    fileoutqc.write("\t".join(header) + "\n")

    # run primer qc
    pool = multiprocessing.Pool(processes = numCores)
    output = [pool.apply_async(analyze, args = (primerId, primers[primerId], fragments[primerId])) for primerId in range(numPrimers)]
    pool.close()
    pool.join()    
    results = [p.get() for p in output]
    
    # save primer qc to disk
    for vec, logtime in results:
        fileoutqc.write("\t".join(vec) + "\n")
        fileoutqclog.write(logtime + "\n")
    
    fileoutqc.close()
    fileoutqclog.close()
    
    #---------- write counts
    pctReasdPassed = round(100 * readPairCounts[index_pass_filter] / numReadPairsTotal, 1) if numReadPairsTotal > 0 else 0.0
    pctReasdPairKept = round(100 * readPairCounts[index_on_target] / numReadPairsTotal, 1) if numReadPairsTotal > 0 else 0.0
    # fileout.write("{:,}".format(numAlignmentsTotal) + "\t" +  "alignments total (R1, R2, or supplementary) in BAM file\n")
    fileoutrc.write("{:,}".format(numReadPairsTotal ) + "\t" +  "alignment pairs (with unique read Ids) in BAM\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_read_supplementary] ) + "\t" +  "alignments supplementary\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_multi_align] ) + "\t" +  "multi-alignments\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_not_proper_pair] ) + "\t" +  "alignment pairs, R1 or R2 not mapped\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_not_same_chrom] ) + "\t" +  "alignment pairs not mapped to same chromosome\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_not_same_loc] ) + "\t" +  "alignment pairs not mapped to same locus (read ends locs > " + str(maxSameLocDist) + "bp)\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_reverse_align] ) + "\t" +  "alignment pairs FF or RR mapping orientation\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_low_mapq] ) + "\t" +  "alignment pairs R1 or R2 not uniquely mapped\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_spe_not_mapped] ) + "\t" +  "alignment pairs dropped, SPE side not mapped\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_pass_filter] ) + "\t" +  "passed filters\n")
    fileoutrc.write(str(pctReasdPassed) + "\t" +  "% passed filters\n")
    fileoutrc.write("{:,}".format(readPairCounts[index_on_target] ) + "\t" +  "on-target alignments, kept\n")
    fileoutrc.write(str(pctReasdPairKept) + "\t" +  "% on-target alignments, kept\n")
    fileoutrc.close()


#------------------------------------------------------------------------------#
def ge(cfg):
    '''
    simple algorithm to estimate gene expression based on marked primers from DB
    '''
    pass


#------------------------------------------------------------------------------#
def combine(cfgs):
    '''
    combine gene expression estimates from all read sets in one table
    '''
    pass


#------------------------------------------------------------------------------#
def summary(cfgs):
    '''
    create final summary of the whole run
    '''
    
    # combine gene expression estimates
    combine(cfgs)
    
    # create summary
    subprocess.check_call("rm -rf summary && mkdir -p summary", shell = True)    
    for readset, cfg in cfgs.items():        
        subprocess.check_call("cp " + readset + "/qc.primer.txt summary/" + readset + ".primer.txt" , shell = True)
        subprocess.check_call("cat " + readset + "/trim.read.txt " + readset + "/qc.read.txt > summary/" + readset + ".read.txt" , shell = True)
        # subprocess.check_call("cp " + readset + "/qc.ge.txt summary/" + readset + ".expression.txt" , shell = True)

#------------------------------------------------------------------------------#
# run all functions in this module
def run(cfgs, numCores):     
    # primer qc
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " filter, clip, and primer qc")
    numThreads = min(3, numCores // len(cfgs))
    joblib.Parallel(n_jobs = numCores)(joblib.delayed(qc)(cfgs[readset], numThreads) for readset in cfgs)
    
    # gene expression
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " gene expression")
    joblib.Parallel(n_jobs = numCores)(joblib.delayed(ge)(cfgs[readset]) for readset in cfgs)
    
    # summary and combine
    summary(cfgs)