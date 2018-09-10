#------------------------------------------------------------------------------#
# general

sectionSeparator = "-" * 50
readIdSeparator = "|"
targetSeparator = "_"
distanceSeparator = "!"
isDebug = True
tagSeparator = "-"
numJobsScale = 16
chunkSize = 1024

#------------------------------------------------------------------------------#
# trim

maxErrorAdapter3 = 0.20
maxErrorAdapter5 = 0.15
minSeqLen = 20
maxAdapterLenDiffRatio = 0.125      # 1/8
polyA = "A" * 15
umiMaxLowqBases = 1
umiMaxNs = 1
readTrimmer = "/home/qiauser/code/read-trimmer/trimmer/run.py"
primerCol = 7   # 0-based

#------------------------------------------------------------------------------#
# align

maxSameLocDist = 2000
minMapqDna = 17
minAlignLen = 25
starIndex = "/home/qiauser/data/star_indices/human"

#------------------------------------------------------------------------------#
# qc

minBpCovered = 20
mtRawMax = 4**12/100
minDistPos = 5
followBasesMin = 25
offSiteTopMax = 10
umiTag = "Mi"
primerIdTag = "Pr"
onTargetTag = "Tt"


#------------------------------------------------------------------------------#
# vc

qimera = "/home/qiauser/code/qiaseq-qimera"
smCounter = "/home/qiauser/code/qiaseq-smcounter-v2/run.py"
exonModels = "/home/qiauser/data/star_indices/human/Ensembl81.human.rna.exonlocs.txt"
refGenome = "/home/qiauser/data/star_indices/human/GRCh38_e92_SIRV.fas"
srBed = "/home/qiauser/data/repeat_files/human/SR_LC_SL.hg38.bed"
repBed = "/home/qiauser/data/repeat_files/human/simpleRepeat.hg38.bed"
bkgErrorDistSimulation = "/home/qiauser/data/varcall_bkg_error_dist/bkg.error.v2.RData"
vcFragLenMax = 300
vcDepthMin = 100