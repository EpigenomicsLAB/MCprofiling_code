library(epistats)
library(optparse)

option_list <- list(
  make_option(c("-b", "--bamfile"),
              help="Bamfile to extract the epialleles from",
              metavar="PATH TO THE BAMFILE"),
  make_option(c("-g", "--Genome"),
              help="Genome used to perform the alignment",
              metavar="PATH TO THE GENOME FASTA"),
  make_option(c("-t", "--threshold"), type = "integer", default = 50,
              help = "Minimum coverage value to select regions [default= %default]",
              metavar = "integer"),
  make_option(c("-s", "--stranded"), type = "logical", default = FALSE,
              help = "Analysis mode [default= %default]",
              metavar = "logical"),
  make_option(c("-m", "--mode"), type = "character", default = "CG",
              help = "Mode string to perform the analysis  [default= %default]",
              metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "path",
              help = "Set the output path [default= %default]",
              metavar = "character"),
  make_option("--modeAnalysis", type = "character", default = "n_cg",
              help = "Option to specify the analysis mode
                      Possible values: 'n_cg' or 'windows'[default= %default]",
              metavar = "Character"),
  make_option(c("-w", "--window"), type = "integer", default = 50,
              help = "The width of the sliding window [default= %default]",
              metavar = "integer"),
  make_option("--step", type = "integer", default = 1,
              help = "Integer indicating the step between two slidind windows [default= %default]",
              metavar = "integer"),
  make_option("--minC", type = "integer", default = 4,
              help = "Integer indicating the minimum number of C in the mode string 
              there should be in the window  [default= %default]",
              metavar = "integer"),
  make_option("--maxC", type = "integer", default = 10,
              help = "Integer indicating the maximum number of C in the mode string 
              there should be in the window [default= %default]",
              metavar = "integer"),
  make_option(c("-n", "--n"), type = "integer", default = 4,
              help = "Integer indicating the maximum number of C in the mode string 
              there should be in the bin if --modeAnalysis is set as 'n_cg' [default= %default]",
              metavar = "integer"),
  make_option("--minBinsize", type = "integer", default = 8,
              help = "Integer indicating the number of C in the mode string 
              there should be in bin [default= %default]",
              metavar = "integer"),
  make_option("--maxBinsize", type = "integer", default = 50,
              help = "Integer indicating the number of C in the mode string
              there should be in the bin [default= %default]",
              metavar = "integer"),
  make_option("--cores", type = "integer", default = 50,
              help = "Integer indicating the number of cores to use [default= %default]",
              metavar = "integer"),
  make_option("--bisu", type = "integer", default = 0,
              help = "Integer indicating the percetage of bisulfite efficiency
              considered to keep the read to analyse[default= %default]",
              metavar = "integer"),
  make_option("--rmAmb", type = "logical", default = TRUE,
              help = "Logical indicating if the reads with ambiguity 
              shoul be discarded or not [default= %default]",
              metavar = "logical"),
  make_option("--retainReads", type = "logical", default = TRUE,
              help = "Logical indicating if the reads without cytosines to 
              perform bisulfite quality check shoul be kept for the analysis
              [default= %default]",
              metavar = "logical"),
  make_option("--getCpos", type = "logical", default = FALSE,
              help = "Option to write a file containing the C coordinates
              analysed by the routine [default= %default]",
              metavar = "logical"),
  make_option("--retainChr", type = "character", default = NULL,
              help = "A character or a vector of characters containing the chr names
              to retain to perform the analysis [default= %default]",
              metavar = "character"),
  make_option("--keepStChr", type = "logical", default = TRUE,
              help = "A logical value indicating if only the standarde chromosomes
              there should be kept for the analysis [default= %default]",
              metavar = "logical"),
  make_option("--keepM", type = "logical", default = FALSE,
              help = "A logical value indicating if the mithocondrial chromosome
              should be kept for the analysis [default= %default]",
              metavar = "logical")
)

opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE, prog = "CMD.R <in.bam> <genome.fa>")
opt = parse_args(opt_parser)

if (is.null(opt$bamfile)){
  print_help(opt_parser)
  stop("Bamfile must be provided", call.=FALSE)
}

if (is.null(opt$Genome)){
  print_help(opt_parser)
  stop("Genome fasta file must be provided", call.=FALSE)
}

bamfile= opt$bamfile
genomepath= opt$Genome 
threshold= opt$threshold
stranded= opt$stranded
mode= opt$mode
mode_regions= opt$modeAnalysis 
## windows
window= opt$window
step= opt$step
min.C= opt$minC
max.C= opt$maxC
## n_cg
n= opt$n
min.binsize= opt$minBinsize
max.binsize= opt$maxBinsize
## optional arguments 
cores= opt$cores
bisu.Thresh= opt$bisu
remove.Amb= opt$rmAmb
retain.reads= opt$retainReads
get.cPos= opt$getCpos
retainChr = opt$retainChr
keepStChr = opt$keepStChr
keepM = opt$keepM

#----------------------Chunk1: load input----------------------------
data=loadInput(bamfile, genomepath)
Genome=data[[2]]
align=data[[1]]

#--------------------Chunk2: filter bam sequences--------------------
#########filter bam sequences
align=filterBAM(algn = align,
                keepStChr = keepStChr, 
                keepM = keepM, 
                retainChr = retainChr
                )

#########filter bam by coverage

if(mode_regions=="windows"){
  min.binsize=window
}
covered=filterByCoverage(align, threshold=threshold, minsize=min.binsize)

if(mode_regions=="windows"){
  #-------------------Chunk3: select windows---------------------------
  windows=makeWindows(gr=covered, 
                      Genome=Genome,
                      window=window, 
                      step=step,
                      mode=mode, 
                      min.C=min.C, 
                      max.C=max.C, 
                      cores=cores)
  cpos=windows$c_pos
  regions=windows$windows
} else if (mode_regions=="n_cg"){
  #------------------Chunk4: select epialleles-------------------------
  epialleles=makeBins(gr=covered, 
                      Genome=Genome, 
                      mode=mode, 
                      n=n, 
                      min.binsize=min.binsize, 
                      max.binsize=max.binsize,
                      cores=cores)
  cpos=epialleles$c_pos
  regions=epialleles$bins
}
#------------------Chunk5:  analyze epialleles
### Analyze epialleles
aname=gsub(tail(unlist(Biostrings::strsplit(bamfile, split="/")), n=1),pattern = "\\.bam", replacement = "")

myfuns=list("dist"=cdist,
            "epi"=epi,
            "singleton"=singleton,
            "maxfreq"=maxfreq,
            "shannon"=shannon,
            "mean_met"=meanMeth,
            "num_cg"=numCG,
            "num_reads"=num_reads)

output= epiAnalysis(align= align,
                    bin= regions,
                    aname=aname,
                    bisu.Thresh = bisu.Thresh,
                    stranded = stranded,
                    mode = mode,
                    remove.Amb = remove.Amb,
                    genome = Genome,
                    pathDir = out,
                    cores = cores,
                    retain.reads = retain.reads,
                    get.cPos = get.cPos,
                    myfuns = myfuns,
                    threshold = threshold)
