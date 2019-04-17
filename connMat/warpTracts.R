#!/data/grossman/pipedream2018/bin/R/R-3.4.3/bin/Rscript

library("optparse")


##
## Gets the directory containing this script, either relative or absolute depending on how script was invoked
##
getBinDir <- function() {
    argv = commandArgs(trailingOnly = FALSE)
    binDir = dirname(substring(argv[grep("--file=", argv)], 8))
    return(binDir)
}

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--input-file", type="character", default=NULL, dest="inputFile",
              help="Input file name or stdin", metavar="character"),
  make_option("--output-root", type="character", default=NULL, dest="outputRoot",   
              help="Output root", metavar="character"),
  make_option("--seed-image", type="character", default=NULL, dest="seedImage",
              help="Seed image", metavar="character"),
  make_option("--reference-image", type="character", default=NULL, dest="referenceImage",
              help="Reference image", metavar="character"),
  make_option("--dist-corr-inverse-warp", type="character", default="", dest="distCorrInverseWarp",
              help="Inverse warp from registration with fixed=reference.", metavar="character"),
  make_option("--dist-corr-affine", type="character", default="", dest="distCorrAffine",
              help="Affine transform from registration with fixed=reference", metavar="character"),
  make_option("--composed-inverse-warp", type="character", default="", dest="composedInverseWarp",
              help="Composed inverse warp from registration with fixed=reference image, if this is specified, it is used in place of the dist-corr warp / affine transforms.", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);

if (length(args)==0) {
  print_help(opt_parser)
  quit(save = "no")
}

opt = parse_args(opt_parser);

source(paste(getBinDir(), "warpTractsHelper.R", sep = "/"));

transformList = c()
invertList = c()

if (file.exists(opt$composedInverseWarp)) {
  transformList = c(opt$composedInverseWarp)
  invertList = c(F)
} else {
  transformList = c(opt$distCorrAffine, opt$distCorrInverseWarp)
  invertList = c(T,F)
}

warpTracts(opt$seedImage, opt$referenceImage, opt$inputFile, paste(opt$outputRoot, "TractsDeformed.Bfloat", sep = ""), transformList, invertList)

print("Done warping tracts")

warpedSeeds = warpSeeds(opt$seedImage, opt$referenceImage, transformList, invertList)

antsImageWrite(warpedSeeds[[2]], paste(opt$outputRoot, "SeedsDeformed.nii.gz", sep = ""))
