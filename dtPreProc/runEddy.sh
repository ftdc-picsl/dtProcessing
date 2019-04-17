#!/bin/bash 

function checkFileExists {

  file=$1

  if [[ -f "$file" ]]; then
    return 0
  else
    echo " Required input files missing or not specified, see usage "
    exit 1
  fi 
}

if [[ $# -eq 0 ]]; then
  echo " 

  `basename $0` -i <data> -b <bvals> -r <bvecs> -m <mask> -a <acqParams> -o <outputRoot> [options]


   Process some data with FSL's eddy (version 5.0.9-eddy-patch). This script is designed for simple data, 
   and doesn't do anything with topup, field mapping, or lsr resampling.


  Required args:

   -i : input data (4D NIFTI)

   -b : bvals

   -r : bvecs
 
   -m : brain mask

   -a : acquisition parameters file

   -o : output root, including path and file root


  Options

   -n : index file. If not specified, it is assumed that all the acquisition parameters are the same, and an index file 
        referencing the first line of the acquisition parameters file is generated. This should work for most data acquired 
        in a single phase encode direction.


   -c : cores, number of qsub cores to request (default = 1). Not used if running locally

   -e : memory to request (default = 4G)

   -p : PEAS 1 (true) or 0 (false) (default = 1) 

   -q : Submit to qsub (1) or run locally (0) (default = 1)

   -t : Replace outliers (default = 0).


   Output is prepended with \${outputRoot}. Main output is 

   .nii.gz               Corrected DWI data
   .eddy_rotated_bvecs   Corrected bvecs

   Other diagnostic information from eddy should be reviewed, see eddy user guide for more information

"

exit 1

fi

binDir=`dirname $0`

acqParams=""
brainMask=""
bvals=""
bvecs=""
data=""
index=""
outputRoot=""

cores=1
doPeas=1
doRepol=0
ram="4G"
submitToQueue=1

while getopts "a:b:c:e:i:m:n:o:p:q:r:t:" OPT; do
    case $OPT in
	a)  # acqParams.txt
	    acqParams=$OPTARG
	    ;;
	b)  # bvals
	    bvals=$OPTARG
	    ;;
	c) # cores
	    cores=$OPTARG
	    ;;
	e)  # memory
	    ram=$OPTARG
	    ;;
	i)  # input image
	    data=$OPTARG
	    ;;
	m) # mask
	    brainMask=$OPTARG
	    ;;
	n) # index
	    index=$OPTARG
	    ;;
	o)  # output
	    outputRoot=$OPTARG
	    ;;
	p)
	    doPeas=$OPTARG 
	    ;;
	q) # Submit to queue
	    submitToQueue=$OPTARG
	    ;;
	r) # bvecs
	    bvecs=$OPTARG
	    ;;
	t)  # replace outliers
	    doRepol=$OPTARG
	    ;;
	\?) # getopts issues an error message
	    exit 1
	    ;;
    esac
done

checkFileExists "$acqParams"
checkFileExists "$brainMask"
checkFileExists "$bvals" 
checkFileExists "$bvecs"
checkFileExists "$data"

outputDir=`dirname $outputRoot`

if [[ ! -d "$outputDir" ]]; then
  mkdir -p "$outputDir"
fi

# If no index file, assume acqParams is one line, and all images have same parameters
if [[ ! -f "$index" ]]; then

    index="${outputDir}/index.txt"

    numMeas=`wc -w $bvals | cut -d ' ' -f 1`
    for i in `seq 1 $numMeas`; do
      echo "1" >> $index 
    done
fi


parallel=""

if [[ $cores -gt 1 ]]; then
  parallel="-pe unihost $cores -binding linear:$cores" 
fi


# Call runEddy_qscript.sh <data> <mask> <acqParams> <index> <bvals> <bvecs> <outputRoot> <doPeas> <doRepol>

exe="${binDir}/runEddy_qscript.sh $data $brainMask $acqParams $index $bvals $bvecs $outputRoot $doPeas $doRepol"

if [[ $submitToQueue -gt 0 ]]; then
    qsub -l h_vmem=${ram},s_vmem=${ram} $parallel -cwd -S /bin/bash -j y -o ${outputDir}/eddyLog.txt $exe
else
  $exe > ${outputDir}/eddyLog.txt

  if [[ $? -gt 0 ]]; then
      echo " eddy returned non-zero exit code $? "
      exit 1
  fi
fi

# Copy the mask and the bvals to the output directory for future processing
cp $bvals ${outputRoot}.bval 
cp $brainMask ${outputRoot}BrainMask.nii.gz
