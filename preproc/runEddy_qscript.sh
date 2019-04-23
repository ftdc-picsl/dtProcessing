#!/bin/bash

if [[ $# -eq 0 ]]; then

echo "
  $0 <data> <mask> <acqParams> <index> <bvals> <bvecs> <outputRoot> <doPeas> <doRepol>

  The full path to the output can be included (create necessary directories before calling this script) 
  
  Uses FSL 5.0.9-eddy-patch.

  This script uses mostly default settings. Some other settings you can try:

  --fwhm=8,0,0,0,0 is the first thing to try if the results don't look good (smooths first iteration)

  --niter=7 more iterations, takes longer. 

  --nvoxhp 2000 sample more voxels for model parameters, also takes longer. Useful for low SNR data

  --slm=linear uses second level model to improve low angular resolution data.

  For speed, you can do

  --flm=linear

  to simplify the eddy correction model. This might be OK on data with eddy reduction built into the sequence.

  To improve QC you may adjust the outlier detection thresholds,

  --ol_nvox=250 --ol_nstd=4

"

exit 1

fi

export FSLDIR=/share/apps/fsl/5.0.9-eddy-patch/

source ${FSLDIR}/etc/fslconf/fsl.sh

export PATH=${FSLDIR}/bin:$PATH

# Necessary to stop FSL's insistence on qsubbing things
unset SGE_ROOT

OMP_NUM_THREADS=1

# This defined on CFN cluster for qsub / qlogin jobs
if [[ ! -z "$NSLOTS" ]]; then
    export OMP_NUM_THREADS=$NSLOTS
fi

# Print path to executable, since version info printed by eddy_openmp is wrong in 5.0.9-eddy-patch
echo "
--- eddy command ---
"
which eddy_openmp
echo "
--------------------
"

echo "
--- Running on $HOSTNAME ---
" 

doPeas=$8
doRepol=$9

eddyOpts="--verbose"

# Default to peas
if [[ $doPeas -eq 0 ]]; then
  eddyOpts="$eddyOpts --dont_peas"
fi

# Default to no repol
if [[ $doRepol -gt 0 ]]; then
  eddyOpts="$eddyOpts --repol"
fi

cmd="eddy_openmp --imain=$1 --mask=$2 --acqp=$3 --index=$4 --bvals=$5 --bvecs=$6 --out=$7 $eddyOpts"

echo "
--- eddy call ---
$cmd
-----------------

"

$cmd
