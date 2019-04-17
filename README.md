# dtProcessing

DTI processing scripts for retrospective data, no topup support.

## Pipeline

Assumes ANTs cortical thickness is run on T1 data

1. dtPreProc/processDTI.pl
2. connMat/dtConnMatSubj.pl


## Longitudinal pipeline

Requires cross-sectional processing and then antsLongitudinalCorticalThickness output

1. dtPreProc/processDTILong.pl
2. connMat/dtConnMatSubj.pl


