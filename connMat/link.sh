# Run from output

for subj in `ls ../../../pipedream2018/crossSectional/dti/`; do 
  for tp in `ls ../../../pipedream2018/crossSectional/dti/$subj`; do 
    mkdir -p $subj/$tp
    cd $subj/$tp 
    for dir in distCorr dt dtNorm; do 
      ln -s /data/grossman/pipedream2018/crossSectional/dti/$subj/$tp/$dir .
    done
    cd -
  done
done
