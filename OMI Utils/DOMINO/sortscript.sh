#!/bin/bash
# Will sort OMNO2 files from the "download_staging" folder to the appropriate
# year/month directories


# Assuming this resides in a staging directory and we want the root directory
# containing the year folders to be one level up (i.e. the location of this 
# script is a folder within the root folder)
cd `dirname $0`
__mydir=`pwd -P`

indir=$__mydir
outdir="$__mydir/.."


while [[ $# > 0 ]]; do
    key="$1"
    case "$key" in
        --indir*)
            indir="${key#*=}"
            ;;
        --outdir*)
            outdir="${key#*=}"
            ;;
        *)
            echo "Option $key not recognized"
            ;;
    esac
    shift
done

cd "$outdir"
for f in $(ls "$indir"/*.he5)
do
    fname=$(basename $f)
    y=${fname:21:4}
    m=${fname:26:2}
    
    if [[ ! -d ${y}/${m} ]]
    then
        mkdir -p ${y}/${m}
    fi
    
    mv -v $f ${y}/${m}/
done

exit 0
