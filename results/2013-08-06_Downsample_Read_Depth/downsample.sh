#!/bin/sh
clear
echo "Run Picard to downsample syntehtic Bam files"
echo "Please input desired downsampling rate"
read dRate
echo $dRate

InputPath=../../data/Synthetic_BAM_files/
OutputPath=../../data/Synthetic/
OutputFolder=$(printf "%.4f" $dRate)
OutputFolder=${OutputFolder/./_}
OutputPath=${OutputPath%%/}/$OutputFolder
echo Output Path
echo $OutputPath
mkdir -p $OutputPath
DownsamplePath=picard-tools-1.96/DownsampleSam.jar

echo ----------------------------------------
echo Start Downsampling

for f in `find $InputPath -maxdepth 1 -name \*.bam -type f -printf "%f\n"`
do
Input=${InputPath%%/}/$f
Output=${OutputPath%%/}/$f
echo $Output

if [ -f $Output ]
then
echo File $f exists already
else
echo File $f does not exist
java -Xmx2g -jar $DownsamplePath INPUT=$Input OUTPUT=$Output RANDOM_SEED=null PROBABILITY=$dRate VALIDATION_STRINGENCY=SILENT
fi

done
