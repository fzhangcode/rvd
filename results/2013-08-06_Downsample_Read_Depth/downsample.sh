#!/bin/sh

dRate=$1

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

for f in $InputPath/*.bam
do
	filename=${f##*/}
	Input=$f
	Output=${OutputPath%%/}/$filename
	echo $Output

	if [ -f $Output ]
		then
			echo File $filename exists already
		else
			echo Downsampling $filename
			java -Xmx2g -jar $DownsamplePath INPUT=$Input OUTPUT=$Output RANDOM_SEED=null PROBABILITY=$dRate VALIDATION_STRINGENCY=SILENT
		fi
done
