#!/bin/bash

workpath=/home/zoulian/Gamos/GamosSim

apppath=$workpath/DICOMPhantom

cd $apppath
pwd

export DICOM_MERGE_ZSLICES=1


echo "Running buildG4DICOM "

buildG4DICOM

echo "THE DICOM",$DICOM_MERGE_ZSLICES

echo @%@%@%@%@%@%@%@%@%@%@
echo End of the gamos run

