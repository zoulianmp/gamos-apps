#!/bin/bash

workpath=/home/zoulian/Gamos/GamosSim

apppath=$workpath/tutorials/Primer

cd $apppath
pwd
appname='exercise0a'
echo "gamos $appname.in 2>&1 | tee out.$appname"

gamos $apppath/$appname.in 2>&1 | tee out.$appname

mv gamos.log ./gamos.$appname.log
mv gamos_error.log ./gamos_error.$appname.log

echo @%@%@%@%@%@%@%@%@%@%@
echo End of the gamos run

echo "**********************************"
echo "If you want to see the root's files, please run the command:"
echo "python /home/zoulian/Gamos/GamosSim/rootbrowser.py"
echo "**********************************"

