#!/bin/sh
starver 'dev'
echo $ROOTSYS
echo $LD_LIBRARY_PATH
cwd=$(pwd)
gdb root4star <<EOF 
run 'lMuDst.C' 'RecoSimGlobal.C+'
EOF
