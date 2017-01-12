# HGCalSimAnalysis

Build analysis:

$source path/to/ilcsoft/init_ilcsoft.sh

$mkdir build
$mkdir bin
$mkdir lib

$cd build
$cmake -C ${ILCSOFT}/ILCSoft.cmake ../
$make install


Run an analysis :

$cd ../
$./bin/astest toto.root
