#!/bin/bash

# Script to set env variables to use updated Geant4Reweight
# Useage: source set_env_variables_for_g4rw.sh
# Author: C Thorpe

# Find the installation of larsim in the tarball
echo "Finding larsim installation"
echo LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH | tr : \\n
LARSIM_DIR=$(echo $LD_LIBRARY_PATH | tr : \\n | grep larsim)
echo $LARSIM_DIR
echo $LARSIM_DIR > LARSIM_DIR.txt

# Remove everything after the word larsim from the address and replace with g4rw dir
TARBALL_PATH=$(sed 's@/larsim.*@@g' LARSIM_DIR.txt)
echo $TARBALL_PATH > TARBALL_PATH.txt

echo "Modifying LD_LIBRARY_PATH"
LIBRARY=$TARBALL_PATH/geant4reweight/v01_08_01/slf7.x86_64.e17.prof/lib
echo $LIBRARY
export LD_LIBRARY_PATH=$LIBRARY:$LD_LIBRARY_PATH
GEANT4RW_PATH=$TARBALL_PATH/geant4reweight/v01_08_01/slf7.x86_64.e17.prof/bin
export PATH=$GEANT4RW_PATH:$PATH

echo "Setting GEANT4RW_DATA_DIR"
export GEANT4RW_DATA_DIR=$TARBALL_PATH/uboonedata/v08_00_00_49/systematics/reint
echo $GEANT4RW_DATA_DIR
