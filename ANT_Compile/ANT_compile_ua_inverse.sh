#!/bin/bash

# USER VARIABLES
# Path to correct version of MATLAB
# Setup below has been tested on ARCHER2 for MATLAB v2021a 
MATLAB_PATH=/usr/local/MATLAB/R2023b
# Path to Ua_ANT repository
REPO_DIR=/mnt/md0/Ua/cases/ANT
# Path to Ua build directory (will be created if it doesn't exist)
UA_BUILD=./UaBuild
# Path to Ua_ANT wrapper files
UA_WRAPPER_FILES=$REPO_DIR/ANT_Inverse
# Path to configuration-specific Ua files to overwrite
UA_CASE_UPDATES=$REPO_DIR/ANT_Inverse/ANT_Inverse_9999
# Path to Ua source directory (default use the one inside UaMITgcm)
UA_SOURCE=/mnt/md0/Ua/UaSource_beta

if [ -e $UA_BUILD ]; then
    # Empty the directory
    rm -rf $UA_BUILD/*
else
    # Create the directory
    mkdir $UA_BUILD
fi

# Copy all Matlab files from UaSource
cp $UA_SOURCE/*.m $UA_BUILD
# Need to collapse a couple of subdirectories for more Matlab files
cp `find $UA_SOURCE/UaUtilities/ -name "*.m"` $UA_BUILD
cp `find $UA_SOURCE/NewestVertexBisection/ -name "*.m"` $UA_BUILD
# Copy mesh2d files
cp `find $UA_SOURCE/Mesh2d/ -name "*.m" ! -name 'inpoly2*'` $UA_BUILD
# Also copy everything from updates folders
cp -r $REPO_DIR/*.m $UA_BUILD
cp -r $UA_WRAPPER_FILES/*.m $UA_BUILD
cp -r $UA_CASE_UPDATES/* $UA_BUILD

# Create the executable
$MATLAB_PATH/bin/mcc -m $UA_BUILD/ANT_MatlabWrapper.m -o Ua -d $UA_BUILD
# Copy just the executable (not the auto-generated run script as we have a custom one) to the current directory
cp $UA_BUILD/Ua ./
echo 'Now copy "Ua" to the Ua executable directory on the server where you will run the model.'
rm -rf $UA_BUILD
