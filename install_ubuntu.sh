#!/bin/bash
#
# NOTE: Please use this only with Ubuntu 22.04 +. It has not been tested on older versions.

spotlight_dir=`pwd`

# Requirements
apt update && apt install -y vim curl wget g++ gcc cmake locate man libomp-dev zip unzip software-properties-common gcc-multilib bash python3 libzmq5-dev git git-lfs libstdc++6
wget "https://drive.usercontent.google.com/download?id=1Psk2yM3NwvUaG2v_XW2w_-t77Yfm5xIe&export=download&confirm=t&uuid=8e86dd1f-e040-45fb-8e1d-f184ed4ff873" -O libtorch.zip && unzip libtorch.zip

# Going to the build folder
cd spotlight_pt_port

# Setting up paths for an auto-build
sed -e '0,/LIBTORCH_LOC/ {/LIBTORCH_LOC/d}' autobuild.sh > ab
echo "LIBTORCH_LOC=/usr/share/libtorch" > autobuild.sh
cat ab >> autobuild.sh
mkdir -p $spotlight_dir/install/models/default

# Putting model data in an accessible location
cp ../spotlight_data/model/*.pt $spotlight_dir/install/models/default
chmod +x setpaths.sh && ./setpaths.sh $spotlight_dir/install/models/default

# Building
chmod +x autobuild.sh && ./autobuild.sh
chmod +x ../scripts/*
cd build

# Post compile
read -p "Add build directory to path? (Y/N) " yn
yn=`echo $yn|cut -c 1`
if [ "$yn" == 'Y' -o "$yn" == 'y' ]; then
	bd=`pwd`
	echo "PATH=\$PATH:$bd" >> $HOME/.bashrc
	cp ../../scripts/* .
fi
