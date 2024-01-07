_term() {
	echo "Interrupt!"
	rm -r build
	exit 1
}

trap _term SIGTERM

LIBTORCH_LOC="$HOME/cpp/libtorch" #TODO: Change this to match the location you downloaded the libtorch folder. Download from: 


# What do you want to call your build directory (Default: build). Pass it on as an argument (./autobuild [dirname]) to give a new name
bd=build
if test -n "$1"
then
	bd=$1
fi

# Setting the path for dnv
loc=`pwd`
cd $loc/../dnv/
DNV_ROOT=`pwd`
cd support
echo $DNV_ROOT
sed -e 's|dnvpath|'"$DNV_ROOT"'|g' local.h > temp
cp local.h local.h.bak
mv temp local.h

echo "Path is set"

mkdir $bd
cd $bd

export LIBRARY_PATH=$LIBRARY_PATH:"$loc/../dnv"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$loc/../dnv"
cmake -DCMAKE_PREFIX_PATH=$LIBTORCH_LOC ..
cmake --build . --config Release #-- -j 4 # (You can parallelize the compilation if needed. -j gives number of parallelly running threads)
cd ..
#cmake . --config Release
