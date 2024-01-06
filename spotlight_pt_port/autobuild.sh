_term() {
	echo "Interrupt!"
	rm -r build
	exit 1
}

trap _term SIGTERM

# What do you want to call your build directory (Default: build). Pass it on as an argument (./autobuild [dirname]) to give a new name
bd=build
if test -n "$1"
then
	bd=$1
fi

mkdir $bd
cd $bd

loc=`pwd`
LIBTORCH_LOC="$HOME/cpp/libtorch" #TODO: Change this to match the location you downloaded the libtorch folder. Download from: 
export LIBRARY_PATH=$LIBRARY_PATH:"$loc/../dnv"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$loc/../dnv"
cmake -DCMAKE_PREFIX_PATH=$LIBTORCH_LOC ..
cmake --build . --config Release
cd ..
#cmake . --config Release
