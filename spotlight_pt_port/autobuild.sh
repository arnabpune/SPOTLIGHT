_term() {
	echo "Interrupt!"
	rm -r build
	exit 1
}

trap _term SIGTERM

bd=build
if test -n "$1"
then
	bd=$1
fi
mkdir $bd
cd $bd
export LIBRARY_PATH=$LIBRARY_PATH:"/home/venkata/dnv"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/venkata/dnv"
cmake -DCMAKE_PREFIX_PATH=$HOME/cpp/libtorch ..
cmake --build . --config Release -- -j 4
cd ..
#cmake . --config Release
