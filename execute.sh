if [ ! -d build ]; then
 mkdir build
fi
if [ ! -d resulats ]; then
 mkdir resultats
fi
cd build
cmake .. && make
