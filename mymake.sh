rm -rf build && mkdir build

cd build
  cmake -Dlibpanmat_DIR=${PWD}/../../pangenome-mat/build ..
  make -j8
