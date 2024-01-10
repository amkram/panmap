rm -rf build && mkdir build

cd build
  cmake -Dlibpanmat_DIR=${PWD}/../../panmat-new/build ..
  make -j8
