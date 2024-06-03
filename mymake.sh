rm -rf build && mkdir build

cd build
  wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
  bzip2 -d htslib-1.19.1.tar.bz2
  tar -xvf htslib-1.19.1.tar
  mkdir htslib_install
  cd htslib-1.19.1
  ./configure --disable-lzma --prefix=${PWD}/../htslib_install
  make
  make install
  cd ..

  cmake -Dlibpanmat_DIR=${PWD}/../../pangenome-mat/build -DTBB_DIR=${PWD}/../../pangenome-mat/build/oneTBB-2019_U9 ..
  make -j8
