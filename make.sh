cd build
wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
rm -rf htslib-1.19.1
bzip2 -d htslib-1.19.1.tar.bz2
tar -xvf htslib-1.19.1.tar
cd htslib-1.19.1
./configure --without-libdeflate --disable-bz2 --disable-largefile --disable-lzma --disable-libcurl
make -j 100
cp libhts.a ../
cd ..