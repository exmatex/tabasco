#!/bin/bash -e

j="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || j=0
((j++))

#module load mpi/openmpi-1.8.4-gcc_4.9.1
#module load mpi/openmpi-1.6.4-gcc.4.7.2
##module load compilers/intel/14.0.2
##module load mpi/intelmpi-4.1.3.045
#
PV=6.6.1
P="charm-${PV}"
URL=http://charm.cs.uiuc.edu/distrib/${P}.tar.gz
A="${URL##*/}"
if [[ ! -f $A ]]
then
    wget -O "$A" "$URL"
    tar -xvf "${A}"
fi
D=${PWD}/install

#cd "${A%%.tar.*}"
cd charm
./build charm++ mpi-linux-amd64 smp --with-numa -j${j} --build-shared

cd mpi-linux-amd64-smp/tests/charm++/simplearrayhello
make
./charmrun hello

cd -
#install {bin,include,lib}/*
#mind the symlinks
mkdir -p ${D}
for i in bin include lib; do
  mkdir -p ${D}/$i
  for j in $i/*; do
    cp $(readlink -f $j) ${D}/$i/
  done
done
