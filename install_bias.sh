#!/bin/bash
mkdir third-party
cd third-party
wget https://www.tuhh.de/ti3/keil/profil/Profil-2.0.8.tgz
tar -xzf Profil-2.0.8.tgz
cd Profil-2.0.8
sed -i '42s/.*/\/\/#define __BIASRAISEDIVIDEBYZERO__/' config/x86-64-Linux-compat-gcc/BiasRnd.h

./configure
make all
make install
