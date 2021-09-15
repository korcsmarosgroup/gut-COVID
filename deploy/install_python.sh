#!/bin/sh
add-apt-repository ppa:deadsnakes/ppa
apt-get update

apt-get -y install python3.6 python3.6-dev python3-setuptools python3-pip
add-apt-repository ppa:cran/libgit2
apt-get update
apt-get -y install libgit2-dev
update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 2
update-alternatives --set python3 /usr/bin/python3.6
pip3 install numpy
pip3 install scipy
pip3 install pandas
pip3 install matplotlib
pip3 install argparse
pip3 install pyfasta
pip3 install networkx
pip3 install pprintpp
