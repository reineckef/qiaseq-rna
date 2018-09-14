set -e

# Wrapper script to compile read-trimmer cython function and QIMERA's folder structure setup

# assuming the code directory is mounted at : /home/qiauser/code on the container

# trimmer
cd /home/qiauser/code/read-trimmer/
python setup.py build_ext --inplace
# qimera
cd /home/qiauser/code/QIMERA
perl Makefile.PL
make
make install
