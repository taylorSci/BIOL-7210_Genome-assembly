#!/bash/bin

#REAPR Installation
cd $toolsDir
mkdir REAPR
cd REAPR
wget ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.18.tar.gz
tar -xzf Reapr_1.0.18.tar.gz -C $toolsDir/REAPR
./configure
make
make install
