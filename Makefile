CXXFLAGS	= -ansi -Wall -pedantic -O3

all: wavelets

wavelets: wavelets.tgz
	tar -xzf wavelets.tgz
	mv wavelets wavelets_src
	cd wavelets_src/src && $(CXX) -o ../../wavelets $(CXXFLAGS) WaveletApp.cpp

wavelets.tgz:
	wget http://staff.washington.edu/dbp/WMTSA/NEPH/wavelets.tgz

clean:
	rm -f wavelets.tgz
	rm -rf wavelets/
	rm -rf wavelets_src
	rm -f wavelets
