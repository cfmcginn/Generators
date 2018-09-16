CXX = g++
CXXFLAGS = -Wall -Werror -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ifndef GENERATORSPATH
$(error GENERATORSPATH is not set at all. Do 'source setGeneratorsPath.sh')
endif


INCLUDE=-I $(PWD) -I $(GENERATORSPATH)
ROOT=$(shell root-config --cflags --glibs)
#PYTHIA8=$(shell pythia8-config --cflags --libs --lhapdf6)
MOREPYTHIA8=-I$(PYTHIA8PATH)/include -O2 -pedantic -W -Wall -Wshadow -fPIC -L$(PYTHIA8PATH)/lib -Wl,-rpath,$(PYTHIA8PATH)/lib -lpythia8 -ldl
#^based on make of examples
LHAPDF6=-I/home/cfmcginn/Packages/Generators/LHAPDF6/LHAPDF6-build/include -L/home/cfmcginn/Packages/Generators/LHAPDF6/LHAPDF6-build/lib -lLHAPDF
FASTJET=$(shell $(FASTJETPATH)/bin/fastjet-config --cxxflags --libs --plugins)


MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdf mkdirOutput bin/forestToGen.exe bin/pythia8CUETP8M2T4.exe bin/statisticalComparison.exe bin/removeDuplicateTNamed.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

bin/pythia8CUETP8M2T4.exe: src/pythia8CUETP8M2T4.C
	$(CXX) src/pythia8CUETP8M2T4.C -o bin/pythia8CUETP8M2T4.exe $(CXXFLAGS) $(ROOT) $(FASTJET) $(MOREPYTHIA8) $(LHAPDF6) $(INCLUDE) 

bin/statisticalComparison.exe: src/statisticalComparison.C
	$(CXX) src/statisticalComparison.C -o bin/statisticalComparison.exe $(CXXFLAGS) $(ROOT) $(INCLUDE) 

bin/removeDuplicateTNamed.exe: src/removeDuplicateTNamed.C
	$(CXX) src/removeDuplicateTNamed.C -o bin/removeDuplicateTNamed.exe $(CXXFLAGS) $(ROOT) $(INCLUDE) 

bin/forestToGen.exe: src/forestToGen.C
	$(CXX) src/forestToGen.C -o bin/forestToGen.exe $(CXXFLAGS) $(ROOT) $(INCLUDE) 

clean:
	rm -f ./*~
	rm -f src/*~
	rm -f include/*~
	rm -f configs/*~
	rm -f ./#*#
	rm -f src/#*#
	rm -f include/#*#
	rm -f configs/#*#
	rm -f bin/*.exe
	rm -rf bin
