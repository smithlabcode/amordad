#    Copyright (C) 2014 University of Southern California,
#                       Andrew D. Smith and Ehsan Behnam
#
#    Authors: Ehsan Behnam and Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

SMITHLAB_CPP = ../smithlab_cpp
CROW = ../crow/include
BOOST = ../boost_1_55_0

ifndef BOOST
$(error Must define BOOST variable)
endif

PROGS = extract_kmer_counts generate_hash_function build_metagenome_hash_table \
	build_graph  build_graph_naively exec_query angle_measurement query_server

#INCLUDEDIRS = $(BOOST)
INCLUDEDIRS += $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

CXX = g++
CXXFLAGS = -Wall -std=c++1y -rdynamic -lboost_date_time -lboost_filesystem -lboost_system -lboost_thread -lpthread 
OPTFLAGS = -O3
DEBUGFLAGS = -DNDEBUG 

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
else
CXXFLAGS += -DHAVE_OMP
LIBS += -fopenmp
endif
endif
ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all:    $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

all: $(PROGS)

extract_kmer_counts : Metagenome.o

generate_hash_function : MGHashFunction.o MGFeatureSet.o Metagenome.o

build_metagenome_hash_table : MGHashTable.o MGHashFunction.o MGFeatureSet.o Metagenome.o

build_graph : MGraph.o Metagenome.o MGHashTable.o

build_graph_naively : MGraph.o Metagenome.o

exec_query : MGraph.o Metagenome.o MGHashTable.o MGHashFunction.o

angle_measurement : Metagenome.o

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o \
	smithlab_os.o smithlab_utils.o OptionParser.o)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
