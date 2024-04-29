
lib_cpsp_HEADERS =	Exception.hh \
					HCore.hh \
					HCoreDatabase.hh \
					HCoreDatabaseFILE.hh \
					HPThreadingOptions.hh \
					HPThreadingHandler.hh \
					PDBsupport.hh

lib_cpspdir = $(includedir)/cpsp

lib_LIBRARIES = libcpsp.a

libcpsp_a_SOURCES =	Exception.cc \
					HCore.cc \
					HCoreDatabaseFILE.cc \
					HPThreadingOptions.cc

libcpsp_a_CXXFLAGS = -I$(top_srcdir)/src

EXTRA_DIST = $(lib_cpsp_HEADERS)

#####################################

SUBDIRS = . gecode SideChain #ilog


