
#####################################
# cpsp gecode library
#####################################

lib_cpspSC_HEADERS =	\
						options/chain_option.h \
						options/protein_like_option.h \
						options/side_chain_optdeg_option.h \
						chain_threading.h \
						chain_threading_all.h

lib_cpspSCdir = $(includedir)/cpsp/gecode

lib_LIBRARIES = libcpspSC.a

libcpspSC_a_SOURCES =	\
						options/chain_option.cpp \
						options/protein_like_option.cpp \
						options/side_chain_optdeg_option.cpp \
						chain_threading.cpp \
						chain_threading_all.cpp

libcpspSC_a_CXXFLAGS = -I$(top_srcdir)/src

EXTRA_DIST = $(lib_cpspSC_HEADERS)