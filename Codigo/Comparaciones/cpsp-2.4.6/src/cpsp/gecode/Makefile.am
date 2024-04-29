
#####################################
# cpsp gecode library
#####################################

lib_cpspgecode_HEADERS =	GC_HPThreading.hh \
							GC_LatticeNeighbored2.hh \
							GC_StlSetRangeIterator.hh \
							GC_ValNearCenter.hh \
							GC_ThreadingSuperSpace.hh \
							GC_ThreadingSpace.hh \
							GC_Search.hh \
							GC_ThreadingSymmBreaker.hh \
							GC_RankViewSel.hh \
							GC_RankViewSel.icc \
							GC_searchExt.hh \
							GC_FirstSolBranching.hh \
							GC_FirstSolBranching.icc \
							symm-break.icc \
							gecodeExtensions.hh

lib_cpspgecodedir = $(includedir)/cpsp/gecode

lib_LIBRARIES = libcpspgecode.a

libcpspgecode_a_SOURCES =	GC_HPThreading.cc \
							GC_LatticeNeighbored2.cc \
							GC_StlSetRangeIterator.cc \
							GC_ValNearCenter.cc \
							GC_ThreadingSpace.cc

libcpspgecode_a_CXXFLAGS = -I$(top_srcdir)/src

EXTRA_DIST = $(lib_cpspgecode_HEADERS)