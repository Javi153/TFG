
#####################################

bin_PROGRAMS =	HPconvert HPcompress HPstruct HPnnet HPdesign \
				HPview HPrand HPdeg HPoptdeg HPsurface \
				HPviewSC HPstructSC HPoptdegSC \
				HPstructCG HPseq HPrep HPrepSC #HPshapes

INCLUDES = -I$(top_srcdir)/src

#####################################
# HPdesign tool
#####################################

HPdesign_SOURCES =	HPdesign.hh HPdesign.cc

HPdesign_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPstruct tool
#####################################

HPstruct_SOURCES =	HPstruct.hh HPstruct.cc

HPstruct_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

HPstructSC_SOURCES = side_chain_threading.cpp

HPstructSC_LDADD =  ../cpsp/SideChain/libcpspSC.a ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPdeg tool
#####################################

HPdeg_SOURCES =	HPdeg.hh HPdeg.cc

HPdeg_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPshapes tool
#####################################

#HPshapes_SOURCES =	HPshapes.hh HPshapes.cc

#HPshapes_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPrand tool
#####################################

HPrand_SOURCES =	HPrand.hh HPrand.cc

HPrand_LDADD =  ../cpsp/libcpsp.a

#####################################
# HPconnect tool
#####################################

HPnnet_SOURCES =	HPnnet.hh HPnnet.cc

HPnnet_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPcompress tool
#####################################

HPcompress_SOURCES =	HPcompress.hh HPcompress.cc

HPcompress_LDADD =  ../cpsp/libcpsp.a

#####################################
# HPconvert tool
#####################################

HPconvert_SOURCES =	HPconvert.hh HPconvert.cc

HPconvert_LDADD =  ../cpsp/libcpsp.a

#####################################
# HPview tool
#####################################

HPview_SOURCES =	HPview.hh HPview.cc

HPview_LDADD =  ../cpsp/libcpsp.a

HPviewSC_SOURCES =	scbb_viewer.h scbb_viewer.cpp

HPviewSC_LDADD =  ../cpsp/libcpsp.a

#####################################
# HPoptdeg tool
#####################################

HPoptdeg_SOURCES =	HPoptdeg.hh HPoptdeg.cc

HPoptdeg_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

HPoptdegSC_SOURCES =	side_chain_optdeg.h side_chain_optdeg.cpp

HPoptdegSC_LDADD =  ../cpsp/SideChain/libcpspSC.a ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
# HPsurface tool
#####################################

HPsurface_SOURCES =	HPsurface.hh HPsurface.cc

HPsurface_LDADD =  ../cpsp/libcpsp.a

#####################################
# HPstructCG tool
#####################################

HPstructCG_SOURCES =	HPstructCG.hh HPstructCG.cc

HPstructCG_LDADD =  ../cpsp/libcpsp.a


#####################################
# HPseq tool
#####################################

HPseq_SOURCES =	HPseq.cc

#####################################
# HPrep tool
#####################################

HPrep_SOURCES =	HPrep.hh HPrep.cc

HPrep_LDADD =  ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

HPrepSC_SOURCES =	HPrepSC.cpp

HPrepSC_LDADD =  ../cpsp/SideChain/libcpspSC.a ../cpsp/gecode/libcpspgecode.a ../cpsp/libcpsp.a

#####################################
