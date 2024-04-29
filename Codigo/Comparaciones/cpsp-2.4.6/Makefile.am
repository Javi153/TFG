
SUBDIRS = src misc tests

# doxygen support
include $(top_srcdir)/aminclude.am

MOSTLYCLEANFILES = $(DX_CLEANFILES)

# new doc rule to ensure that the picture is copied
doc-doxygen.hh : 
	./config.status --file=src/doxygen.hh
doc: doc-doxygen.hh src/doxygen.hh doxygen-doc 
if DX_COND_html
	cp $(top_srcdir)/misc/doxygen/CPSP-logo-100.png $(top_srcdir)/$(DX_DOCDIR)/html/.
endif

# for exporting
DOXYGENSUPPORTFILES = acinclude.m4 aminclude.am
COREDBFILES =	\
				CoreDB/Cubic/1/cub*txt \
				CoreDB/Cubic/2/cub*txt \
				CoreDB/Cubic/3/cub*txt \
				CoreDB/Cubic/4/cub*txt \
				CoreDB/Cubic/5/cub*txt \
				CoreDB/Cubic/6/cub*txt \
				CoreDB/Cubic/7/cub*txt \
				CoreDB/Cubic/8/cub*txt \
				CoreDB/Cubic/9/cub*txt \
				CoreDB/Cubic/10/cub*txt 

EXTRA_DIST = $(DOXYGENSUPPORTFILES) $(COREDBFILES) LICENSE

# test rule
TSTDIR = tests

test:	all
	@echo "-------------------------------------------------------------------"
	@echo "Run test"
	@make test -s -C $(TSTDIR)
	@echo "-------------------------------------------------------------------"
