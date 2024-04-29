

TESTSOURCE = \
				testHCore.cc \
				testHCoreDatabaseFILE.cc

TESTOBJS=$(TESTSOURCE:%.cc=%.o)

TESTEXECS = $(TESTSOURCE:%.cc=%.test)


DIFFS=$(TESTSOURCE:.cc=.diff)

TESTDEPENDFILE = .depend

EXTRA_DIST = $(TESTSOURCE) $(TESTSOURCE:.cc=.verified)

LDFLAGS =	$(GLOBALLDFLAGS) -L$(top_srcdir)/src/cpsp -L$(top_srcdir)/src/cpsp/gecode/ \
			-lcpspgecode -lcpsp -lbiu
CPPFLAGS = $(GLOBALCPPFLAGS) -I$(top_srcdir)/src


#################################################################

test: $(TESTEXECS) $(DIFFS)

$(TESTEXECS): $(TESTDEPENDFILE) $(TESTOBJS) 
	@$(CXX) $*.o $(LDFLAGS) -o $@


.test.diff: 
	@./$< > $*.lastout
	@name=$<; \
	name=$${name/test/}  ;\
	if [ -e $*.verified ] ; then \
		if diff -w $*.verified $*.lastout >$@ ; then \
			printf " Testing %-40s OK\n" $$name; \
		else \
			printf " Testing %-40s FAILED\n" $$name;\
		    echo "----------------------------------------" ;\
		    echo "Diff for $$name (verified '<', current '>')" ;\
	  		cat $@ ;\
			echo "----------------------------------------" ;\
			echo "[ If new output is correct: mv $*.lastout $*.verified ]" ;\
		fi ;\
		$(RM) $@ ;\
	else \
	    echo "----------------------------------------" ;\
	    echo "Output of $$name" ;\
		cat $*.lastout ;\
	    echo "----------------------------------------" ;\
		echo "WARNING: verified output of $$name is missing." ;\
	    echo "Please generate a file $*.verified that" ;\
	    echo "contains the verfied output of $$name."; \
		echo "[ If output is correct: mv $*.lastout $*.verified ]" ;\
	fi

$(TESTDEPENDFILE):
		@$(CXX) $(INCLUDES) $(CPPFLAGS) -MM   $(TESTSOURCE) > $(TESTDEPENDFILE)

-include $(TESTDEPENDFILE)

.cc.o:  
	@$(CXX) $(INCLUDES) $(CPPFLAGS) -c $< 



clean-local:
	$(RM) -rf $(TESTOBJS) $(TESTEXECS) $(TESTDEPENDFILE) $(TESTEXECS:=.lastout)\
		$(EXECS:=.diff) *~ $(TESTEXECS:=.exe)
