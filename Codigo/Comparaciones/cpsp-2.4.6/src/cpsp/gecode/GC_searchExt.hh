/*
 *  Main authors:
 *     Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *
 *  Contributing authors:
 *     Sebastian Will http://www.bioinf.uni-freiburg.de/~will/
 *
 *  Copyright:
 *     Martin Mann, 2007
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef GC_SEARCHEXT_HH_
#define GC_SEARCHEXT_HH_

/**
 * Symmetry Breaking Search
 *
 * implements symmetry breaking ala Backofen, Will (Backofen, Will, CP, 1999)
 *
 * contributed by Sebastian Will
 *
 */

#include "gecode/minimodel.hh"
#include "gecode/support/dynamic-array.hh"

namespace Gecode {


  template <class Int, class T>
  Int
  countDFS(T* initialSpace, unsigned int c_d, unsigned int a_d, Search::Statistics* stat, Gecode::Search::Stop* st, Int maxSols) {

        // init depth first search engine
    Gecode::DFS<T> dfsEngine(   initialSpace, c_d, a_d, st);

    Int sols = 0;

        // get first solution if available
    T* nextSol = dfsEngine.next();
        // while solution found
    while( nextSol != NULL && sols < maxSols) {

        sols++; // count for statistics

        delete nextSol;             // cleaning
        nextSol = dfsEngine.next(); // get next solution
    }
    if (stat != NULL)
       *stat = dfsEngine.statistics();
    return sols;
  }
  template <class Int, class T>
  Int
  countDDS(T* initialSpace, unsigned int c_d, unsigned int a_d, Search::Statistics* stat, Gecode::Search::Stop* st, Int maxSols) {
  	// DDS in current Gecode 1.3.1 not supported --> redirect to use DFS
  	return countDDS<Int,T>(initialSpace, c_d, a_d, stat, st, maxSols);
  }


   namespace Search {





       /**
	  \brief Symmetry Breaker

	  abstract class that is to be specialized for breaking symmetries of a certain problem
	  
	  used as argument to SBViewValBranching
	  
       */
       
       template<class View, class Val>
       class SymmBreaker {
       protected:
	   /// store of cp variables
	   BoolVarArray cps;
	   /// indexing from cps[index] to symmetry index
	   Support::DynamicArray<int> symOf;
	   
	   /// tell constraint of x,pos,val under symmetry
	   virtual
	   void tellSymm(Space *home,int symmetry, ViewArray<View> &x, int pos, Val val) const {
	       BoolVar b = reifySymm(home,symmetry,x,pos,val);
	       post(home,b==1);
	   }

	   /// reify constraint of x[pos] with val under symmetry
	   virtual
	   BoolVar reifySymm(Space *home,int symmetry, ViewArray<View> &x, int pos, Val val) const=0;

	   /// update the Cp variables
	   void
	   updateCps(Space *home, ViewArray<View> &x, int pos, Val val);
	   
	   /// post the implications that break symmetries
	   void
	   postImplications(Space *home, ViewArray<View> &x, int pos, Val val) const;
	   
       public:
	   /// constructor with number of symmetries
	   SymmBreaker(Space *home, int sym_num)
	     : cps(BoolVarArray(home,sym_num,1,1)), // all symms are valid to break
	       symOf(sym_num)
	   {
	   	   // init of the symOf array
	   	 for (int i=0; i< sym_num; i++)
	   	   symOf[i] = i;
	   }
	   
	   /// constructor with number of symmetries and initially valid symmetries
	   // generate Cp variables in space *home
	   // with initial cp for deactivating the handling
	   // of certain symmetries
	   SymmBreaker(Space *home, int sym_num, BoolVarArray &init_cp)
	       : cps(), symOf(sym_num) 
	   {
	      assert(sym_num == init_cp.size());
	         // copy init_cp
	      cps.update(home, false, init_cp);
	   	     // init of the symOf array
	   	  for (int i=0; i< sym_num; i++)
	   	    symOf[i] = i;
	   }
	   
	   /// copy constructor
	   SymmBreaker(Space *home, bool share,const SymmBreaker<View,Val> &symm) {
	   	 assert(symm.cps.size()>0);
	   	 int newSize = 1; // first will be copied for sure (== identity)
	   	 for (int i=1; i<symm.cps.size(); i++)
	   	   if (!symm.cps[i].assigned() || symm.cps[i].val() != false)
	   	     newSize++;
	   	 cps = BoolVarArray(home, newSize);
	   	 symOf = Support::DynamicArray<int>(newSize);
	   	 cps[0].update(home,share,(BoolVar&)symm.cps[0]);
	   	 symOf[0] = symm.symOf[0];
	   	 int cpI = 1;  // fill remaining cp vars and symOf entries
	   	 for (int i=1; i<symm.cps.size(); i++) {
	   	   if (!symm.cps[i].assigned() || symm.cps[i].val() != false) {
	         cps[cpI].update(home,share,(BoolVar&)symm.cps[i]);
	         symOf[cpI] = symm.symOf[i];
	         cpI++;
	   	   }
	   	 }
//	       cps.update(home,share,(BoolVarArray&)symm.cps);
	   }

	   virtual
	   ~SymmBreaker() {}

	   /// control of branching
	   /// returns true as long as SymmBreaker wants to branch
	   virtual bool
	   status(const Space *home) const {return true;}
	   
	   /// post constraints for the symmetry breaking branching in alternative a in {0,1}
	   void
	   tell(Space *home, int a, ViewArray<View> &x, int pos, Val val);

	   /// print (for debugging)
	   virtual
	   void
	   print(void) const;
       };
       

       /**
	* \brief Generic branching with symmetry breaking
	*
	*/
       template <class View, class Val, class ViewSel, class ValSel, class SymmBreaker>
       class SBViewValBranching : public ViewValBranching<View,Val,ViewSel,ValSel> {
       protected:
	   SymmBreaker symmetries; // handler for symmetry breaking
	   
	   SBViewValBranching(Space* home, bool share, SBViewValBranching& b);
       public:
	   /// Constructor for creation
	   SBViewValBranching(Space* home, ViewArray<View>& x);
	   
	   /// Constructor for creation with pre-initialized symmetries object
	   SBViewValBranching(Space* home, ViewArray<View>& x, SymmBreaker &sym);
	   
	   /// Check status of branching, return true if alternatives left
	   virtual bool status(const Space* home) const;
	   
	   /// Perform commit for branching description \a d and alternative \a a
	   virtual ExecStatus commit(Space* home, const BranchingDesc* d,
				     unsigned int a);
	   /// Perform cloning
	   virtual Actor* copy(Space* home, bool share);
       };
   }
}
/* end of Symmetry Breaking */

#include "symm-break.icc"
#endif /*GC_SEARCHEXT_HH_*/
