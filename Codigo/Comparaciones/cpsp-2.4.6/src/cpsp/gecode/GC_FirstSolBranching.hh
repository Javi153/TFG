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

#ifndef GC_FIRSTSOLBRANCHING_HH_
#define GC_FIRSTSOLBRANCHING_HH_

#include <gecode/kernel.hh>


namespace cpsp
{
  namespace gecode
  {
  	
  	int solutionLvl;

	/**
	 * A ViewValBranching that stops search after the first solution have been
	 * found.
	 */
   template <class View, class Val, class ViewSel, class ValSel>
   class FirstSolBranching : public Gecode::ViewValBranching<View,Val,ViewSel,ValSel> {
   protected:
   		 /// tag if the first solution was already found
	   mutable bool firstSolFound; 
	   
	     /// branching level where a solution was found
//	   static int solutionLvl;
	     /// level of this branching in the search tree
	   mutable int myLvl;
	   
	   using Gecode::ViewValBranching<View,Val,ViewSel,ValSel>::x;
	   
	   FirstSolBranching(Gecode::Space* home, bool share, FirstSolBranching& b);
   public:
	   /// Constructor for creation
	   FirstSolBranching(Gecode::Space* home, Gecode::ViewArray<View>& x);
	   
	   /// Check status of branching, return true if alternatives left
	   virtual bool 
	   status(const Gecode::Space* home) const;
	   
	   /// Perform cloning
	   virtual Gecode::Actor* 
	   copy(Gecode::Space* home, bool share);
	   
	   /// Perform branching
	   virtual Gecode::ExecStatus 
	   commit(Gecode::Space* home, const Gecode::BranchingDesc* d, unsigned int a);
   };
   
#include "cpsp/gecode/GC_FirstSolBranching.icc"

  } // namespace gecode
} // namespace cpsp


#endif /*GC_FIRSTSOLBRANCHING_HH_*/
