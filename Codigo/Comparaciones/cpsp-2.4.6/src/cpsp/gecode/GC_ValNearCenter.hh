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

#ifndef GC_VALNEARCENTER_HH_
#define GC_VALNEARCENTER_HH_

#include "gecode/int.hh"
#include "biu/LatticeFrame.hh"

namespace cpsp
{
  namespace gecode
  {

	class GC_ValNearCenter
	{
	public:
			//! The lattice frame the indexing of the variables
			//! is based on, has to be static due to the generic template
			//! structure of the Branchings.
		static const biu::LatticeFrame* latFrame;

			//! returns the value of x next to the center of latFrame
		int val(const Gecode::Space* home, Gecode::Int::IntView x) const;
		
			//! apply the value n to the view x
		Gecode::ModEvent 
		  tell(	Gecode::Space* home, unsigned int a, 
				Gecode::Int::IntView x, int n);
	};

  } // namespace gecode
} // namespace cpsp

#endif /*GC_VALNEARCENTER_HH_*/
