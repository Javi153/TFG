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

#include "cpsp/gecode/GC_ValNearCenter.hh"

#include <gecode/iter.hh>

namespace cpsp
{
  namespace gecode
  {


	const biu::LatticeFrame* GC_ValNearCenter::latFrame = NULL;

	int	
	GC_ValNearCenter::val(const Gecode::Space* home, Gecode::Int::IntView x) const
	{
		using namespace Gecode::Int;
		using namespace Gecode::Iter::Ranges;
		
		int retVal = x.min();
		
			// init iterator 
		ToValues< ViewRanges< IntView > > it;
		ViewRanges< IntView >  r(x);
		it.init( r );

			// init helper variables		
		biu::IntPoint center = latFrame->getCenter();
		double actDist = DBL_MAX, minDist = actDist;
			// check values
		unsigned int c = 0;
		while ( it() ) {	// still valid
			c++;
				// get current core-center distance
			actDist = latFrame->getPoint(it.val()).distance(center);
				// check if new minimum
			if (actDist < minDist) {
				minDist = actDist;
				retVal = it.val();
			}
				// increase iterator 
			++it;
		}
		assert(c == x.size());
		
		return retVal;
	}
	
	Gecode::ModEvent
	GC_ValNearCenter::tell(	Gecode::Space* home, unsigned int a, 
							Gecode::Int::IntView x, int n)
	{
		assert((a == 0) || (a == 1));
		return (a == 0) ? x.eq(home,n) : x.nq(home,n);
	}

  } //namespace gecode
} // namespace cpsp
