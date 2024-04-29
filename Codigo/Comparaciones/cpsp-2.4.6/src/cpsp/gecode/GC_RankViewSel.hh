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

#ifndef GC_RANKVIEWSEL_HH_
#define GC_RANKVIEWSEL_HH_


namespace cpsp {
  namespace gecode {

	/**
	 * A view selection for IdxIntView that selects hierachically first by the
	 * rank of the view and if equal as best so far by the given second view 
	 * selector.
	 */
  	template <class VS>
  	class GC_RankViewSel {
	protected:
		  //! best rank seen so var
		unsigned int rank;
		  //! second view selector in case of equal rank of two views
		VS vs;
	public:
		  //! inits the selector
		Gecode::ViewSelStatus
		init(const Gecode::Space* home, const IdxIntView&x);

		  //! does the selection and returns the selection status for the 
		  //! given view
		Gecode::ViewSelStatus
		select(const Gecode::Space* home, IdxIntView x);
  	};

  }  // namespace gecode
}  // namespace gecode


#include "GC_RankViewSel.icc"

#endif /*GC_RANKVIEWSEL_HH_*/
