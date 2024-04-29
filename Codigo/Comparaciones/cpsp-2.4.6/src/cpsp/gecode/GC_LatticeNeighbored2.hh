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

#ifndef GC_LATTICENEIGHBORED2_HH_
#define GC_LATTICENEIGHBORED2_HH_

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>

#include <biu/LatticeFrame.hh>

namespace cpsp
{
  namespace gecode
  {

  	typedef Gecode::BinaryPropagator<	Gecode::Int::IntView,
										Gecode::Int::PC_INT_DOM> GC_BinProp;

	class GC_LatticeNeighbored2 : public GC_BinProp

	{
	private:

		//TODO rausnehmen wenn moeglich
		biu::LatticeFrame::index_set& getIndexedNeighbors(
			const biu::LatticeFrame::index_type center, biu::LatticeFrame::index_set& toFill) const ;

			//! Removes all elements from x that are not neighbored to center.
		Gecode::ModEvent removeNonNeighbors(	Gecode::Space* home,
												Gecode::Int::IntView& x,
												int center);
			//! Removes all elements from x0/x1 that are not neighbored to
			//! an element in x1/x0
		Gecode::ModEvent removeNonNeighbors(Gecode::Space* home);

			//! the initial maxPropSize value
		static const unsigned int MAXPROPSIZEINIT;


	protected:
			// variablen von superklasse uebernehmen
		using GC_BinProp::x0;
		using GC_BinProp::x1;

			//! the indexed neighborhood of lattice for fast access
		const biu::LatticeFrame::index_set	*			neighborhood;

			//! Consistency level for neighbor propagation.
			//! ICL_VAL = only if one domain is bound
			//! ICL_BND = only propagating if one domain size is lower than maxPropSize
			//! ICL_DOM = full propagation
		Gecode::IntConLevel 			conLevel;

			//! the maximal size for the minimal domain so that the neighbor
			//! constraint is propagated
		unsigned int 					maxPropSize;

			//! if true, two successive domains have to be smaller or equal
			//! than maxPropSize so that the neighbor constraint is propagated
		bool							maxPropSizeBin;

			//! the minimal size for the maximal domain so that the neighbor
			//! constraint is propagated
		unsigned int					minPropSize;

			//! cancels propagation during next propagate() call
		bool							cancelProp;

			/// Constructor for cloning \a p
		GC_LatticeNeighbored2(Gecode::Space* home, bool share, GC_LatticeNeighbored2& p);

			/// Constructor for posting \a p
		GC_LatticeNeighbored2(	Gecode::Space* home,
						Gecode::Int::IntView x0,
						Gecode::Int::IntView x1,
  						const biu::LatticeFrame::index_set * indexedNeighVecs,
						Gecode::IntConLevel conLvl,
						unsigned int maxPrSize,
						bool maxPrSizeBin,
						unsigned int minPrSize);

	public:
			/// post a binary neighbor constraint
		static Gecode::ExecStatus post(	Gecode::Space* home,
									Gecode::Int::IntView x0,
									Gecode::Int::IntView x1,
			  						const biu::LatticeFrame::index_set * indexedNeighVecs,
			  						Gecode::IntConLevel conLvl =Gecode::ICL_BND,
			  						unsigned int maxPrSize = MAXPROPSIZEINIT,
									bool maxPrSizeBin = false,
									unsigned int minPrSize = 0);

	    	/// Copy propagator during cloning
	    virtual Gecode::Actor* copy(Gecode::Space* home, bool share);

		    /// Cost function
	    virtual Gecode::PropCost cost(void) const;

		    /// Perform propagation
	    virtual Gecode::ExecStatus propagate(Gecode::Space* home);

	};

  } // namespace gecode
} // namespace cpsp

#endif /*GC_LATTICENEIGHBORED2_HH_*/

