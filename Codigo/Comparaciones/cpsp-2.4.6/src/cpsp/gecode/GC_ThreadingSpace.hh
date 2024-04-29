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

#ifndef GC_THREADINGSPACE_HH_
#define GC_THREADINGSPACE_HH_

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/minimodel.hh>

#include <string>
#include <vector>

#include "cpsp/HCore.hh"
#include "cpsp/HPThreadingOptions.hh"
#include "cpsp/gecode/GC_ThreadingSymmBreaker.hh"
#include "cpsp/gecode/GC_StlSetRangeIterator.hh"

#include "cpsp/gecode/gecodeExtensions.hh"

#include "cpsp/gecode/GC_ThreadingSuperSpace.hh"

namespace cpsp
{
  namespace gecode
  {


  class IntDomSize {
  public:
	  unsigned int operator()(const Gecode::Int::IntView& i) const {
	    return i.size();
	  }
  };


	//! BranchingType is used to guide the branchings stated in the 
	//! ThreadingSpace constructor.
	enum BranchingType { 
			BR_NONE = 0,	//!< no branching is posted
			BR_DFS = 1,		//!< DFS ViewValbranching on all variables is posted
	//  BR_DDS  ==>  currently not supported  ==> no branching posted 
			BR_DDS = 2,		//!< DFS ViewValbranching on all singlet or core border elements, DDS ViewValbranching is posted on all remaining vars ==>  currently not supported  ==> no branching posted
			BR_SYM = 4 		//!< symmetrie breaking is posted as the first branching (not sufficient to solve the problem alone) 
		};


	//! CSP model to predict optimal structures in the HP model following the
	//! CPSP approach.
	class GC_ThreadingSpace : public SuperSpace
	{
	public:
			//! sequence features
		enum SeqFeatures{ P_SINGLET, H_SINGLET };

		typedef std::map< SeqFeatures, std::vector <unsigned int> > SeqFeatureMap;

		typedef std::vector<unsigned int> HullLevel;
		
			//! neighborhood propagation level
		enum NeighPropLvl{ 
				CUSTOM_2,	//!< CUSTOM_2 : using an advanced binary propagator
				NO_PROP 	//!< NO_PROP  : post no neighboring propagator 
			};
		
		
	protected:

			//! the position domains for the sequence elements
		Gecode::IntVarArray domains;
		
			//! ranks for the view selection ordering (returned by getRank(..))
		std::vector<unsigned int> ranks;

	public:
		GC_ThreadingSpace(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType = BR_DFS|BR_SYM,
			NeighPropLvl neighPropLvl = CUSTOM_2);
		GC_ThreadingSpace(bool share, GC_ThreadingSpace& toCopy);
		virtual ~GC_ThreadingSpace();

		virtual SuperSpace* copy(bool share);

			//! Returns a vector of indexed points thats a structure for
			//! the current CPSP.
		virtual biu::LatticeFrame::index_vec getSolution() const ;
 		virtual void print() const;

		virtual void handleSolution(GC_ThreadingSpace* solution)
		{}
		
		 /*! returns the rank of the variable with the given index
		  * 
		  * @param index of the variable of interest
		  * @return its rank
		  */ 
		virtual
		unsigned int
		getRank(const int index) const;

/*
	*/

    virtual
    int getIndex(const Gecode::Int::IntView &v)  const	 {
    	for (int i=0;i <domains.size(); i++) {
    		Gecode::Int::IntView vd(domains[i]);
    		if (vd.variable() == v.variable())
    			return i;
    	}
    	return -1;
	}
	virtual
    int getIndex(const Gecode::VarBase *vb) const
    {
    	for (int i=0;i <domains.size(); i++) {
    		if (domains[i].variable() == vb)
    			return i;
    	}
    	return -1;
	}
	};


	

	
	//! CSP model that extends the CPSP approach so that resulting structures
	//! differ in a given number of absolute positions
	class GC_ThreadingSpacePosDist : public GC_ThreadingSpace {
	public:
	
		typedef std::vector<Gecode::IntArgs> SolutionVec;

			//! position information of the found solutions
		SolutionVec* solPositions;

	protected:
			//! the maximal number of solution positions that may be equal to
			//! corresponding positions in already found solutions
		unsigned int maxEqual;

	public:
	
		GC_ThreadingSpacePosDist(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType = BR_DFS|BR_SYM,
			SolutionVec* solPositions_ = NULL,
			unsigned int maxEqual_ = 0)
		  : GC_ThreadingSpace( sequence, latFrame, neighVecs, seqFeatureMap,
		  		hullLvl, hCore, shiftVec, branchingType),
		  	solPositions(solPositions_), maxEqual(maxEqual_) 
		{
			// apply addition constraints of already found solutions on 
			// different cores
			constrain(this);
		}
			
		GC_ThreadingSpacePosDist(bool share, GC_ThreadingSpacePosDist& toCopy)
		  :	GC_ThreadingSpace(share, toCopy), 
			solPositions(toCopy.solPositions),
			maxEqual(toCopy.maxEqual)
		{}
		
		virtual ~GC_ThreadingSpacePosDist() 
		{}

		virtual SuperSpace* copy(bool share)
		{ return new GC_ThreadingSpacePosDist(share,*this); }

			//! adds additional constraints if BAB search is used to ensure
			//! a minimal position difference between the next and all other
		virtual void constrain(GC_ThreadingSpacePosDist* lastSolution)
		{
				// add one distance constrain for each solution found so far
			for(SolutionVec::iterator sol = solPositions->begin(); 
					sol != solPositions->end(); sol++) 
			{
				Gecode::atmost(this, domains, (*sol), maxEqual, Gecode::ICL_DOM);
			}
		}
			//! handles the information of a solution to ensure the minimal
			//! distance of the next solution to all found
		virtual void handleSolution(GC_ThreadingSpacePosDist* solution)
		{
//@TODO : add all symmetric solutions !?!?!?
			assert(solPositions != NULL);
			solPositions->push_back(Gecode::IntArgs(domains.size()));	// new sol data
			Gecode::IntArgs& y = *(solPositions->rbegin()); // access to last element
			for (int i = y.size(); i--; )			// daten kopieren 
				y[i] = solution->domains[i].val();
		}

	};
	
	
	

	
	/*! CSP model that extends the CPSP approach so that resulting structures
	 * differ in at least one H-monomer position
	 */
	class GC_ThreadingSpaceHcoreDist : public GC_ThreadingSpace {
		
	protected:
		
		const HullLevel *hullLvl;
		const size_t hCoreSize;

	public:
	
		GC_ThreadingSpaceHcoreDist(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl_,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType = BR_DFS|BR_SYM
			);
		
		GC_ThreadingSpaceHcoreDist(bool share, GC_ThreadingSpaceHcoreDist& toCopy)
		  :	GC_ThreadingSpace(share, toCopy),
	  		hullLvl(toCopy.hullLvl),
	  		hCoreSize(toCopy.hCoreSize)
		{}
		
		virtual ~GC_ThreadingSpaceHcoreDist() 
		{}

		virtual SuperSpace* copy(bool share)
		{ return new GC_ThreadingSpaceHcoreDist(share,*this); }

			//! adds additional constraints if BAB search is used to ensure
			//! a minimal Hcore difference between the next and all other
		virtual void constrain(GC_ThreadingSpaceHcoreDist* lastSolution)
		{
			
			  // temporary variables to post the new constraint on
			Gecode::IntVarArgs hVars(hCoreSize);
			Gecode::IntArgs hPos(hCoreSize);
			size_t curPos = 0;
			
			  // collect Hcore variables and Hcore assignments of current solution
			for (size_t i=0; i<hullLvl->size(); i++) {
				if (hullLvl->at(i) == 0) {
					hVars[curPos] = domains[i];
					hPos[curPos]  = lastSolution->domains[i].val();
					curPos++;
				}
			}
			
			  // ensure that at least one Hcore variable is differently assigned
			Gecode::atmost(this, hVars, hPos, (hCoreSize-1), Gecode::ICL_DOM);
			
		}
			//! handles the information of a solution to ensure the minimal
			//! distance of the next solution to all found
		virtual void handleSolution(GC_ThreadingSpaceHcoreDist* solution)
		{
		}

	};
	
	
	//! CSP model that extends the CPSP approach so that resulting structures
	//! differ in a given number of moves
	class GC_ThreadingSpaceMoveDist : public GC_ThreadingSpace {
	public:
	
		typedef std::vector<Gecode::IntArgs> SolutionVec;

			//! move information of the found solutions
		SolutionVec* solMoves;

	protected:
			//! the maximal number of solution moves that may be equal to
			//! corresponding positions in already found solutions
		unsigned int maxEqual;
		
			//! absolute moves between two variables
		Gecode::IntVarArray absMoves;
		
			//! position information for each monomer of superclass
		using GC_ThreadingSpace::domains; 
		
	public:
	
		GC_ThreadingSpaceMoveDist(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType = BR_DFS|BR_SYM,
			SolutionVec* solMoves_ = NULL,
			unsigned int maxEqual_ = 0)
		  : GC_ThreadingSpace( sequence, latFrame, neighVecs, seqFeatureMap,
		  		hullLvl, hCore, shiftVec, branchingType, NO_PROP),
		  	solMoves(solMoves_), maxEqual(maxEqual_), 
		  	absMoves(this,domains.size()-1)
		  	 
		{
			std::set<int> tmp;
			assert(neighVecs != NULL);
			for(biu::LatticeFrame::index_set::const_iterator it = neighVecs->begin();
					it != neighVecs->end(); it++)
				tmp.insert(*it);
			int neighMin = *(tmp.begin()), neighMax = *(tmp.rbegin());
			
			using namespace Gecode;
			// adding redundant neighboring constraints for absMoves propagation
			for (int i=absMoves.size(); i--; ) {
					// init move var
				absMoves[i] = Gecode::IntVar(this, neighMin, neighMax);
					// shrink move var to allowed moves
				GC_StlSetRangeIterator it(&tmp);
				Gecode::Int::IntView(absMoves[i]).narrow<GC_StlSetRangeIterator>(this, it);
					// add propagation on move var and position vars
				post(this, (domains[i+1]-domains[i]) == absMoves[i], Gecode::ICL_DOM);
			}
			// apply addition constraints of already found solutions on 
			// different cores
			constrain(this);
		}
			
		GC_ThreadingSpaceMoveDist(bool share, GC_ThreadingSpaceMoveDist& toCopy)
		  :	GC_ThreadingSpace(share, toCopy), 
			solMoves(toCopy.solMoves),
			maxEqual(toCopy.maxEqual), absMoves()
		{
			absMoves.update(this, share, toCopy.absMoves);
		}
		
		virtual ~GC_ThreadingSpaceMoveDist() 
		{}

		virtual SuperSpace* copy(bool share)
		{ return new GC_ThreadingSpaceMoveDist(share,*this); }

			//! adds additional constraints if BAB search is used to ensure
			//! a minimal position difference between the next and all other
		virtual void constrain(GC_ThreadingSpaceMoveDist* lastSolution)
		{
				// add one distance constrain for each solution found so far
			for(SolutionVec::iterator sol = solMoves->begin(); 
					sol != solMoves->end(); sol++) 
			{
				Gecode::atmost(this, absMoves, (*sol), maxEqual, Gecode::ICL_DOM);
			}
		}
			//! handles the information of a solution to ensure the minimal
			//! distance of the next solution to all found
		virtual void handleSolution(GC_ThreadingSpaceMoveDist* solution)
		{
//@TODO : add all symmetric solutions !?!?!?

//@TODO : POS_DIST : was nach kernwechsel? alte löschen? irgendwie umrechnen? normieren?
			assert(solMoves != NULL );
			// original data
			solMoves->push_back(Gecode::IntArgs(absMoves.size()));	// new sol data
			Gecode::IntArgs& y = *(solMoves->rbegin()); // access to last element
			for (int i = y.size(); i--; ) {			// daten kopieren 
				y[i] = solution->absMoves[i].val();
			}
			
		}

	};
	
	
	
	//! CSP model that extends the CPSP approach so that only structure shapes
	//! are enumerated. These shapes differ in the way they cross an H-core.
	//! For eache shape one representing structure is calculated
	class GC_ThreadingSpaceShapes : public GC_ThreadingSpaceHcoreDist {
	
	protected:
	
		using GC_ThreadingSpace::domains;
		using GC_ThreadingSpace::ranks;

	public:
	
		GC_ThreadingSpaceShapes(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType = BR_DFS|BR_SYM);
			
		GC_ThreadingSpaceShapes(bool share, GC_ThreadingSpaceShapes& toCopy)
		  :	GC_ThreadingSpaceHcoreDist(share, toCopy)
		{}
		
		virtual ~GC_ThreadingSpaceShapes() 
		{}

		virtual SuperSpace* copy(bool share)
		{ return new GC_ThreadingSpaceShapes(share,*this); }

	};
	
	

  } // namespace gecode
} // namespace cpsp

#endif /*GC_THREADINGSPACE_HH_*/
