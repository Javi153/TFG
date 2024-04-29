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

#ifndef GC_HPTHREADING_HH_
#define GC_HPTHREADING_HH_

#include "biu/LatticeDescriptorCUB.hh"

#include "cpsp/HPThreadingHandler.hh"
#include "cpsp/gecode/GC_ThreadingSpace.hh"

namespace cpsp
{

  /**
   * The "gecode" namespace collects all classes neccessary for an
   * implementation of CPSP with the Gecode(TM) constraint programming
   * library.
   *
   * http://www.gecode.org
   */
  namespace gecode
  {

		/**
		 * This is a HPThreadingHandler implementation using the constraint
		 * propagation library Gecode(TM)
		 */
	class GC_HPThreading : public cpsp::HPThreadingHandler
	{

	public:

			//! constants that define the core selection mode
		enum CORESELECTS { FIRST, ALL_BEST, ALL };	// steuert kernauswahl
			//! possible run modes of the threading
		enum RUNMODE { COUNT, ENUMERATE };
			//! possible solution distance measures
		enum DISTANCE { DIST_POSITIONS, DIST_ABSMOVES, DIST_SHAPES, DIST_HCORE, DIST_NONE }; 
		


	protected:
			//! the lattice descriptor of the 3D-cubic lattice for conversion issues
		const biu::LatticeDescriptorCUB cubicLatDescr;

			//! whether or not to do a verbose output
		bool verboseOutput;
			//! the lattice descriptor of the lattice to do the threading in
		const biu::LatticeDescriptor* latDescr;
			//! the lattice frame to do the threading in
		biu::LatticeFrame* latFrame;
			//! the opened H-core database
		HCoreDatabase* coreDB;
			//! the HP-sequence to find optimal structure for
		std::string sequence;
			//! the hull levels for the sequence positions
		GC_ThreadingSpace::HullLevel hullLvl;
			//! the neccessary core size
		unsigned int coreSize;
			//! minimal even/odd distribution of core positions in the cubic#
			//! lattice
		unsigned int minEvenOddHs;
			//! contacts inside the sequence
		unsigned int contactsInSeq;
			//! maximal distance of a P to an H in the sequence
		unsigned int maxPCoreDistance;
			//! maps the sequence indices to their features
		GC_ThreadingSpace::SeqFeatureMap seqFeatureMap;
			//! whether or not to take symmetric structures into account
		bool useSymmetryBreaking;
			//! whether or not to use DDS
		bool useDDS;
			//! maximal structures to find
		unsigned int maxStructures;
			//! which cores to take into account
		CORESELECTS coreSelection;
			//! the run mode for generateStructures()
		RUNMODE runMode;
			//! whether or not generateStructures() does any output or not
		bool noOutput;
			//! whether ot not to print search statistics
		bool statisticsOut;
			//! wether or not to run the search engine with or without batch
			//! recomputation during searhc (with == less propagations and more
			//! space consumption, see Gecode help)
		bool noRecomputation;
			//! the distance measure between solutions
		DISTANCE distanceMode;
			//! the minimal distance between two solutions
		unsigned int distanceMin;
			//! whether to print normalized move strings or not
		bool normalizeStructures;
		
		
		template<class ENGINE, class SPACE>
		unsigned int
		enumerateStructures(SPACE *space, unsigned int maxNum, 
			Gecode::Search::Statistics *stat, std::ostream &out);
		
		template<class Int, class SPACE>
		Int
		countBAB(SPACE *space, unsigned int c_d, unsigned int a_d, 
			Gecode::Search::Statistics* stat, Gecode::Search::Stop* st, 
			Int maxNum); 

	public:
			//! construction
		GC_HPThreading(const biu::LatticeDescriptor* latDescr_ = NULL,
			HCoreDatabase* db = NULL, const std::string seq = "")
			throw(cpsp::Exception);
		virtual ~GC_HPThreading();

			//! enumerates maxNum structures of space and prints them to std::cout
		unsigned int
		enumerateStructures(SuperSpace * space, unsigned int maxNum,
		Gecode::Search::Statistics *globalStat = NULL, 
		std::ostream &out = std::cout);

			//! counts maxNum structures of space and returns their number
		unsigned int
		countStructures(SuperSpace * space,  unsigned int maxNum,
		bool useDDS = true, Gecode::Search::Statistics *statistic = NULL);

			// abstract function implementation
		virtual
		unsigned int
		generateStructures();

			//! enumerates structures and writes their string representation 
			//! into structs
		virtual
		unsigned int
		enumerateStructures(StringVec & structs);

	public:

			//! setter for verboseOutput
		void
		setVerboseOutput( const bool val ) { verboseOutput = val; }
			//! getter for verboseOutput
		const bool
		getVerboseOutput( void ) const { return verboseOutput; }

			//! setter for latDescr (creates a new LatticeFrame)
		void
		setLatticeDescriptor( const biu::LatticeDescriptor* val );
			//! getter for latDescr
		const biu::LatticeDescriptor*
		getLatticeDescriptor( void ) const { return latDescr; }

			//! setter for coreDB
		void
		setCoreDB( HCoreDatabase* db ) { coreDB = db; }
			//! getter for coreDB
		HCoreDatabase*
		getCoreDB( void ) const { return coreDB; }

			//! setter for the HP-sequence
		void
		setSequence( const std::string& seq ) throw(cpsp::Exception);
			//! getter for sequence
		const std::string&
		getSequence( void ) { return sequence; }

			//! getter for coreSize
		const unsigned int
		getCoreSize( void ) const { return coreSize; }

			//! getter for contactsInSeq
		const unsigned int
		getContactsInSeq( void ) const { return contactsInSeq; }

			//! setter for useSymmetryBreaking
		void
		setSymmetryBreaking( const bool val ) { useSymmetryBreaking = val; }
			//! getter for useSymmetryBreaking
		const bool
		getSymmetryBreaking( void ) const { return useSymmetryBreaking; }

		//! setter for useDDS
		void
		setUsingDDS( const bool val ) { useDDS = val; }
			//! getter for useDDS
		const bool
		getUsingDDS( void ) const { return useDDS; }

		//! setter for maxStructures
		void
		setMaxStructures( const unsigned int val ) { maxStructures = val; }
			//! getter for maxStructures
		const unsigned int
		getMaxStructures( void ) const { return maxStructures; }

		//! setter for coreSelection
		void
		setCoreSelection( const CORESELECTS val ) { coreSelection = val; }
			//! getter for coreSelection
		const CORESELECTS
		getCoreSelection( void ) const { return coreSelection; }

		//! setter for runMode
		void
		setRunMode( const RUNMODE val ) { runMode = val; }
			//! getter for runMode
		const RUNMODE
		getRunMode( void ) const { return runMode; }

		//! setter for noOutput
		void
		setNoOutput( const bool val ) { noOutput = val; if(noOutput) setVerboseOutput(false);}
			//! getter for noOutput
		const bool
		getNoOutput( void ) const { return noOutput; }

		//! setter for statisticsOut
		void
		setStatOutput( const bool val ) { statisticsOut = val; }
			//! getter for statisticsOut
		const bool
		getStatOutput( void ) const { return statisticsOut; }

		//! setter for noRecomputation
		void
		setNoRecomputation( const bool val ) { noRecomputation = val; }
			//! getter for statisticsOut
		const bool
		getNoRecomputation( void ) const { return noRecomputation; }

		//! setter for distanceMode
		void 
		setDistanceMode( const DISTANCE val ) { distanceMode = val; }
			//! getter for distanceMode
		const DISTANCE
		getDistanceMode( void ) const { return distanceMode; }
		
		//! setter for distanceMin
		void 
		setDistanceMin( const unsigned int val ) { distanceMin = val; }
			//! getter for distanceMin
		const unsigned int
		getDistanceMin( void ) const { return distanceMin; }
		
		//! setter for normalizeStructures
		void 
		setNormalizeStructures( const bool val ) { normalizeStructures = val; }
			//! getter for normalizeStructures
		const bool
		getNormalizeStructures( void ) const { return normalizeStructures; }
		

			//! extracts the equivalence class string encoding from a given
			//! structure
		std::string
		getEquivClassString( const std::string& sequence
					, const biu::MoveSequence & structure ) const;


			//! creates a unique string encoding for a given vector in 3D integer space
		std::string
		getXYZstring( const biu::IntPoint& vec ) const;

	};
	
	template<class ENGINE, class SPACE>
	inline
	unsigned int
	GC_HPThreading::enumerateStructures(SPACE *space,
		unsigned int maxNum, Gecode::Search::Statistics *stat, 
		std::ostream &out)
	{
		unsigned int sols = 0;
		SPACE* s = dynamic_cast<SPACE*>(space);
			// e = searchengine
		ENGINE dfs_e( s,
					(noRecomputation?1:Gecode::Search::Config::c_d),
					(noRecomputation?1:Gecode::Search::Config::a_d),
					NULL );
			// first solution if any
		SPACE* solution = dfs_e.next();
			// enumerate solutions
		while (solution != NULL && ++sols <= maxNum) {
			// ausgabe der struktur
			if (normalizeStructures) {
				out	<<latFrame->getString(
						latFrame->getDescriptor()->normalizeSequence(
							latFrame->indicesToAbsMoves(
								solution->getSolution())))
					<<std::endl;
			} else {
				out	<<latFrame->getString(
							latFrame->indicesToAbsMoves(
								solution->getSolution()))
					<<std::endl;
			}
			// handle solution if neccessary
			solution->handleSolution(solution);
			delete solution;		// clear memory
			solution = dfs_e.next();
		} ;
		
		if (solution != NULL)
			delete solution;		// clear memory

		if (stat != NULL) {
				// get statistics
			*stat = dfs_e.statistics();
		}
		return (sols > maxNum ? maxNum : sols);
	}

	template<class Int, class SPACE>
	inline
	Int
	GC_HPThreading::countBAB(SPACE *space, unsigned int c_d, unsigned int a_d, 
		Gecode::Search::Statistics* stat, Gecode::Search::Stop* st,
		Int maxNum) 
	{
	        // init branch and bound search engine
	    Gecode::BAB<SPACE> engine( space, c_d, a_d, st);
			// counter
	    Int sols = 0;
	
	        // get first solution if available
	    SPACE* nextSol = engine.next();
	        // while solution found
	    while( nextSol != NULL && sols <= maxNum) {
	
	        sols++; // count for statistics
	        	// handle solution
			nextSol->handleSolution(nextSol);
	        delete nextSol;             // clean up
	        nextSol = engine.next();	// get next solution
	    }
	    if (nextSol != NULL)
	    	delete nextSol;				// clean up
	    if (stat != NULL)
	       *stat = engine.statistics();
	    return (sols > maxNum ? maxNum : sols);
	}
	

  } // namespace gecode
} // namespace cpsp

#endif /*GC_HPTHREADING_HH_*/
