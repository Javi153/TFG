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

#include "cpsp/HCore.hh"
#include <biu/assertbiu.hh>

namespace cpsp
{
	
	HCore::HCore() :
		contacts(0), size(0), numOfPSinglets(-1), numOfHSinglets(-1),
		evenPositions(0)
	{
		coreFrame.reset();	// init coreFrame 
	}

	HCore::HCore(const unsigned int contacts_, const unsigned int size_) :
		contacts(contacts_), size(size_), numOfPSinglets(-1), 
		numOfHSinglets(-1), evenPositions(0)
	{
		coreFrame.reset();	// init coreFrame 
	}
	
	HCore::~HCore()
	{
	}
	
	void 
	HCore::reset(const unsigned int contacts_, const unsigned int size_) {
		contacts = contacts_;
		size = size_;
		points.clear();
		coreFrame.reset();
		actFrameSize = 0;
		numOfPSinglets = -1;
		numOfHSinglets = -1;
		evenPositions = 0;
	}
	
	void 
	HCore::addPoint(const CorePoint& toAdd) {
		points.insert(toAdd);			// insertion
		coreFrame.resize(toAdd);		// resize core frame if neccessary
		if (toAdd.isEven())				// check if new position is even
			evenPositions++;
	}	
	
	void 
	HCore::indexHulls(	const biu::LatticeFrame* latFrame, 
						const unsigned int maxHullLevel) {
assertbiu(latFrame != NULL, "tried to index the hulls without LatticeFrame.");
		
		if ( actFrameSize != latFrame->getFrameSize() ) 
		{	
				// new indexing neccessary
			actFrameSize = latFrame->getFrameSize();
			hull.clear();			// triggers new indexing 
			numOfPSinglets = -1;	// triggers new indexing for singlets 
			numOfHSinglets = -1;
		}
			// abort if nothing to do
		if ( hull.size() > maxHullLevel)
			return;
			// neccessary vector to shift the core center into the frame center
		biu::IntPoint shift = latFrame->getCenter() - coreFrame.getCenter();
			// is new core indexing neccessay?
		if (hull.size() == 0) {
			hull.push_back(biu::LatticeFrame::index_set());
			for (	CPointSet::const_iterator it = points.begin(); 
					it != points.end(); it++) {
				// shift point and add index
				hull[0].insert(latFrame->getIndex((*it) + shift));	
			}
		}
			// get indexed neighbors for faster hull calculation
		biu::LatticeFrame::index_set neighInd = latFrame->getIndexedNeighborhood();
			// index hulls on demand if not calculated yet - 
			// hull[i+1] contains all points of hull i+1 
			// INCLUSIVE all points of hull i that are not part of the core
		for (unsigned int lev = hull.size(); lev <= maxHullLevel; lev++) {
			if (lev>1) 	// insert previous hull if not core
				hull.push_back(biu::LatticeFrame::index_set(hull[lev-1]));
			else		// new empty hull
				hull.push_back(biu::LatticeFrame::index_set());
				// for all indices from the previous hull (lev-1)
			for (	biu::LatticeFrame::index_set::const_iterator it = hull[lev-1].begin();
					it != hull[lev-1].end(); it++) {
					// for all neighborvectors
				for (biu::LatticeFrame::index_set::const_iterator neigh = neighInd.begin(); neigh!=neighInd.end(); neigh++) {
						// if not core element
					if (hull[0].find((*it) + *neigh) == hull[0].end()) {
							// add corepos+neighborVec to hull
						hull[lev].insert((*it) + *neigh);
					}
				}
			}
		}
	} // indexHulls(..)	
	
	const biu::LatticeFrame::index_set&
	HCore::getIndexedHull(	const biu::LatticeFrame* latFrame, 
							const unsigned int hullLevel) {
		indexHulls(latFrame, hullLevel);
		return hull[hullLevel];
	} // getIndexedHull(..)
	
	int
	HCore::getNumOfHSinglets(const biu::LatticeFrame* latFrame) {
		if (!isValid())
			numOfHSinglets = -1;
		else {
			if (	actFrameSize != latFrame->getFrameSize() 
					|| numOfHSinglets == -1) 	// recalculation neccessary
			{
				biu::LatticeFrame::index_set	core = getIndexedHull(latFrame, 0),
								hull = getIndexedHull(latFrame, 1);
				biu::LatticeFrame::index_set neighInd = latFrame->getIndexedNeighborhood();
				biu::LatticeFrame::index_set::const_iterator neigh;
				unsigned int n=0;
					// init
				numOfHSinglets = 0;
				hSinglets.clear();
					// go through all possible H-singlet positions of the core
				for (	biu::LatticeFrame::index_set::const_iterator it = core.begin();
						it!= core.end();
						it++) {
						// init
					n = 0;
						// check if there is a 2nd neighbor in first hull
					for (neigh=neighInd.begin(); neigh != neighInd.end(); neigh++) {
						if ( hull.find((*it)+(*neigh)) != hull.end())
							n++;
						if (n>1) {	// at least 2 neighbors in hull 1
							numOfHSinglets++;	// found an H-singlet
							hSinglets.insert(*it);	// store
							break;				// go to next 
						}
					}
				}
			}
		}
		return numOfHSinglets;
	}
	
	int
	HCore::getNumOfPSinglets(const biu::LatticeFrame* latFrame) {
		if (!isValid())
			numOfPSinglets = -1;
		else {
			if (	actFrameSize != latFrame->getFrameSize() 
					|| numOfPSinglets == -1) 	// recalculation neccessary
			{
				biu::LatticeFrame::index_set	core = getIndexedHull(latFrame, 0),
								hull = getIndexedHull(latFrame, 1);
				biu::LatticeFrame::index_set neighInd = latFrame->getIndexedNeighborhood();
				biu::LatticeFrame::index_set::const_iterator neigh;
				unsigned int n=0;
					// init
				numOfPSinglets = 0;
				pSinglets.clear();
					// go through all possible P-singlet positions of the hull 1
				for (	biu::LatticeFrame::index_set::const_iterator it = hull.begin();
						it!= hull.end();
						it++) {
						// init
					n = 0;
						// check if there is a 2nd neighbor in core
					for (neigh=neighInd.begin(); neigh != neighInd.end(); neigh++) {
						if ( core.find((*it)+(*neigh)) != core.end())
							n++;
						if (n>1) {	// min 2 neighbors in core
							numOfPSinglets++;	// found p-singlet position
							pSinglets.insert(*it);	// store
							break;				// go to next
						}
					}
				}
			}
		}
		return numOfPSinglets;
	}
	
		//! Returns the indexed P-singlet positions of the first hull,
		//! shifted to the center of the lattice frame.
	const biu::LatticeFrame::index_set&
	HCore::getIndexedPSinglets(const biu::LatticeFrame* latFrame) {
		getNumOfPSinglets(latFrame);
		return pSinglets;
	}

	const biu::LatticeFrame::index_set&
	HCore::getIndexedHSinglets(const biu::LatticeFrame* latFrame) {
		getNumOfHSinglets(latFrame);
		return hSinglets;
	}
	
	unsigned int 
	HCore::getMinEvenOdd(void) {
		return ( (points.size()-evenPositions) > evenPositions) ?
					evenPositions :
					points.size()-evenPositions;
	}

} // namespace cpsp
