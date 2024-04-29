/*
 *  Main authors:
 *     Martin Mann 
 *
 *  Contributing authors:
 *     Sebastian Will 
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

#ifndef HCORE_HH_
#define HCORE_HH_


#include <algorithm>			// fuer min / max
#include <set>
#include <vector>
#include <limits.h>

#include <biu/Point.hh>
#include <biu/LatticeFrame.hh>

namespace cpsp
{
	
		/**
		 * An H-core represents the possible coordinates
		 * for the H-monomers of a lattice protein in the 
		 * HP- or HPNX-model.
		 * 
		 */
	class HCore
	{
	public:
			//! an H-core position
		typedef biu::IntPoint CorePoint;
			//! a std::set of CorePoint objects
		typedef std::set<CorePoint> CPointSet;
	
	protected:
	
			//! The CoreFrame handles the maximal dimensions of the H-core.
		class CoreFrame {
		public:
				//! minimal coordinates
			int minX, minY, minZ;
				//! maximal coordinates
			int maxX, maxY, maxZ;
			
				//! Resets the core frames data.
			void reset() { 
				minX = minY = minZ = INT_MAX; 
				maxX = maxY = maxZ = 0; 
			}
				//! Resizes the core frame if necessary.
			void resize(const biu::IntPoint& p) {
				minX = std::min( minX, p.getX() );
				minY = std::min( minY, p.getY() );
				minZ = std::min( minZ, p.getZ() );
				maxX = std::max( maxX, p.getX() );
				maxY = std::max( maxY, p.getY() );
				maxZ = std::max( maxZ, p.getZ() );
			}
				//! Returns the center of the core frame
			biu::IntPoint getCenter() const { 
				return biu::IntPoint(	minX + ((maxX-minX)/2), 
										minY + ((maxY-minY)/2), 
										minZ + ((maxZ-minZ)/2)); 
			}
			
				//! Returns the minimal corner of the core frame.
			biu::IntPoint getMinCorner() const {
				return biu::IntPoint( minX, minY, minZ );
			}
			
				//! Returns the maximal corner of the core frame.
			biu::IntPoint getMaxCorner() const {
				return biu::IntPoint( maxX, maxY, maxZ );
			}
			
				//! Returns the maximal dimension of the core frame.
			int getMaxDimension() const {
				return	std::max(	maxX-minX, 
									std::max(	maxY-minY,
												maxZ-minZ))
						+ 1 ;
			}
		};
		
		typedef std::vector<biu::LatticeFrame::index_set> IndSetVec;
		
			//! the number of H-H-contacts in the h-core
		unsigned int	contacts;
			//! the number of CorePoints of the h-core
		unsigned int	size;	 
		
		CPointSet		points;			//!< the set of core elements.
		
		int				numOfPSinglets;	//!< the number of P-singlets
		int				numOfHSinglets;	//!< the number of H-singlets

		CoreFrame		coreFrame;		//!< the core frame of this H-core

			//! the last frame size, that was used for indexing
		unsigned int actFrameSize;
		
			//! the indexed hull data, shifted to the center of latFrame
		IndSetVec		hull;
		
			//! the indexed P-singlet positions of the first hull
		biu::LatticeFrame::index_set			pSinglets;
			//! the indexed H-singlet positions of the core
		biu::LatticeFrame::index_set			hSinglets;
	
			//! Calculacte the indexed hulls of the H-core, shifted
			//! to the center of the lattice frame.
		void indexHulls(const biu::LatticeFrame* latFrame, 
						const unsigned int maxHullLevel);
						
			//! The number of even positions of the core.
		unsigned int	evenPositions;

	public:
		
		HCore();
		HCore(const unsigned int contacts, const unsigned int size);
		virtual ~HCore();
		
			//! Returns the surrounding core frame.
		const CoreFrame * const
		getCoreFrame() const { return &coreFrame; }
		
			//! Returns the size of the H-core.
			//! The return value might differ from the available
			//! CorePoints of the H-core, if this core is NOT VALID.
		unsigned int 
		getSize() const { return size; }
		
			//! Returns the number of H-H-contacts in the H-core.
		unsigned int 
		getContacts() const { return contacts; }
		
			//! Returns if this core is ready to be used. This is
			//! not the case, if not enough points are added.
		bool 
		isValid() const { return size == points.size(); }
		
			//! Returns the unshifted H-core coordinates
		const biu::IPointSet& 
		getPoints() const { return points; }
		
			//! Resizes the core to the given size and deletes the
			//! current point data. The H-core data has to be refilled
			//! afterwards to make it valid.
		void 
		reset(const unsigned int contacts, const unsigned int size);
		
			//! Adds a CorePoint object to the H-core.
		void 
		addPoint(const CorePoint& toAdd);
		
			//! Returns the center of the H-core.
		biu::IntPoint 
		getCenter() const { return coreFrame.getCenter(); }
		
			//! Returns the indexed hull of the H-core, shifted to the
			//! center of the lattice frame.
		const biu::LatticeFrame::index_set&
		getIndexedHull( const biu::LatticeFrame* latFrame,
										const unsigned int hullLevel);
										
			//! Returns the maximal dimension of the surrounding core frame.
		int 
		getMaxDimension() const { return coreFrame.getMaxDimension(); }
		
			//! Returns the number of P-singlets, if the core is valid, 
			//! otherwise -1 is returned.
		int 
		getNumOfPSinglets(const biu::LatticeFrame* latFrame);

			//! Returns the number of H-singlets, if the core is valid, 
			//! otherwise -1 is returned.
		int 
		getNumOfHSinglets(const biu::LatticeFrame* latFrame);
		
			//! Returns the indexed P-singlet positions of the first hull,
			//! shifted to the center of the lattice frame.
		const biu::LatticeFrame::index_set&
		getIndexedPSinglets(const biu::LatticeFrame* latFrame);

			//! Returns the indexed H-singlet positions of the core,
			//! shifted to the center of the lattice frame.
		const biu::LatticeFrame::index_set&
		getIndexedHSinglets(const biu::LatticeFrame* latFrame);
		
			//! Returns the smaller number of even vs. odd positions.
		unsigned int 
		getMinEvenOdd(void);
	};

} // namespace cpsp

#endif /*HCORE_HH_*/
