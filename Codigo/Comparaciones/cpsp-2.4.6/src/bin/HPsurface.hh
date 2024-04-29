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

#ifndef HPSURFACE_HH_
#define HPSURFACE_HH_

#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>

	/** HPsurface offers several properties of a lattice proteins surface.
	 * - number of 'open' H-contacts
	 */
	class HPsurface
	{
	protected:
	
			// the descriptor of the lattice to use
		biu::LatticeDescriptor* latDescr;

		  //! the sequence of the current protein
		std::string seq;
		
		  //! the point representation of the structure
		biu::IPointVec points;
		
		  //! the sorted set of structure points
		biu::IPointSet pointSet;
	
			//! Inits the allowed parameter setting.
		void initAllowedArguments(biu::OptionMap & allowedArgs,
									std::string &infoText ) const;
		
		  //! calculates the surface positions and contacts of the given monomers
		  //! @param surfaceIdx indices of points of monomers on the surface
		  //! @param hull points of the first hull in contacts with the given monomers
		  //! @param monomerType the type of the monomers in focus (H or P) 
		unsigned int getSurface( std::set<size_t> & surfaceIdx, biu::IPointSet & hull, const char monomerType);
		
		  //! calculates the number of internal contacts between monomers of 
		  //! type1 and type2.
		  //! @param type1 the sequence character of monomer type 1
		  //! @param type2 the sequence character of monomer type 2
		  //! @param lat the lattice for neighboring check
		  //! @return the number of internal contacts
		unsigned int getInternalContacts( char type1, char type2, const biu::LatticeModel& lat );		
	
	public:
		HPsurface();
		virtual ~HPsurface();
		
		  // starts surface investigation with the given program arguments
		int run( int argc, char** argv );
	};

#endif /*HPSURFACE_HH_*/
