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

#ifndef HCOREDATABASEXML_HH_
#define HCOREDATABASEXML_HH_


#include "cpsp/HCoreDatabase.hh"
#include <string>
#include <fstream>


namespace cpsp
{

	typedef std::vector<int> intVec;

		/**
		 * The HCoreDatabaseFILE provides an interface to a textfile
		 * based H-core database.
		 * 
		 * It handles the access to HCore objects from the database.
		 */
	class HCoreDatabaseFILE : public HCoreDatabase
	{
	protected:
	
		unsigned int coreSize;			//!< the size of the currently selected HCores
		unsigned int minHH, maxHH;		//!< HH-contact limits for HCore selection
		const std::string DBROOTPATH;	//!< the root path of the database files
	
		intVec actHHcontacts;			//!< handles available contacts for the current coreSize
		std::ifstream *actFile;			//!< current file stream the data is read from
		std::string actFileNamePref;	//!< the filename of the current file streams without energy und postfix
		int actFileHHcontInd;			//!< contact index of the current file stream in actEnergies
		
		// Exception Klassen
		class ExcEOF {};				//!< endOfFile exception class / no further data available
	
	public:
	
		HCoreDatabaseFILE(const std::string& dbRootPath);
		HCoreDatabaseFILE(const HCoreDatabaseFILE& toCopy);
		virtual ~HCoreDatabaseFILE();

			//! Initializes the database access to control the HCore objects
			//! returned by getNextCore()
			//! @return == true if the initialization was successfull
		virtual bool initCoreAccess(	const biu::LatticeDescriptor& latDescr, 
										const unsigned int size, 
										const unsigned int minContacts=0, 
										const unsigned int maxContacts=UINT_MAX-2);
										
			//! Fills succesively all HCore objects in the ranges given by 
			//! initAccess(...) in descending order according to the number of
			//! contacts.
			//!
			//! @return returns whether or not filling was successfull and if 
			//!			the filled HCore toFill is valid or not.
		virtual bool getNextCore(HCore& toFill);

			//! Returns the core size of the HCore objects returned by 
			//! getNextCore().
		virtual unsigned int getActCoreSize() { return coreSize; }
		
			//! Returns the minimal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual unsigned int getActMinHHcontacts() { return minHH; }
		
			//! Sets the minimal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual void setActMinHHcontacts(const unsigned int minContacts) { 
			minHH = minContacts; 
		}
		
			//! Returns the maximal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual unsigned int getActMaxHHcontacts() { return maxHH; }
		
			//! Returns whether or not the database is connected and open 
			//! for access.
		virtual bool isConnected() const { return true; }
		
	};

} // namespace cpsp

#endif /*HCoreDatabaseFILE_HH_*/
