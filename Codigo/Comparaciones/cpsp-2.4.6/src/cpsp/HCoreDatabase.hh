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

#ifndef HCOREDATABASE_HH_
#define HCOREDATABASE_HH_


#include "cpsp/HCore.hh"
#include <biu/LatticeDescriptor.hh>

namespace cpsp
{
		/**
		 * An HCoreDatabase handles the access to stored HCore objects.
		 */
	class HCoreDatabase
	{
	public:
		HCoreDatabase() {};
		virtual ~HCoreDatabase() {};

		// memberfunktionen
		
			//! Initializes the database access to control the HCore objects
			//! returned by getNextCore()
			//! @return == true if the initialization was successfull
		virtual bool initCoreAccess(	const biu::LatticeDescriptor& latDescr, 
										const unsigned int size, 
										const unsigned int minContacts=0, 
										const unsigned int maxContacts=UINT_MAX-2
										) = 0;
										
			//! Fills succesively all HCore objects in the ranges given by 
			//! initAccess(...) in descending order according to the number of
			//! contacts.
			//!
			//! @return returns whether or not filling was successfull and if 
			//!			the filled HCore toFill is valid or not.
		virtual bool getNextCore(HCore& toFill) = 0;
		
			//! Returns the core size of the HCore objects returned by getNextCore().
		virtual unsigned int getActCoreSize() = 0;
		
			//! Returns the minimal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual unsigned int getActMinHHcontacts() = 0;
		
			//! Sets the minimal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual void setActMinHHcontacts(const unsigned int minContacts) = 0;
		
			//! Returns the maximal number of contacts of the HCore objects
			//! returned by getNextCore().
		virtual unsigned int getActMaxHHcontacts() = 0;
		
			//! Returns whether or not the database is connected and open 
			//! for access.
		virtual bool isConnected() const  = 0;
		
	};

} // namespace cpsp

#endif /*HCOREDATABASE_HH_*/
