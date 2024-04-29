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

#ifndef HPTHREADINGHANDLER_HH_
#define HPTHREADINGHANDLER_HH_

#include <vector>
#include <string>

namespace cpsp
{
	typedef std::vector<std::string> StringVec;
		/**
		 * This abstract class handles the structure generation in the 
		 * HP-lattice-modell
		 */
	class HPThreadingHandler
	{
	public:


		HPThreadingHandler()
		{}
		virtual ~HPThreadingHandler()
		{}
		
			//! starts structure generation with the internal parameters
			//! @return the number of generated structures
		virtual 
		unsigned int
		generateStructures() = 0;
		
			//! starts structure enumeration with the internal parameters and
			//! stores them in structs
			//! @param structs the container to store the found structures in
			//! @return the number of generated structures
		virtual
		unsigned int
		enumerateStructures(StringVec & structs) = 0;
		
	};

} // namespace cpsp

#endif /*HPTHREADINGHANDLER_HH_*/
