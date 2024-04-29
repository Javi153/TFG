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

#ifndef HPSTRUCT_HH_
#define HPSTRUCT_HH_

#include <string>
#include <map>
#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptor.hh>
#include <cpsp/Exception.hh>
#include <cpsp/HCoreDatabase.hh>
#include <cpsp/gecode/GC_HPThreading.hh>



	class HPstructParameters {

	protected:
	
			//! initialisation of the allowed program arguments and the
			//! corresponding additinal informations
		virtual
		void 
		initAllowedArguments(biu::OptionMap & allowedArgs,
								std::string &infoText ) const;
		
			//! the lattice descriptor of the lattice to predict in
		biu::LatticeDescriptor * latDescr;
			//! the H-core database the threading will work with 						
		cpsp::HCoreDatabase * coreDB;
			//! minimal output
		bool minimalOutput;
		
	public:
			//! Construction and initialisation of the threading object.
			//!
			//! @param argc the argument number given by the "main" function
			//! @param argv the string arguments given by the "main" function
			//! @param threading the GC_HPThreading object to initialise
		HPstructParameters(int argc, char** argv, cpsp::gecode::GC_HPThreading & threading) 
			throw(cpsp::Exception);
		
			//! destruction
		virtual ~HPstructParameters();
		
			//! whether or not to do minimal output
		bool
		silent(void) const { return minimalOutput; }
	
	};


#endif /*HPSTRUCT_HH_*/
