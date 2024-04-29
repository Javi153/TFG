/*
 *  Main authors:
 *     Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *
 *  Copyright:
 *     Martin Mann, 2009
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef HPREP_HH_
#define HPREP_HH_

#include <string>
#include <map>
#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptor.hh>
#include <cpsp/Exception.hh>
#include <cpsp/HCoreDatabase.hh>
#include <cpsp/gecode/GC_HPThreading.hh>



	class HPrepParameters {

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
			//! move sequence of the structure; will have length 0 if none was given
		biu::MoveSequence givenStructure;
		
	public:
			//! Construction and initialisation of the threading object.
			//!
			//! @param argc the argument number given by the "main" function
			//! @param argv the string arguments given by the "main" function
			//! @param threading the GC_HPThreading object to initialise
		HPrepParameters(int argc, char** argv, cpsp::gecode::GC_HPThreading & threading) 
			throw(cpsp::Exception);
		
			//! destruction
		virtual ~HPrepParameters();
		
			//! whether or not to do minimal output
		bool
		silent(void) const { return minimalOutput; }

			//! whether or not a structure was given; i.e. class has to be extracted
		bool
		moveStringGiven(void) const { return givenStructure.size()!=0; }

			//! access to the given structure
		const biu::MoveSequence&
		getGivenMoveSequence(void) const {return givenStructure;}

	};


#endif /*HPREP_HH_*/
