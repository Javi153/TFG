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

#ifndef HPRAND_HH_
#define HPRAND_HH_


#include <string>
#include <biu/OptionParser.hh>
#include <cpsp/Exception.hh>

	class HPrand;

	class HPrandParameters {
	protected:
	
			//! initialisation of the allowed program arguments and the
			//! corresponding additinal informations
		virtual
		void 
		initAllowedArguments(biu::OptionMap & allowedArgs,
								std::string &infoText ) const;
				
			//! whether or not no error occured during parameter parsing				
		bool ok;
		
	public:
			//! Construction and initialisation of the HPrand object.
			//!
			//! @param argc the argument number given by the "main" function
			//! @param argv the string arguments given by the "main" function
			//! @param threading the HPrand object to initialise
		HPrandParameters(int argc, char** argv, HPrand & randHP) 
			throw(cpsp::Exception);
		
			//! destruction
		virtual ~HPrandParameters();
		
			//! whether or not no error occured during parameter parsing  
		bool
		wasOK(void) const { return ok; }
	};



	class HPrand {
	protected:
			//! whether or not to do minimal output
		bool minOutput;
			//! number of sequences to generate
		unsigned int numToGen;
			//! sequence length to generate
		unsigned int seqLength;
			//! number of Hs in seq
		int numOfHs;
		
	public:
	
		HPrand();
		
			//! generates random sequences depending on the given HPrand settings
		void
		generateSeqs(void) const;
		
		
		// getter and setter functions
	public:
			//! setter for minOutput
		void 
		setMinOutput( const bool val ) { minOutput = val; }
			//! getter for minOutput
		const bool
		isMinOutput( void ) const { return minOutput; }
		
		//! setter for numToGen
		void 
		setSeqNumber( const unsigned int val ) { numToGen = val; }
			//! getter for numToGen
		const unsigned int
		getSeqNumber( void ) const { return numToGen; }
		
		//! setter for seqLength
		void 
		setSeqLength( const unsigned int val ) { seqLength = val; }
			//! getter for seqLength
		const unsigned int
		getSeqLength( void ) const { return seqLength; }
		
		//! setter for numOfHs
		void 
		setNumOfHs( const int val ) { numOfHs = val; }
			//! getter for numOfHs
		const int
		getNumOfHs( void ) const { return numOfHs; }
		
	};
	
	

#endif /*HPRAND_HH_*/
