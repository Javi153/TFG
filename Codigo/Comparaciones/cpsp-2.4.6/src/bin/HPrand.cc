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



#include "HPrand.hh"
#include "version.hh"

#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

//////////////////////////////////////////////////////////////////
// THE PARAMETER PARSER FOR HPrand
//////////////////////////////////////////////////////////////////

	HPrandParameters::HPrandParameters(int argc, char** argv, 
		HPrand & randHP) throw(cpsp::Exception) 
	{
		biu::OptionMap allowedArgs;
		std::string infoText;
		initAllowedArguments(allowedArgs,infoText);	// init
		
			// parse programm arguments	
		biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
														(char**)argv, infoText);
			// arguments parseable and all mandatory arguments given
		ok = opts.noErrors();
		if ( ok ) {
				// help output
			if (opts.getBoolVal("help")) {
				opts.coutUsage();
				throw cpsp::Exception("",0);
			}
			if (opts.getBoolVal("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
			
				// minimal output
			randHP.setMinOutput(opts.getBoolVal("s"));
			
			int tmp = 0;
				// sequence length
			if ((tmp = opts.getIntVal("len")) > 0) {
				randHP.setSeqLength((unsigned int)tmp);
			} else {
				throw cpsp::Exception("Error in arguments : len <= 0", -2);
			}
				// number of sequences
			if ((tmp = opts.getIntVal("n")) > 0) {
				randHP.setSeqNumber((unsigned int)tmp);
			} else {
				throw cpsp::Exception("Error in arguments : n <= 0", -2);
			}
				// number of Hs in sequence if given
			if (opts.argExist("H") ) {
				if ((tmp = opts.getIntVal("H")) >= 0) {
					if ((unsigned int)tmp <= randHP.getSeqLength())
						randHP.setNumOfHs(tmp);
					else 
						throw cpsp::Exception("Error in arguments : H > len", -2);
				} else {
					throw cpsp::Exception("Error in arguments : H < 0", -2);
				}
			}

			// init randomizer
			if (opts.argExist("seed")) {
				int tmp = opts.getIntVal("seed");
				if (tmp < 0) 
					throw cpsp::Exception("Error in arguments : seed < zero", -2);
				srand ( tmp );
			} else {
				srand( time(0) );
			}

		}
	}
	
	HPrandParameters::~HPrandParameters()
	{
	}

	void
	HPrandParameters::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(	
								"len", false, biu::COption::INT, 
								"the length of the sequences (value >= 1)"));
		allowedArgs.push_back(biu::COption(	
								"n", true, biu::COption::INT, 
								"the number of sequences to generate (value >= 1)",
								"1"));
		allowedArgs.push_back(biu::COption(	
								"H", true, biu::COption::INT, 
								"number of Hs in the sequence (value >= 0)"));
		allowedArgs.push_back(biu::COption(	
								"s", true, biu::COption::BOOL, 
								"silent - minimal output"));
		allowedArgs.push_back(biu::COption(	
								"seed", true, biu::COption::INT, 
								"seed value of the randomizer"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPrand generates random HP-sequences of a given length ")
			+ std::string("and prints them to STDOUT.");
	
	} // initArguments



//////////////////////////////////////////////////////////////////
// HPrand IMPLEMENTATIONS
//////////////////////////////////////////////////////////////////


	HPrand::HPrand() :
		minOutput(false),
		numToGen(1),
		seqLength(0),
		numOfHs(-1)	 
	{
	}
	
	
	void
	HPrand::generateSeqs() const {
		if (seqLength == 0 || numToGen == 0)
			return;
			
			// the sequence temp string
		std::string seq(seqLength, 'P');
		
		if (numOfHs >= 0) {
			for (int i=numOfHs; i>0; i--)
				seq[i] = 'H';
			// randomize sequence and output
			for (unsigned int n = 0; n<numToGen; n++) {
				std::random_shuffle(seq.begin(), seq.end());
				std::cout << seq <<std::endl;
			} 
		} else {
			// randomize sequence and output
			std::string::size_type i=0;
			for (unsigned int n = 0; n<numToGen; n++) {
				for ( i=0; i<seqLength; i++) {
					seq[i] = ((rand()%2 == 0)?'H':'P');
				}				
				std::cout << seq <<std::endl;
			}			
		}
	}
	
	
//////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////
	

	int 
	main( int argc, char **  argv ) {
		try {
			
			// generate threading object with user parameter settings
			HPrand randHP;

			// check user parameters and init threading object
			HPrandParameters params(argc,argv, randHP);
			if (!params.wasOK())
				return -1;
			
			if (!randHP.isMinOutput()) {
				std::cout	<<"\n=================================="
							<<"\n   HPrand - CPSP-tool-library"
							<<"\n=================================="
							<<"\n\n";
			}
						
			// generate sequences 
			randHP.generateSeqs();
			
			if (!randHP.isMinOutput()) {
				std::cout	<<"\n==================================\n"
							<<std::endl;
			}
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}

		return 0;
	} 

