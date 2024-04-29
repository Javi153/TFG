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


#include "HPcompress.hh"
#include "version.hh"
#include <biu/util/Util_String.h>
#include <cpsp/Exception.hh>


	HPcompress::HPcompress() 
	{
	}

	HPcompress::~HPcompress()
	{
	}
	
	int
	HPcompress::convert(int argc, char ** argv)
	{
		using namespace biu;
		using namespace util;
		try {
			

			// check user parameters 
			biu::OptionMap allowedArgs;
			std::string infoText;
			initAllowedArguments(allowedArgs,infoText);	// init
			
				// parse programm arguments	
			biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
															(char**)argv, infoText);
				// arguments not parseable or not all mandatory arguments given
			if (!opts.noErrors()) {
				throw cpsp::Exception("",-1);
			}
				// help output
			if (opts.argExist("help")) {
				opts.coutUsage();
				throw cpsp::Exception("",0);
			}
				// help output
			if (opts.argExist("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
				// run mode
			bool compressing = !opts.argExist("d");
				// get sequence
			std::string seq = Util_String::str2upperCase(opts.getStrVal("seq"));
			if ( seq.size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no sequence given",-2);
				// init allowed alphabet
			std::string alphN = "1234567890", alphHP = "HP";
				// check if sequence is valid alphabet sequence
			if ( seq.find_first_not_of(compressing?alphHP:alphHP+alphN) < seq.size() ) 
				throw cpsp::Exception("given sequence contains non-supported characters",-2);
				// check if last is a character
			if ( seq.find_first_not_of(alphHP,seq.size()-1) < seq.size() )
				throw cpsp::Exception("last char in sequence is not in 'HPhp'",-2);
				
			bool verboseOut = opts.argExist("v");
			
			if (verboseOut) {
				std::cout	<<"\n===================================="
							<<"\n   HPcompress - CPSP-tool-library"
							<<"\n===================================="
							<<"\n\n";
			}
				// the output string
			std::string output = "";
			std::string::size_type i = 0, length = 1;
				// generate sequence 
			if (compressing) {
				i = 0;
				while( i < seq.size() ) {
					length = std::min(seq.find_first_not_of(seq[i],i+1),seq.size()) - i;
					if (length > 2) {	// compress
						output += Util_String::int2str(length) + seq[i]; 
					} else { // just add one or two chars
						output.resize(output.size()+length, seq[i]);
					}
					i += length;
				}
			} else {
				i = seq.find_first_of(alphN,0);	// init i
				std::string::size_type lastUnCompr = 0;	// == 
				int number = 0;
				while( i < seq.size() ) {
						// append uncompressed sub-sequence
					output += seq.substr(lastUnCompr,i-lastUnCompr);
						// read number of letters to expand
					length = std::min(seq.find_first_not_of(alphN,i+1),seq.size()) - i;
					number = Util_String::str2int(seq.substr(i,length));
					i += length;
					lastUnCompr = i+1;
						// append char sequence of length 'number'
					output.resize(output.size()+number,seq[i]);
						// search next number
					i = seq.find_first_of(alphN,i);	
				}
					// append uncompressed sub-sequence
				output += seq.substr(lastUnCompr,seq.size()-lastUnCompr);
			}
			std::cout <<output <<std::endl;
			
			if (verboseOut) {
				std::cout	<<"\n====================================\n"
							<<std::endl;
			}
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}

		return 0;
	}

	void	
	HPcompress::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const 
	{
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"the HP-sequences to (de)compress"));
		allowedArgs.push_back(biu::COption(	
								"d", true, biu::COption::BOOL, 
								"decompress given sequence"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPcompress (de)compresses HP-sequences")
			+ std::string(" and prints them to STDOUT.")
			+ std::string("\n\n  e.g. HHHHPPPPPH <--> 4H5PH");
	}


	int main(int argc, char** argv) {
			// create (de)compression object
		HPcompress hpc;
			// run (de)compression
		return hpc.convert(argc,argv);
	}

