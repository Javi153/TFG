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

#ifndef HPCOMPRESS_HH_
#define HPCOMPRESS_HH_

#include <biu/OptionParser.hh>


	/**
	 * HPcompress is a tool for conversion of HP-sequence from normal/expanded
	 * form to a compressed/condensed representation and vice versa.
	 */
	class HPcompress {
	protected:

			//! inits the allowed parameters for HPcompress
		void initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const; 
		
	public:
			//! construction
		HPcompress();
			//! destruction
		~HPcompress();
		
			//! converts the given sequence to the compressed/decompressed form
			//! @return execution status of parameter parsing and the conversion
		int convert(int argc, char**argv);
	};



#endif /*HPCOMPRESS_HH_*/
