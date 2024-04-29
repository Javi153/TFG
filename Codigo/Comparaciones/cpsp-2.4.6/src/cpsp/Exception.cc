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

#include "cpsp/Exception.hh"

namespace cpsp
{

	unsigned int Exception::errorCount = 0;


	Exception::Exception(): 
		outString("\nException raised\n"), retValue(-1)  
	{ 
		errorCount++; 
	}
	
	Exception::Exception(const std::string &outString) : 
		outString(outString), retValue(-1) 
	{ 
		errorCount++; 
	}
	
	Exception::Exception(const std::string &outString, const int retValue) : 
		outString(outString), retValue(retValue) 
	{ 
		errorCount++; 
	}
	
	Exception::~Exception() {
	}

	std::ostream& operator<< (std::ostream& os, const Exception& p) {
		if (&p == NULL)
			os <<"NULL";
		else
			os <<p.toString();
		return os;
	}

} // namespace cpsp
