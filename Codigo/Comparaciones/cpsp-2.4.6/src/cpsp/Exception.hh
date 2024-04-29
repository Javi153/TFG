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

#ifndef EXCEPTION_HH_
#define EXCEPTION_HH_

#include <string>
#include <iostream>

/**
 * This namespace includes all classes and functionalities to predict optimal
 * and suboptimal structures for lattice proteins in the HP-modell using the 
 * approach introduced by Rolf Backofen and Sebastian Will.
 * The method is known as "Constraint-based Protein Structure Prediction"
 * (CPSP) and published as
 * 
 * "A constraint-based approach to fast and exact structure prediction in
 * three-dimensional protein models", R. Backofen and S. Will, Journal of
 * Constraints, 2005.
 * 
 */
namespace cpsp
{

	class Exception
	{
	private:
		std::string outString;
	public:
		Exception();
		Exception(const std::string &outString);
		Exception(const std::string &outString, const int retValue);
		virtual ~Exception();

		std::string toString() const { return outString; }
	
		int retValue;
		static unsigned int errorCount;
	};
	
	std::ostream& operator<< (std::ostream& os, const Exception &p);
	

} // namespace cpsp

#endif /*EXCEPTION_HH_*/
