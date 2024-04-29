/*
 *  Main authors:
 *     Mohamad Rabbath <rabbath@informatik.uni-freiburg.de>
 *
 *  Contributing authors:
 *     Martin Mann <mmann@informatik.uni-freiburg.de>
 *     Sebastian Will <will@informatik.uni-freiburg.de>
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef SIDE_CHAIN_OPTDEG_OPTION_HH
#define SIDE_CHAIN_OPTDEG_OPTION_HH
#include <iostream>
#include "chain_option.h"
namespace cpsp{
   //side chain threading
   namespace scth{
	class SideChainDegOptions {
	public:
		SideChainDegOptions();
		///Side chain option
		SideChainOptions sideChainOptions;
		
		///Seeds number
		int Seed;

		///The length of the sequence
		int Length;

		///The threshold of suitable seeds
		int Threshold;

		///Set the accepted solutions one more than the threshold
		void setThreshold(int);

		///The first sequence to start with
		std::string startSeq;

		/// temperature for metropolis local search
		float TEMPERATURE; 

		///The cutoff value that is necessary for the MonteCarlo Algorithm 
		double Cutoff;

		//!< break after such many steps
    		int maxSteps;
    		
		//!< break after this time limit
		int timeLimit; 
		
		///Target degneracy
		unsigned int tgtDeg; 

		/// Give statistics, default is true
		bool Stat;

		///Verbose output, default is false
		bool verboseOut;

		/// The path directory of the HPrand binray in order to generate a random sequence
		std::string RandomizerPath;
	};
   }
}
#endif
