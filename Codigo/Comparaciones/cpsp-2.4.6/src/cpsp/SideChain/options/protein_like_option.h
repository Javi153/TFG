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

#ifndef PROTEIN_LIKE_OPTION_HH
#define PROTEIN_LIIKE_OPTION_HH
#include "chain_option.h"
using namespace std;
namespace cpsp{
   //side chain threading
   namespace scth{
	class ProteinLikeOptions {
	   public:
		/// The path directory of the HPrand binray in order to generate a random sequence
		string RandomizerPath;

		/// The threshold ( maximum degeneracy to be accepted)
		int Threshold;
	
		///The length of the string to be generated
		string SeqLen;

		///To give statistics or not, default is false
		bool statistics;

		///The constructor
		ProteinLikeOptions();

		///The necessary options (The ProteinLikeOptions needs an instance of the SideChainOptions to simplify the implementation)
		SideChainOptions sideChainOptions;	
	};
  }
}
#endif
