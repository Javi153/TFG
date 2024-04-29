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

#include "protein_like_option.h"

namespace cpsp{
	//side chain threading
	namespace scth{
	ProteinLikeOptions::ProteinLikeOptions(){
	   RandomizerPath="";
	   sideChainOptions.withOutput=false;
	   statistics=false;
	}
	
   }
}
