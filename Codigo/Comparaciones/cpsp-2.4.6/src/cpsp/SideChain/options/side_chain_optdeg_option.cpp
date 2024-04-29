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

#include"side_chain_optdeg_option.h"
namespace cpsp{
   //side chain threading
   namespace scth{
	SideChainDegOptions::SideChainDegOptions(){
	   sideChainOptions.withOutput=false;
	   Length=30;
	   Threshold=10000;
	   sideChainOptions.SolutionsNum=Threshold+1;
	   startSeq="";
	   Seed=-1;
	   tgtDeg=100; //?
	   timeLimit=-1;
	   TEMPERATURE=0.4;
	   Cutoff=0.01; //?
	   maxSteps=100000;
	   Stat=false;
	   verboseOut=true;
	   RandomizerPath="";
	}
	void SideChainDegOptions::setThreshold(int th){
	   Threshold=th;
	   sideChainOptions.SolutionsNum=Threshold+1;
	}
   }
}
