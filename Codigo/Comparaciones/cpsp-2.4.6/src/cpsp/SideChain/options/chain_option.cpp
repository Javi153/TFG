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

#include"chain_option.h"
namespace cpsp{
   //side chain threading
   namespace scth{
	SideChainOptions::SideChainOptions()
	 :	Output(NORMAL)
	 	, Draw(false)
	 	, BreakSym(true)
	 	, withOutput(true)
	 	, verbose(false)
	 	, SolutionsNum(1)
	 	, coreDB(NULL)
	 	, Sequence(NULL)
	 	, Description(NULL)
	 	, JmlHome(NULL)
	 	, ViewerHome(NULL)
	 	, onlyHcoreRepresentatives(false)
	 	, countOnly(false)
	 	, bestOnly(true)
	{
	   //reading the JmlHome from the environment variables
	    char * jPath;
	   jPath=getenv ("JMOL_HOME");
	   if(jPath!=NULL){
		JmlHome=new std::string(jPath);
	   }
	}


	SideChainOptions::SideChainOptions(const SideChainOptions &rhs)
	 :	Output(rhs.Output)
	 	, Draw(rhs.Draw)
	 	, BreakSym(rhs.BreakSym)
	 	, withOutput(rhs.withOutput)
	 	, verbose(rhs.verbose)
	 	, SolutionsNum(rhs.SolutionsNum)
	 	, coreDB(rhs.coreDB)
	 	, Sequence(NULL)
	 	, Description(NULL)
	 	, JmlHome(NULL)
	 	, ViewerHome(NULL)
	 	, onlyHcoreRepresentatives(rhs.onlyHcoreRepresentatives)
	 	, countOnly(rhs.countOnly)
	 	, bestOnly(rhs.bestOnly)
 	{

	   if(rhs.Sequence!=NULL)
	   	this->Sequence=new std::string(*(rhs.Sequence));
	   if(rhs.Description!=NULL)
	   	this->Description=new std::string(*(rhs.Description));
	   if(rhs.JmlHome!=NULL)
	   	this->Description=new std::string(*(rhs.JmlHome));
	   if(rhs.ViewerHome!=NULL)
	   	this->Description=new std::string(*(rhs.ViewerHome));
	}
	
	SideChainOptions::~SideChainOptions(){
	   if (Description != NULL) delete(Description);
	   Description=NULL;
	   if (Sequence != NULL) delete(Sequence);
	   Sequence=NULL;
	   if (JmlHome != NULL) delete(JmlHome);
	   JmlHome=NULL;
	   if (ViewerHome != NULL) delete(ViewerHome);
	   ViewerHome=NULL;
	}
	
   }
}
