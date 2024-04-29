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

#ifndef CHAIN_THREADING_OPTIONS
#define CHAIN_THREADING_OPTIONS
#include<iostream>

#include <cpsp/HCoreDatabase.hh>

namespace cpsp{
   //side chain threading
   namespace scth{
	enum OPTION_OUTPUT { NORMAL,MOVES,TEST};
	class SideChainOptions{
	   public:
		/** The output style 
		* 0 if the output wanted to be ordinary
		* 1 if it should be moves
		* Other if it should be moves for drawing
		*/
		int Output;

		/**
		* A boolean variable to decide whether to draw or not
		*/
		bool Draw;
		
		/**
		* A boolean variable to introduce breaking the symmetry or not, it is true by default
		*/
		bool BreakSym;
		
		///This option is useful for generating the protein like sequences binary, it means the program will be in the silent mode
		bool withOutput;

		/// do verbose output
		bool verbose;

		///The number of desired solutions
		int SolutionsNum;
		
		//! the H-core database the threading will work with 						
		cpsp::HCoreDatabase * coreDB;
//		///The root of the Cores database
//		std::string* DBRoot;
	
		///The input sequence of amino acids
		std::string* Sequence;

		///The description of the lattice
		std::string* Description;

		///The directory of the jmlHome
		std::string* JmlHome;

		///The directory of the side chain viewer binary, in order to view directly					
		std::string* ViewerHome;
		
		//! if true: only one representive structure for each H-core assignment is produced, otherwise all structures are enumerated
		bool onlyHcoreRepresentatives;
		
		//! if true: structures are only counted; otherwise they are printed
		bool countOnly;
		
		//! if true: only best optimal structures are enumerated; otherwise all up to the given maximal number
		bool bestOnly;

		///The constructor
		SideChainOptions();
		SideChainOptions(const SideChainOptions &rhs);
//		//This constructor is so ugly so do not use it :)
//		SideChainOptions(int,bool,bool,const std::string& ,const std::string& ,const std::string& ,int,const std::string&  ,const std::string& );

		///The destructor
		~SideChainOptions();
	};
   }
}
#endif
