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

#include <biu/OptionParser.hh>
#include <biu/util/Util_String.h>
#include <cpsp/Exception.hh>

#include "cpsp/SideChain/chain_threading_all.h"
#include <cpsp/HCoreDatabaseFILE.hh>

#include <algorithm>

#include "version.hh"



////////////////////////////////////////////////////////////////////////////////

void 
initAllowedArguments(	biu::OptionMap & allowedArgs
						, std::string &infoText );


////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
	
	int exitVal = 0;
	cpsp::HCoreDatabase * coreDB;

	try {
		biu::OptionMap allowedArgs;
		std::string infoText;
		initAllowedArguments(allowedArgs,infoText);	// init

		cpsp::scth::SideChainOptions option;
		
		option.Output = cpsp::scth::MOVES;
		option.Draw = false;

		// parse programm arguments
		biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
				(char**)argv, infoText);
		
		// arguments parseable and all mandatory arguments given
		if (opts.noErrors()) {
			// help output
			if (opts.getBoolVal("help")) {
				opts.coutUsage();
			} else if (opts.argExist("version")) {
				giveVersion();
			} else {
				
				// enforce that only H-core representatives are generated
				option.onlyHcoreRepresentatives = true;
				
				// check if only optimal structures have to be calculated
				option.bestOnly = !opts.argExist("all");

				// check if only structure counting is needed
				option.countOnly = opts.argExist("count");

				// using symmetryBreaking
				option.BreakSym = ! opts.getBoolVal("noSymmBreak");

				// print found structures 
				option.withOutput= !opts.getBoolVal("s");
				
				option.verbose = opts.getBoolVal("v");
				
				if (option.verbose && ! option.withOutput) {
					throw cpsp::Exception("Error in arguments : cannot do verbose and silent output together ! ;)", -2);
				}

				// maximal structure number to find
				int tmp = opts.getIntVal("maxSol");
				if (tmp < 1)
					throw cpsp::Exception("Error in arguments : maxSol <= 0", -2);
				option.SolutionsNum = tmp;


				// H-core database
				char * tmpStr = NULL;
				if (opts.argExist("dbPath")) {
					if (opts.getStrVal("dbPath").size() == 0) {
						throw cpsp::Exception(
								"H-core database Error: dbPath parameter present but no database given",-1);
					}
					coreDB = new cpsp::HCoreDatabaseFILE(opts.getStrVal("dbPath"));
				} else if ((tmpStr = getenv("CPSP_COREDB")) != NULL) {
					coreDB = new cpsp::HCoreDatabaseFILE(std::string(tmpStr));
				} else {
					throw cpsp::Exception(
							"H-core database Error: no database available (FILE)",-1);
				}
				// testen
				if (!coreDB->isConnected()) {
					throw cpsp::Exception(
							"H-core database Error: can't open database",-1);
				}
				option.coreDB = coreDB;

				// lattice modell
				if (opts.getStrVal("lat") == "CUB")
					option.Description = new std::string("CUB");
				else if (opts.getStrVal("lat") == "FCC")
					option.Description = new std::string("FCC");
				else
					throw cpsp::Exception(
							"Error in arguments: lattice is not supported",-2);

				// sequence
				std::string s = opts.getStrVal("seq");
				if ( s.size() == 0 ) 
					throw cpsp::Exception("Error in arguments: no sequence given",-2);
				for (std::string::iterator si = s.begin(); si!=s.end(); si++) {
					if (*si == 'h') 		{ *si = 'H';
					} else if (*si == 'p')	{ *si = 'P';
					} else if (*si != 'H' && *si != 'P') {
						throw cpsp::Exception(
								"Error in arguments: sequence is no HP-sequence : "+s,-2);
					}
				}
				option.Sequence = new std::string(s);
			}

		} else {
			throw cpsp::Exception(
					"Parameter Error: not parseable",-1);
		}
		
		if (option.withOutput) {
			std::cout	<<"\n===================================="
						<<"\n   HPrepSC - CPSP-tool-library"
						<<"\n===================================="
						<<"\n\n";
		}

		cpsp::scth::SideChainThreadingHandler handler(option);
		
		if (option.withOutput) {
			std::cout	<<"\n====================================\n"
						<<std::endl;
		}

	
	} catch (cpsp::Exception& ex) {
		std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
		exitVal = ex.retValue;
	}
	
	if (coreDB != NULL) delete coreDB;
	
	return exitVal;
}

////////////////////////////////////////////////////////////////////////////////


void 
initAllowedArguments(	biu::OptionMap & allowedArgs
						, std::string &infoText )
{
	infoText = "HPrepSC calculates a representive subset of optimal structures of a given "
				"sequence in the 3D-lattice HP-model including side chains."
				" All structures differ in at least one H-side-chain position."
				" The structures are given in "
				"absolute move representation whereby we use the following "
				"encoding : F/B +-x, L/R +- y, U/D +-z."
				" Side chain positions are given in brackets."
				"\nTo find them an"
				" H-core database is used with optimal and suboptimal cores."
				" Beginning with the biggest possible core iteratively all"
				" cores up to a given suboptimality level are taken into account"
				" to design the sequences."			
				"\n\n - the option <all> overwrites the <allbest> option";
	
	allowedArgs.push_back(biu::COption(
							"seq", false, biu::COption::STRING,
							"HP-sequence to fold"));
	allowedArgs.push_back(biu::COption(
							"dbPath", true, biu::COption::STRING,
							"file based database root path (or use $CPSP_COREDB)"));
	allowedArgs.push_back(biu::COption(
							"lat", true, biu::COption::STRING,
							"lattice to use : CUB (cubic), FCC (face centered cubic)", "CUB"));
	allowedArgs.push_back(biu::COption(
							"noSymmBreak", true, biu::COption::BOOL,
							"do not break symmetries"));
	allowedArgs.push_back(biu::COption(
							"maxSol", true, biu::COption::INT,
							"maximum of solutions to search for (value > 0)",
							"10"));
	allowedArgs.push_back(biu::COption(
							"all", true, biu::COption::BOOL,
							"generates optimal and suboptimal structures up to the given maximal number (maxSol)"));
	allowedArgs.push_back(biu::COption(
							"count", true, biu::COption::BOOL,
							"only counting structures"));
	allowedArgs.push_back(biu::COption(
							"s", true, biu::COption::BOOL,
							"minimal output"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"verbose output"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(	
							"version", true, biu::COption::BOOL, 
							"version information of this program"));

}


////////////////////////////////////////////////////////////////////////////////

