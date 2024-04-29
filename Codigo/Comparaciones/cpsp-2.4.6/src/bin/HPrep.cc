/*
 *  Main authors:
 *     Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *
 *  Copyright:
 *     Martin Mann, 2009
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */


#include "HPrep.hh"
#include "version.hh"

#include <iostream>
#include <biu/util/Util_String.h>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <cpsp/HCoreDatabaseFILE.hh>
#include <cpsp/gecode/GC_HPThreading.hh>

//////////////////////////////////////////////////////////////////
// THE PARAMETER PARSER FOR GC_HPTHREADING
//////////////////////////////////////////////////////////////////

	HPrepParameters::HPrepParameters(int argc, char** argv,
		cpsp::gecode::GC_HPThreading & threading) throw(cpsp::Exception) :
		latDescr(NULL), coreDB(NULL), minimalOutput(false)
	{
		biu::OptionMap allowedArgs;
		std::string infoText;
		initAllowedArguments(allowedArgs,infoText);	// init
			// parse programm arguments
		biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
														(char**)argv, infoText);
			// arguments parseable and all mandatory arguments given
		if (opts.noErrors()) {
				// help output
			if (opts.getBoolVal("help")) {
				opts.coutUsage();
				throw cpsp::Exception("",0);
			}
			if (opts.getBoolVal("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
				// set run mode
			if (opts.getBoolVal("count")) {
				threading.setRunMode(cpsp::gecode::GC_HPThreading::COUNT);
			} else {
				threading.setRunMode(cpsp::gecode::GC_HPThreading::ENUMERATE);
			}

				// verbose output
			threading.setVerboseOutput(opts.getBoolVal("v"));
			if (threading.getVerboseOutput() || opts.getBoolVal("stat"))
				threading.setStatOutput(true);
			minimalOutput = opts.getBoolVal("s");
			if (minimalOutput)
				threading.setVerboseOutput(false);
			threading.setNoOutput(minimalOutput);
				// using symmetryBreaking
			threading.setSymmetryBreaking(! opts.getBoolVal("noSymmBreak"));
				// using DDS for counting only
//			threading.setUsingDDS(threading.getRunMode() == cpsp::gecode::GC_HPThreading::COUNT && !opts.getBoolVal("noDDS"));
// no DDS currently in Gecode 1.3.1
			threading.setUsingDDS(false);

			threading.setNoRecomputation(opts.getBoolVal("noRecomp"));
				// core selection mode
				// <all> overwrites <allbest>
			if (opts.getBoolVal("all")) {
				threading.setCoreSelection(cpsp::gecode::GC_HPThreading::ALL);
			} else if (opts.getBoolVal("allbest")) {
				threading.setCoreSelection(cpsp::gecode::GC_HPThreading::ALL_BEST);
			}
			if (threading.getCoreSelection() == cpsp::gecode::GC_HPThreading::FIRST) {
				threading.setMaxStructures(1);
			} else {
				// maximal structure number to take into find
				int tmp = opts.getIntVal("maxSol");
				if (tmp < 1)
					throw cpsp::Exception("Error in arguments : maxSol <= 0", -2);
				threading.setMaxStructures((unsigned int) tmp);
			}

				// lattice model
			if (opts.getStrVal("lat") == "CUB")
				latDescr = new biu::LatticeDescriptorCUB();
			else if (opts.getStrVal("lat") == "FCC")
				latDescr = new biu::LatticeDescriptorFCC();
			else
				throw cpsp::Exception(
					"Error in arguments: lattice is not supported",-2);
			threading.setLatticeDescriptor(latDescr);

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
			threading.setSequence(s);

				// check if structure was given --> extract equivalence class
			if (opts.argExist("abs")) {
				  // check if move string is valid
				if (!latDescr->getAlphabet()->isAlphabetString(opts.getStrVal("abs"))) {
					throw cpsp::Exception("Error in arguments: absolute move string given is not valid in the current lattice",-2);
				}
				  // get move sequence of the move string
				givenStructure = latDescr->getSequence(opts.getStrVal("abs"));
				  // check if of correct length
				if ((givenStructure.size()+1) != s.size()) {
					throw cpsp::Exception("Error in arguments: absolute move string structure and sequence differ in length",-2);
				}
			}

				// energy limits
			int minE = INT_MIN/2+2;
			int maxE = 0;
// @TODO option wieder reinnehmen
//			minE = opts.getIntVal("minE");
//			maxE = opts.getIntVal("maxE");
				// semantische fehler pruefen
			if (minE > maxE)
				throw cpsp::Exception("Error in arguments : minE > maxE", -2);

				// contact limits for database init
			unsigned int minContacts = -maxE + threading.getContactsInSeq();
			unsigned int maxContacts = -minE + threading.getContactsInSeq();

				// H-core database
			char * tmpStr = NULL;
			if (opts.argExist("dbPath")) {
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
			if (!coreDB->initCoreAccess(	*latDescr,
											threading.getCoreSize(),
											minContacts,
											maxContacts)) {
				throw cpsp::Exception(
					"H-core database Error: access init failed",-1);
			}
			threading.setCoreDB(coreDB);
			
				// whether or not to do normalization of move strings
			threading.setNormalizeStructures(opts.getBoolVal("normalize"));

				// set distance measure mode			
			threading.setDistanceMode(cpsp::gecode::GC_HPThreading::DIST_HCORE);

			if (threading.getDistanceMode() != cpsp::gecode::GC_HPThreading::DIST_NONE
				&&  threading.getUsingDDS()) 
			{
				throw cpsp::Exception(
					"Error in arguments: distance mode != 'NONE' and using DDS is currently not supported (use -noDDS)",-2);
			}
				// set minimal distance
			if (opts.getIntVal("distMin") < 0) {
				throw cpsp::Exception(
					"Error in arguments: minimal distance = '"
					+opts.getStrVal("distMin")+"' is not >= 0",-2);
			} else {
				threading.setDistanceMin((unsigned int)opts.getIntVal("distMin"));
			}
			
		} else {
			throw cpsp::Exception(
				"Parameter Error: not parseable",-1);
		}
	}


	HPrepParameters::~HPrepParameters()
	{
			// aufraeumen
		if( coreDB != NULL) {
			delete coreDB; coreDB = NULL;
		}
		if( latDescr != NULL ) {
			delete latDescr; latDescr = NULL;
		}
	}

	void
	HPrepParameters::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(
								"seq", false, biu::COption::STRING,
								"HP-sequence to fold"));
		allowedArgs.push_back(biu::COption(
								"abs", true, biu::COption::STRING,
								"structure in absolute moves; if given, a normalized string representation of the according equivalence class is produced"));
		allowedArgs.push_back(biu::COption(
								"dbPath", true, biu::COption::STRING,
								"file based database root path (or use $CPSP_COREDB)"));
		allowedArgs.push_back(biu::COption(
								"lat", true, biu::COption::STRING,
								"lattice (CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(
								"allbest", true, biu::COption::BOOL,
								"search all solutions with best energy"));
		allowedArgs.push_back(biu::COption(
								"all", true, biu::COption::BOOL,
								"search all solutions"));
		allowedArgs.push_back(biu::COption(
								"noSymmBreak", true, biu::COption::BOOL,
								"do not break symmetries"));
		allowedArgs.push_back(biu::COption(
								"normalize", true, biu::COption::BOOL,
								"print structures as normalized move strings"));
		allowedArgs.push_back(biu::COption(
								"maxSol", true, biu::COption::INT,
							"maximum of solutions to search for (value > 0)",
								biu::util::Util_String::int2str(INT_MAX-2)));
		allowedArgs.push_back(biu::COption(
								"noRecomp", true, biu::COption::BOOL,
								"no space recomputation"));
		allowedArgs.push_back(biu::COption(
								"count", true, biu::COption::BOOL,
								"only counting structures"));
		allowedArgs.push_back(biu::COption(
								"s", true, biu::COption::BOOL,
								"silent - minimal output"));
		allowedArgs.push_back(biu::COption(
								"v", true, biu::COption::BOOL,
								"verbose output"));
		allowedArgs.push_back(biu::COption(
								"stat", true, biu::COption::BOOL,
								"print search tree statistics"));
		allowedArgs.push_back(biu::COption(
								"help", true, biu::COption::BOOL,
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPrep calculates a representive subset of optimal structures of a given ")
			+ std::string("sequence in the 3D-lattice HP-model.")
			+ std::string(" All structures differ in at least one H-monomer position.")
			+ std::string(" The structures are given in ")
			+ std::string("absolute move representation whereby we use the following ")
			+ std::string("encoding : F/B +-x, L/R +- y, U/D +-z.")
			+ std::string("\nTo find them an")
			+ std::string(" H-core database is used with optimal and suboptimal cores.")
			+ std::string(" Beginning with the biggest possible core iteratively all")
			+ std::string(" cores up to a given suboptimality level are taken into account")
			+ std::string(" to design the sequences.")			
			+ std::string("\n\n - the option <all> overwrites the <allbest> option");

	} // initArguments

//////////////////////////////////////////////////////////////////
//  MAIN FUNCTION
//////////////////////////////////////////////////////////////////
	int main(int argc, char**argv) {
		using namespace cpsp;
		using namespace cpsp::gecode;

		try {

			// generate threading object with user parameter settings
			cpsp::gecode::GC_HPThreading threading;

			// check user parameters and init threading object
			HPrepParameters params(argc,argv, threading);

			if (!params.silent()) {
				std::cout	<<"\n=================================="
							<<"\n   HPrep - CPSP-tool-library"
							<<"\n=================================="
							<<"\n\n";
			}

			  // check if we have to convert or to generate
			if (params.moveStringGiven()) {

				// convert to class
				std::cout <<threading.getEquivClassString( threading.getSequence(), params.getGivenMoveSequence() );

			} else {

				// generate classes
				unsigned int numOfSeqs = threading.generateStructures();

				if (params.silent() && threading.getRunMode() == cpsp::gecode::GC_HPThreading::COUNT)
					std::cout <<numOfSeqs <<std::endl;

			}

			if (!params.silent()) {
				std::cout	<<"\n==================================\n"
							<<std::endl;
			}

		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}


		return 0;
	}
