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


#include "HPdeg.hh"
#include "version.hh"
#include <cstdlib>
#include <iostream>
#include <biu/util/Util_String.h>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <cpsp/HCoreDatabaseFILE.hh>
#include <cpsp/gecode/GC_HPThreading.hh>

//////////////////////////////////////////////////////////////////
// THE PARAMETER PARSER FOR GC_HPTHREADING
//////////////////////////////////////////////////////////////////

	HPdegParameters::HPdegParameters(int argc, char** argv, 
		cpsp::gecode::GC_HPThreading & threading) throw(cpsp::Exception) :
		latDescr(NULL), coreDB(NULL), silent(true)	
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
			threading.setRunMode(cpsp::gecode::GC_HPThreading::COUNT);
			
				// verbose output
			threading.setVerboseOutput(opts.getBoolVal("v"));
			silent = opts.getBoolVal("s");
			threading.setNoOutput(true);
				// using symmetryBreaking
			threading.setSymmetryBreaking(! opts.getBoolVal("noSymmBreak"));
				// using DDS for counting only
//			threading.setUsingDDS(!opts.getBoolVal("noDDS"));
// no DDS in current Gecode 1.3.1
			threading.setUsingDDS(false);

				// core selection mode
			threading.setCoreSelection(cpsp::gecode::GC_HPThreading::ALL_BEST);
				// maximal structure number to find
			threading.setMaxStructures(UINT_MAX-2);
			
				// lattice modell
			if (opts.getStrVal("lat") == "CUB")
				latDescr = new biu::LatticeDescriptorCUB();
			else if (opts.getStrVal("lat") == "FCC")
				latDescr = new biu::LatticeDescriptorFCC();
			else 
				throw cpsp::Exception(
					"Error in arguments: lattice is not supported",-2);
			threading.setLatticeDescriptor(latDescr);
			
				// sequence
			if ( opts.getStrVal("seq").size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no sequence given",-2);
			
			threading.setSequence(opts.getStrVal("seq"));
			
				// energie limits
			int minE = INT_MIN/2+2;
			int maxE = 0;
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
				// test
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

		} else {
			throw cpsp::Exception(
				"Parameter Error: not parseable",-1);
		}

	}


	HPdegParameters::~HPdegParameters()
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
	HPdegParameters::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"HP-sequence to fold"));
		allowedArgs.push_back(biu::COption(	
								"dbPath", true, biu::COption::STRING, 
								"file based database root path (or use $CPSP_COREDB)"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"noSymmBreak", true, biu::COption::BOOL, 
								"do not break symmetries"));
/*
		allowedArgs.push_back(biu::COption(	
								"noDDS", true, biu::COption::BOOL, 
								"do not use Dynamic Decomposition Search"));
*/		
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"s", true, biu::COption::BOOL, 
								"silent - minimal output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPdeg calculates the degeneracy of a given")
			+ std::string(" HP-sequence, that is the number of optimal structures")
			+ std::string(" in the 3D-lattice HP-model the sequence can adopt.")
			+ std::string("\nTo find them an")
			+ std::string(" H-core database is used with optimal and suboptimal cores.")
			+ std::string(" Beginning with the biggest possible core iteratively all")
			+ std::string(" cores up to a given suboptimality level are taken into account")
			+ std::string(" to design the sequences.");			
	
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
			HPdegParameters params(argc,argv, threading);
			
			if (!params.isSilent()) {
				std::cout	<<"\n=================================="
							<<"\n   HPdeg - CPSP-tool-library"
							<<"\n=================================="
							<<"\n\n";
			}

			// generate structures to get degeneracy 
			unsigned int deg = threading.generateStructures();
			
			// get energy
			int energy = 0;
			threading.setMaxStructures(2);
			HCoreDatabase *coreDB = threading.getCoreDB();
			coreDB->initCoreAccess(	*threading.getLatticeDescriptor(), 
									coreDB->getActCoreSize(),
									coreDB->getActMinHHcontacts(),
									coreDB->getActMinHHcontacts());
			std::vector<std::string> abs;
			threading.enumerateStructures(abs);
			if (abs.size() > 0) {
				biu::LatticeModel lattice(threading.getLatticeDescriptor());
				biu::IPointVec points = lattice.absMovesToPoints(lattice.parseMoveString(abs[0]));
				std::string seq = threading.getSequence();
				for (std::string::size_type i=0; i<seq.size(); i++) {
					if (seq[i] != 'H')
						continue;
					for (std::string::size_type j=i+2; j<seq.size(); j++) {
						if (seq[j] != 'H')
							continue;
						energy -= ( lattice.areNeighbored(points[i],points[j]) ? 1 : 0 );  
					}
				}
			}
			
			if (params.isSilent()) {
				std::cout <<deg <<std::endl;
			} else {
				std::cout	<<" sequence   : " <<threading.getSequence() <<"\n"
							<<" lattice    : " <<threading.getLatticeDescriptor()->getName() <<"\n"
							<<" energy     : " <<energy <<"\n"
							<<" degeneracy : " <<deg <<"\n";
				std::cout.flush();
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

