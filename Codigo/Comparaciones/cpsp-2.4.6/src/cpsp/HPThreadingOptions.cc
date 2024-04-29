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

#include "cpsp/HPThreadingOptions.hh"

#include <limits.h>

#include <biu/util/Util_String.h>

#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>

#include "cpsp/HCoreDatabaseFILE.hh"
//#include "cpsp/HCoreDatabasePGSQL.hh"


namespace cpsp
{

	HPThreadingOptions::HPThreadingOptions()  :
		sequence(NULL), latDescr(NULL), lattice(NULL), 
		coreSelection(FIRST), minE(INT_MIN/2+2), maxE(0), maxSols(0),
		hullLvl(), maxPCoreDistance(0), coreDB(NULL),
		onlyCountingStructures(false), verboseOutput(false),
		minEvenOddHs(0),
		contactsInSeq(0),
		useSymmetryBreaking(true),
		useNoDDS(false),
		normalizeStructures(false)
		
	{
	}

	HPThreadingOptions::HPThreadingOptions(int argc, char** argv) 
		throw(cpsp::Exception) :
		sequence(NULL), latDescr(NULL), lattice(NULL), 
		coreSelection(FIRST), minE(INT_MIN/2+2), maxE(0), maxSols(0),
		hullLvl(), maxPCoreDistance(0), coreDB(NULL),
		onlyCountingStructures(false), verboseOutput(false),
		minEvenOddHs(0),
		contactsInSeq(0),
		useSymmetryBreaking(true),
		useNoDDS(false),
		normalizeStructures(false)
		
	{
			// init allowed arguments
		biu::OptionMap allowedArgs;
		std::string infoText;
		initAllowedArguments(allowedArgs,infoText);
		
			// parse programm arguments	
		biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
														(char**)argv, infoText);
		
		if (opts.noErrors()) {
	
			if (opts.getBoolVal("help")) {
				throw cpsp::Exception(
					"",0);
			}
			
				// mode flags 
			verboseOutput = opts.getBoolVal("v");
			onlyCountingStructures = opts.getBoolVal("count");
			
				// sequence
			std::string tmpSeq = opts.getStrVal("seq");
			unsigned int coreSize = setSequence(tmpSeq);
						
				// lattice modell
			if (opts.getStrVal("lat") == "CUB")
				latDescr = new biu::LatticeDescriptorCUB();
			else if (opts.getStrVal("lat") == "FCC")
				latDescr = new biu::LatticeDescriptorFCC();
			else 
				throw cpsp::Exception(
					"Error in arguments: lattice is not supported",-2);
			
				
			lattice = new biu::LatticeFrame(latDescr, 2*sequence->size()+4);
			
				// core selection mode
				// <all> overwrites <allbest> if both given
			if (opts.getBoolVal("all")) {
				coreSelection = ALL;
			} else if (opts.getBoolVal("allbest")) {
				coreSelection = ALL_BEST;
			}
			if (coreSelection == FIRST) {
				maxSols = 1;
			}
				// energie bounds
// @TODO option wieder reinnehmen				
//			minE = opts.getIntVal("minE");
//			maxE = opts.getIntVal("maxE");
			
				// contact bounds for database access
			unsigned int minContacts = -maxE + contactsInSeq;
			unsigned int maxContacts = -minE + contactsInSeq;
			
				// check for semantic errors
			if (minE > maxE)
				throw cpsp::Exception("Error in arguments : minE > maxE", -2);
				
				// maximal structure number to enumerate
			if (coreSelection != FIRST) {
				int tmp = opts.getIntVal("maxSol");
				if (tmp < 1)
					throw cpsp::Exception("Error in arguments : maxSol <= 0", -2);
				maxSols = (unsigned int) tmp;
			}

			// use symmetry breaking
			useSymmetryBreaking = ! opts.getBoolVal("noSymmBreak");
			
			// use no dynamic decomposition search (DDS)
			// use DDS only in counting mode
			useNoDDS = !onlyCountingStructures || opts.getBoolVal("noDDS");
			
			
			normalizeStructures = opts.getBoolVal("normalize");

				// H-core database

			if (opts.argExist("dbPath")) {
				coreDB = new HCoreDatabaseFILE(opts.getStrVal("dbPath"));
			} else {
				throw cpsp::Exception(
					"H-core database Error: no database available (FILE or SQL)",-1);
			}
				// test connection
			if (!coreDB->isConnected()) {
				throw cpsp::Exception(
					"H-core database Error: can't open database",-1);
			}
/*
			std::cout <<"coreSize = " <<coreSize <<std::endl;
			std::cout <<"minContacts = " <<minContacts <<std::endl;
			std::cout <<"maxContacts = " <<maxContacts <<std::endl;
			std::cout <<"latDescr = " <<latDescr->getName() <<std::endl;
*/			
			if (!coreDB->initCoreAccess(	*latDescr,
											coreSize,
											minContacts,
											maxContacts)) {
				throw cpsp::Exception(
					"H-core database Error: access init failed",-1);
			}
		} else {
			throw cpsp::Exception(
				"Parameter Error: not parseable",-1);
		}
	}
	
	HPThreadingOptions::~HPThreadingOptions() 
	{
			// aufraeumen
		if( coreDB != NULL) {
			delete coreDB; coreDB = NULL;
		}
		if( lattice != NULL ) {
			delete lattice; lattice = NULL;
		}
		if( latDescr != NULL ) {
			delete latDescr; latDescr = NULL;
		}
		if( sequence != NULL ) {
			delete sequence; sequence = NULL;
		}
	}
	
	unsigned int 
	HPThreadingOptions::setSequence(const std::string &seq) throw(cpsp::Exception)
	{
		sequence = new std::string(biu::util::Util_String::str2upperCase(seq));
			// check if valid
		biu::Alphabet hpAlph("HP",1);
		if (!hpAlph.isAlphabetString(*sequence) || sequence->size() < 3) {
			throw cpsp::Exception(
				"Error in arguments: sequence is no valid HP-sequence",-2);
		}
		
		size_t i=0, pHead=0, pTail=0, pSubSeq=0, actLvl=0;

			// get P-head und P-tail
		pHead = sequence->find_first_not_of('P');
		pTail = sequence->find_last_of('P');
		pTail = (pTail == std::string::npos)? 0 : sequence->size() - pTail;
		
			// core size init
		unsigned int coreSize=(pHead<sequence->size())?1:0;
		if (coreSize == 1 && pHead%2==0)	// first H-pos is even
			minEvenOddHs = 1;
			// hull level
		hullLvl = std::vector<unsigned int>(sequence->size(),0);	// init
			// set P-head
		for (i=0; i<pHead; i++) {
			hullLvl[i] = pHead-i;
		}
			// init first hullLvl 
		for (i=pHead+1; i<sequence->size(); i++) {
			if (sequence->at(i) == 'H') {				// end of p-subseq
				pSubSeq = std::max(pSubSeq, actLvl);	// set p-subseq 
				actLvl = 0;								// reinit
				coreSize++;								// kern size ++
				if (i%2==0) minEvenOddHs++;				// even H position
			} else {
				actLvl++;
			}
			hullLvl[i] = actLvl;			// set hull level 
		}
			// adjust backwards
		for (i=sequence->size()-pTail-1; i>pHead+1; i--) {
			if (hullLvl[i-1] > (hullLvl[i]+1)) {
				hullLvl[i-1] = hullLvl[i]+1;
			}
		}
			// get min of even and odd H positions
		minEvenOddHs = std::min( minEvenOddHs, coreSize-minEvenOddHs);
		
			// get singlet position
			// P-singlets
		i = sequence->find("PHP");
		seqFeatureMap[H_SINGLET].clear();
		while ( i != std::string::npos ) {
			seqFeatureMap[H_SINGLET].push_back(i+1);
			i = sequence->find("PHP",i+2);
		}
			// h-singlets
		i = sequence->find("HPH");
		seqFeatureMap[P_SINGLET].clear();
		while ( i !=  std::string::npos ) {
			seqFeatureMap[P_SINGLET].push_back(i+1);
			i = sequence->find("HPH",i+2);
		}
		
			// sequence internal HH-contacts
		contactsInSeq = 0;
		i = sequence->find("HH");
		while ( i !=  std::string::npos ) {
			contactsInSeq++;
			i = sequence->find("HH", i+1);
		}
		
			// set lattice frame size 
		maxPCoreDistance = 2* (	std::max(	(pSubSeq/2) +1, 
											std::max( pTail, pHead))
								+2);
		
		return coreSize;
	}
	
	
	void
	HPThreadingOptions::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(	
								"dbPath", false, biu::COption::STRING, 
								"file based database root path"));
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"HP-sequence to fold"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"all", true, biu::COption::BOOL, 
								"search all solutions"));
		allowedArgs.push_back(biu::COption(	
								"noSymmBreak", true, biu::COption::BOOL, 
								"do not break symmetries"));
		allowedArgs.push_back(biu::COption(	
								"noDDS", true, biu::COption::BOOL, 
								"do not use Dynamic Decomposition Search"));
		allowedArgs.push_back(biu::COption(	
								"allbest", true, biu::COption::BOOL, 
								"search all solutions with best energy"));
/*		allowedArgs.push_back(biu::COption(	
								"minE", true, biu::COption::INT, 
					"lower energy bound for the search (value <= maxE <= 0)",
								biu::util::Util_String::int2str(INT_MIN/2+2)));
		allowedArgs.push_back(biu::COption(	
								"maxE", true, biu::COption::INT, 
							"upper energy bound for the search (value <= 0)",
								"0"));
*/		allowedArgs.push_back(biu::COption(	
								"maxSol", true, biu::COption::INT, 
							"maximum of solutions to search for (value > 0)",
								biu::util::Util_String::int2str(INT_MAX-2)));
		allowedArgs.push_back(biu::COption(	
								"count", true, biu::COption::BOOL, 
								"only counting structures"));
		allowedArgs.push_back(biu::COption(	
								"normalize", true, biu::COption::BOOL, 
								"print structures as normalized move strings"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));

		infoText = " - the option <all> overwrites the <allbest> option";
	
	} // initArguments

} // namespace cpsp
