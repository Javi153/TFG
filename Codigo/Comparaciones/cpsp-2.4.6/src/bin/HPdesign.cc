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

#include "HPdesign.hh"
#include "version.hh"

#include <iomanip>
#include <biu/Timer.hh>
#include <biu/LatticeFrame.hh>
#include <cpsp/HCoreDatabaseFILE.hh>
#include <cpsp/gecode/GC_HPThreading.hh>

namespace cpsp
{
	
	/////////////////////////////////////////////////
	// HPdesignOptions implementations
	/////////////////////////////////////////////////


	HPdesignOptions::HPdesignOptions(int argc, char** argv) 
		throw(cpsp::Exception) :
		opts(NULL), lattice(NULL), structure(NULL), coreDB(NULL), subOptLvl(0),
		maxDegeneracy(10), verboseOut(false), minHs(2), maxHs(99999),
		all(false), maxSeqs(1), minimalOut(false), structsRel()
	{
			// init allowed arguments
		biu::OptionMap allowedArgs;
		std::string infoText;
		initAllowedArguments(allowedArgs,infoText);
		
			// parse programm arguments	
		opts = new biu::COptionParser(	allowedArgs, argc, 
										(char**)argv, infoText);
		
		if (opts->noErrors()) {
	
			if (opts->getBoolVal("help")) {
				opts->coutUsage();
				throw cpsp::Exception("",0);
			}
			if (opts->getBoolVal("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
			
			verboseOut = opts->getBoolVal("v");
			minimalOut = opts->getBoolVal("s");
			if (minimalOut)
				verboseOut = false;
			all = opts->getBoolVal("all");
			if (all) {
				maxSeqs = UINT_MAX-2;
			}

				// lattice modell
			if (opts->getStrVal("lat") == "CUB")
				latDescr = new biu::LatticeDescriptorCUB();
			else if (opts->getStrVal("lat") == "FCC")
				latDescr = new biu::LatticeDescriptorFCC();
			else { 
				throw cpsp::Exception(
					"Error in arguments: lattice is not supported",-2);
			}
			lattice = new biu::LatticeModel(latDescr);
			
				// structure
			if (opts->argExist("relMoves") ) {
				structure = new biu::IPointVec(
									lattice->relMovesToPoints(
										lattice->parseMoveString(
											opts->getStrVal("struct")
									)));
			} else {
				structure = new biu::IPointVec(
									lattice->absMovesToPoints(
										lattice->parseMoveString(
											opts->getStrVal("struct")
									)));
			}
			
				// the set of relative move strings to structure
			const biu::AutomorphismVec autM = lattice->getDescriptor()->getAutomorphisms();
			biu::AutomorphismVec::const_iterator a;
			biu::IPointVec act(*structure);
			for (a = autM.begin(); a!= autM.end(); a++) {
				for (biu::IPointVec::size_type pos=0; pos < structure->size(); pos++) {
					act[pos] = (*a) * (structure->at(pos));
				}
				structsRel.insert(lattice->getString(lattice->pointsToRelMoves(act)));
			}
			
			

				// H-core database
			char * tmpStr = NULL;
			if (opts->argExist("dbPath")) {
				coreDB = new cpsp::HCoreDatabaseFILE(opts->getStrVal("dbPath"));
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
											structure->size())) {
				throw cpsp::Exception(
					"H-core database Error: access init failed",-1);
			}
			
				// suboptimality level
			int tmp = opts->getIntVal("subOpt");
			if (tmp < 0) {
				throw cpsp::Exception(
					"Parameter Error: subOpt < 0",-1);
			} else {
				subOptLvl = (unsigned int)tmp;
			}

				// maximal degeneracy
			tmp = opts->getIntVal("maxDeg");
			if (tmp < 1) {
				throw cpsp::Exception(
					"Parameter Error: maxDeg < 1",-1);
			} else {
				maxDegeneracy = (unsigned int)tmp;
			}
			
			tmp = opts->getIntVal("minH");
			if (tmp < 2) {
				throw cpsp::Exception(
					"Parameter Error: minH < 2",-1);
			} else {
				minHs = (unsigned int)tmp;
			}
			
			tmp = opts->getIntVal("maxH");
			if (tmp < 0) {
				throw cpsp::Exception(
					"Parameter Error: maxH < 0",-1);
			} else {
				maxHs = std::min((unsigned int)tmp, (unsigned int)structure->size());
				  // check if useful limits
				if (minHs > maxHs) {
					throw cpsp::Exception(
						"Parameter Error: minH > maxH",-1);
				}
			}
			
			if (opts->argExist("maxSeq")) {
				tmp = opts->getIntVal("maxSeq");
				if (tmp < 0) {
					throw cpsp::Exception(
						"Parameter Error: maxSeq < 0",-1);
				} else {
					maxSeqs = (unsigned int)tmp;
				}
			}
			
			statistics = opts->argExist("stat");
			

		} else {
			throw cpsp::Exception(
				"Parameter Error: not parseable",-1);
		}
	}
	
	HPdesignOptions::~HPdesignOptions()
	{
		if (coreDB != NULL) delete coreDB;
		if (structure != NULL) delete structure;
		if (lattice != NULL) delete lattice;
		if (opts != NULL) delete opts;
	}

	
	bool
	HPdesignOptions::structFound(const cpsp::StringVec & structs) const {
		assert(lattice != NULL);
		bool found = false;
		for (cpsp::StringVec::const_iterator i = structs.begin(); 
				!found && i!=structs.end();
				i++) 
		{
			 // empty strings are skipped (may occure by GC_HPThreading::enumerateStructures(..) )
			if (i->compare("") == 0) continue;			
			 // get relative move string
			std::string str = lattice->getString(
								lattice->absMovesToRelMoves(
									lattice->getDescriptor()->getSequence(*i)));
			 // search if one of the symmetrics to the initial one
			found = std::find(structsRel.begin(),structsRel.end(),str) != structsRel.end();
		}
		return found;
	}
	

	void
	HPdesignOptions::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(	
								"struct", false, biu::COption::STRING, 
								"the structure as a move string"));
		allowedArgs.push_back(biu::COption(	
								"dbPath", true, biu::COption::STRING, 
								"file based database root path (or use $CPSP_COREDB)"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"relMoves", true, biu::COption::BOOL, 
								"structure given in relative moves"));
		allowedArgs.push_back(biu::COption(	
								"subOpt", true, biu::COption::INT, 
							"suboptimal core levels taken into account (value >= 0)",
								"0"));
		allowedArgs.push_back(biu::COption(	
								"maxDeg", true, biu::COption::INT, 
							"maximal degeneracy allowed for a sequence (value >= 1)",
								"1"));
		allowedArgs.push_back(biu::COption(	
								"maxSeq", true, biu::COption::INT, 
							"maximal number of sequences to design (value >= 0)"));
		allowedArgs.push_back(biu::COption(	
								"minH", true, biu::COption::INT, 
							"minimal number of Hs in sequence (value >= 2)",
								"2"));
		allowedArgs.push_back(biu::COption(	
								"maxH", true, biu::COption::INT, 
							"maximal number of Hs in sequence (value >= 0)",
								"99999"));
		allowedArgs.push_back(biu::COption(	
								"all", true, biu::COption::BOOL, 
								"find all sequences"));
		allowedArgs.push_back(biu::COption(	
								"s", true, biu::COption::BOOL, 
								"silent - minimal output"));
		allowedArgs.push_back(biu::COption(	
								"stat", true, biu::COption::BOOL, 
								"print statistics of the run"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPdesign designs for a given structure sequences that will")
				+ std::string(" maximally fold into a given number of optimal structures")
				+ std::string(" (degeneracy), the given structure within.\n\nTo find them an")
				+ std::string(" H-core database is used with optimal and suboptimal cores.")
				+ std::string(" Beginning with the biggest possible core iteratively all")
				+ std::string(" cores up to a given suboptimality level are taken into account")
				+ std::string(" to design the sequences.");
	
	} // initArguments


	/////////////////////////////////////////////////
	// HPdesign implementations
	/////////////////////////////////////////////////

	
	
	
	HPdesign::HPdesign(HPdesignOptions *opts_) :
		opts(opts_), stat(opts->doStatistic()?new Statistic():NULL)
	{
		initStructPoints();
	}
	
	HPdesign::~HPdesign()
	{
		if (stat != NULL) { delete stat; stat = NULL; }
	}
	
	void
	HPdesign::generateSeqs()
	{
			// timer for statistics
		biu::Timer timer; timer.start();

		biu::IPointVec::size_type maxHH = opts->getMaxHs(), hCount = maxHH;
		const biu::LatticeDescriptor *latD = opts->getLattice()->getDescriptor();
		unsigned int i = 0, lastMaxHH = 0, designedSeqs = 0;
			// for each number of Hs allowed (and if still seqs to design)
		while ((opts->getMaxSeqs() > designedSeqs) && hCount >= (biu::IPointVec::size_type)opts->getMinHs()) {
			lastMaxHH = INT_MAX;
				// for each suboptimality level for this core size
			for(i=opts->getSubOptLvl()+1; lastMaxHH>0 && i>0; i--) {
				lastMaxHH--;
				if (opts->getCoreDB()->initCoreAccess(*latD, hCount, 0, lastMaxHH) ) {
						// check for all best in initialised coreDB
					designedSeqs += checkCores( opts->getCoreDB(), lastMaxHH, (opts->getMaxSeqs()-designedSeqs) );
				} else {
					lastMaxHH = 0; // loop breaking
				}
			}
			hCount--;	// reduce core size
		}
		if (!opts->silent()) {
			if (designedSeqs == 0) {	
				std::cout <<"\n ==> no sequence found\n"
						<<std::endl;
			} else {
				std::cout <<"\n ==> designed sequences = " <<designedSeqs <<std::endl;
			} 
		}
		
		if (stat != NULL) {
			stat->time = timer.stop();
			std::cout	<<"\n - statistics :\n"
						<<"\n    runtime      : " <<stat->time <<" ms"
						<<"\n    seqs. tested : " <<stat->generatedSeqs
						<<"\n    cores tested : " <<stat->coresTested
						<<"\n" <<std::endl; 
		}
	}

	unsigned int
	HPdesign::checkCores(HCoreDatabase * coreDB, unsigned int & contacts, unsigned int maxSeqs) {
		HCore core;
		contacts = 0;
		unsigned int numOfSeqs = 0;
		std::set< std::string > foundSeqs;
		if (! coreDB->getNextCore(core))
			return false;
		else {
			contacts = core.getContacts();
			if (opts->verboseOutput())
				std::cout <<"\n --> next cores:  size = "<< core.getSize()
					<<", HH-contacts = " <<contacts
					<<std::endl;
			coreDB->setActMinHHcontacts(contacts);
			std::string actSeq;
			std::vector< PosSet >::const_iterator s;
			biu::IntPoint shift;
			do {
				if (stat != NULL) { stat->coresTested++; }
				for (s = structPoints.begin(); s != structPoints.end(); s++) {
					if (isSubset(*s, core.getPoints(),shift)) {
//						std::cerr <<" was gefunden ... nun noch seq ausrechnen!" <<std::endl;
						actSeq = std::string(s->size(),'P');
//						std::cout <<"shift : "<<shift <<std::endl ;
						for(biu::IPointSet::const_iterator it = core.getPoints().begin(); it != core.getPoints().end(); it++) {
							StructPos act = *(s->find(*it+shift));
							actSeq[act.getSeqPos()] = 'H';
						}
						foundSeqs.insert(actSeq);
					}
				}
			} while (coreDB->getNextCore(core));
			for(std::set< std::string >::const_iterator seqs = foundSeqs.begin(); 
				numOfSeqs < maxSeqs && seqs != foundSeqs.end(); 
				seqs++) 
			{
				if(checkSeq(*seqs))
					numOfSeqs++;
			}
		}
		return numOfSeqs;
	}


	void 
	HPdesign::initStructPoints(void) {
		const biu::AutomorphismVec autM = opts->getLattice()->getDescriptor()->getAutomorphisms();
		biu::AutomorphismVec::const_iterator a;
		structPoints.resize(autM.size());
		std::vector< PosSet >::iterator sp = structPoints.begin();
		const biu::IPointVec* orig = opts->getStructure();
		biu::IPointVec::const_iterator op;
		std::string::size_type pos = 0;
		for (a = autM.begin(); a!= autM.end(); a++) {
			pos = 0;
			for(op = orig->begin(); op != orig->end(); op++) {
				sp->insert( StructPos((*a) * (*op),pos));	// apply automorphism a to point op
				pos++;
			}
			sp++; // next automorph point set to fill
		}
	}	
	
	
	bool
	HPdesign::isSubset(const PosSet &s, const biu::IPointSet &sub, biu::IntPoint & shift) {
		if (s.size() < sub.size()) // sub has to be smaller
			return false;
		if (s.size() == 0)	// is true for empty sets
			return true;
		
//		std::cout <<"\nstruct : ";
//		{for (PosSet::const_iterator k = s.begin(); k!= s.end(); k++) std::cout <<*k <<", ";} 
//		std::cout <<"\nsub    : ";
//		{for (biu::IPointSet::const_iterator k = sub.begin(); k!= sub.end(); k++) std::cout <<*k <<", ";}
//		std::cout <<std::endl<<std::endl; 
		PosSet::size_type maxTests = s.size()-sub.size()+1;
		PosSet::const_iterator sIt = s.begin();
		biu::IPointSet::const_iterator pIt = sub.begin();
		for (;maxTests>0; maxTests--) {
			pIt = sub.begin();
			shift = *sIt - *pIt;
			for (pIt++; pIt != sub.end(); pIt++) {
				if (  s.find(*pIt + shift) == s.end() ) 
					break;
			}
			if (pIt == sub.end()) // than all points found in s
				return true;
			sIt++;
		}
			
		return false;
	}


	bool 
	HPdesign::checkSeq(const std::string& seq) 
	{
		if (stat != NULL) {
			stat->generatedSeqs++;
		}
		using namespace cpsp;
		using namespace cpsp::gecode;
		GC_HPThreading threading;
		threading.setLatticeDescriptor(opts->getLattice()->getDescriptor());
		threading.setSequence(seq);
			// testen
		if (!opts->getCoreDB()->isConnected()) {
			throw cpsp::Exception(
				"H-core database Error: can't open database",-1);
		}
		if (!opts->getCoreDB()->initCoreAccess(	*(opts->getLattice()->getDescriptor()),
										threading.getCoreSize(),
										
										threading.getContactsInSeq(),
										UINT_MAX-2)) {
			throw cpsp::Exception(
				"H-core database Error: access init failed",-1);
		}
		threading.setCoreDB(opts->getCoreDB());
		threading.setCoreSelection(GC_HPThreading::ALL_BEST);
		threading.setVerboseOutput(false);
		threading.setRunMode(GC_HPThreading::COUNT);
		threading.setUsingDDS(false);
		threading.setSymmetryBreaking(true);
		threading.setMaxStructures(opts->getMaxDegeneracy()+1);
		threading.setNoOutput(true);
		
		cpsp::StringVec structs;
//		unsigned int deg = threading.generateStructures();
		unsigned int deg = threading.enumerateStructures(structs);
		
		
		  // check if degeneracy limit is maintained
		  // check if structure is among the optimal ones
		if (	deg > 0 && deg <= opts->getMaxDegeneracy()
				&& opts->structFound(structs))
		{
			if (!opts->silent()) {
				std::cout <<"  " <<seq <<" : degeneracy = " <<std::setw(6)  <<deg
					<<",  energy = -"  <<std::setw(5) <<(opts->getCoreDB()->getActMinHHcontacts()-threading.getContactsInSeq()) 
					<<std::endl;
			} else {
				std::cout <<seq <<" " <<std::setw(6) <<deg <<std::endl;
			}
			return true;
		}
		return false; 
	}	

}

	/////////////////////////////////////////////////
	// main
	/////////////////////////////////////////////////

	int
	main(int argc, char **  argv) {
		using namespace cpsp;
		try {
			
			// check user parameters
			HPdesignOptions opts(argc,argv);
			
			if (!opts.silent()) {
				std::cout	<<"\n=================================="
							<<"\n   HPdesign - CPSP-tool-library"
							<<"\n=================================="
							<<"\n" <<std::endl;
			}
			// generate deisgn object with user parameter settings
			HPdesign hpDesign(&opts);
			
			// generate sequences 
			hpDesign.generateSeqs();
			
			if (!opts.silent()) {
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
