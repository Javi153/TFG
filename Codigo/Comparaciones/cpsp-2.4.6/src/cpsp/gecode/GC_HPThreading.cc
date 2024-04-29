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



#include "GC_HPThreading.hh"
#include <biu/Timer.hh>
#include <biu/util/Util_String.h>
#include "cpsp/HCoreDatabase.hh"
#include "GC_ThreadingSpace.hh"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>

namespace cpsp
{
  namespace gecode
  {
	GC_HPThreading::GC_HPThreading(const biu::LatticeDescriptor* latDescr_,
			HCoreDatabase* db, const std::string seq) throw(cpsp::Exception)
	  :
		cubicLatDescr(),
		verboseOutput(false),
		latDescr(NULL),
		latFrame(NULL),
		coreDB(NULL),
		sequence(""),
		hullLvl(),
		coreSize(0),
		minEvenOddHs(0),
		contactsInSeq(0),
		maxPCoreDistance(0),
		seqFeatureMap(),
		useSymmetryBreaking(true),
		useDDS(false),
		maxStructures(1),
		coreSelection(FIRST),
		runMode(ENUMERATE),
		noOutput(false),
		statisticsOut(false),
		noRecomputation(false),
		distanceMode(DIST_NONE),
		normalizeStructures(false)
	{
		setLatticeDescriptor(latDescr_);
		setSequence(seq);
		setCoreDB(db);
	}

	GC_HPThreading::~GC_HPThreading()
	{
		if(latFrame != NULL) { delete latFrame; latFrame = NULL; }

	}

	void
	GC_HPThreading::setLatticeDescriptor( const biu::LatticeDescriptor* val )
	{
		latDescr = val;
		if(latFrame != NULL) { delete latFrame; latFrame = NULL; }
		if (latDescr != NULL)
			latFrame = new biu::LatticeFrame(latDescr, 1);
	}

	void
	GC_HPThreading::setSequence( const std::string& seq ) throw(cpsp::Exception)
	{
		sequence = biu::util::Util_String::str2upperCase(seq);
		coreSize = 0;
		minEvenOddHs = 0;
		contactsInSeq = 0;
		maxPCoreDistance = 0;

		if (sequence.size() == 0)
			return;
			// gueltigkeit pruefen
		biu::Alphabet hpAlph("HP",1);
//		if (!hpAlph.isAlphabetString(sequence) || sequence.size() < 3) {
		if (!hpAlph.isAlphabetString(sequence)) {
			throw cpsp::Exception(
				"Error in arguments: sequence is no valid HP-sequence",-2);
		}

		size_t i=0, pHead=0, pTail=0, pSubSeq=0, actLvl=0;

			// p-kopf und p-schwanz bestimmen
		pHead = sequence.find_first_not_of('P');
		if (pHead == std::string::npos) pHead = 0;
		pTail = sequence.find_last_not_of('P');
		pTail = (pTail == std::string::npos || pTail <= pHead)? 0 : sequence.size() - pTail;

			// core size init
		coreSize=(pHead<sequence.size())?1:0;
		if (coreSize == 1 && pHead%2==0)	// dann erste H-pos gerade
			minEvenOddHs = 1;
		else
			minEvenOddHs = 0;
			// huellen level
		hullLvl = GC_ThreadingSpace::HullLevel(sequence.size(),0);	// init
			// p-kopf setzen
		for (i=0; i<pHead; i++) {
			hullLvl[i] = pHead-i;
		}
			// erster hullLvl init
		for (i=pHead+1; i<sequence.size(); i++) {
			if (sequence.at(i) == 'H') {				// ende von p-subseq
				pSubSeq = std::max(pSubSeq, actLvl);	// p-subseq setzen
				actLvl = 0;								// reinit
				coreSize++;								// kerngroesse ++
				if (i%2==0) minEvenOddHs++;				// gerade H position
			} else {
				actLvl++;
			}
			hullLvl[i] = actLvl;			// level setzen
		}
			// rueckwaerts angleichen
		for (i=sequence.size()-pTail-1; i>pHead+1; i--) {
			if (hullLvl[i-1] > (hullLvl[i]+1)) {
				hullLvl[i-1] = hullLvl[i]+1;
			}
		}
			// min von even und odd H positionen bestimmen
		minEvenOddHs = std::min( minEvenOddHs, coreSize-minEvenOddHs);

			// singlet positionen bestimmen
			// P-singlets
		i = sequence.find("PHP");
		seqFeatureMap[GC_ThreadingSpace::H_SINGLET].clear();
		while ( i != std::string::npos ) {
			seqFeatureMap[GC_ThreadingSpace::H_SINGLET].push_back(i+1);
			i = sequence.find("PHP",i+2);
		}
			// h-singlets
		i = sequence.find("HPH");
		seqFeatureMap[GC_ThreadingSpace::P_SINGLET].clear();
		while ( i !=  std::string::npos ) {
			seqFeatureMap[GC_ThreadingSpace::P_SINGLET].push_back(i+1);
			i = sequence.find("HPH",i+2);
		}

			// sequence internal HH-contacts
		contactsInSeq = 0;
		i = sequence.find("HH");
		while ( i !=  std::string::npos ) {
			contactsInSeq++;
			i = sequence.find("HH", i+1);
		}

		maxPCoreDistance = 2* (	std::max(	(pSubSeq/2) +1,
											std::max( pTail, pHead))
								+2);

	}
	

		//! enumerates structures and returns their number
	unsigned int
	GC_HPThreading::enumerateStructures(SuperSpace * space,
		unsigned int maxNum, Gecode::Search::Statistics *stat, 
		std::ostream &out)
	{
			// no enumeration using DDS currently
		if (useDDS) 
			return 0;
			// solution number
		unsigned int sols = 0;
		switch (distanceMode) {
			case DIST_NONE : {				
					GC_ThreadingSpace* s = dynamic_cast<GC_ThreadingSpace*>(space);
					assert( s != NULL );
					return enumerateStructures<Gecode::DFS<GC_ThreadingSpace>,GC_ThreadingSpace> 
						( s, maxNum, stat, out); 
					break;
				}
			case DIST_SHAPES : {
					GC_ThreadingSpaceShapes* s = dynamic_cast<GC_ThreadingSpaceShapes*>(space);
					assert( s != NULL );
					return enumerateStructures<Gecode::BAB<GC_ThreadingSpaceShapes>,GC_ThreadingSpaceShapes> 
						( s, maxNum, stat, out); 
					break;
				}
			case DIST_POSITIONS : {
					GC_ThreadingSpacePosDist* s = dynamic_cast<GC_ThreadingSpacePosDist*>(space);
					assert( s != NULL );
					return enumerateStructures<Gecode::BAB<GC_ThreadingSpacePosDist>,GC_ThreadingSpacePosDist> 
						( s, maxNum, stat, out); 
					break;
				}
			case DIST_ABSMOVES : {
					GC_ThreadingSpaceMoveDist* s = dynamic_cast<GC_ThreadingSpaceMoveDist*>(space);
					assert( s != NULL );
					return enumerateStructures<Gecode::BAB<GC_ThreadingSpaceMoveDist>,GC_ThreadingSpaceMoveDist> 
						( s, maxNum, stat, out); 
					break;
				}
			case DIST_HCORE : {
					GC_ThreadingSpaceHcoreDist* s = dynamic_cast<GC_ThreadingSpaceHcoreDist*>(space);
					assert( s != NULL );
					return enumerateStructures<Gecode::BAB<GC_ThreadingSpaceHcoreDist>,GC_ThreadingSpaceHcoreDist> 
						( s, maxNum, stat, out); 
					break;
				}
			default : 
					assertbiu(false,"distance mode not supported");
					break;
		}
		return sols;
	}

	unsigned int
	GC_HPThreading::countStructures(SuperSpace * space,
		unsigned int maxNum, bool useDDS, Gecode::Search::Statistics *actStat)
	{
		switch(distanceMode) {

				// normal counting
			case DIST_NONE : {
				GC_ThreadingSpace* s = dynamic_cast<GC_ThreadingSpace*>(space);
				assert(s != NULL);
				if (useDDS) {	// using DDS
					assertbiu(false,"DDS not supported");
				} else {		// using DFS
					return Gecode::countDFS<unsigned int, GC_ThreadingSpace> (
								s,
    							(noRecomputation?1:Gecode::Search::Config::c_d),
    							(noRecomputation?1:Gecode::Search::Config::a_d),
								actStat, NULL,
								maxNum
							);
				}
				break;
			}

				// normal counting
			case DIST_SHAPES : {
				GC_ThreadingSpaceShapes* s = dynamic_cast<GC_ThreadingSpaceShapes*>(space);
				assert(s != NULL);
				if (useDDS) {	// using DDS
					assertbiu(false,"DDS not supported");
				} else {		// using BAB
					return countBAB<unsigned int, GC_ThreadingSpaceShapes> (
								s,
    							(noRecomputation?1:Gecode::Search::Config::c_d),
    							(noRecomputation?1:Gecode::Search::Config::a_d),
								actStat, NULL,
								maxNum
							);
				}
				break;
			}

				// BAB counting introducing distance constraints to already 
				// found solutions on positions		
			case DIST_POSITIONS : {
				GC_ThreadingSpacePosDist* s = dynamic_cast<GC_ThreadingSpacePosDist*>(space);
				assert(s != NULL);
				if (useDDS) {
					assertbiu(false,"DDS not supported");
				} else {
						// using branch and bound
					return countBAB<unsigned int, GC_ThreadingSpacePosDist> ( 
								s,
    							(noRecomputation?1:Gecode::Search::Config::c_d),
    							(noRecomputation?1:Gecode::Search::Config::a_d),
								actStat, NULL,
								maxNum
							);
				} 
				break;
			}
				// BAB counting introducing distance constraints to already 
				// found solutions on absolute moves				
			case DIST_ABSMOVES : {
				GC_ThreadingSpaceMoveDist* s = dynamic_cast<GC_ThreadingSpaceMoveDist*>(space);
				assert(s != NULL);
				if (useDDS) {
					assertbiu(false,"DDS not supported");
				} else {
						// using branch and bound
					return countBAB<unsigned int, GC_ThreadingSpaceMoveDist> ( 
								s,
    							(noRecomputation?1:Gecode::Search::Config::c_d),
    							(noRecomputation?1:Gecode::Search::Config::a_d),
								actStat, NULL,
								maxNum
							);
				} 
				break;
			}
				// BAB counting introducing distance constraints to already 
				// found solutions on absolute moves				
			case DIST_HCORE : {
				GC_ThreadingSpaceHcoreDist* s = dynamic_cast<GC_ThreadingSpaceHcoreDist*>(space);
				assert(s != NULL);
				if (useDDS) {
					assertbiu(false,"DDS not supported");
				} else {
						// using branch and bound
					return countBAB<unsigned int, GC_ThreadingSpaceHcoreDist> ( 
								s,
    							(noRecomputation?1:Gecode::Search::Config::c_d),
    							(noRecomputation?1:Gecode::Search::Config::a_d),
								actStat, NULL,
								maxNum
							);
				} 
				break;
			}
			default : { 
					assertbiu(false,"distance mode not supported");
				break;
			}
		}
		return 0;
	}

	unsigned int
	GC_HPThreading::generateStructures(void) {

		assert(latFrame != NULL);
		assert(coreDB != NULL);
		
		if (sequence.size() == 0)
			return 0;
		
		HCore core;
		unsigned int sols = 0, oldSols = 0, oldContacts = 0,
			oldContactSols = 0, maxContactsFound = 0;

		if (!noOutput)
			std::cout <<"==> HP-sequence = " <<sequence <<std::endl;

		Gecode::Search::Statistics *globStat= NULL, *actStat=NULL;
		if (verboseOutput) {
			std::cout << "\ngenerating structures via Gecode\n" << std::endl;
		}
		if (statisticsOut) {
			globStat = new Gecode::Search::Statistics;
			actStat = new Gecode::Search::Statistics;
		}
			// solution storage for enumeration with minimal distance
		GC_ThreadingSpacePosDist::SolutionVec* solPos = NULL, *solMove = NULL;

			// timer setzen
		biu::Timer timer; timer.start();
		biu::Timer sTimer;
		double searchTime= 0.0;

		GC_ThreadingSymmBreaker::GlobalShiftVec* symmetryShiftVec = NULL;
			// neighborhood vectors for neighboring constraints
		biu::LatticeFrame::index_set indexedNeighVecs;
		 
			// set maximal frame size needed
		latFrame->setFrameSize(	maxPCoreDistance
								+ coreDB->getActCoreSize());

		while (oldSols < maxStructures && coreDB->getNextCore(core)) {
			try {
					// guide output and energy level control
				if (core.getContacts() > oldContacts) {
					oldContacts = core.getContacts();
					if (!noOutput)
						std::cout	<<"\n==> energy = -"
									<<(oldContacts - contactsInSeq)
									<<std::endl;
				} 
				
				if (core.getContacts() < oldContacts) {
					if (!noOutput)
						std::cout	<<"\n ==> number generated structures = "
									<<(oldSols - oldContactSols)
									<<std::endl;
						// number of solutions since last energy level change
					oldContactSols = oldSols;
					oldContacts = core.getContacts();
					if (!noOutput)
						std::cout	<<"\n==> energy = -"
									<<(oldContacts - contactsInSeq)
									<<std::endl;
				}
					// set minimal frame size
				if (distanceMode == DIST_NONE ) {	
						// not possible if distance of solution is forced 
						// (different positions indices etc.)  
					latFrame->setFrameSize(	maxPCoreDistance
											+ core.getMaxDimension());
				}

				if (verboseOutput) {
					std::cout <<" => next core ("<<core.getContacts()
						<<" contacts, "
						<<core.getNumOfPSinglets(latFrame)
						<<" P-singlets) "
						<<std::endl;
				}
					// naechster, wenn nicht genug singlet-positionen vorhanden
				if ( core.getNumOfPSinglets(latFrame)
							< (int)  seqFeatureMap[GC_ThreadingSpace::P_SINGLET].size()
					|| core.getNumOfHSinglets(latFrame)
							< (int)  seqFeatureMap[GC_ThreadingSpace::H_SINGLET].size())
				{
					if (verboseOutput) {
						std::cout	<<"     --> skipped (singlet positions)"
									<<std::endl;
					}
					continue;
				}
					// naechster, wenn cubisches gitter und falsches verhaeltnis
					// zwischen gerade und ungeraden kernpositionen
				if ( latFrame->getDescriptor()->getName() == "cub"
					&& core.getMinEvenOdd() != minEvenOddHs )
				{
					if (verboseOutput) {
						std::cout	<<"     --> skipped (cubic even/odd ratio) "
									<<std::endl;
					}
					continue;
				}

					// generate global shift vectors for symmetry breaking
				if (useSymmetryBreaking) {
					symmetryShiftVec = GC_ThreadingSymmBreaker::generateGlobalShiftVec(
										&core, latFrame, symmetryShiftVec);
				}
				
				indexedNeighVecs.clear();
				indexedNeighVecs = latFrame->getIndexedNeighborhood();

					// new space
				GC_ThreadingSpace* s = NULL;
				int branching = BR_NONE;
				switch (distanceMode) {
					case DIST_NONE : 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpace(&sequence,
						latFrame, &indexedNeighVecs, &seqFeatureMap, &hullLvl, &core,
						symmetryShiftVec, branching);
						break;
					case DIST_SHAPES : 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceShapes(&sequence,
						latFrame, &indexedNeighVecs, &seqFeatureMap, &hullLvl, &core,
						symmetryShiftVec, branching);
						break;
					case DIST_POSITIONS :
						if (solPos == NULL)
							solPos = new GC_ThreadingSpacePosDist::SolutionVec();
						else
							solPos->clear(); // removed found sequences because core changed 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpacePosDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching,
							solPos, 
							(std::max(sequence.size()-distanceMin,(std::string::size_type)0)) );
						break;
					case DIST_ABSMOVES :
						if (solMove == NULL)
							solMove = new GC_ThreadingSpaceMoveDist::SolutionVec(); 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceMoveDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching,
							solMove, 
							(std::max(sequence.size()-1-distanceMin,(std::string::size_type)0)) );
						break;
					case DIST_HCORE :
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceHcoreDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching);
						break;
					default :
						break;
				}

				sTimer.start();
				if (runMode == ENUMERATE) {
					sols = enumerateStructures(s, maxStructures-oldSols, actStat);
				} else {
					sols = countStructures(s, maxStructures-oldSols, useDDS, actStat);
				}
				searchTime += sTimer.stop();
				if (globStat != NULL) {
					globStat->clone		+= actStat->clone;
					globStat->fail		+= actStat->fail;
					globStat->commit	+= actStat->commit;
//					globStat->decomp	+= actStat->decomp;
					globStat->propagate	+= actStat->propagate;
					globStat->memory	+= actStat->memory;
				}
					// increase maximal found contacts
				if (sols > 0) {
					maxContactsFound = std::max( maxContactsFound, 
												(oldContacts - contactsInSeq));
				}

				oldSols += sols;
				delete s;
					// wenn nur beste gesucht werden, dann u.u. energieschranke
					// anpassen
				if (	coreSelection == ALL_BEST
						&& sols > 0) {
					if (!noOutput)
						std::cout	<<"   + " <<(sols)
									<<" structures" <<std::endl;
					coreDB->setActMinHHcontacts(core.getContacts());
				}
				if (symmetryShiftVec != NULL) { delete symmetryShiftVec; symmetryShiftVec = NULL; }
			} catch (Gecode::Exception e) {
				std::cout << "Exception in Gecode: "
				 << e.what() << "." << std::endl
				<< "Stopping..." << std::endl;
			}
		}
		indexedNeighVecs.clear();
		if (symmetryShiftVec != NULL) { delete symmetryShiftVec; symmetryShiftVec = NULL; }
		if (!noOutput)
			std::cout	<<"\n ==> number generated structures = "
						<<(oldSols - oldContactSols)
						<<"\n\n- total number of generated structures = "
						<<std::setw(15) <<oldSols
						<<"\n- maximal contacts reachable           = "
						<<std::setw(15) << maxContactsFound
						<<std::endl;

		if (globStat != NULL) {
				// statistik ausgabe
			std::cout << "- statistics"
			<< "\n\tsolutions:    " << oldSols
			<< "\n\tpropagations: " << globStat->propagate 
			<< "\n\tfailures:     " << globStat->fail
			<< "\n\tclones:       " << globStat->clone 
//			<< "\n\tdecompos.:    " << globStat->decomp 
			<< "\n\tcommits:      " << globStat->commit 
			<< "\n\tpeak memory:  "
				<< static_cast<int>((globStat->memory+1023) / 1024) << " KB"
			<< "\n\tsearchtime:   " << searchTime <<" ms" 
			<< "\n\truntime:      " << timer.stop() <<" ms"
			<< std::endl;
			
				// garbage collection
			delete globStat; globStat = NULL;
			delete actStat;  actStat = NULL;
		}
		
		if (solPos != NULL) { solPos->clear(); delete solPos; solPos = NULL; }
		if (solMove != NULL) { solMove->clear(); delete solMove; solMove = NULL; }

		return oldSols;
	}

		// generates structures and stores them in structs
	unsigned int
	GC_HPThreading::enumerateStructures(std::vector<std::string> & structs) {

		assert(latFrame != NULL);
		assert(coreDB != NULL);
		
		structs.clear();
		
		if (sequence.size() == 0)
			return 0;

		HCore core;
		unsigned int sols = 0, oldSols = 0, oldContacts = INT_MAX,
			oldContactSols = 0;

		GC_ThreadingSymmBreaker::GlobalShiftVec* symmetryShiftVec = NULL;
			// neighborhood vectors for neighboring constraints
		biu::LatticeFrame::index_set indexedNeighVecs;

			// solution storage for enumeration with minimal distance
		GC_ThreadingSpacePosDist::SolutionVec* solPos = NULL, *solMove = NULL;

			// set maximal frame size needed
		latFrame->setFrameSize(	maxPCoreDistance
								+ coreDB->getActCoreSize());

		while (oldSols < maxStructures && coreDB->getNextCore(core)) {
			try {
					// set minimal frame size
					// set minimal frame size
				if (distanceMode == DIST_NONE || distanceMode == DIST_HCORE ) {	
						// not possible if distance of solution is forced 
						// (different positions indices etc.)  
					latFrame->setFrameSize(	maxPCoreDistance
											+ core.getMaxDimension());
				}

					// naechster, wenn nicht genug singlet-positionen vorhanden
				if ( core.getNumOfPSinglets(latFrame)
							< (int)  seqFeatureMap[GC_ThreadingSpace::P_SINGLET].size()
					|| core.getNumOfHSinglets(latFrame)
							< (int)  seqFeatureMap[GC_ThreadingSpace::H_SINGLET].size())
				{
					continue;
				}
					// naechster, wenn cubisches gitter und falsches verhaeltnis
					// zwischen gerade und ungeraden kernpositionen
				if ( latFrame->getDescriptor()->getName() == "cub"
					&& core.getMinEvenOdd() != minEvenOddHs )
				{
					continue;
				}

					// generate global shift vectors for symmetry breaking
				if (useSymmetryBreaking) {
					symmetryShiftVec = GC_ThreadingSymmBreaker::generateGlobalShiftVec(
										&core, latFrame, symmetryShiftVec);
				}
				
				indexedNeighVecs.clear();
				indexedNeighVecs = latFrame->getIndexedNeighborhood();

					// new space
					// new space
				SuperSpace* s = NULL;
				int branching = BR_NONE;
				switch (distanceMode) {
					case DIST_NONE : 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpace(&sequence,
						latFrame, &indexedNeighVecs, &seqFeatureMap, &hullLvl, &core,
						symmetryShiftVec, branching);
						break;
					case DIST_SHAPES : 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceShapes(&sequence,
						latFrame, &indexedNeighVecs, &seqFeatureMap, &hullLvl, &core,
						symmetryShiftVec, branching);
						break;
					case DIST_POSITIONS :
						if (solPos == NULL)
							solPos = new GC_ThreadingSpacePosDist::SolutionVec();
						else
							solPos->clear(); // clear found sols due to new core  
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpacePosDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching,
							solPos, 
							(std::max(sequence.size()-distanceMin,(std::string::size_type)0)) );
						break;
					case DIST_ABSMOVES :
						if (solMove == NULL)
							solMove = new GC_ThreadingSpacePosDist::SolutionVec(); 
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceMoveDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching,
							solMove, 
							(std::max(sequence.size()-1-distanceMin,(std::string::size_type)0)) );
						break;
					case DIST_HCORE :
						branching = (useDDS?BR_DDS:BR_DFS)|(useSymmetryBreaking?BR_SYM:BR_NONE);
						s = new GC_ThreadingSpaceHcoreDist(&sequence,
							latFrame, &indexedNeighVecs, &seqFeatureMap, 
							&hullLvl, &core, symmetryShiftVec, 
							branching );
						break;
					default :
						break;
				}
					// ausgabe der energieniveaus
				if (core.getContacts() < oldContacts) {
						// number of solutions since last energy level change
					oldContactSols = oldSols;
					oldContacts = core.getContacts();
				}
				std::ostringstream str;
				sols = enumerateStructures(s, maxStructures-oldSols, NULL, str);
				std::string retStr = str.str();
				size_t b=0, e=retStr.find('\n');
				do {
					structs.push_back(retStr.substr(b,(e==retStr.npos?retStr.size()-b:e-b)));
					b = e+1;
					e = retStr.find('\n',b);
				} while(e != retStr.npos); 
				
				oldSols += sols;
				delete s;
					// wenn nur beste gesucht werden, dann u.u. energieschranke
					// anpassen
				if (	coreSelection == ALL_BEST && sols > 0) {
					coreDB->setActMinHHcontacts(core.getContacts());
				}
				if (symmetryShiftVec != NULL) { delete symmetryShiftVec; symmetryShiftVec = NULL; }
			} catch (Gecode::Exception e) {
				std::cout << "Exception in Gecode: "
				 << e.what() << "." << std::endl
				<< "Stopping..." << std::endl;
			}
		}
		indexedNeighVecs.clear();
		if (symmetryShiftVec != NULL) { delete symmetryShiftVec; symmetryShiftVec = NULL; }
		if (solPos != NULL) { solPos->clear(); delete solPos; solPos = NULL; }
		if (solMove != NULL) { solMove->clear(); delete solMove; solMove = NULL; }

		return oldSols;
	}

	std::string
	GC_HPThreading
	::getXYZstring( const biu::IntPoint& p ) const {
		const biu::LatticeNeighborhood& cubNeigh = cubicLatDescr.getNeighborhood();
		std::multiset< biu::Move > moves;
		int xDist = p.getX();
		int yDist = p.getY();
		int zDist = p.getZ();
		biu::IntPoint d;
		int direction = (xDist > 0) ? +1 : -1;
		for (int i=(xDist/direction); i>0; i--) {
			d = biu::IntPoint( direction, 0, 0 );
			moves.insert(cubNeigh.getElement( d ).getMove());
		}
		direction = (yDist > 0) ? +1 : -1;
		for (int i=(yDist/direction); i>0; i--) {
			d = biu::IntPoint( 0, direction, 0 );
			moves.insert(cubNeigh.getElement( d ).getMove());
		}
		direction = (zDist > 0) ? +1 : -1;
		for (int i=(zDist/direction); i>0; i--) {
			d = biu::IntPoint( 0, 0, direction );
			moves.insert(cubNeigh.getElement( d ).getMove());
		}

		std::string sortedMoves = "";
		for (std::multiset<biu::Move>::const_iterator m=moves.begin(); m!=moves.end();m++)
			sortedMoves += cubicLatDescr.getAlphabet()->getString(*m);

		return sortedMoves;
	}


	std::string
	GC_HPThreading
	::getEquivClassString( 	const std::string& sequence
							, const biu::MoveSequence & structure ) const
	{
		std::string normalizedString = "";

		  // get all symmetric move sequences
		std::set< biu::MoveSequence > symmStructures
			= latDescr->getAllSymmetricSequences( structure );

		  // identify subsequence of interest (between first and last H-monomer)
		const size_t firstHpos = sequence.find_first_of('H');
		const size_t lastHpos = sequence.find_last_of('H');

		std::vector< biu::IntPoint > equivClass, normEquivClass;
		biu::IntPoint curDist;

		for (std::set< biu::MoveSequence >::const_iterator curStructure = symmStructures.begin();
				curStructure != symmStructures.end(); curStructure++)
		{
			equivClass.clear();
			curDist = biu::IntPoint(0,0,0);

			  // open first H-segment description
			std::string classString = "[";

			  // extract the distance description between all H-monomers
			for (size_t i=(firstHpos+1); i<=lastHpos; i++) {
				  // extend current H-H distance description
				curDist += latDescr->getNeighborhood().getElement(curStructure->at(i-1));
				  // if next H is reached
				if (sequence.at(i)=='H') {
					  // if we are within a segment
					if (sequence.at(i-1) == 'H') {
						  // add monomer
						if (*classString.rbegin()!='[')
							classString +=",";
						classString += getXYZstring(curDist);
//						classString += latDescr->getAlphabet()->getString( curStructure->at(i-1) );
					} else {
						  // close last segment that is still open
						classString += "]";
						  // add segment distance description
						classString += getXYZstring(curDist);
						  // open new segment starting with current H
						classString += "[";
					}
					  // clear current H-H distance encoding for next filling
					curDist = biu::IntPoint(0,0,0);
				}
			}
			  // close last H-segment
			classString += "]";

			  // store smallest string representation
			if (curStructure == symmStructures.begin() || classString.compare( normalizedString ) < 0)
				normalizedString = classString;

			  // debug output
//			std::cout <<latDescr->getString(*curStructure) << " : " <<classString <<"\n";

		}



		return normalizedString;
	}


  } // namespace gecode
} // namespace cpsp

