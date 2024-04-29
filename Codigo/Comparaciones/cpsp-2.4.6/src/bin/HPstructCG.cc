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


#include "HPstructCG.hh"
#include "version.hh"

#include <sstream>
#include <fstream>
#include <map>
#include <iomanip>
#include <algorithm>

#include <biu/util/Util_String.h>
#include <cpsp/Exception.hh>

#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>




///########################################################################
//   SymmFreeStructStore class implementations
///########################################################################


		//! construction
	SymmFreeStructStore::SymmFreeStructStore(const biu::LatticeDescriptor* const latDescr_,
			const size_t seqLengthToStore_)
	 : store(), latDescr(latDescr_), seqLengthToStore(seqLengthToStore_)
	{
		assertbiu(latDescr != NULL, "SymmFreeStructStore construction : lattice descriptor not given!");
	}
	
		//! destruction
	SymmFreeStructStore::~SymmFreeStructStore()
	{}
	
		/** Adds the given move sequence normalized to the store.
		 */
	void 
	SymmFreeStructStore::addStructure(biu::MoveSequence & moves, int energy)
	{
		assertbiu(moves.size() == seqLengthToStore, "Sequence length differs the expected length!");
		store.insert(StructEpair(latDescr->getAlphabet()->compress(latDescr->normalizeSequence(moves)), energy));
	}
	
		/** Access to first element stored.
		 */
	SymmFreeStructStore::const_iterator 
	SymmFreeStructStore::begin() const
	{ 
		return store.begin();
	}
	
		/** Iterator to the end of the storage.
		 */
	SymmFreeStructStore::const_iterator
	SymmFreeStructStore::end() const
	{
		return store.end();
	}
	
	SymmFreeStructStore::size_type 
	SymmFreeStructStore::size() const
	{
		return store.size();
	}
	
	void 
	SymmFreeStructStore::clear()
	{
		store.clear();
	}
	
	SymmFreeStructStore::const_iterator 
	SymmFreeStructStore::removeStructure(const const_iterator & toDel)
	{
			// get iterator to successor
		const_iterator nextIt = toDel;
		nextIt++;
			// delete	
		store.erase(toDel);
			// return
		return nextIt;	
	}
	
	size_t 
	SymmFreeStructStore::getElementNumber(const int energy)
	{
		size_t number = 0;
		for (const_iterator it = begin(); it!=end(); it++) {
			if (it->second == energy)
				number++;
		}
		return number;
	}
	
	biu::MoveSequence 
	SymmFreeStructStore::getStructure(const SymmFreeStructStore::const_iterator & it) const
	{
		assertbiu( it != end() , "cant access end iterator");
		return (biu::MoveSequence)latDescr->getAlphabet()->decompress(it->first, seqLengthToStore);
	}
	
	int 
	SymmFreeStructStore::getEnergy(const SymmFreeStructStore::const_iterator & it) const
	{
		assertbiu( it != end() , "cant access end iterator");
		return it->second;
	}
	
	void
	SymmFreeStructStore::ensureMaxE(const int maxE) 
	{
			// find first with E <= maxE
		StoreClass::iterator it = store.begin();
		while(it != store.end() && it->second > maxE) {
			it++;
		}
			// remove all leading positions in front of it (only that have E > maxE)
		store.erase(store.begin(), it);
	}


///########################################################################
//   HPstructCG class implementations
///########################################################################





	HPstructCG::HPstructCG()
	 :	seq(""), lattice(NULL), structStore(NULL)
	{
	}
	
	HPstructCG::~HPstructCG()
	{
			// clean up data structures
		if (structStore != NULL) { delete structStore; structStore = NULL;}
		if (lattice != NULL) { delete lattice; lattice = NULL;}
	}
	
	int 
	HPstructCG::appendMonomer( const size_t nextPos )
	{
		assertbiu(nextPos >= 0 && nextPos < seq.size(), "index to add is not within sequence length");
		
		SymmFreeStructStore* newStore = new SymmFreeStructStore(lattice->getDescriptor(), nextPos);
		
		const biu::LatticeNeighborhood& neighs = lattice->getNeighborhood();
		biu::LatticeNeighborhood::const_iterator actNeigh = neighs.begin();
			// create a set of neighbor vectors
		biu::IPointSet neighSet;
		for (actNeigh = neighs.begin(); actNeigh != neighs.end(); actNeigh++) 
			neighSet.insert(*actNeigh);
		
		int minEseen = INT_MAX; // minimal energy seen so far

			// extend all stored structures
		for (SymmFreeStructStore::const_iterator actMvs = structStore->begin();
				actMvs != structStore->end(); actMvs++) 
		{
				// skip structure if minimal energy not reachable
			if (seq[nextPos] == 'P' 
				&& minEseen != INT_MAX 	// at least one structure seen so far
				&& (structStore->getEnergy(actMvs) > (minEseen + deltaE)) )
				continue;
				// get current structure information
			biu::MoveSequence actStr = structStore->getStructure(actMvs);
			int actE = structStore->getEnergy(actMvs);
				// convert move sequence to 3D points
			biu::IPointVec actPnt = lattice->absMovesToPoints( actStr );
				// generate all self avoiding extensions			 
			for (actNeigh = neighs.begin(); actNeigh != neighs.end(); actNeigh++) {
					// get new possible leading position
				bool isSelfavoiding = true;
				biu::IntPoint tailPoint = actPnt[nextPos-1]+(*actNeigh);
				int energyGain = 0;
				if (seq[nextPos] == 'H') {	// check selfavoidingness and calculate new contacts
					for (biu::IPointVec::size_type i = 0; i < actPnt.size(); i++) {
						if (actPnt[i] == tailPoint) {  // selfavoidingness check
							isSelfavoiding = false;
							break;	// not selfavoiding --> go to next possibility
						}
						if ( (i+1) < nextPos	// not successive 
							&& seq[i] == 'H'	// possible H-H contact
							&& neighSet.find(tailPoint - actPnt[i]) != neighSet.end()) // neighbored
						{
							energyGain--;
						}
					}
				} else {	// check for selfavoidingness only
					if (std::find(actPnt.begin(), actPnt.end(), tailPoint) != actPnt.end()) {
						isSelfavoiding = false; // not selfavoiding --> go to next possibility
					}
				}
					// store new structure if energy is small enough
				if (isSelfavoiding && (minEseen == INT_MAX || (actE + energyGain) <= (minEseen + deltaE)) ) {
						// generate new move sequence by appending new move
					biu::MoveSequence newStr = actStr;
					newStr.push_back(lattice->getDescriptor()->getAlphabet()->getIndex(actNeigh->getMove()));
						// add to storage
					newStore->addStructure( newStr, (actE + energyGain) );
						// get new minimal energy seen
					if (minEseen > (actE + energyGain)) {
						minEseen = (actE + energyGain);
							// remove all structures from new store with too high energy
						newStore->ensureMaxE(minEseen+deltaE);
					}
				}
			}
		}
			// remove all stored elements that are not within the energy band
		newStore->ensureMaxE(minEseen+deltaE);
		
			// discard old store		
		SymmFreeStructStore* tmp = structStore;
		structStore = newStore;
		newStore = NULL;
		tmp->clear();
		delete tmp; tmp = NULL;
	
		return ( structStore->size() == 0 ? INT_MAX : minEseen );	
	}
	
	void
	HPstructCG::printStorage(std::ostream & out, const int maxE,
		const int max2print ) const 
	{
		int printed = 0;
		for (SymmFreeStructStore::const_iterator it = structStore->begin();
				printed<max2print && it != structStore->end(); it++) 
		{
			if (it->second <= maxE) {
				out	<<lattice->getString(structStore->getStructure(it))		// structure (normalized) 
					<<" " <<std::setw(8)  <<it->second 	// energy
					<<"\n";
				printed++;
			}
		}
	}
	
	
	int
	HPstructCG::run(int argc, char** argv)
	{
		using namespace biu;
		using namespace util;
		
			// the lattice descriptor
		LatticeDescriptor* descr = NULL;
		
		try {
			

			// check user parameters 
			biu::OptionMap allowedArgs;
			std::string infoText;
			initAllowedArguments(allowedArgs,infoText);	// init
			
				// parse programm arguments	
			COptionParser opts = COptionParser(	allowedArgs, argc, 
												(char**)argv, infoText);
				// arguments not parseable or not all mandatory arguments given
			if (!opts.noErrors()) {
				throw cpsp::Exception("",-1);
			}
				// help output
			if (opts.argExist("help")) {
				opts.coutUsage();
				throw cpsp::Exception("",0);
			}
			if (opts.argExist("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
				// set lattice
			std::string latStr = opts.getStrVal("lat");
			if (latStr.compare("SQR") == 0)
				descr = new LatticeDescriptorSQR();
			else if (latStr.compare("CUB") == 0)
				descr = new LatticeDescriptorCUB();
			else if (latStr.compare("FCC") == 0)
				descr = new LatticeDescriptorFCC();
			else
				throw cpsp::Exception("Unknown lattice type '"+latStr+"'",-2);
			lattice = new biu::LatticeModel(descr);

				// check sequence input
			seq = opts.getStrVal("seq");
				// to upper case 
			for (std::string::iterator p = seq.begin(); p != seq.end(); p++) {
				if (*p == 'p') *p = 'P';
				else if (*p == 'h') *p = 'H';
			}
			if ( seq.size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no sequence given",-2);
			if ( seq.size() < 3 ) 
				throw cpsp::Exception("Nothing to do for such a short sequence",-2);
			biu::Alphabet hpAlph("HP",1);
			if (!hpAlph.isAlphabetString(seq))
				throw cpsp::Exception("Error in arguments: given sequence "+seq+" is no HP-sequence",-2);
				
				// energy band to store
			deltaE = opts.getIntVal("deltaE");
			if (deltaE < 0)
				throw cpsp::Exception("Error in arguments: given deltaE parameter is smaller than zero",-2);
				
			bool printBestOnly = opts.getBoolVal("bestOnly");
				
				// init structure store
			structStore = new SymmFreeStructStore(descr,1);

				// check for verbose output
			bool verboseOut = opts.argExist("v");
			
				// check for counting mode
			bool countingOnly = opts.argExist("count");
			
				// maximal number of minimal energy structures to print
			int max2print = opts.argExist("max2print") 
								? opts.getIntVal("max2print") 
								: INT_MAX;
			
				// the outstream to write to
			std::ostream& out = std::cout;
			
			if (verboseOut) {
				out	<<"\n==================================="
					<<"\n   HPstructCG - CPSP-tool-library"
					<<"\n==================================="
					<<"\n\n"
					<<"\n HP-sequence      : " <<seq
					<<"\n lattice type     : " <<lattice->getDescriptor()->getName()
					<<"\n delta E          : " <<deltaE
					<<std::endl <<std::endl;
			}
			
				// init chain growth
			biu::MoveSequence initialMove = descr->getSequence("F"); 
			structStore->addStructure( initialMove, 0 );
			std::string::size_type curStructLength = 2;
			int minEfound = INT_MAX;
				// elongate chain
			for (; curStructLength < seq.size(); curStructLength++) {
				if (verboseOut) {
					out <<std::setw(5)<<(curStructLength+1)
						<<". monomer : appending to " <<std::setw(12) <<structStore->size()
						<<" structures with minimal E = "
						<<std::setw(8) <<(minEfound==INT_MAX ? 0 : minEfound) 
						<<std::endl;
				}
					// append next monomer
				minEfound = appendMonomer( curStructLength );
					// check if appending failed
				if (minEfound == INT_MAX) {	// appending failed
					out	<<"\n --> Appending of monomer " <<(curStructLength+1)
						<<" failed! Chain growth stopped.\n" <<std::endl;
					throw cpsp::Exception("",-3);
					
				}
			}
			
			if (verboseOut) {
				out <<"\n\n Folded structures and energies : \n";
			}
				// handle results
			if (!countingOnly) {
				out <<std::endl;
				printStorage(out, (printBestOnly ? minEfound : INT_MAX), max2print );
			}
			out <<"\n number of structures = " <<std::setw(12) <<structStore->size()
				<<"\n minimal energy found = " <<std::setw(12) <<minEfound
				<<"\n min. E structures    = " <<std::setw(12) <<structStore->getElementNumber(minEfound)
				<<std::endl;
			
			if (verboseOut) {
				out	<<"\n====================================\n"
					<<std::endl;
			}
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}
		
		if (lattice != NULL) { delete lattice; lattice = NULL;}
		if (descr != NULL) delete descr;

		return 0;
	}
	

		
	void
	HPstructCG::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const
	{
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"the HP-sequence to fold"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (SQR, CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"deltaE", true, biu::COption::INT, 
								"energy band of structure to keep (>=0)", "0"));
		allowedArgs.push_back(biu::COption(	
								"count", true, biu::COption::BOOL, 
								"only count structures"));
		allowedArgs.push_back(biu::COption(
								"max2print", true, biu::COption::INT,
							"maximal number of optimal structures to print (value > 0)",
								biu::util::Util_String::int2str(INT_MAX-2)));
		allowedArgs.push_back(biu::COption(	
								"bestOnly", true, biu::COption::BOOL, 
								"print only structures with minimal energy"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPstructCG implements a chain growth folding")
			+ std::string(" algorithm. The structures are generated with increasing")
			+ std::string(" length, whereby only structures within a given energy")
			+ std::string(" band to the best found so far are kept.")
			+ std::string("\n The resulting structures are printed as normalized")
			+ std::string(" absolute move strings to STDOUT.")
			+ std::string("\n\nThe moves are encoded using:")
			+ std::string("\n F/B : +-x\n L/R : +-y\n U/D : +-z");
	}
	



	int main(int argc, char** argv) {
		HPstructCG hpcg;
			// run conversion
		return hpcg.run(argc, argv);
	}

	/*

##########################################################################
HPPHPPPHPPPPPHHPPHHP


===================================
   HPstructCG - CPSP-tool-library
===================================


 HP-sequence      : HPPHPPPHPPPPPHHPPHHP
 lattice type     : sqr
 delta E          : 0

    3. monomer appending to            1 structures with minimal E =        0
    4. monomer appending to            2 structures with minimal E =        0
    5. monomer appending to            1 structures with minimal E =       -1
    6. monomer appending to            2 structures with minimal E =       -1
    7. monomer appending to            6 structures with minimal E =       -1
    8. monomer appending to           16 structures with minimal E =       -1
    9. monomer appending to            3 structures with minimal E =       -2
   10. monomer appending to            5 structures with minimal E =       -2
   11. monomer appending to           14 structures with minimal E =       -2
   12. monomer appending to           37 structures with minimal E =       -2
   13. monomer appending to          102 structures with minimal E =       -2
   14. monomer appending to          277 structures with minimal E =       -2
   15. monomer appending to           10 structures with minimal E =       -3
   16. monomer appending to            1 structures with minimal E =       -4

 --> Appending of monomer 16 failed! Chain growth stopped.



 
##########################################################################

	 * */
	
