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

#ifndef THREADINGOPTIONS_HH_
#define THREADINGOPTIONS_HH_

#include "cpsp/Exception.hh"
#include <biu/assertbiu.hh>
#include <biu/OptionParser.hh>
#include <string>
#include <map>
#include <biu/LatticeDescriptor.hh>
#include <biu/LatticeFrame.hh>
#include "cpsp/HCoreDatabase.hh"

namespace cpsp
{

	class HPThreadingOptions
	{
	public:
		
			//! constants that define the core selection mode
		enum CORESELECTS { FIRST, ALL_BEST, ALL };			// steuert kernauswahl
	
	private:
	
			//! initialisation of the allowed program arguments and the
			//! corresponding additinal informations
		void initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const;
			//! the HP-sequence the structures are generated for
		std::string*			sequence;
		
			//! the descriptor for the underlying lattice model
		const biu::LatticeDescriptor*	latDescr;
		
			//! the underlying lattice model
		biu::LatticeFrame*		lattice;
		
			//! the core selection mode
		CORESELECTS				coreSelection;
			
			//! upper and lower energy bound for structure prediction
		int 					minE, maxE;
		
			//! maximal number of structure to generate
		unsigned int			maxSols;
		
			//! The hull level for each position of the sequence.
			//! The core has hull level 0, all non-core but core-neighbored
			//! positions have hull level 1 and so on.
		std::vector<unsigned int>	hullLvl;
		
			//! The maximal distance of a P-element of the sequence 
			//! to an H-element.
		unsigned int			maxPCoreDistance;
		
			//! the H-core database, the cores are coming from
		HCoreDatabase*			coreDB;
		
			//! sequence features
		enum SeqFeatures{ P_SINGLET, H_SINGLET };
		
		typedef std::map< int, std::vector <unsigned int> > FeatureMap;
		
			//! maps the sequence indices to their features
		FeatureMap	seqFeatureMap;
		
			//! tag if counting in mode
		bool		onlyCountingStructures;
		
			//! tag if verbose output is wanted
		bool 		verboseOutput;
		
			//! The minimal number of even vs. odd H positions in sequence.
		unsigned int	minEvenOddHs;
		
			//! The number of sequence internal HH-contacts.
		unsigned int	contactsInSeq;

		//! turn on symmetry breaking
		bool useSymmetryBreaking;
		
		//! turn off dynamic decomposing search
		bool useNoDDS;
		
			//! whether or not to print normalized move string representations
		bool normalizeStructures;
		
	public:
			//! construction
			//! @param argc the argument number given by the "main" function
			//! @param argv the string arguments given by the "main" function
		HPThreadingOptions(int argc, char** argv) 
			throw(cpsp::Exception);
		HPThreadingOptions();
		
		virtual ~HPThreadingOptions();
		
		


		const std::string* const 
		getSequence() const { return sequence; }
		unsigned int setSequence(const std::string &seq) throw(cpsp::Exception);

		const bool
		getUseSymmetryBreaking() const { return useSymmetryBreaking; }
		void setUseSymmetryBreaking(bool val) { useSymmetryBreaking = val; }
		
		const bool
		getUseNoDDS() const { return useNoDDS; }
		void setUseNoDDS(bool val) { useNoDDS = val; }
		
		const biu::LatticeFrame* const 
		getLattice() const { return lattice; }
		void setLattice(biu::LatticeFrame* lat) { lattice = lat; latDescr = lat->getDescriptor();}
		
		const CORESELECTS 
		getCoreSelection() const { return coreSelection; }
		void setCoreSelection(CORESELECTS cs) { coreSelection = cs; }
		
		const int 
		getMinimalEnergy() const { return minE; }
		void setMinimalEnergy(int val) { minE = val; }
		
		const int 
		getMaximalEnergy() const { return maxE; }
		void setMaximalEnergy(int val) { maxE = val;}
		
		const unsigned int 
		getMaximalSolutions() const { return maxSols; }
		void setMaximalSolutions(unsigned int val) { maxSols = val; }
		
		const std::vector<unsigned int> *const 
		getHullLevel() const { 
			return &hullLvl; 
		}
		
		const unsigned int 
		getMaxPCoreDistance() const { 
			return maxPCoreDistance; 
		}
		
		void 
		setLatticeFrameSize(const unsigned int newFrameSize) {
			assertbiu(lattice != NULL,"no lattice present");
			lattice->setFrameSize(newFrameSize);
		}
		
			//! returns the current H-core database
		HCoreDatabase* 
		getHCoreDB(void) { return coreDB; }
		void setHCoreDB(HCoreDatabase* db) {assertbiu(db!=NULL,"no core data base present"); coreDB = db; }
		
			//! returns the number of P-singlet positions in the sequence
		unsigned int 
		getNumOfPSinglets(void) { 
			return seqFeatureMap[P_SINGLET].size(); 
		}

			//! returns the number of H-singlet positions in the sequence
		unsigned int 
		getNumOfHSinglets(void) { 
			return seqFeatureMap[H_SINGLET].size(); 
		}
		
			//! returns the P-singlet corresponding indices 
		std::vector <unsigned int>& 
		getPSinglets(void) { 
			return seqFeatureMap[P_SINGLET];
		}

			//! returns the H-singlet corresponding indices 
		std::vector <unsigned int>& 
		getHSinglets(void) { 
			return seqFeatureMap[H_SINGLET];
		}
		
			//! returns whether or not in counting mode
		bool 
		countingMode(void) { return onlyCountingStructures; }
		bool countingMode(bool newVal) { onlyCountingStructures = newVal; return onlyCountingStructures; }
		
			//! returns whether or not in verbose output mode
		bool 
		verboseOut(void) { return verboseOutput; }
		bool verboseOut(bool newVal) { verboseOutput = newVal; return verboseOutput; }
		
			//! Returns the minimum of even vs. odd H positions of the sequence.
		unsigned int 
		getMinEvenOddHs(void) { return minEvenOddHs; }
		
			//! Returns the number of sequence internal HH-contacts.
		unsigned int
		getContactsInSequence(void) { return contactsInSeq; }
		
		//! setter for normalizeStructures
		void 
		setNormalizeStructures( const bool val ) { normalizeStructures = val; }
			//! getter for normalizeStructures
		const bool
		getNormalizeStructures( void ) const { return normalizeStructures; }
		
	};

} // namespace cpsp

#endif /*THREADINGOPTIONS_HH_*/
