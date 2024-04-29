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

#ifndef HPSTRUCTCG_HH_
#define HPSTRUCTCG_HH_

#include <limits.h>

#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>

	/** 
	 * SymmFreeStructStore is a storage container for move string 
	 * representations of equal length. Only the normalized sequences are stored 
	 * that ensures the exclusion of symmetric structures in the stored set.
	 * 
	 */
	class SymmFreeStructStore {
	protected :

			//! the element type for the internal storage
		typedef std::pair<biu::Alphabet::CSequence,int> StructEpair;
		
	 		/*! Compares the energy values and returns
			 *  true if the energy of the first CStateEnergyPair is greater
			 *  than the energy of the second. This ensures, that the 
			 *  elements in a sorted container are in energy decreasing order.
			 *  Compares the compressed structure move sequences lexicographically
			 *  if they have the same energy value.
			 */ 
		struct	StructEpairCompare {
			bool operator()(const StructEpair& a, const StructEpair& b) const 
			{
				if (a.second != b.second) 
					return a.second > b.second;
				assertbiu(a.first.size()==b.first.size(), "Elements to compare differ in sequence length.");
				for (size_t i=0; i< a.first.size(); i++) {
					if (a.first[i] < b.first[i])
						return true;
					if (a.first[i] > b.first[i])
						return false;
				}
				return false;
			}
		};

			//! the internal storage type		 	
		typedef std::set< StructEpair, StructEpairCompare > StoreClass;
		
			//! the internal data store
		StoreClass store;
		
			//! the lattice descriptor to normalize structures
		biu::LatticeDescriptor const * const latDescr;
	 	
	 		//! length of the sequences in this storage
	 	const size_t seqLengthToStore;
	 	
	 	
	public :
	
			//! a constant iterator to access the elements of the neighborhood
		typedef StoreClass::const_iterator const_iterator;
	
			//! a constant iterator to access the elements of the neighborhood
		typedef StoreClass::size_type size_type;
	
	
	public :
		
			//! construction
		SymmFreeStructStore(const biu::LatticeDescriptor* const latDescr,
			const size_t seqLengthToStore);
		
			//! destruction
		virtual ~SymmFreeStructStore();
		
			/** Adds the given move sequence normalized to the store.
			 */
		virtual void addStructure(biu::MoveSequence & moves, int energy);
		
			/** Access to first element stored.
			 */
		virtual const_iterator begin() const;
		
			/** Iterator to the end of the storage.
			 */
		virtual const_iterator end() const;
		
			/** Returns the number of elements in the storage.
			 */
		virtual size_type size() const;
		
			/** Removes all elements from the storage.
			 */
		virtual void clear();
		
			/** Removes the element given by the iterator and returns the 
			 * element after in the storage.
			 */
		virtual const_iterator removeStructure(const const_iterator & toDel);
		
			/** Returns the number of elements with the given energy-
			 */
		virtual size_t getElementNumber(const int energy);
		
			/** Returns the decompressed sequence the iterator points to.
			 */
		virtual biu::MoveSequence getStructure(const const_iterator & it) const;
		
			/** Returns the energy of the structure the iterator points to.
			 */
		virtual int getEnergy(const const_iterator & it) const;
		
			/** Removes all elements with energy above maxE.
			 * @param maxE the maximal energy allowed for stored elements
			 */
		virtual void ensureMaxE( const int maxE );
		
	};

	/**
	 * HPstructCG is the implementation of a chain growth folding algorithm 
	 * for the structure prediction in the HP-model. Starting from the tailing
	 * 2 monomer sequence, all structures within a given energy threshold to the
	 * minimum reachable so far are kept. These are extended by appending 
	 * successively the remaining sequence monomers (from tail to front) to the 
	 * structures, discarding structures that fall below the threshold.
	 */
	 class HPstructCG {
	 	
	 protected:
	 
	 		//! the HP sequence to fold
	 	std::string seq;
	 	
	 		//! the lattice to work on
	 	biu::LatticeModel* lattice;
	 	
			//! holds the current structures											
		SymmFreeStructStore * structStore;
		
			//! the energy band to store
	 	int deltaE;
	 	
	 protected:
	 
			//! Inits the allowed parameters for HPstructCG.
		void initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const;
											
			//! Generates all structures based on structStore that are extended
			//! with the given monomer and have an energy within the energy band
			//! @param nextMonomerPos the index of the next sequence monomer to 
			//!                       add
			//! @return the minimal energy reached by appending or INT_MAX if
			//!         no appending was possible 
		int appendMonomer( const size_t nextMonomerPos );
		
			//! Prints all structures in the current structStore content to 
			//! stream that have an energy less or equal to the given energy
			//! @param out the stream to write to
			//! @param maxE the maximal energy a printed structure should have
			//! @param max2print the maximal number of structures to print 
		void printStorage(std::ostream & out, const int maxE = INT_MAX,
				int max2print = INT_MAX) const;
					
	 public:
	 		//! construction
	 	HPstructCG();
	 	
	 		//! destruction
	 	~HPstructCG();
	 	
	 		//! runs the chain growth algorithm according to the given program 
	 		//! parameters
	 		//! @return the execution status of parameter parsing and conversion
	 	int run(int argc, char** argv);
	 	
	 }; // HPstructCG

#endif /*HPSTRUCTCG_HH_*/
