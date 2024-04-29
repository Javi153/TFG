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

#ifndef HPDESIGN_H_
#define HPDESIGN_H_


#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <cpsp/Exception.hh>
#include <cpsp/HCoreDatabase.hh>
#include <cpsp/HPThreadingHandler.hh>
#include <cpsp/HPThreadingOptions.hh>
#include <string>

namespace cpsp
{
	
	typedef std::set<std::string> StringSet;
	
		//! a lattice point that belongs to a structure in a given sequential 
		//! position
	class StructPos : public biu::IntPoint {
	protected:
			//! the position of point in the structure
		std::string::size_type seqPos;
	public:
			//! construction
		StructPos(const biu::IntPoint& p, std::string::size_type pos) :
			biu::IntPoint(p), seqPos(pos)
		{}
			//! copy construction
		StructPos(const StructPos& p) :
			biu::IntPoint(p), seqPos(p.seqPos)
		{}
			//! copy construction
		StructPos(const biu::IntPoint& p) :
			biu::IntPoint(p), seqPos(0)
		{}
			//! destruction
		virtual ~StructPos() 
		{}
			//! returns the sequential position in the structure of this point
		std::string::size_type 
		getSeqPos(void) const
		{ return seqPos; }
			//! sets the sequential position in the structure of this point
		void
		setSeqPos(std::string::size_type pos)
		{ seqPos = pos; }
	};
	
		/**
		 * The options for HPdesign are parsed, validated and the basic
		 * objects for the further process are build.
		 */
	class HPdesignOptions {
	protected:

		biu::COptionParser* opts;
		 
			//! initialisation of the allowed program arguments and the
			//! corresponding additinal informations
		void 
		initAllowedArguments(	biu::OptionMap & allowedArgs,
								std::string &infoText ) const;
		
			//! the lattice descriptor
		biu::LatticeDescriptor* latDescr;
			//! the lattice model
		biu::LatticeModel * lattice;
			//! the structure to find a sequence for
		biu::IPointVec * structure;
			//! the H-core database
		HCoreDatabase * coreDB;
			//! suboptimality levels taken into account
		unsigned int subOptLvl; 
			//! maximal degeneracy a designed structure might have
		unsigned int maxDegeneracy;
			//! whether to do a verbose output or not
		bool verboseOut;
			//! the minimal number of Hs in the sequence
		unsigned int minHs;
			//! the maximal number of Hs in the sequence
		unsigned int maxHs;
			//! whether or not to find all sequences that fullfill the params
		bool all;
			//! the maximal number of sequences to find
		unsigned int maxSeqs; 
			//! whether or not to do only minimal output
		bool minimalOut;
			//! the structure all have to fold in coded in relative moves
		StringSet structsRel;
			//! whether or not statistics are necessary
		bool statistics;
		

	public:
			//! construction
		HPdesignOptions(int argc, char**argv) throw(cpsp::Exception);
			//! destruction
		virtual ~HPdesignOptions();
		
		const biu::LatticeModel*
		getLattice(void) const { return lattice; }
		
		const biu::IPointVec*
		getStructure(void) const { return structure; }
		
		HCoreDatabase*
		getCoreDB(void) { return coreDB; }
		
		unsigned int 
		getSubOptLvl(void) const { return subOptLvl; }
		
		unsigned int
		getMaxDegeneracy(void) const { return maxDegeneracy; }
		
		bool
		verboseOutput(void) const { return verboseOut; }
		
		unsigned int
		getMinHs(void) const { return minHs; }	
		
		unsigned int
		getMaxHs(void) const { return maxHs; }	
		
		bool
		findAll(void) const { return all; }
		
		unsigned int
		getMaxSeqs(void) const { return maxSeqs; }
		
		bool
		silent(void) const { return minimalOut; }
		
		const StringSet &
		getRelStructures(void) const { return structsRel; }
		
			//! returns wether or not structs includes the global structure to 
			//! fullfill
		bool
		structFound(const cpsp::StringVec & structs) const;
		
		bool
		doStatistic(void) const { return statistics; }
	};
	
	
	
	class HPdesign
	{
	protected:
	
		class Statistic {
		public:
			double time;
			unsigned int generatedSeqs;
			unsigned int coresTested;
			Statistic() : time(0.0), generatedSeqs(0), coresTested(0)
			{}
		};
	
			//! the given options
		HPdesignOptions *opts;
		
		typedef std::set< StructPos > PosSet;
		
			//! the sets of points for the structure in all automorphisms for
			//! the lattice
		std::vector< PosSet > structPoints;
		
		Statistic * stat;
		
		
			//! generates all automorphisms of the structures
		void initStructPoints(void);
		bool
		isSubset(const PosSet &s, const biu::IPointSet &sub, biu::IntPoint & shift);
		
			//! checks whether or not the sequence fullfills the requirements
		bool 
		checkSeq(const std::string& seq);	
	public:
		HPdesign(HPdesignOptions *opts);
		virtual ~HPdesign();
		
		virtual 
		void 
		generateSeqs(void);
		
			//! Checks if there are sequences designable with the given 
			//! initialised core database. The contacts yields the number of 
			//! contacts of the cores available and maximally checked.
			//! @param coreDB the initilised core database 
			//! @param contacts the number of contacts of the available cores
			//! @param maxSeqs the maximal number of sequences to design
			//! @return the number of sequences designed 
		virtual
		unsigned int
		checkCores(HCoreDatabase * coreDB, unsigned int & contacts,
		unsigned int maxSeqs);
	};

}

#endif /*HPDESIGN_H_*/
