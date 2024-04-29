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

#ifndef HPCONNECT_HH_
#define HPCONNECT_HH_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <queue>

#include <biu/OptionParser.hh>

#include <biu/LatticeDescriptor.hh>

#include <cpsp/Exception.hh>
#include <cpsp/HCoreDatabase.hh>
#include <cpsp/HPThreadingHandler.hh>

#include <biu/Graph_UD.hh>
//#include <boost/graph/adjacency_list.hpp>

		//! compressed HP-sequence representation 
	typedef std::vector<unsigned char> CSequence;
	typedef std::set<std::string> StringSet;
	
	typedef std::queue<CSequence> SeqQueue;
	typedef std::map<CSequence,unsigned int> CSeqDegHash;

	
	class HPconnect {
	public:
	
		typedef biu::Graph_UD Net;
		typedef biu::Graph_UD::CompLabel ComponentLabel;
		
		
	protected:
	
			//! Manages statistical data of HPconnect runs.
		class Statistic {
		public:
				//! the runtime
			double time;
				//! number sequences hashed during the run
			CSeqDegHash::size_type hashEntries;
				//! number of sequences in the neutral net
			cpsp::StringVec::size_type netSize;
				//! number of connected components of the net
			int conComponents;
				//! construction
			Statistic() : time(0.0), hashEntries(0), conComponents(0) 
			{}
		};
	
			//! the lattice descriptor of the lattice to predict in
		biu::LatticeDescriptor * latDescr;
			//! the lattice model to predict in
		biu::LatticeModel * lattice;
			//! the H-core database the threading will work with 						
		cpsp::HCoreDatabase * coreDB;
			//! minimal output
		bool minimalOutput;
		
			//! the structure all seqs have to adopt (coded in relative moves)
		StringSet structsRel;
		
			//! true if a structure was given
		bool structGiven;

			//! the HP-sequences to connect
		cpsp::StringVec seqs;
		
			//! the connection network
		Net net;
		
			//! maps compressed sequences to their degeneracy 
		CSeqDegHash seq2deg;
		
		unsigned int maxDeg;
		
		biu::Alphabet alphabetHP;
		
		bool doLexOrdering;
		
		Statistic * stat;
		
		bool verboseOut;
		
		std::string::size_type seqSize;
		
		
		
			//! test whether or not two sequences are neighbored in the net
		bool
		areNeighbored( const std::string & s1, const std::string & s2) const;
		
		
		void
		initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const;

			//! initialises the threading object for the next run with the 
			//! given sequence (i.e. reinit of H-core database access)
		bool
		initNextThreading(const std::string & seq, 
			cpsp::HPThreadingHandler & threading);
												
		CSequence
		compressSeq(const std::string& str) const;									

		std::string
		uncompressSeq(const CSequence& cs, size_t seqLength) const;
			
			//! returns wether or not sols includes the global structure to 
			//! fullfill
		bool
		structFound(const cpsp::StringVec & sols) const;

		cpsp::StringVec
		processQueue(cpsp::HPThreadingHandler & t, SeqQueue & q);
		
			//! checks whether or not a sequence is part of the neutral net
			//! using the CPSP approach
		bool
		isNetSequence(const std::string & seq, 
			cpsp::HPThreadingHandler & t, cpsp::StringVec &sols, 
			bool errorOutput = false);

			//! checks whether or not a sequence is part of the neutral net
			//! using the CPSP approach
		bool
		isNetSequence(const std::string & sequence, cpsp::HPThreadingHandler & t,
			bool errorOutput = false);
		
	public:
		
		bool
		initThreading(int argc, char** argv, 
			cpsp::HPThreadingHandler & threading, SeqQueue & toCheck);
			//! construction
		HPconnect();
			//! destruction
		virtual 
		~HPconnect();
		
			//! reads all available HP-sequences from input and appends them
			//! to seqVec
		void
		readData(std::istream & input, cpsp::StringVec & seqVec);
		
			//! connects all sequences in the net that differ only in one 
			//! position
		void
		directConnect( void );
		
			//! prints the current network to the given out-stream
	    void 
	    printDot(std::ostream & output) const;

			//! prints the current sequences and the associated indices in the 
			//! graph if withInd == true
		void
		printSeqs(std::ostream & output, bool withInd = true) const;
		
			//! calculates the connected components of the current network,
			//! the component assignment will be given in label
			//! @return number of connected components
		int
		connectedComponents( ComponentLabel & label );
		
			//! calculates the connected components of the current network
			//! @return number of connected components
		int
		connectedComponents();
		
		int
		run(int argc, char** argv);
			    
	    	
	public:
	
			//! removes all sequences
		void
		clearSequences( void ) { seqs.clear(); net = Net(); }
			//! setter for seqs
		void 
		setSequences( const cpsp::StringVec& val ) { clearSequences(); seqs = val; }
			//! getter for seqs
		const cpsp::StringVec&
		getSequences( void ) const { return seqs; }
			//! adds the sequences
		void
		addSequences( const cpsp::StringVec& toAdd ) { 
			seqs.insert(seqs.end(), toAdd.begin(), toAdd.end()); }		
			//! adds a sequence
		void
		addSequences( const std::string& toAdd ) { seqs.push_back(toAdd); }		
		
		
	};

#endif /*HPCONNECT_HH_*/
