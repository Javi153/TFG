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


#include "HPnnet.hh"
#include "version.hh"

#include <iomanip>
#include <fstream>
#include <algorithm>
#include <queue>

#include <biu/Timer.hh>
#include <biu/util/Util_String.h>

#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeModel.hh>
#include <cpsp/Exception.hh>
#include <cpsp/HCoreDatabaseFILE.hh>
#include <cpsp/gecode/GC_HPThreading.hh>



///////////////////////////////////////////////////////////////////
// HPconnect CLASS IMPLEMENTATIONS
///////////////////////////////////////////////////////////////////

	HPconnect::HPconnect()
	 : latDescr(NULL), lattice(NULL), coreDB(NULL), minimalOutput(true),
	   alphabetHP("HP",1), stat(NULL), seqSize(0)
	{
		clearSequences();
	}
	
	void
	HPconnect::readData(std::istream & input, cpsp::StringVec & seqVec) 
	{
		std::string tmp="", lineBuffer="";
		while (std::getline(input, tmp, '\n')) {
			if (!tmp.empty())
				seqVec.push_back(biu::util::Util_String::chompStr(tmp));
		}
	}

	HPconnect::~HPconnect() 
	{
		if (coreDB		!= NULL) { delete coreDB; }
		if (lattice		!= NULL) { delete lattice; }
		if (latDescr	!= NULL) { delete latDescr; }
		if (stat		!= NULL) { delete stat; }
	}
	
	void
	HPconnect::printDot( std::ostream & output ) const {
		net.printDOT(output);
	}	

	void
	HPconnect::printSeqs(std::ostream & output, bool withInd) const {
		for (cpsp::StringVec::size_type s = 0; s < seqs.size(); s++) {
			if (withInd)
				output <<std::setw(8) <<s <<"  ";
			output <<seqs[s] <<std::endl;
		}
	}	    
	
	bool
	HPconnect::areNeighbored( const std::string & s1, const std::string & s2) const {
		if (s1.size() == s2.size()) {
			unsigned int differences = 0;
			for (std::string::size_type i = 0; differences < 2 && i<s1.size(); i++)
				if (s1[i]!=s2[i])
					differences++;
			return differences == 1;
		}
		return false;
	}
	
	void
	HPconnect::directConnect( void ) {
		cpsp::StringVec::size_type i,j;
		net = Net(seqs.size());	// clear network and initialize with all node indices
		  // add all edges
		for (i=0; i < seqs.size() ; i++) {
			for (j = i+1; j < seqs.size(); j++) {
				if (areNeighbored(seqs[i], seqs[j])) {
					net.addEdge(i,j);
				}
			}
		}
	}

	int
	HPconnect::connectedComponents( ComponentLabel & label ) {
		return (int)net.connectedComponents(label);
	}
	
	int
	HPconnect::connectedComponents() {
		HPconnect::ComponentLabel label;
		return connectedComponents(label);
	}
	
	CSequence
	HPconnect::compressSeq(const std::string& str)
	 const {
	 	using namespace biu;
	 	const Alphabet::Sequence& s = alphabetHP.getSequence(str);
			// base of the number system used for the compression
		size_t base = alphabetHP.getAlphabetSize();
			// the number of positions that can be saved in a byte
		size_t numOfPos = (size_t) (std::log(256.) / std::log((double) base));
			// ceil (s->size() / numOfPos) bytes are necessary to save the
			// sequence
		std::vector<unsigned char> cs(
		 (size_t) std::ceil((double) s.size() / (double) numOfPos), 0);

			// 1 byte = s[0] * base^0 + ... + s[numOfPos] * base^numOfPos
		for (size_t i = 0; i < s.size(); i++) {
			cs[i/numOfPos] +=
			 alphabetHP.getIndex(s.at(i)) * (int) pow(base, i%numOfPos);
		}
		return cs;
	}

	std::string
	HPconnect::uncompressSeq(const std::vector<unsigned char>& cs,
		size_t seqLength) const 
	{
			// base of the number system used for the compression
		size_t base = alphabetHP.getAlphabetSize();
			// the number of positions that can be saved in a byte
		size_t numOfPos = (size_t) (std::log(256.) / std::log((double) base));
		std::string uncompressedSeq = "";
			// element of cs which is currently processed
		unsigned char c = cs[0];

		assertbiu(seqLength <= cs.size()*numOfPos,
		 "The compressed sequence can't hold <seqLength> many elements.");

			// 1 byte = s[0] * base^0 + ... + s[numOfPos-1] * base^(numOfPos-1)
		for (size_t i = 0; i < seqLength; i++) {
				// update c if new element of cs is processed
			if (i%numOfPos == 0) {
				c = cs[i/numOfPos];
			}
			uncompressedSeq += alphabetHP.getString(alphabetHP.getElement(
								c % base));
			c /= base;
		}

		return uncompressedSeq;
	}
	
	bool
	HPconnect::isNetSequence(const std::string & seq, 
		cpsp::HPThreadingHandler & t, cpsp::StringVec &sols, bool errorOutput) 
	{
			// init CPSP threading object
		initNextThreading( seq, t);
			// conversion to gecode implementation reference
		cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)t;
		
		sols.clear();
		unsigned int deg = threading.enumerateStructures(sols);
			// mark as seen
		CSequence cseq = compressSeq(seq);
		seq2deg[cseq] = deg;
			// checking for errors
		if (deg > maxDeg) {
			if (errorOutput) {
				std::cerr <<"\n\tError: " 
					<<"Degeneracy of the given sequence higher than allowed (deg. > "
					<< maxDeg <<" ).\n\n";
			}
			return false;
		}
		if (deg == 0) {
			if (errorOutput) {
				std::cerr <<"\n\tError: " 
					<<"Degeneracy of the given sequence is zero (no optimal "
					<<"structure found).\n\n";
			}
			return false;
		}
		if (!structFound(sols)) {
			if (errorOutput) {
				std::cerr <<"\n\tError: " 
					<<"Given structure is not optimal for checked sequence.\n\n";
//				for (size_t x = 0; x<sols.size(); x++)  {
//					std::cout <<" - "<<sols[x] <<"\n";
//				}
			}
			return false;
		}
		
		return true;		
	}

	bool
	HPconnect::isNetSequence(const std::string & seq, 
		cpsp::HPThreadingHandler & t, bool errorOutput) 
	{
		cpsp::StringVec sols;
		return isNetSequence( seq, t, sols, errorOutput);
	}
	
	bool
	HPconnect::initThreading(int argc, char** argv, 
		cpsp::HPThreadingHandler & t, SeqQueue & toCheck) 
	{
		
		cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)t;
		
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
				return false;
			}
			if (opts.getBoolVal("version")) {
				giveVersion();
				return false;
			}
			
			if ( !opts.argExist("seq") && !opts.argExist("file"))  {
				std::cerr <<"Error in arguments: neither sequence nor a file given\n";
				return false;
			}
			
				// set run mode
			threading.setRunMode(cpsp::gecode::GC_HPThreading::ENUMERATE);

				// verbose output
			threading.setStatOutput(false);
			minimalOutput = opts.getBoolVal("s");
			threading.setVerboseOutput(false);
			threading.setNoOutput(minimalOutput);
				// using symmetryBreaking
			threading.setSymmetryBreaking(! opts.getBoolVal("noSymmBreak"));
				// using DDS for counting only
			threading.setUsingDDS(false); // no DDS so far for enumeration

			threading.setNoRecomputation(opts.getBoolVal("noRecomp"));
				// core selection mode
			threading.setCoreSelection(cpsp::gecode::GC_HPThreading::ALL_BEST);
			// maximal structure number to take into find
			int tmp = opts.getIntVal("maxDeg");
			if (tmp < 1) {
				std::cerr <<"Error in arguments : maxSol <= 0\n";
				return false;
			}
			threading.setMaxStructures((unsigned int) tmp+1);
			maxDeg = (unsigned int)tmp;

				// sequence
			bool seqInitialized = false;
			std::string initialSeq = "";
			if (opts.argExist("seq")) {
				initialSeq = opts.getStrVal("seq");
				if ( initialSeq.size() == 0 )  {
					std::cerr <<"Error in arguments: given sequence is empty";
					return false;
				}
				seqSize = initialSeq.size();
				for (std::string::iterator si = initialSeq.begin(); si!=initialSeq.end(); si++) {
					if (*si == 'h') 		{ *si = 'H';
					} else if (*si == 'p')	{ *si = 'P';
					} else if (*si != 'H' && *si != 'P') {
						std::cerr <<"Error in arguments: sequence is no HP-sequence : " 
									<<initialSeq <<"\n";
						return false;
					}
				}
				if (doLexOrdering) {
					std::string tmp = initialSeq;
					std::reverse(tmp.begin(),tmp.end());
					if ( tmp < initialSeq ) {
						std::cerr <<"\n\tError: " 
							<<"Please give reverted sequence and corresponding structure.\n\n";
						return false;
					}
				}
				seqInitialized = true;
			}
				// read file content
			cpsp::StringVec initSeqVec;
			if (opts.argExist("file")){
				std::ifstream input;
				input.open(opts.getStrVal("file").c_str());
				if (input.fail()) {
					std::cerr <<"Error in arguments: could not open file '"
					<<opts.getStrVal("file") <<"'\n";
					return false;
				}
				readData(input, initSeqVec);
				input.close();
				if (!seqInitialized) {
					if (initSeqVec.size() == 0 ) {
						std::cerr <<"Error in arguments: no sequence given and"
							<<" sequence file '" <<opts.getStrVal("file") 
							<<"' contains no sequence.\n";
						return false;
					}
					seqSize = initSeqVec.begin()->size();
				}
			}

				// verbose output if not silent mode choosen
			verboseOut = !opts.argExist("s");

				// lattice modell
			if (lattice != NULL) {
				delete lattice;
				lattice = NULL;
			}
			if (latDescr != NULL) { 
				delete latDescr; 
				latDescr = NULL;
			}
				
			if (opts.getStrVal("lat") == "CUB")
				latDescr = new biu::LatticeDescriptorCUB();
			else if (opts.getStrVal("lat") == "FCC")
				latDescr = new biu::LatticeDescriptorFCC();
			else {
				std::cerr <<"Error in arguments: lattice is not supported\n";
				return false;
			}
			threading.setLatticeDescriptor(latDescr);
			lattice = new biu::LatticeModel(latDescr);
			
				// structure
			structGiven = opts.argExist("struct");
			structsRel.clear();
			if (structGiven) {
				std::string str = opts.getStrVal("struct");
				if ( str.size() == 0 ) {
					std::cerr <<"Error in arguments: no sequence given";
					return false;
				}
				if (str.size() != (seqSize-1)) {
					std::cerr <<"Error in arguments: structure and sequence are of"
							<<" different length\n";
					return false;
				}
				if (!latDescr->getAlphabet()->isAlphabetString(str)) {
					std::cerr <<"Error in arguments: structure is no allowed"
							<<" Absative move string in this lattice\n";
					return false;
				}
				
				biu::IPointVec strPoints(lattice->absMovesToPoints(
										lattice->parseMoveString(
											str )));

				const biu::AutomorphismVec autM = lattice->getDescriptor()->getAutomorphisms();
				biu::AutomorphismVec::const_iterator a;
				biu::IPointVec act = strPoints;
				for (a = autM.begin(); a!= autM.end(); a++) {
					for (biu::IPointVec::size_type pos=0; pos < strPoints.size(); pos++) {
						act[pos] = (*a) * strPoints[pos];
					}
					structsRel.insert(lattice->getString(lattice->pointsToRelMoves(act)));
				}
				
				if(verboseOut) {
					std::cout <<"\n symmetric structures\n";
					for (StringSet::const_iterator it = structsRel.begin(); it != structsRel.end(); it++)
						std::cout <<" rel = " <<(*it) <<"\n";
				}
			}

				// energie limits
			int minE = INT_MIN/2+2;
			int maxE = 0;
// @TODO option wieder reinnehmen
//			minE = opts.getIntVal("minE");
//			maxE = opts.getIntVal("maxE");
				// semantische fehler pruefen
			if (minE > maxE) {
				std::cerr <<"Error in arguments : minE > maxE";
				return false;
			}

				// contact limits for database init
//			unsigned int minContacts = -maxE + threading.getContactsInSeq();
//			unsigned int maxContacts = -minE + threading.getContactsInSeq();

				// H-core database
			char * tmpStr = NULL;
			if (coreDB != NULL) { delete coreDB; coreDB = NULL;}
			if (opts.argExist("dbPath")) {
				coreDB = new cpsp::HCoreDatabaseFILE(opts.getStrVal("dbPath"));
			} else if ((tmpStr = getenv("CPSP_COREDB")) != NULL) {
				coreDB = new cpsp::HCoreDatabaseFILE(std::string(tmpStr));
			} else {
				std::cerr << "H-core database Error: no database available (FILE)\n";
				return false;
			}
				// testen
			if (!coreDB->isConnected()) {
				std::cerr << "H-core database Error: can't open database\n";
				return false;
			}
			threading.setCoreDB(coreDB);
			
			doLexOrdering = opts.argExist("lexOrd");
			
			if (opts.argExist("stat")) {
				if (stat != NULL) delete stat;
				stat = new Statistic();
			}
			
				// check initial sequence if adopting the structure as optimum
			if (initialSeq != "") {
				seqInitialized = isNetSequence( initialSeq, threading, true);
				if (seqInitialized) {
						// add to queue of neutral net HP-sequences to expand
					CSequence cseq = compressSeq(initialSeq);
					toCheck.push(cseq);
				} else { // sequence is no neutral net element
					seqSize = 0;  // reinit sequence length
					return false;
				}
			}
			
				// check sequences from file if adopting the structure as optimum
			if (initSeqVec.size() > 0) {
				for (cpsp::StringVec::const_iterator it = initSeqVec.begin();
					it != initSeqVec.end(); it++) 
				{
					if (it->size() != seqSize) {
						std::cerr  <<"Error in arguments: sequences in file '"
							<<opts.getStrVal("file") <<"' differ in length\n";
						return false;
					}
					if (!isNetSequence( *it, threading, true)) {
						return false;
					}
						// add to queue of neutral net HP-sequences to expand
					CSequence cseq = compressSeq(*it);
					toCheck.push(cseq);
				}					
			}
			
			return true;
		} 
		std::cerr << "Parameter Error: not parseable\n";
		return false;
	}
	
	bool
	HPconnect::initNextThreading(const std::string & seq,
		cpsp::HPThreadingHandler & t) 
	{
		assert(coreDB != NULL);
		assert(latDescr != NULL);
		
		cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)t;
		
			// next sequence
		threading.setSequence(seq);
		
			// reinit core db 
		if (!coreDB->initCoreAccess(	*latDescr,
										threading.getCoreSize(),
										threading.getContactsInSeq(),
										UINT_MAX/2)) {
			std::cerr <<"H-core database Error: access init failed\n";
			return false;
		}
		threading.setCoreDB(coreDB);
		
		return  true;
	}


	void
	HPconnect::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const {
		allowedArgs.push_back(biu::COption(
								"seq", true, biu::COption::STRING,
								"HP-sequence to expand"));
		allowedArgs.push_back(biu::COption(
								"file", true, biu::COption::STRING,
								"file with HP-sequences to expand"));
		allowedArgs.push_back(biu::COption(
								"struct", true, biu::COption::STRING,
								"structure all sequences have to fold in (abs. Moves)"));
		allowedArgs.push_back(biu::COption(
								"dbPath", true, biu::COption::STRING,
								"file based database root path (or use $CPSP_COREDB)"));
		allowedArgs.push_back(biu::COption(
								"lat", true, biu::COption::STRING,
								"lattice (CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(
								"noSymmBreak", true, biu::COption::BOOL,
								"do not break symmetries"));
		allowedArgs.push_back(biu::COption(
								"maxDeg", true, biu::COption::INT,
							"maximal degeneration of the sequences (value > 0)",
								biu::util::Util_String::int2str(1)));
		allowedArgs.push_back(biu::COption(
								"noRecomp", true, biu::COption::BOOL,
								"no space recomputation"));
		allowedArgs.push_back(biu::COption(
								"s", true, biu::COption::BOOL,
								"silent - minimal output"));
		allowedArgs.push_back(biu::COption(
								"stat", true, biu::COption::BOOL,
								"print run statistics"));
		allowedArgs.push_back(biu::COption(
								"lexOrd", true, biu::COption::BOOL,
								"use always lexicographically smaller of sequence and reversion"));
		allowedArgs.push_back(biu::COption(
								"help", true, biu::COption::BOOL,
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPnnet calculates the neutral net for a given ")
				+ std::string("sequence and structure. All sequences have to have ")
				+ std::string("a degeneracy below a given level, have to form ")
				+ std::string("the given structure as an optimal, and to be a ")
				+ std::string("direct or indirect neighbor of the initial sequence ")
				+ std::string("by mutating sequence positions.");

	} // initArguments
	
	
	bool
	HPconnect::structFound(const cpsp::StringVec & sols) const {
		assert(lattice != NULL);
		if (!structGiven)
			return true;
		bool found = false;
		for (cpsp::StringVec::const_iterator i = sols.begin(); !found && i!=sols.end();i++) {
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
	
	cpsp::StringVec
	HPconnect::processQueue(cpsp::HPThreadingHandler & t, SeqQueue & q) {

		cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)t;

		std::string actSeq = "", tmp="";
		std::string * next = &actSeq;
		std::string::size_type i = 0;
		CSequence cseq;
		cpsp::StringVec sols, otherNets;
		while(!q.empty()) {
			actSeq = uncompressSeq(q.front(),seqSize);
			std::cout <<" "<<std::setw(4) <<seqs.size() <<"  "<<actSeq <<" " <<seq2deg[q.front()] <<std::endl;
			seqs.push_back(actSeq);
			q.pop();
			for (i = 0; i< actSeq.size(); i++) {
					// flip position i
				if (actSeq[i] == 'H') actSeq[i] = 'P';
				else actSeq[i] = 'H';
				
				if (doLexOrdering) {
					tmp = actSeq;
						// tmp == reverted sequence of actSeq
					std::reverse_copy(actSeq.begin(),actSeq.end(),tmp.begin());
	
						// reference to lexicographically smaller sequence				
					next = (tmp < actSeq) ? &tmp : &actSeq;
				} else {
					next = &actSeq;
				}

				cseq = compressSeq(*next);
					// if not already processed
				if (seq2deg.find(cseq)==seq2deg.end()) {
					
//					std::cout <<"\n new " <<(*next) <<" with deg = " <<deg <<"";
					 // if degeneracy low enough and able to build searched structure
					if ( isNetSequence(*next, threading) ) {
						if (verboseOut) {
							std::cout <<std::setw(5) <<" " <<"  --> " <<(*next) <<" added"<<std::endl;
						}
						q.push(cseq);
					} else if (	seq2deg[cseq] <= maxDeg && seq2deg[cseq] > 0) {
						otherNets.push_back(*next);
					}

//				} else {
//					std::cout <<"\n    " <<(*next) <<" already seen";
				}
				
								
					// flip position i back
				if (actSeq[i] == 'H') actSeq[i] = 'P';
				else actSeq[i] = 'H';
			}
			
		}
		return otherNets;
	}

	int
	HPconnect::run(int argc, char** argv) {
		seq2deg.clear();		// clear hash of seen sequences
		cpsp::StringVec sols;	// holds generated structures

		cpsp::gecode::GC_HPThreading threading;	// object to determine degeneracy

			// timer for statistics
		biu::Timer timer; timer.start();
		
		SeqQueue toCheck;	// queue of sequences of the neutral net (to expand)
		if (initThreading(argc,argv,threading, toCheck)) {	// init was ok

			seqs.clear();
			
				// process queue until empty
			std::cout	<<"\n###################################################"
						<<"\n##  the neutral net (ID,SEQUENCE,DEGENERACY) :  ###"
						<<"\n###################################################"
						<<"\n" <<std::endl;
			cpsp::StringVec otherNets = processQueue(threading, toCheck);
			
			if (stat != NULL) {
				stat->hashEntries = seq2deg.size();
			}
			
				// print network
			std::cout	<<"\n###################################################"
						<<"\n##  the neutral net in DOT format :  ##############"
						<<"\n###################################################"
						<<"\n" <<std::endl;
			directConnect();
			printDot(std::cout);
			
				// give remaining nets with different structure
			if (otherNets.size() > 0) {
				std::cout	<<"\n###################################################"
							<<"\n##  sequences of other neutral nets :  ############"
							<<"\n###################################################"
							<<"\n" <<std::endl;
				for (cpsp::StringVec::const_iterator it = otherNets.begin(); it != otherNets.end(); it++) {
					std::cout <<std::setw(7) <<" "<<(*it)<<" " <<seq2deg[compressSeq(*it)]<<"\n";
				}
				std::cout.flush();
			}
			
		} else {
			return -1;
		}
		
		if (stat != NULL) {
			stat->time = timer.stop();
			stat->conComponents = connectedComponents();
			stat->netSize = seqs.size();

			std::cout	<<"\n###################################################"
						<<"\n##  statistics :  #################################"
						<<"\n###################################################"
						<<"\n"
						<<"\n    runtime           : " <<stat->time <<" ms"
						<<"\n    checked sequences : " <<stat->hashEntries
						<<"\n    sequences of net  : " <<stat->netSize
						<<"\n    independent nets  : " <<stat->conComponents
						<<"\n" <<std::endl; 
		}
		
		std::cout	<<"\n###################################################"
					<<"\n" <<std::endl; 
		return 0;
	}


///////////////////////////////////////////////////////////////////
// THE MAIN FUNCTION
///////////////////////////////////////////////////////////////////

	int main(int argc, char** argv) {
		HPconnect hpConnect;
		
/*		std::ifstream input;
		input.open("test.seq");
		if (input.good()) {
			hpConnect.readData(input);
			input.close();
		} else {
			std::cerr <<"could not open file 'test.seq'\n";
		}
		
		hpConnect.directConnect();
			std::cout <<std::endl;
		hpConnect.printSeqs(std::cout);
			std::cout <<std::endl;
		hpConnect.printDot(std::cout);
			std::cout <<std::endl;

		HPconnect::ComponentLabel label;
		std::cout	<<"connected components = " 
					<<hpConnect.connectedComponents(label)
					<<std::endl;
		for (HPconnect::ComponentLabel::size_type i = 0; i< label.size(); i++) {
			std::cout <<std::setw(8) <<i <<" : " <<std::setw(4) <<label[i] <<std::endl;
		}
		std::cout <<std::endl;
*/

		return hpConnect.run(argc,argv);
	}
