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


#include "HPoptdeg.hh"
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
// HPoptdeg CLASS IMPLEMENTATIONS
///////////////////////////////////////////////////////////////////

HPoptdeg::HPoptdeg()
:	latDescr(NULL)
, lattice(NULL)
, coreDB(NULL)
, minimalOutput(true)
, seq2deg()
, limitDeg(0)
, tgtDeg(0)
, alphabetHP("HP",1)
, stat(NULL)
, verboseOut(false)
, hpnx(false)
, seqLength(0)
, startSeq()
, TEMP(0.0)
, maxSteps(0)
, timeLimit(0)
, cutoff(0.0)
, threadingHandler(cpsp::gecode::GC_HPThreading())
, steps(0)
{
}

HPoptdeg::~HPoptdeg() 
{
	if (coreDB		!= NULL) { delete coreDB; }
	if (lattice		!= NULL) { delete lattice; }
	if (latDescr	!= NULL) { delete latDescr; }
	if (stat		!= NULL) { delete stat; }
}

CSequence
HPoptdeg::compressSeq(const std::string& str)
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
HPoptdeg::uncompressSeq(const std::vector<unsigned char>& cs,
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

unsigned int 
HPoptdeg::read_dev_urandom() const {
	static const int uintsize = sizeof(unsigned int);
	static const char *random_dev="/dev/urandom";

	std::ifstream in;
	in.open (random_dev, std::ios::binary );
	char buf[uintsize];
	in.read(buf,uintsize);
	unsigned int random_number = *((unsigned int *)&buf);
	in.close();
	return random_number;
}


bool
HPoptdeg::initThreading(int argc, char** argv) 
{

	cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)threadingHandler;

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

		// verbose output if not silent mode choosen
		verboseOut = !opts.argExist("s");

		unsigned int seed;
		if (opts.argExist("seed")) {
			seed=opts.getIntVal("seed");
		} else {
			// seed random number generator
			seed=read_dev_urandom();
		}
		if (verboseOut) {
			std::cout << "seed: "<<seed<<std::endl;
		}
		srandom(seed);

		// handle sequence length and start sequence
		if ( !opts.argExist("seq") && !opts.argExist("seqlen"))  {
			std::cerr <<"Error in arguments: neither start sequence nor a sequence length given\n";
			return false;
		}

		if ( opts.argExist("seq") && opts.argExist("seqlen")
				&& opts.getIntVal("seqlen") != (int)opts.getStrVal("seq").length()
		)  {
			std::cerr <<"Warning: start sequence length differs from given sequence length\n";
			return false;
		}

		if ( opts.argExist("seq") ) {
			startSeq=opts.getStrVal("seq");
			seqLength=startSeq.length();
		} else {
			seqLength=opts.getIntVal("seqlen");
		}

		// set run mode
		threading.setRunMode(cpsp::gecode::GC_HPThreading::COUNT);

		// verbose output
		threading.setStatOutput(false);
		threading.setVerboseOutput(false);
		minimalOutput = opts.getBoolVal("s");
		threading.setNoOutput(true);
		// using symmetryBreaking
		threading.setSymmetryBreaking(! opts.getBoolVal("noSymmBreak"));
		// using DDS for counting only
		threading.setUsingDDS(false); // no DDS so far for enumeration

		threading.setNoRecomputation(opts.getBoolVal("noRecomp"));
		// core selection mode
		threading.setCoreSelection(cpsp::gecode::GC_HPThreading::ALL_BEST);
		// maximal structure number to take into find
		threading.setMaxStructures((unsigned int) limitDeg);


		if (opts.getIntVal("limitDeg") < 0) {
			std::cerr << "Error: limitDeg<0"<<std::endl;
		}
		limitDeg = (unsigned int) opts.getIntVal("limitDeg");

		if (opts.getIntVal("tgtDeg") < 0) {
			std::cerr << "Error: tgtDeg<0"<<std::endl;
		}
		tgtDeg = (unsigned int) opts.getIntVal("tgtDeg");

		TEMP = opts.getDoubleVal("temp");
		cutoff = opts.getDoubleVal("cutoff");

		maxSteps  = opts.getIntVal("maxSteps");
		if (opts.argExist("timeLimit")) {
			timeLimit = opts.getIntVal("timeLimit");
		} else {
			timeLimit = INT_MAX;
		}

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

		// energie limits
		int minE = INT_MIN/2+2;
		int maxE = 0;

		// contact limits for database init
		int minContacts = -maxE + threading.getContactsInSeq();
		int maxContacts = -minE + threading.getContactsInSeq();

		assertbiu(minContacts >= 0, "minimal contacts smaller than ZERO");
		assertbiu(maxContacts >= 0, "maximal contacts smaller than ZERO");

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

		if (opts.argExist("stat")) {
			if (stat != NULL) delete stat;
			stat = new Statistic();
		}		

		return true;
	} 
	std::cerr << "Parameter Error: not parseable\n";
	return false;
}

bool
HPoptdeg::initNextThreading(const std::string & seq, int limitDeg) 
{
	assert(coreDB != NULL);
	assert(latDescr != NULL);

	cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)threadingHandler;

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

	threading.setMaxStructures((unsigned int) limitDeg);

	return  true;
}


void
HPoptdeg::initAllowedArguments(biu::OptionMap & allowedArgs,
		std::string &infoText ) const {
	allowedArgs.push_back(biu::COption(
			"seq", true, biu::COption::STRING,
			"HP-sequence as start"));
	allowedArgs.push_back(biu::COption(
			"seqlen", true, biu::COption::STRING,
			"sequence length"
			, biu::util::Util_String::int2str(27)));
	allowedArgs.push_back(biu::COption(
			"dbPath", true, biu::COption::STRING,
			"file based database root path (or use $CPSP_COREDB)"));
	allowedArgs.push_back(biu::COption(
			"lat", true, biu::COption::STRING,
			"lattice (CUB, FCC)"
			, "CUB"));
	allowedArgs.push_back(biu::COption(
			"noSymmBreak", true, biu::COption::BOOL,
			"do not break symmetries"));
	allowedArgs.push_back(biu::COption(
			"maxSteps", true, biu::COption::INT,
			"maximal number of optimization steps (value > 0)",
			biu::util::Util_String::int2str(10000)));
	allowedArgs.push_back(biu::COption(
			"timeLimit", true, biu::COption::INT,
			"time limit for optimization run (value > 0)"));
	allowedArgs.push_back(biu::COption(
			"seed", true, biu::COption::INT,
			"seed for random number generator"));
	allowedArgs.push_back(biu::COption(
			"limitDeg", true, biu::COption::INT,
			"maximal degeneracy for start sequences (value > 0)",
			biu::util::Util_String::int2str(5000)));
	allowedArgs.push_back(biu::COption(
			"tgtDeg", true, biu::COption::INT,
			"target degeneracy (value > 0)",
			biu::util::Util_String::int2str(1)));
	allowedArgs.push_back(biu::COption(
			"temp", true, biu::COption::DOUBLE,
			"temperature for metropolis local search"
			, "0.4"));
	allowedArgs.push_back(biu::COption(
			"cutoff", true, biu::COption::DOUBLE,
			"cutoff for approximation in metropolis local search"
			, "0.01"));
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
			"help", true, biu::COption::BOOL,
			"program parameters and help"));
	allowedArgs.push_back(biu::COption(	
			"version", true, biu::COption::BOOL, 
			"version information of this program"));

	infoText = std::string("HPoptdeg finds sequences with low degeneracy using local search.");

} // initArguments

double HPoptdeg::logmetropolis(int deg,int next_deg) const {
	return exp(-(log(next_deg/(double)deg))/TEMP);
}

int
HPoptdeg::cpsp_deg(std::string sequence, int limitDeg) {

	if ( initNextThreading(sequence, limitDeg) ) {
		// conversion to gecode implementation reference
		cpsp::gecode::GC_HPThreading & threading = (cpsp::gecode::GC_HPThreading &)threadingHandler;
		int deg = threading.generateStructures();
		return (deg<limitDeg)?deg:-deg;
	}
	else {
		std::cerr << std::endl << "init of cpsp threading failed."<<std::endl;
		return 0;
	}
}

int 
HPoptdeg::get_degeneracy(std::string sequence, unsigned int limitDeg) {
	CSequence cseq = compressSeq(sequence);

	CSeqDegHash::iterator it = seq2deg.find(cseq);
	if (it == seq2deg.end()
			|| (it->second<0 && - it->second < limitDeg)
	) {
		seq2deg[cseq] = cpsp_deg(sequence,limitDeg);
		if (verboseOut) { std::cout <<"N "; std::cout.flush(); }
	} else {
		if (verboseOut) { std::cout <<"C "; std::cout.flush(); }
	}
	return seq2deg[cseq];
}

std::string 
HPoptdeg::point_mutate(const std::string &sequence, int pos) const {
	std::string result=sequence;
	result[pos]=result[pos]=='H'?'P':'H'; // switch position pos
	return result;
}

std::string 
HPoptdeg::mutate(const std::string &sequence,bool hpnx) const {
	int pos=random()%sequence.length();
	if (hpnx) {
		const char *tab;
		switch(sequence[pos]) {
		case 'H': tab="PNX"; break;
		case 'P': tab="HNX"; break;
		case 'N': tab="HPX"; break;
		default: tab="HPN"; break;
		}
		std::string sequence2=sequence;
		sequence2[pos] = tab[random()%3];
		return sequence2;
	} else {
		return point_mutate(sequence,pos);
	}
}

std::string
HPoptdeg::random_hp(int len) const {
	std::string sequence="";
	for (int i=0; i<len; i++)
		sequence += (random() > RAND_MAX/2 )?'H':'P';
	return sequence;
}

// propability H:P:N:X = 6:1:1:1
std::string 
HPoptdeg::random_hpnx(int len) const {
	static std::string tab = "HHHHHHPNX";
	std::string sequence="";
	for (int i=0; i<len; i++)
		sequence += tab[random() % tab.length()];
	return sequence;
}

double 
HPoptdeg::random_double() const {
	return (random()/(RAND_MAX+1.0));
}


void
HPoptdeg::optimizeDegeneracy() {
	CSequence cseq;
	std::string sequence;

	steps=0;

	int start_time=time(0L);

	int deg=-1;

	if (startSeq != "") {
		steps++;
		sequence=startSeq;
		if (verboseOut) { std::cout << "G "<<steps<<" "<<sequence<<" "; std::cout.flush(); }
		deg=get_degeneracy(sequence,limitDeg);
		if (verboseOut) { std::cout <<deg<<std::endl; }
	}

	// predict structures for random sequences to find a low deg sequence
	while ( deg<=0 ) {
		steps++;
		if (hpnx)
			sequence = random_hpnx(seqLength);
		else
			sequence = random_hp(seqLength);

		if (verboseOut) { std::cout << "G "<<steps<<" "<<sequence<<" "; std::cout.flush(); }
		deg=get_degeneracy(sequence,limitDeg);
		if (verboseOut) { std::cout <<deg<<std::endl; }
	}

	if (verboseOut) {
		std::cout <<"0 "<<sequence<<" N "<<deg<<" A "<<time(0L)-start_time<<std::endl;
	}

	// now we have a sequence with deg>0 (and 'deg<limitDeg')

	std::string new_sequence;
	int new_deg;

	std::string best_sequence=sequence;
	int best_deg=deg;

	steps=0;
	// run monte carlo
	while(	(deg > (int)tgtDeg) 
			&& (steps < maxSteps) 
			&& ((time(0L)-start_time < timeLimit)) ) 
	{

		steps++;
		new_sequence = mutate(sequence,hpnx);  

		//int approx_limitDeg = - (int)(log(0.03)*TEMP) + deg;
		int approx_limitDeg = (int)(deg*exp(-TEMP*log(cutoff))); // for logmetropolis

		if (verboseOut) { 
			std::cout <<steps<<" "<<new_sequence
			<<" "
			// << "("<<approx_limitDeg<<")"
			;
			std::cout.flush();
		}
		new_deg=get_degeneracy(new_sequence,approx_limitDeg);
		if (verboseOut) { 
			std::cout<<new_deg<<" ";
		}

		if (new_deg>0) {
			if (new_deg<deg || random_double() < logmetropolis(deg,new_deg) ) {
				if (verboseOut) { std::cout << "A "; }
				//if (new_deg>=deg) std::cout <<" "<<logmetropolis(deg,new_deg);
				deg=new_deg;
				sequence=new_sequence;
				if (deg<best_deg) {
					best_deg=deg;
					best_sequence=sequence;
				}
			} else {
				if (verboseOut) { std::cout << "D "; }
				//std::cout <<logmetropolis(deg,new_deg)
			}
		}
		else 
			if (verboseOut) { std::cout <<"I "; }

		if (verboseOut) { 
			std::cout <<time(0L)-start_time;
			std::cout << std::endl;
		}
	}

	if (verboseOut) {
		std::cout <<std::endl;
		std::cout <<"Best ";
	}
	std::cout <<best_sequence<<" "<<best_deg<<std::endl;

}


int
HPoptdeg::run(int argc, char** argv) {

	seq2deg.clear();		// clear hash of seen sequences

	// timer for statistics
	biu::Timer timer; timer.start();

	if (initThreading(argc,argv)) {	// init was ok

		if (verboseOut) {
			std::cout	<<"\n###################################################"
			<<"\n##  Optimize Degeneracy :                       ###"
			<<"\n###################################################"
			<<"\n" <<std::endl;
		}

		optimizeDegeneracy();

		if (stat != NULL) {
			stat->hashEntries = seq2deg.size();
		}

	} else {
		//std::cerr << "ERROR: init failed."<<std::endl;
		return -1;
	}

	if (stat != NULL) {
		stat->time = timer.stop();

		if (verboseOut) {
			std::cout	<<"\n###################################################"
			<<"\n##  statistics :  #################################"
			<<"\n###################################################"
			<<"\n";
		}
		std::cout
		<<  "    runtime           : " <<stat->time <<" ms"
		<<"\n    steps             : " <<steps
		<<"\n    cache size        : " <<stat->hashEntries<<" entries"
		<<std::endl; 
	}

	if (verboseOut) {		
		std::cout	<<"\n###################################################"
		<<"\n" <<std::endl;
	} 
	return 0;
}


///////////////////////////////////////////////////////////////////
// THE MAIN FUNCTION
///////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	HPoptdeg hpOptdeg;

	return hpOptdeg.run(argc,argv);
}
