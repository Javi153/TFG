/*
 *  Main authors:
 *     Mohamad Rabbath <rabbath@informatik.uni-freiburg.de>
 *
 *  Contributing authors:
 *     Martin Mann <mmann@informatik.uni-freiburg.de>
 *     Sebastian Will <will@informatik.uni-freiburg.de>
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#include "side_chain_optdeg.h"
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
#include <math.h>

#include <ctime>
#include <cmath>
#include <stdio.h>
using namespace std;

//Generate random sequence
std::string SidechainOptdeg::generateSequence(){
	int lower=0,upper=1,j=0;
	std::string str("");
	srand((unsigned)time(NULL));
	for(int i=0;i<option.Length;i++){
		j = lower + rand() % upper;
		if(j>0)
			str+="H";
		else str+="P";
	}
	return str;
	/*std::stringstream ss;
     std::string str;
     ss << option.Length;
     ss >> str;
    std::string* s=new std::string();
   //Randomizing the sting knowing the ranomizer path
   string binary=option.RandomizerPath+" -s -len="+str+" > proteinliketemp.out";

   char buff[255];
   strcpy(buff, binary.c_str());
   int st1=system(buff);
   //Getting the generated string
   ifstream ifs("proteinliketemp.out");
   std::getline(ifs, *s);

   ifs.close();
   int st2=system("rm -f proteinliketemp.out");
   if(st1<0 || st2<0)
	std::cout<<"Error in randomizing, check the permission to write"<<std::endl;
  std::string path=*s;
  delete s;
  return path;*/
}

SidechainOptdeg::SidechainOptdeg() 
 :	latDescr(NULL)
 	, lattice(NULL)
 	, coreDB(NULL)
 	, minimalOutput(true)
 	, seq2deg()
 	, limitDeg(0)
 	, tgtDeg(0)
 	, alphabetHP("HP",1)
 	, stat(NULL)
 	, option()
 	, verboseOut(false)
 	, hpnx(false)
 	, seqLength(0)
 	, startSeq()
 	, TEMPERATURE(0.0)
 	, maxSteps(0)
 	, timeLimit(0)
 	, cutoff(0.0)
 	, steps(0)
{

}

SidechainOptdeg::~SidechainOptdeg() 
{
	if (coreDB		!= NULL) { delete coreDB; }
	if (lattice		!= NULL) { delete lattice; }
	if (latDescr	!= NULL) { delete latDescr; }
	if (stat		!= NULL) { delete stat; }
}

/*CSequence SidechainOptdeg::compressSeq(const std::string& str)
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
	    alphabetHP.getIndex(s.at(i)) * (int) pow((double)base,(double) (i%numOfPos));
    }
    return cs;
}

std::string SidechainOptdeg::uncompressSeq(const std::vector<unsigned char>& cs,
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
}*/

unsigned int 
SidechainOptdeg::read_dev_urandom() const {
	/*static const int uintsize = sizeof(unsigned int);
    static const char *random_dev="/dev/urandom";

    std::ifstream in;
    in.open (random_dev, std::ios::binary );
    char buf[uintsize];
    in.read(buf,uintsize);
    unsigned int random_number = *((unsigned int *)&buf);
    in.close();*/
	srand((unsigned)time(NULL));
	unsigned int random_number=rand();
	return random_number;
}

int
SidechainOptdeg::initThreading(int arg, char** args) 
{
	if(arg>=2){
		std::string firstOpt=args[1];
		if(firstOpt=="--help" || firstOpt=="-help"){

			std::cout<<"-dbPath=.. to set the DB Root (or use ENV variable $CPSP_COREDB)"<<std::endl;
			std::cout<<"\t type=str"<<std::endl;
			std::cout<<"[-lat=..] to set the description of the lattice (CUB for cubic or FCC for faced centered cubic), default is CUB"<<std::endl;
			std::cout<<"\t type=str"<<std::endl;
			std::cout<<"[-limitDeg]=.. the threshold of the degeneracy, default is 10000"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"[-seed=..] seed random generator"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"-seqlen=.. the length of the generated sequences"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"[-start=..] the sequence to start with.  If you do not have any specific starting sequence, please randomizing option."<<std::endl;
			std::cout<<"\t type=str"<<std::endl;;
			std::cout<<"[-cutoff=..] the cuttoff"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"-temp the temprature"<<std::endl;
			std::cout<<"\t type=double"<<std::endl;;
			std::cout<<"[-maxSteps=..] max steps"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"[-timeLimit=..] the limit for the running time in milli seconds"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"-tgtDeg=.. the target degneracy"<<std::endl;
			std::cout<<"\t type=int"<<std::endl;;
			std::cout<<"[-stat] for statistical information (true or false)"<<std::endl;
			std::cout<<"\t type=bool"<<std::endl;;
			std::cout<<"[-s] in case minimum output (true or false)"<<std::endl;
			std::cout<<"\t type=bool"<<std::endl;;
			
			giveVersion();
			return 1;
		}
		else{
			for(int i=1;i<arg;i++){
				std::string optionStr=args[i];
				string::size_type loc = optionStr.find( "=", 0 );
				if( loc != string::npos ) {
					string optionKind=optionStr.substr(0,loc);
					if(loc+1<optionStr.size()){
						string optionVal=optionStr.substr(loc+1);
						if(optionKind=="-dbPath"){
							if (optionVal.size() == 0) {
								cout<< "ERROR : H-code database : parameter dbPath present but no path given"<<std::endl;
								return -1;
							}
							//				   option.sideChainOptions.DBRoot=new string(optionVal);
							option.sideChainOptions.coreDB = new cpsp::HCoreDatabaseFILE(optionVal);
						}
						else if(optionKind=="-lat"){
							option.sideChainOptions.Description=new string(optionVal);
						}

						else if(optionKind=="-threshold"){
							int thr=atoi(optionVal.c_str());
							option.setThreshold(thr);
						}
						else if(optionKind=="-seed"){
							int sNum=atoi(optionVal.c_str());
							option.Seed=sNum;
						}
						else if(optionKind=="-seqlen"){

							int len=atoi(optionVal.c_str());
							option.Length=len;
						}

						else if(optionKind=="-start"){
							option.startSeq = string(optionVal);
						}
						else if(optionKind=="-maxSteps"){
							int msteps=atoi(optionVal.c_str());
							option.maxSteps=msteps;
						}
						else if(optionKind=="-cutoff"){
							double cut=atof(optionVal.c_str());
							option.Cutoff=cut;
						}
						else if(optionKind=="-temp"){
							double tmp=atof(optionVal.c_str());
							option.TEMPERATURE=tmp;
						}
						else if(optionKind=="-stat"){
							if(optionVal=="true")
								option.Stat=true;
							else option.Stat=false;
						}
						else if(optionKind=="-s"){
							if(optionVal=="true")
								option.verboseOut=true;
							else option.verboseOut=false;
						}
						else if(optionKind=="timeLimit"){
							int tlimit=atoi(optionVal.c_str());
							option.timeLimit=tlimit;
						}

					}
					else{

					}
				} 

				else if(optionStr=="-stat"){
					option.Stat=true;
				}

				else if(optionStr=="-s"){
					option.verboseOut=false;
				}
				else {
					cout<< "Bad command, please see --help"<<std::endl;
					return -1;
				}
			}//End for
		}

		// set default for dbPath via env var
		if (option.sideChainOptions.coreDB == NULL) {
			char * tmpStr = NULL;
			if ((tmpStr = getenv("CPSP_COREDB")) != NULL) {

				//		option.sideChainOptions.DBRoot=new string(tmpStr);
				option.sideChainOptions.coreDB = new cpsp::HCoreDatabaseFILE(tmpStr);
			} else {
				cout <<"\n ERROR : no core database given nor CPSP_COREDB environment variable specified!\n" <<std::endl;
				return -1;
			}
		}


	}//end if


	if(option.startSeq==""){
		option.startSeq=generateSequence();
	}

	if((int)option.startSeq.size()!=option.Length){
		std::cerr <<"Warning: start sequence length differs from given sequence length\n";
		option.Length=option.startSeq.size();
	}
	unsigned int seed;
	if (option.Seed!=-1) {
		seed=option.Seed;
	} else {
		// seed random number generator
		seed=read_dev_urandom();
	}

	srandom(seed);

	if (option.timeLimit!=-1) {
		timeLimit = option.timeLimit;
	} else {
		timeLimit = INT_MAX;
	}

	TEMPERATURE = option.TEMPERATURE;
	cutoff = option.Cutoff;
	maxSteps  = option.maxSteps;
	seqLength = option.Length;
	minimalOutput=!option.verboseOut;
	verboseOut = option.verboseOut;
	if(option.startSeq!="")
		startSeq = option.startSeq;
	tgtDeg=option.tgtDeg;
	limitDeg=option.Threshold;
	if (option.Stat == true) {
		if (stat != NULL) delete stat;
		stat = new Statistic();
	}
	return 0;
}

double SidechainOptdeg::logmetropolis(int deg,int next_deg) const {
	return exp(-(log(next_deg/(double)deg))/TEMPERATURE);
}



std::string 
SidechainOptdeg::point_mutate(const std::string &sequence, int pos) const {
	std::string result=sequence;
	result[pos]=result[pos]=='H'?'P':'H'; // switch position pos
	return result;
}

std::string 
SidechainOptdeg::mutate(const std::string &sequence,bool hpnx) const {
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
SidechainOptdeg::random_hp(int len) const {
	std::string sequence="";
	for (int i=0; i<len; i++)
		sequence += (random() > RAND_MAX/2 )?'H':'P';
	return sequence;
}

// propability H:P:N:X = 6:1:1:1
std::string 
SidechainOptdeg::random_hpnx(int len) const {
	static std::string tab = "HHHHHHPNX";
	std::string sequence="";
	for (int i=0; i<len; i++)
		sequence += tab[random() % tab.length()];
	return sequence;
}

double 
SidechainOptdeg::random_double() const {
	return (random()/(RAND_MAX+1.0));
}

bool
SidechainOptdeg::initNextThreading(const std::string & seq, int limitDeg) 
{
	assert(option.sideChainOptions.Sequence!=NULL);
	*(option.sideChainOptions.Sequence)=seq;
	option.setThreshold(limitDeg);
	return  true;
}

int
SidechainOptdeg::cpsp_deg(std::string sequence, int limitDeg) {
	if ( initNextThreading(sequence, limitDeg) ) {
		cpsp::scth::SideChainThreadingHandler handler(option.sideChainOptions);
		int deg = handler.getSolutionsNum();
		return (deg<limitDeg)?deg:-deg;
	}
	else {
		std::cerr << std::endl << "init of cpsp threading failed."<<std::endl;
		return 0;
	}
}

int 
SidechainOptdeg::get_degeneracy(std::string sequence, unsigned int limitDeg) {
	const biu::Alphabet::Sequence& s = alphabetHP.getSequence(sequence);
	CSequence cseq = alphabetHP.compress(s);

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

void
SidechainOptdeg::optimizeDegeneracy() {
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
	while(deg>(int)tgtDeg && steps<maxSteps && (time(0L)-start_time < timeLimit) ) {

		steps++;
		new_sequence = mutate(sequence,hpnx);  

		//int approx_limitDeg = - (int)(log(0.03)*TEMP) + deg;
		int approx_limitDeg = (int)(deg*exp(-TEMPERATURE*log(cutoff))); // for logmetropolis

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
SidechainOptdeg::run(int argc, char** argv) {

	seq2deg.clear();		// clear hash of seen sequences

	// timer for statistics
	biu::Timer timer; timer.start();
	int status=initThreading(argc,argv);			
	if (status==0) {	// init was ok

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
		optdeg_helper_view();
		return status;
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

int main(int argc, char**argv) {
	SidechainOptdeg shpOptdeg;

	return shpOptdeg.run(argc,argv);

}
