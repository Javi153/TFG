#ifndef CHAIN_THREADING_ALL
#define CHAIN_THREADING_ALL
#include <biu/LatticeFrame.hh>
#include <cpsp/HCore.hh>
#include "chain_threading.h"
#include <iostream>

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>
#include "options/chain_option.h"
#include <cpsp/gecode/GC_ThreadingSymmBreaker.hh>
#include <algorithm>

using namespace Gecode;
using namespace std;
namespace {
/// Stop opbject for controlling stoppage based on both time and
/// failures.
class FailTimeStop : public Search::Stop {
private:
	Search::TimeStop *ts;
	Search::FailStop *fs;
	FailTimeStop(int fails, int time) {
		ts = new Search::TimeStop(time);
		fs = new Search::FailStop(fails);
	}
public:
	bool stop(const Search::Statistics& s) {
		return fs->stop(s) || ts->stop(s);
	}
	/// Create appropriate stop-object
	static Search::Stop* create(int fails, int time) {
		if (fails < 0 && time < 0) return NULL;
		if (fails < 0) return new Search::TimeStop( time);
		if (time  < 0) return new Search::FailStop(fails);
		return new FailTimeStop(fails, time);
	}
};
}

namespace cpsp{
//side chain threading
namespace scth{

///This class is for side chain handler, it applies the constraint threading for the most optimal Core till the less until a solution is found
class SideChainThreadingHandler{
private:
	/** A simple function returns the number of Hs in the sequence
	 *@param seq The sequence string
	 *@return the number of Hs
	 */
	int calculateCoreSize(std::string* seq);

	/** This function returns as a pair the number of odd and even Hs in the sequence
	 *@return pair representing the number of odd (first) and even (second) number of Hs
	 */
	pair<int,int> checkParity(const std::string* seq);

	/** This function returns as a pair the number of odd and even points of the HCore
	 *@return pair representing the number of odd (first) and even (second) points of the HCore
	 */
	pair<int,int> checkHCoreParity(const HCore* hCore);

	/** In the case of Cubic lattice only (not the FCC), this function will apply the parity check of the sequence against the Core
	 *@return true if the parity check of the sequence against the hCore holds otherwise false
	 */		 
	bool checkSequenceHCoreParity(const std::string *seq, const HCore* hCore);

	/// number of found solutions
	int solutionsNum;
public:
	///The constructor gets the Side chain options as the only input
	SideChainThreadingHandler(const SideChainOptions& );
	int getSolutionsNum();
};

///This class to handle running the side chain threading
class RunSideChainThreading{

	//The same input for the ChainThreading contructor,
	//with an option to decide how to print
public:

	///Just to delete brackets
	static void deleteBrackets(std::string& absMoves){
		for(size_t i=0;i<absMoves.size();){
			if(absMoves[i]=='(' || absMoves[i]==')' ){
				absMoves.erase(i,1);
			}
			i=i+1;
		}
	}

	/** Draw one solution
	 *@param seq The sequence string
	 *@param solution The solution as set of moves
	 *@param jmlHome The path of the jmlHome 
	 *@param viewerHome The path of the side chain viewer binary
	 *@param Description The lattice description as CUB (Cubic) or FCC
	 */
	static void drawElement(const std::string* seq,std::string solution,const std::string* jmlHome, const std::string* viewerHome,const std::string* Description ){
		assert(viewerHome!=NULL);
		assert(Description!=NULL);
		//Just make it upper case
		std::string modifiedSeq= *seq;
//		modifiedSeq=RunSideChainThreading::StringToUpper(modifiedSeq);
		//Remove the brackets from the solution
		std::string solutionWithoutBrackets=solution;
		RunSideChainThreading::deleteBrackets(solutionWithoutBrackets);
		string lat;
		//Draw
		if(*(Description)=="fcc" || *(Description)=="FCC")
			lat="FCC";
		else lat="CUB";
		string binary=*viewerHome+" -seq="+modifiedSeq+" -abs="+solutionWithoutBrackets+" -jmolHome="+*jmlHome+" -jmol -lat="+lat+"";
//		char buff[255];
//		strcpy(buff, binary.c_str());
//		int st1=system(buff);
	}

	/** Run the side chain threading
	 *@param sequence The string sequence
	 *@param latFram The lattice frame
	 *@param neighVecs The neighbor vectors (for effeciency, though it is redundant)
	 *@param hCore The HCore
	 *@param option The side chains options
	 */
	template <template<class> class Engine>
	static int run(const std::string *sequence,const biu::LatticeFrame* latFrame, const biu::LatticeFrame::index_set* neighVecs,cpsp::HCore* hCore,const cpsp::scth::SideChainOptions& option, const int maxSol){
		assert(sequence!=NULL);
		assert(latFrame!=NULL);
		assert(neighVecs!=NULL);
		assert(hCore!=NULL);
		//Obtaining the shifting vector
		biu::IPointVec* symmetryShiftVec=NULL;
		if(option.BreakSym)
			symmetryShiftVec = cpsp::gecode::GC_ThreadingSymmBreaker::generateGlobalShiftVec(
					hCore, latFrame, symmetryShiftVec);
		int solutions=0; //The number of solutions
		SideChainThreading* sideChainInstance=new SideChainThreading(sequence,latFrame,neighVecs,hCore,symmetryShiftVec,option.BreakSym,option.onlyHcoreRepresentatives);

		unsigned int n_p = 0;
		unsigned int n_b = 0;
		if (sideChainInstance->status() != SS_FAILED) {
			n_p = sideChainInstance->propagators();
			n_b = sideChainInstance->branchings();
		}

//		Search::Stop* stop = FailTimeStop::create(-1,-1);
		SideChainThreading* s = dynamic_cast<SideChainThreading*>(sideChainInstance);
		assert(s!=NULL);
		Engine<SideChainThreading> e(s,Search::Config::c_d,Search::Config::a_d,NULL);
		delete s;
		while(solutions< maxSol  ){
			SideChainThreading* ex = dynamic_cast<SideChainThreading*>(e.next());
			if (ex == NULL)
				break;
			
			if (!option.countOnly) {
				 // write move string to stream
				std::string sl = ex->print(sequence,latFrame,option.Output);
			}
			
//			if(option.withOutput){
//				//Selecting draw option
//				if(option.Draw && option.Output==cpsp::scth::MOVES){
//					cpsp::scth::RunSideChainThreading::drawElement(sequence,sl,option.JmlHome,option.ViewerHome,option.Description);
//					std::cout<<"Enter to 'c' continue, or any char to stop drawing"<<std::endl;
//					char drawChar;
//					std::cin>>drawChar;
//					if(drawChar!='c')
//						option.Draw=false;
//				}
//			}

			solutions++;
			delete ex;
		} 
		return solutions;
	}
};
}
}
#endif
