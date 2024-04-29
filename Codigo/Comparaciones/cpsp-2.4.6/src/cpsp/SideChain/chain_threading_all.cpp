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

#include "chain_threading_all.h"

#include <cpsp/HCoreDatabaseFILE.hh>
#include <biu/LatticeDescriptor.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <cpsp/HCoreDatabase.hh>
#include <biu/Point.hh>
#include <stdlib.h>
#include <iomanip>

using namespace cpsp;
using namespace std;
namespace cpsp{
//side chain threading
namespace scth{
//Global functionto convert to lower cases


int SideChainThreadingHandler::calculateCoreSize(std::string* seq){
	int coreSize=0;
	for(size_t i=0;i<seq->size();i++){
		if((*seq)[i]=='H'){
			coreSize++;
		}
	}
	return coreSize;
}

// Just to check the parity for the sequence
pair<int,int> SideChainThreadingHandler::checkParity(const std::string* seq){
	int odd=0;
	int even=0;
	for(size_t i=0;i<seq->size();i++){
		if((*seq)[i]=='H'){
			if((i % 2) == 0)
				even+=1;
			else odd+=1;
		}
	}
	return pair<int,int>(odd,even);
}
//Just to check the parity of the H-Core
pair<int,int> SideChainThreadingHandler::checkHCoreParity(const HCore* hCore){
	assert(hCore!=NULL);
	biu::IPointSet points=hCore->getPoints();
	int odd=0;
	int even=0;
	for (biu::IPointSet::const_iterator t=points.begin(); t!= points.end();t++){
		biu::IntPoint point=*t;
		int sum=point.getX()+point.getY()+point.getZ();
		if( (sum %2) == 0)
			even+=1;
		else odd+=1;
	}
	return pair<int,int>(odd,even);
}

//to check if the HCore is suitable for the Lattice
bool SideChainThreadingHandler::checkSequenceHCoreParity(const std::string *seq, const HCore* hCore){
	pair<int,int> parity=checkParity(seq);
	pair<int,int> coreParity=checkHCoreParity(hCore);
	if(parity.first==coreParity.first || parity.second==coreParity.first)
		return true;
	else return false;
}


//Draw one element
/*static void drawElement(const std::string* seq,const std::string& jmlHome, const std::string& viewerHome ){
	}*/
//Do Threading main operation
SideChainThreadingHandler::SideChainThreadingHandler(const SideChainOptions& option){
	solutionsNum=0;
//	std::string* dbRootPath=option.DBRoot;
	cpsp::HCoreDatabase* db = option.coreDB; //new cpsp::HCoreDatabaseFILE (*dbRootPath);
	assert(db!=NULL);
	assert(db->isConnected());
	std::string* seq = option.Sequence;
	assert(seq!=NULL);
	
	if (option.withOutput)
		std::cout <<"==> HP-sequence = " <<*seq <<std::endl;


	//Choose the descriptor, Descriptor factory
	biu::LatticeDescriptor* latDescr;
	bool useParity=false;
	if(*(option.Description)=="CUB"){
		latDescr=new biu::LatticeDescriptorCUB();
		useParity=true;
	}
	else{
		latDescr=new biu::LatticeDescriptorFCC();
	}
	//Getting the size of the core and inializing the parameters of the side chain threading process
	int size=this->calculateCoreSize(seq);
	biu::LatticeFrame* latFrame= new biu::LatticeFrame(latDescr,2*(seq->size()+2));
	db->initCoreAccess(*latFrame->getDescriptor() ,size);
	biu::LatticeFrame::index_set* neighVecs=new biu::LatticeFrame::index_set(latFrame->getIndexedNeighborhood ());
	HCore* hCore=new HCore();
	int allSolutions=0;
	int curSolutions=0;

	int oldContacts=0;//Just initializing
	int maxContacts=0;
	int coreNumber=0;
	while(db->getNextCore(*hCore)){
		//for(Iterator it=)
		int newContacts=hCore->getContacts();
//		if(allSolutions>0 && newContacts!=oldContacts && oldContacts!=0)
//			break;
		if (newContacts!=oldContacts) {
			coreNumber=0;
			if (allSolutions >= option.SolutionsNum) {
				break;
			}
			if (oldContacts!=0) {
				if(option.verbose) {
					std::cout<<"No solutions found for HCore with "<< oldContacts<<" number of contacts"<<std::endl;
				}
				if (option.withOutput)
					std::cout	<<"\n ==> number generated structures = "
								<<(curSolutions)
								<<std::endl;
			}
			curSolutions = 0;
			if (option.withOutput) {
				std::cout	<<"\n==> energy = -"
							<<newContacts
							<<std::endl;
			}
			oldContacts=newContacts;
		} else {
			coreNumber++;
		}
		 // store maximum
		maxContacts = std::max(maxContacts,newContacts);
		
		//Make sure that the number of solutions match that of the options
//		if(allSolutions>0 && option.SolutionsNum>-1){
//			option.SolutionsNum=option.SolutionsNum-allSolutions;
//			//std::cout<<allSolutions<<std::endl;
//			//std::cout<<option.SolutionsNum<<std::endl;
//		}
//		if(option.SolutionsNum == 0) break;
		//Using parity check only for cubic lattice
		if(useParity && !checkSequenceHCoreParity(seq,hCore)) {
			if(option.verbose) {
				std::cout<<"No solution found for HCore number: "<<coreNumber<<" with "<<oldContacts<<" number of contacts"<<std::endl;
			}
			continue;
		}
		
		int sol=0;
		if (option.onlyHcoreRepresentatives) {
			sol=scth::RunSideChainThreading::run<BAB>(seq,latFrame,neighVecs,hCore,option, option.SolutionsNum-allSolutions);
		} else {
			sol=scth::RunSideChainThreading::run<DFS>(seq,latFrame,neighVecs,hCore,option, option.SolutionsNum-allSolutions);
		}
		
		if (sol > 0) {
			if (option.verbose) {
				std::cout	<<"   + " <<(sol)
							<<" structures" <<std::endl;
			}
			if (option.bestOnly) {
				db->setActMinHHcontacts(hCore->getContacts());
			}
		}
		
		if(sol==0 && option.verbose) {
			std::cout<<"No solution found for HCore number: "<<coreNumber<<" with "<<oldContacts<<" number of contacts"<<std::endl;
		}
		//Calculating the solutions
		curSolutions += sol;
		allSolutions += sol;
		/*if(sol>0)
		   break;*/
	}
	if (option.withOutput) {
		std::cout	<<"\n ==> number generated structures = "
					<<(curSolutions)
					<<"\n\n- total number of generated structures = "
					<<std::setw(15) <<allSolutions
					<<"\n- maximal contacts reachable           = "
					<<std::setw(15) << maxContacts
					<<std::endl;
	}
	solutionsNum=allSolutions;
	delete(neighVecs);
	delete(latFrame);
	delete(latDescr);
	/*delete(seq);
	   delete(dbRootPath);*/
	delete(hCore);
//	delete(db);
}

//Just to retrun the number of solutions
int SideChainThreadingHandler::getSolutionsNum(){
	return solutionsNum;
}


}
}

