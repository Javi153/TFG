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

#include "chain_threading.h"
#include <biu/Point.hh>
#include <iostream>
#include <cpsp/gecode/GC_Search.hh>
#include <cpsp/gecode/GC_RankViewSel.hh>
#include <cpsp/gecode/GC_ThreadingSymmBreaker.hh>
#include "cpsp/gecode/GC_FirstSolBranching.hh"
#include "options/chain_option.h"
#include <biu/LatticeDescriptorCUB.hh>
namespace cpsp{
//side chain threading
namespace scth{

SideChainThreading::SideChainThreading(const std::string *sequence,
		const biu::LatticeFrame* latFrame, 
		const biu::LatticeFrame::index_set* neighVecs,cpsp::HCore* hCore,
		const biu::IPointVec* shiftVec,
		bool withBreakingSymmetry,
		bool branchHonly)
	: cpsp::gecode::SuperSpace()
	, ranks( 2*sequence->size(), 100)
{
	typedef cpsp::gecode::MyViewType MyViewType;
	
	//Set all the domains
	assert(sequence!=NULL);

	//Settting the domain of the Hs-Array
	setHDomain(latFrame,hCore);
	//Setting the domain of the BHs-Array
	setBHDomain(latFrame,hCore);
	//Setting the domain of the Ps-Array
	setOutCoreDomain(sequence,latFrame,hCore);
	//Now applying all constrains
	//All H's are in the HCore
	inCore(hCore);

	//All the back bones of H's are 1 distance unit far from the H-Core
	attachedToCore(hCore );

	//All none H are out of the HCore
	outCore(sequence,hCore );

	//Just for test
	/*biu::LatticeFrame::index_set test= hCore->getIndexedHull(latFrame, 0);
	   for (biu::LatticeFrame::index_set::const_iterator t=test.begin(); t!= test.end();t++){
	      std::cout<<(*t)<<std::endl;
	      biu::IntPoint tpoint=latFrame->getPoint(*t);
	      std::cout<<" " <<tpoint.getX()<<" "<<tpoint.getY()<<" "<<tpoint.getZ()<<std::endl;
	   }
	  biu::IntPoint ph1(7,8,9);
	  int phi1=latFrame->getIndex(ph1);
	  Hs[3]=IntVar(this,phi1,phi1);
	  biu::IntPoint ph2(8,8,8);
	  int phi2=latFrame->getIndex(ph2);
	  Hs[0]=IntVar(this,phi2,phi2);
	  biu::IntPoint ph3(8,8,9);
	  int phi3=latFrame->getIndex(ph3);
	  Hs[2]=IntVar(this,phi3,phi3);
	  biu::IntPoint ph4(9,8,9);
	  int phi4=latFrame->getIndex(ph4);
	  Hs[1]=IntVar(this,phi4,phi4);*/
	/*biu::IntPoint p1(7,9,9);
	  int pi1=latFrame->getIndex(p1);
	  BHs[3]=IntVar(this,pi1,pi1);
	  biu::IntPoint p2(8,9,8);
	  int pi2=latFrame->getIndex(p2);
	  BHs[0]=IntVar(this,pi2,pi2);
	  biu::IntPoint p3(8,9,9);
	  int pi3=latFrame->getIndex(p3);
	  BHs[2]=IntVar(this,pi3,pi3);
	  biu::IntPoint p4(9,9,9);
	  int pi4=latFrame->getIndex(p4);
	  BHs[1]=IntVar(this,pi4,pi4); 
	  biu::IntPoint p5(8,8,7);
	  int pi5=latFrame->getIndex(p5);
	  Ps[0]=IntVar(this,pi5,pi5);
	  biu::IntPoint p6(9,8,8);
	  int pi6=latFrame->getIndex(p6);
	  Ps[1]=IntVar(this,pi6,pi6);
	  biu::IntPoint p7(8,9,7);
	  int pi7=latFrame->getIndex(p7);
	  BPs[0]=IntVar(this,pi7,pi7);
	  biu::IntPoint p8(9,9,8);
	  int pi8=latFrame->getIndex(p8);
	  BPs[1]=IntVar(this,pi8,pi8);*/
	//End test

	//Neighbouring constrains for the back bones
	BsNeighbours(sequence,neighVecs);

	//Neighbourig constrains between back bones and side chains P
	B_PsNeighbours(neighVecs);

	//Neighbourig constrains between back bones and side chains P
	B_HsNeighbours(neighVecs);

	//Self avoiding
	selfAvoid();
	//Branching
	Gecode::ViewArray<MyViewType> allViews(this,Hs.size()+BHs.size()+BPs.size()+Ps.size());
	Gecode::ViewArray<MyViewType> HViews(this,Hs.size());
	int qindx=0;

	for(int indx=0;indx<Hs.size();indx++){
		allViews[qindx] = MyViewType(Hs[indx],qindx);
		HViews[qindx] = MyViewType(Hs[indx],qindx);
		ranks[qindx] -= 20;
		qindx++;
	}

	for(int indx=0;indx<BHs.size();indx++){
		allViews[qindx] = MyViewType(BHs[indx],qindx);
		ranks[qindx] -= 10;
		qindx++;
	}

	for(int indx=0;indx<BPs.size();indx++){
		allViews[qindx] = MyViewType(BPs[indx],qindx);
		qindx++;
	}

	for(int indx=0;indx<Ps.size();indx++){
		allViews[qindx] = MyViewType(Ps[indx],qindx);
		qindx++;
	}

	cpsp::gecode::GC_ValNearCenter::latFrame = latFrame;
	if(withBreakingSymmetry){
		Gecode::BoolVarArray cpInit = cpsp::gecode::GC_ThreadingSymmBreaker::initCpByCore(this, hCore, latFrame, shiftVec);
		cpsp::gecode::GC_ThreadingSymmBreaker symmBreaker(this,latFrame, shiftVec, cpInit);
		new (this)
		Gecode::Search::SBViewValBranching< MyViewType,
				int,
//				Gecode::Int::Branch::ByDegreeMax,
				cpsp::gecode::GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>,
				cpsp::gecode::GC_ValNearCenter,
				cpsp::gecode::GC_ThreadingSymmBreaker >
			(this, allViews, symmBreaker);
	}
	
	if (branchHonly) {
		//This is branching without breaking the symmetry
		new (this) Gecode::ViewValBranching< MyViewType
							, int
//							, Gecode::Int::Branch::ByDegreeMax
							, cpsp::gecode::GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
							, cpsp::gecode::GC_ValNearCenter> 
					(this, HViews);
		//This is branching without breaking the symmetry
		new (this) cpsp::gecode::FirstSolBranching< MyViewType
							, int
//							, Gecode::Int::Branch::ByDegreeMax
							, cpsp::gecode::GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
							, cpsp::gecode::GC_ValNearCenter> 
					(this, allViews);
	} else {
		//This is branching without breaking the symmetry
		new (this) Gecode::ViewValBranching< MyViewType
							, int
//							, Gecode::Int::Branch::ByDegreeMax
							, cpsp::gecode::GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
							, cpsp::gecode::GC_ValNearCenter> 
					(this, allViews);
	}
}

SideChainThreading::SideChainThreading(bool share, SideChainThreading& old)
  :
	cpsp::gecode::SuperSpace(share, old)
	, ranks( old.ranks)
{
	Hs.update(this, share, old.Hs);
	Ps.update(this, share, old.Ps);
	BHs.update(this, share, old.BHs);
	BPs.update(this, share, old.BPs);
}

Gecode::Space* SideChainThreading::copy(bool share)
{
	return new SideChainThreading(share, *this);
}



/* returns the rank of the variable with the given index
 * 
 * @param index of the variable of interest
 * @return its rank
 */ 
unsigned int
SideChainThreading::
getRank(const int index) const
{
	assert( index >= 0 && index < (int)ranks.size() /*index out of range*/);
	return ranks[index];
}



SideChainThreading::~SideChainThreading(){}

void SideChainThreading::inCore(cpsp::HCore* hCore){
	Hs=IntVarArray(this,hCore->getSize(),HDomain);
	//Hs=IntVarArray(this,hCore->getSize(),1000,3000);	
}

void SideChainThreading::attachedToCore(cpsp::HCore* hCore){
	BHs=IntVarArray(this,hCore->getSize(),BHDomain);
	//BHs=IntVarArray(this,hCore->getSize(),1000,3000);	
}
//passing the sequence as parameter for the effeciency
void SideChainThreading::printSol(const std::string* seq,const biu::LatticeFrame* latFrame){
	int hIndx=0;
	int pIndx=0;
	//	   std::cout<<"#Cubic Protein Solution:"<<std::endl;
	for(size_t i=0;i<seq->size();i++){
		if( (*seq)[i]=='H' ){
			biu::IntPoint point=latFrame->getPoint(Hs[hIndx].val());
			std::cout<<"h: " <<point.getX()<<" "<<point.getY()<<" "<<point.getZ()<<std::endl;

			point=latFrame->getPoint(BHs[hIndx].val());
			std::cout<<"bh: " <<point.getX()<<" "<<point.getY()<<" "<<point.getZ()<<std::endl;
			hIndx++;
		}
		else{
			biu::IntPoint point=latFrame->getPoint(Ps[pIndx].val());
			std::cout<<"p: " <<point.getX()<<" "<<point.getY()<<" "<<point.getZ()<<std::endl;

			point= latFrame->getPoint(BPs[pIndx].val());
			std::cout<<"bp: " <<point.getX()<<" "<<point.getY()<<" "<<point.getZ()<<std::endl;
			pIndx++;
		}
	}
}

std::string SideChainThreading::printAsMoves(const std::string* seq,const biu::LatticeFrame* latFrame){
	std::string latName=latFrame->getDescriptor()->getName();
	const biu::LatticeDescriptorCUB* cubicDescInstance=
		dynamic_cast<const biu::LatticeDescriptorCUB*>(latFrame->getDescriptor());
		if(cubicDescInstance){
			int hIndx=0;
			int pIndx=0;
			std::string s="";

			//		std::cout<<"#"<< latName<<" Protein Solution:"<<std::endl;
			for(size_t i=0;i<seq->size();i++){
				biu::IntPoint point;
				biu::IntPoint bpoint;
				if( (*seq)[i]=='H' ){
					point=latFrame->getPoint(Hs[hIndx].val());
					bpoint=latFrame->getPoint(BHs[hIndx].val());
					hIndx++;
				}
				else{
					point=latFrame->getPoint(Ps[pIndx].val());
					bpoint= latFrame->getPoint(BPs[pIndx].val());
					pIndx++;
				}
				if( bpoint.getX()!=point.getX() ){
					if( bpoint.getX() > point.getX() )
						s+="(L)";
					else s+="(R)";
				}
				else if( bpoint.getY()!=point.getY() ){
					if( bpoint.getY() > point.getY() )
						s+="(B)";
					else s+="(F)";

				}
				else{
					if( bpoint.getZ() > point.getZ() )
						s+="(D)";
					else s+="(U)";
				}
				//This is for the move to next backbone
				if(i<seq->size()-1){
					biu::IntPoint nextbpoint;
					if( (*seq)[i+1]=='H' ){
						nextbpoint=latFrame->getPoint(BHs[hIndx].val());
					}
					else{
						nextbpoint=latFrame->getPoint(BPs[pIndx].val());
					}
					if( bpoint.getX()!=nextbpoint.getX() ){
						if( bpoint.getX() > nextbpoint.getX() )
							s+="L";
						else s+="R";
					}
					else if( bpoint.getY()!=nextbpoint.getY() ){
						if( bpoint.getY() > nextbpoint.getY() )
							s+="B";
						else s+="F";

					}
					else{
						if( bpoint.getZ() > nextbpoint.getZ() )
							s+="D";
						else s+="U";
					}

				}//End if
			}//End fore
			std::cout<<s<<std::endl;
			return s;
		}

		else{
			int hIndx=0;
			int pIndx=0;
			std::string s="";
			//		std::cout<<"#"<< latName<<" Protein Solution:"<<std::endl;
			for(size_t i=0;i<seq->size();i++){
				biu::IntPoint point;
				biu::IntPoint bpoint;
				if( seq->at(i)=='H' ){
					point=latFrame->getPoint(Hs[hIndx].val());
					bpoint=latFrame->getPoint(BHs[hIndx].val());
					hIndx++;
				}
				else{
					point=latFrame->getPoint(Ps[pIndx].val());
					bpoint= latFrame->getPoint(BPs[pIndx].val());
					pIndx++;
				}
				biu::IntPoint absVec = point - bpoint;
				//Get the moves for the FCC
				biu::Move absMove =
					latFrame->getDescriptor()->getNeighborhood().getElement(absVec).getMove();

				std::string moveStr =
					latFrame->getDescriptor()->getAlphabet()->getString( absMove );
				//
				s+="("+moveStr+")";
				//This is for the move to next backbone
				if(i<seq->size()-1){
					biu::IntPoint nextbpoint;
					if( (*seq)[i+1]=='H' ){
						nextbpoint=latFrame->getPoint(BHs[hIndx].val());
					}
					else{
						nextbpoint=latFrame->getPoint(BPs[pIndx].val());
					}
					biu::IntPoint nextAbsVec = nextbpoint - bpoint;
					//Get the moves for the FCC
					biu::Move nextAbsMove =
						latFrame->getDescriptor()->getNeighborhood().getElement(nextAbsVec).getMove();

					std::string nextMoveStr =
						latFrame->getDescriptor()->getAlphabet()->getString( nextAbsMove );
					//
					s+=nextMoveStr;

				}//End if

			}
			std::cout<<s<<std::endl;
			return s;
		}
}

//Printing for drawing
void SideChainThreading::printAsDrowMoves(const std::string* seq,const biu::LatticeFrame* latFrame){
	int hIndx=0;
	int pIndx=0;
	std::string s="";
	//	   std::cout<<"#"<< latFrame->getDescriptor()->getName()<<" Protein Solution:"<<std::endl;
	for(size_t i=0;i<seq->size();i++){
		biu::IntPoint point;
		biu::IntPoint bpoint;
		if( (*seq)[i]=='H'){
			point=latFrame->getPoint(Hs[hIndx].val());
			bpoint=latFrame->getPoint(BHs[hIndx].val());
			hIndx++;
		}
		else{
			point=latFrame->getPoint(Ps[pIndx].val());
			bpoint= latFrame->getPoint(BPs[pIndx].val());
			pIndx++;
		}
		if( bpoint.getX()!=point.getX() ){
			if( bpoint.getX() > point.getX() )
				s+="LR";
			else s+="RL";
		}
		else if( bpoint.getY()!=point.getY() ){
			if( bpoint.getY() > point.getY() )
				s+="BF";
			else s+="FB";

		}
		else{
			if( bpoint.getZ() > point.getZ() )
				s+="DU";
			else s+="UD";
		}
		//This is for the move to next backbone
		if(i<seq->size()-1){
			biu::IntPoint nextbpoint;
			if( (*seq)[i+1]=='H' ){
				nextbpoint=latFrame->getPoint(BHs[hIndx].val());
			}
			else{
				nextbpoint=latFrame->getPoint(BPs[pIndx].val());
			}
			if( bpoint.getX()!=nextbpoint.getX() ){
				if( bpoint.getX() > nextbpoint.getX() )
					s+="L";
				else s+="R";
			}
			else if( bpoint.getY()!=nextbpoint.getY() ){
				if( bpoint.getY() > nextbpoint.getY() )
					s+="B";
				else s+="F";

			}
			else{
				if( bpoint.getZ() > nextbpoint.getZ() )
					s+="D";
				else s+="U";

			}

		}//End if

	}//End fore
	std::cout<<s<<std::endl;
}

std::string SideChainThreading::print(const std::string* seq,const biu::LatticeFrame* latFrame,int option){
	assert(latFrame!=NULL);
	assert(seq!=NULL);
	//Option 0 means print the coordination
	if(option==cpsp::scth::NORMAL){
		printSol(seq,latFrame);
		return "";
	}
	else if(option==cpsp::scth::MOVES){
		return printAsMoves(seq,latFrame);
	}
	else{
		printAsDrowMoves(seq,latFrame);
		return "";
	}
}
//All none H are out of the HCore
void  SideChainThreading::outCore(const std::string *sequence,cpsp::HCore* hCore){
	assert(sequence!=NULL);
	Ps=IntVarArray(this,sequence->size()-hCore->getSize(),OutCoreDomain);
	BPs=IntVarArray(this,sequence->size()-hCore->getSize(),OutCoreDomain);
	//Ps=IntVarArray(this,sequence->size()-hCore->getSize(),1000,3000);
	//BPs=IntVarArray(this,sequence->size()-hCore->getSize(),1000,3000);
}

void SideChainThreading::setHDomain(const biu::LatticeFrame* latFrame,cpsp::HCore* hCore){
	biu::LatticeFrame::index_set hCoreDomain= hCore->getIndexedHull(latFrame, 0);
	biu::LatticeFrame::index_type HDomainArr[hCoreDomain.size()];
	int i=0;
	for (biu::LatticeFrame::index_set::const_iterator it=hCoreDomain.begin(); it!= hCoreDomain.end();it++){
		HDomainArr[i]=*it;
		i++;
	}
	HDomain=IntSet(HDomainArr,hCoreDomain.size());
}

void SideChainThreading::setBHDomain(const biu::LatticeFrame* latFrame,cpsp::HCore* hCore){
	biu::LatticeFrame::index_set bhCoreDomain= hCore->getIndexedHull(latFrame, 1);
	biu::LatticeFrame::index_type BHDomainArr[bhCoreDomain.size()];
	int i=0;
	for (biu::LatticeFrame::index_set::const_iterator it=bhCoreDomain.begin(); it!= bhCoreDomain.end();it++){
		BHDomainArr[i]=*it;
		i++;
	}
	BHDomain=IntSet(BHDomainArr,bhCoreDomain.size());
}

//Setting the out-core domain
void SideChainThreading::setOutCoreDomain(const std::string *sequence,const biu::LatticeFrame* latFrame,cpsp::HCore* hCore){
	int maxHull=getMaxHull(sequence);
	biu::LatticeFrame::index_set outCoreDomain= hCore->getIndexedHull(latFrame, maxHull);
	biu::LatticeFrame::index_type outDomainArr[outCoreDomain.size()];
	int i=0;
	for (biu::LatticeFrame::index_set::const_iterator it=outCoreDomain.begin(); it!= outCoreDomain.end();it++){
		outDomainArr[i]=*it;
		i++;
	}
	OutCoreDomain=IntSet(outDomainArr,outCoreDomain.size());
}

//Neighbourig constrains between back bones and side chains P
void SideChainThreading::B_PsNeighbours(const biu::LatticeFrame::index_set* neighVecs){
	assert(Ps.size()==BPs.size());
	for(int i=0; i<Ps.size(); i++) {
		GECODE_ES_FAIL(this,cpsp::gecode::GC_LatticeNeighbored2::post(this,Ps[i],BPs[i],neighVecs,Gecode::ICL_DOM,150));
	}
}

//Neighbourig constrains between back bones and side chains H
void SideChainThreading::B_HsNeighbours(const biu::LatticeFrame::index_set* neighVecs){
	assert(Hs.size()==BHs.size());
	for(int i=0; i<Hs.size(); i++) {
		GECODE_ES_FAIL(this,cpsp::gecode::GC_LatticeNeighbored2::post(this,Hs[i],BHs[i],neighVecs,Gecode::ICL_DOM,150));
	}
}

//Neighbouring constrains for the back bones
void SideChainThreading::BsNeighbours(const std::string *sequence,const biu::LatticeFrame::index_set* neighVecs){
	assert(sequence!=NULL);
	//The view of all back bones
	ViewArray<Int::IntView> backBones=ViewArray<Int::IntView>(this,sequence->size());
	int hIndex=0;
	int pIndex=0;
	for(size_t i=0;i<sequence->size();i++){
		if((*sequence)[i]=='H'){
			backBones[i]=BHs[hIndex];
			hIndex++;
		}
		else{
			backBones[i]=BPs[pIndex];
			pIndex++;
		}
	}

	for(int i=1; i<backBones.size(); i++) {
		GECODE_ES_FAIL(this,cpsp::gecode::GC_LatticeNeighbored2::post(this,backBones[i],backBones[i-1],neighVecs,Gecode::ICL_DOM,150));
	}
}

//Self avoiding
void SideChainThreading::selfAvoid(){
	distinct(this,Hs, Gecode::ICL_DOM);
	//Everything except the HCore
	IntVarArgs B_PsVars(Ps.size()+BPs.size()+BHs.size());
	int k=0;
	for(int i=0;i<Ps.size();i++){
		B_PsVars[k]=Ps[i];
		k++;
	}
	for(int i=0;i<BHs.size();i++){
		B_PsVars[k]=BHs[i];
		k++;
	}

	for(int i=0;i<BPs.size();i++){
		B_PsVars[k]=BPs[i];
		k++;
	}
	distinct(this,B_PsVars,Gecode::ICL_DOM);
}

//Get the maximum hull possible hull
int SideChainThreading::getMaxHull(const std::string *sequence){
	assert(sequence!=NULL);
	int maxHull=0;
	bool first=true;
	for(size_t i=0;i<sequence->size();){
		//Escaping the hs
		while((*sequence)[i]=='H' ){

			i++;
			if(i>=sequence->size())
				break;
		}

		//Break when finishing the sequence
		if(i>=sequence->size())
			break;

		//Count the ps
		int pNumber=0;
		int newDistance=0;
		while((*sequence)[i]=='P' ){
			pNumber++;
			i++;
			if(i>=sequence->size())
				break;
		}//End while

		//If this is the first sequence of ps then take the distance according to the most left p
		if(first){
			first=false;
			newDistance=pNumber;
		}
		//If this is the last sequence of ps then take the distance according to the most right p
		else if(i>=sequence->size()){
			newDistance=pNumber;
		}
		//Take the distance according to the middle p
		else{
			if(pNumber % 2==0)
				newDistance=pNumber/2;
			else newDistance=(pNumber+1)/2;
		}
		maxHull=std::max(maxHull,newDistance);
		if(i>=sequence->size())
			break;

	}//end for
	//Because of the side chain
	maxHull=maxHull+2;
	return maxHull;
}

}

}

/*int main(int arg, char** args){

   return 0;
}*/
