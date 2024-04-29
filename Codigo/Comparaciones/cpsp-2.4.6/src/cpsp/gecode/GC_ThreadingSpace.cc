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



//	#define TMPOUT

#include <set>

#include "cpsp/gecode/GC_ThreadingSpace.hh"
#include "cpsp/gecode/GC_Search.hh"
#include "cpsp/gecode/GC_StlSetRangeIterator.hh"
#include "cpsp/gecode/GC_LatticeNeighbored2.hh"
#include "cpsp/gecode/GC_FirstSolBranching.hh"
#include "cpsp/gecode/GC_ValNearCenter.hh"
#include "cpsp/gecode/GC_RankViewSel.hh"

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/minimodel.hh>
#include <gecode/int/branch.hh>

#include <biu/LatticeFrame.hh>

#ifdef TMPOUT
	#include <iostream>
#endif

namespace cpsp
{
  namespace gecode
  {

  		// just a helper function ... TODO : delete ;)
	void print(biu::LatticeFrame::index_set& x, const biu::LatticeFrame* latFrame) {
		std::cout<<"( ";
		for (biu::LatticeFrame::index_set::const_iterator it = x.begin(); it != x.end(); it++)
			std::cout <<latFrame->getPoint(*it) <<"| ";
		std::cout <<")";
		std::cout.flush();
	}

  		// just a helper function ... TODO : delete ;)
	void GC_ThreadingSpace::print() const {
      std::cerr<<"\n";
      for (int i=0; i< domains.size(); i++)  {
    	const Gecode::Int::IntView vd(domains[i]);
    	Gecode::Int::ViewRanges< Gecode::Int::IntView > v(vd);
      	std::cerr <<" " <<i <<"(";
      	while(v()) {
      		if (v.min() != v.max())
      			std::cerr <<v.min() <<".." <<v.max()<<"," ;
      		else
      			std::cerr <<v.min() <<",";
      		++v;
      	}
      	std::cerr <<")\n";
      }
      std::cerr<<"\n";
		}

	GC_ThreadingSpace::GC_ThreadingSpace(const std::string *sequence,
			const biu::LatticeFrame* latFrame,
			const biu::LatticeFrame::index_set* neighVecs,
			SeqFeatureMap* seqFeatureMap,
			const HullLevel *hullLvl,
			HCore* hCore,
			const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
			int branchingType,
			NeighPropLvl neighPropLvl ) 
	:
		domains(	this, sequence->size(),
					0, latFrame->getMaxIndex())
		, ranks( sequence->size(), 100)
	{

		assert(hCore != NULL);
			// temp variables
		biu::LatticeFrame::index_set data;
		GC_StlSetRangeIterator dataIt;
		Gecode::IntSet dataSet;
		Gecode::IntVarArgs hVars(hCore->getSize()), pVars(domains.size() - hCore->getSize());
		int hVarI = 0, pVarI = 0;

			// domain init

		int i=0;

		for (i=0; i<domains.size(); i++) {
				// get domain data
			data = hCore->getIndexedHull(latFrame, hullLvl->at(i));
				// create iterator 
			dataIt.init(&data);
				// shrink domains to given data in dataIt via IntView
			Gecode::Int::IntView actDom(domains[i]);
				// use macro for security
			GECODE_ME_FAIL(this, actDom.inter(this,dataIt));
				// init distinct viewvecs
			if (hullLvl->at(i) == 0) {
				hVars[hVarI] = domains[i]; hVarI++;
			} else {
				pVars[pVarI] = domains[i]; pVarI++;
			}
		}

		
			// p-singlet init
		std::vector<unsigned int>::const_iterator it;
			// get domain data
		data = hCore->getIndexedPSinglets(latFrame);
			// create iterator 
		dataIt.init(&data);
			// shrink domains 
		for (	it = (*seqFeatureMap)[P_SINGLET].begin();
				it != (*seqFeatureMap)[P_SINGLET].end(); it++)
		{
				// view to singlet domain
			Gecode::Int::IntView actDom(domains[*it]);
				// use macro for security
			GECODE_ME_FAIL(this, actDom.inter(this,dataIt));
				// decrease rank of P singlets
			ranks[*it]   -= 10;
		}

			// h-singlet init
			// get domain data
		data = hCore->getIndexedHSinglets(latFrame);
			// create iterator 
		dataIt.init(&data);
			// shrink domains 
		for (	it = (*seqFeatureMap)[H_SINGLET].begin();
				it != (*seqFeatureMap)[H_SINGLET].end(); it++)
		{
				// view to singlet domain
			Gecode::Int::IntView actDom(domains[*it]);
				// use macro for security
			GECODE_ME_FAIL(this, actDom.inter(this,dataIt));
		}

			// handle ranks of HP borders
			// NOTE : singlets are handled twice --> will have lower rank !
		std::string seqPattern = "HP";
		std::string::size_type pos = 0;
		while( (pos = sequence->find(seqPattern,pos)) != std::string::npos) {
			ranks[pos]   -= 10;
			ranks[pos+1] -= 10;
			pos++;
		}
		seqPattern = "PH"; pos = 0;
		while( (pos = sequence->find(seqPattern,pos)) != std::string::npos) {
			ranks[pos]   -= 10;
			ranks[pos+1] -= 10;
			pos++;
		}

		
			// init constraints 

			// ALL DIFFERENT
			// Supports value (ICL_VAL, default), bounds (ICL_BND),
			// and domain-consistency (ICL_DOM)
		Gecode::distinct(this, hVars, Gecode::ICL_DOM);
		Gecode::distinct(this, pVars, Gecode::ICL_DOM);

			// binary lattice neighbor constraints along the sequence
		switch(neighPropLvl) {
		
		case CUSTOM_2 : {	
						for(i=1; i<domains.size(); i++) {
							GECODE_ES_FAIL(this,
							GC_LatticeNeighbored2::post(	this,
													domains[i-1], domains[i],
													neighVecs,
													Gecode::ICL_DOM,
													INT_MAX
												));
						}
						break; }
		case NO_PROP : 
		default : 	//  NO NEIGHBOR PROPAGATION INSERTED
					break;
		}

			// GC_ValNearCenter lattice frame init
		GC_ValNearCenter::latFrame = latFrame;

		// branching
		
		if ( branchingType != BR_NONE ) {

			assertbiu( !(branchingType & BR_DDS) , "DDS branching currently not supported"); 

			Gecode::ViewArray<MyViewType> all(this,domains.size());
			for (int i=0; i<domains.size(); i++) {
				all[i] = MyViewType(domains[i],i);
			}
			
				// symmetry breaking
			if ( branchingType & BR_SYM ) {
	
				Gecode::BoolVarArray cpInit = GC_ThreadingSymmBreaker::initCpByCore(this, hCore, latFrame, shiftVec);
				GC_ThreadingSymmBreaker symmBreaker(this,latFrame, shiftVec, cpInit);
				
				typedef Gecode::Search::SBViewValBranching < MyViewType, int
//		        			, Gecode::Int::Branch::ByDegreeMax
							, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
							, GC_ValNearCenter
							, GC_ThreadingSymmBreaker 
						> MySymmBranching;
				
				// final branching on all variables
				new (this) MySymmBranching(this, all, symmBreaker);
			}
	
				// DFS branching
			if ( branchingType & BR_DFS ) {
		        typedef Gecode::ViewValBranching< MyViewType, int
//		        			, Gecode::Int::Branch::ByDegreeMax
							, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
		        			, GC_ValNearCenter
		        		> MyBranching;
				
				// final branching on all variables
				new (this) MyBranching(this, all);
			}
		}
	}



	GC_ThreadingSpace::GC_ThreadingSpace(bool share, GC_ThreadingSpace& toCopy):
		SuperSpace(share, toCopy), domains(), ranks(toCopy.ranks)
	{
		domains.update(this, share, toCopy.domains);
	}


	GC_ThreadingSpace::~GC_ThreadingSpace()
	{
	}


	SuperSpace*
	GC_ThreadingSpace::copy(bool share) {
		return new GC_ThreadingSpace(share, *this);
	}

	biu::LatticeFrame::index_vec
	GC_ThreadingSpace::getSolution() const {
		biu::LatticeFrame::index_vec retVec;
		for (int i=0; i<domains.size();i++) {
			if (domains[i].assigned()) {
				retVec.push_back((biu::LatticeFrame::index_type)domains[i].val());
			} else {
				std::cerr <<"Error: Variable "<<i<<" is not assigned !!! size = "<<domains[i].size()<<"\n"; 
				std::cerr <<"       Are you using DDS ?\n";
			}
		}
		return retVec;
	}


	
	 /* returns the rank of the variable with the given index
	  * 
	  * @param index of the variable of interest
	  * @return its rank
	  */ 
	unsigned int
	GC_ThreadingSpace::
	getRank(const int index) const
	{
		assert( index >= 0 && index < (int)ranks.size() /*index out of range*/);
		return ranks[index];
	}



///////////////////////////////////////////////////////////////////7
// GC_ThreadingSpaceShapes
///////////////////////////////////////////////////////////////////7



	GC_ThreadingSpaceShapes::GC_ThreadingSpaceShapes(
		const std::string *sequence,
		const biu::LatticeFrame* latFrame,
		const biu::LatticeFrame::index_set* neighVecs,
		SeqFeatureMap* seqFeatureMap,
		const HullLevel *hullLvl,
		HCore* hCore,
		const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
		int branchingType)
	  : GC_ThreadingSpaceHcoreDist( sequence, latFrame, neighVecs, seqFeatureMap,
	  		hullLvl, hCore, shiftVec, BR_NONE) 
	{

			// GC_ValNearCenter lattice frame init
		GC_ValNearCenter::latFrame = latFrame;

		// same branching but use old ranks for symm breaking
		
		if ( branchingType != BR_NONE ) {

			assertbiu( !(branchingType & BR_DDS) , "DDS branching currently not supported"); 

				// viewarray for branching
			Gecode::ViewArray<MyViewType> all(this,domains.size());
			Gecode::ViewArray<MyViewType> Hvars(this,hCore->getSize());
			Gecode::ViewArray<MyViewType> Pvars(this,domains.size() - hCore->getSize());
			int iH=0, iP=0;
			for (int i=0; i<domains.size(); i++) {
				all[i] = MyViewType(domains[i],i);
				if (sequence->at(i) == 'H') {
					Hvars[iH] = MyViewType(domains[i],i); iH++;
				} else {
					Pvars[iP] = MyViewType(domains[i],i); iP++;
				}
			}
				// add P singlets to H var array 
			Gecode::ViewArray<MyViewType> FirstVars(this,hCore->getSize()+(*seqFeatureMap)[P_SINGLET].size());
			std::vector<unsigned int>::const_iterator it;
			for (iH=0; iH < Hvars.size(); iH++ ) {
				FirstVars[iH] = MyViewType(Hvars[iH]);
			}
			for (	it = (*seqFeatureMap)[P_SINGLET].begin();
					it != (*seqFeatureMap)[P_SINGLET].end(); it++)
			{
				FirstVars[iH] = MyViewType(domains[*it],*it); iH++;
			}

				
				// symmetry breaking
			if ( branchingType & BR_SYM ) {
	
				Gecode::BoolVarArray cpInit = GC_ThreadingSymmBreaker::initCpByCore(this, hCore, latFrame, shiftVec);
				GC_ThreadingSymmBreaker symmBreaker(this,latFrame, shiftVec, cpInit);

				typedef Gecode::Search::SBViewValBranching < MyViewType,int
//								, Gecode::Int::Branch::ByDegreeMax
							, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
							, GC_ValNearCenter
							, GC_ThreadingSymmBreaker 
						> MySymmBranching;
				
				// final branching on all variables
				new (this) MySymmBranching(this, all, symmBreaker);
			}
	
				// DFS branching
			if ( branchingType & BR_DFS ) {
		        typedef Gecode::ViewValBranching< MyViewType, int
//		        				, Gecode::Int::Branch::ByDegreeMax
	        				, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
	        				, GC_ValNearCenter
                        > MyBranching;
		                                        
				// first branch exhaustively on all H variables
//				new (this) MyBranching(this, all);
				new (this) MyBranching(this, FirstVars);
//				new (this) MyBranching(this, Hvars);
				
		        typedef FirstSolBranching< MyViewType, int
//		        				, Gecode::Int::Branch::ByDegreeMax
	        				, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
	        				, GC_ValNearCenter
                        > MyOneSolBranching;
		                                        
				// final branching for one solution on P variables
				new (this) MyOneSolBranching(this, Pvars);
			}
			
		}
	}



	GC_ThreadingSpaceHcoreDist::
	GC_ThreadingSpaceHcoreDist(const std::string *sequence,
		const biu::LatticeFrame* latFrame,
		const biu::LatticeFrame::index_set* neighVecs,
		SeqFeatureMap* seqFeatureMap,
		const HullLevel *hullLvl_,
		HCore* hCore,
		const GC_ThreadingSymmBreaker::GlobalShiftVec* shiftVec,
		int branchingType
		)
	  : GC_ThreadingSpace( sequence, latFrame, neighVecs, seqFeatureMap,
	  		hullLvl_, hCore, shiftVec, BR_NONE),
  		hullLvl(hullLvl_),
  		hCoreSize(hCore->getSize())
	{
		
		//####### CALCULATE NEW RANKS WHERE ranks(H) < ranks(P)  #############//
		if (branchingType != BR_NONE) {
			
			  // initialize ranks of H and P variables differently ranks(H) < ranks(P)
			for (size_t i=0; i<hullLvl->size(); i++) {
				if (hullLvl->at(i) == 0) {
					ranks[i] = 100;
				} else {
					ranks[i] = 200;
				}
			}
			
				// handle ranks of HP borders
				// NOTE : singlets are handled twice --> will have lower rank !
			std::vector<unsigned int>::const_iterator it;
			std::string seqPattern = "HP";
			std::string::size_type pos = 0;
			while( (pos = sequence->find(seqPattern,pos)) != std::string::npos) {
				ranks[pos]   -= 10;
				ranks[pos+1] -= 10;
				pos++;
			}
			seqPattern = "PH"; pos = 0;
			while( (pos = sequence->find(seqPattern,pos)) != std::string::npos) {
				ranks[pos]   -= 10;
				ranks[pos+1] -= 10;
				pos++;
			}
			

				// GC_ValNearCenter lattice frame init
			GC_ValNearCenter::latFrame = latFrame;

			// branching
			
			if ( branchingType != BR_NONE ) {

				assertbiu( !(branchingType & BR_DDS) , "DDS branching currently not supported"); 

					// viewarray for branching
				Gecode::ViewArray<MyViewType> all(this,domains.size());
				Gecode::ViewArray<MyViewType> Hvars(this,hCore->getSize());
				Gecode::ViewArray<MyViewType> Pvars(this,domains.size() - hCore->getSize());
				int iH=0, iP=0;
				for (int i=0; i<domains.size(); i++) {
					all[i] = MyViewType(domains[i],i);
					if (sequence->at(i) == 'H') {
						Hvars[iH] = MyViewType(domains[i],i); iH++;
					} else {
						Pvars[iP] = MyViewType(domains[i],i); iP++;
					}
				}
					
					// symmetry breaking
				if ( branchingType & BR_SYM ) {
		
					Gecode::BoolVarArray cpInit = GC_ThreadingSymmBreaker::initCpByCore(this, hCore, latFrame, shiftVec);
					GC_ThreadingSymmBreaker symmBreaker(this,latFrame, shiftVec, cpInit);

					typedef Gecode::Search::SBViewValBranching < MyViewType,int
//								, Gecode::Int::Branch::ByDegreeMax
								, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
								, GC_ValNearCenter
								, GC_ThreadingSymmBreaker 
							> MySymmBranching;
					
					// final branching on all variables
					new (this) MySymmBranching(this, all, symmBreaker);
				}
		
					// DFS branching
				if ( branchingType & BR_DFS ) {
			        typedef Gecode::ViewValBranching< MyViewType, int
//		        				, Gecode::Int::Branch::ByDegreeMax
		        				, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
		        				, GC_ValNearCenter
                            > MyBranching;
			                                        
    				// first branch exhaustively on all H variables
    				new (this) MyBranching(this, Hvars);
    				
			        typedef FirstSolBranching< MyViewType, int
//		        				, Gecode::Int::Branch::ByDegreeMax
		        				, GC_RankViewSel<Gecode::Int::Branch::ByDegreeMax>
		        				, GC_ValNearCenter
                            > MyOneSolBranching;
			                                        
    				// final branching for one solution on P variables
    				new (this) MyOneSolBranching(this, Pvars);
				}
			}


		}
	}


  } // namespace gecode
}  // namespace cpsp


