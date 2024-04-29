/*
 *  Main authors:
 *     Sebastian Will 
 *
 *  Contributing authors:
 *     Martin Mann
 *
 *  Copyright:
 *     Sebastian Will, 2007
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef GC_THREADINGSYMMBREAKER_HH_
#define GC_THREADINGSYMMBREAKER_HH_


#include <gecode/kernel.hh>
#include <gecode/search.hh>
#include <gecode/int.hh>
#include <biu/LatticeFrame.hh>
#include "cpsp/gecode/GC_Search.hh"
#include "GC_searchExt.hh" // == local copy of symmbreaking out of dds-gecode-branch

#ifdef DEBUG
  #include <iostream>
#endif

// ------------------------------------------------------------

namespace cpsp { namespace gecode { 
    
/**
 * @brief SymmBreaker for Graph Coloring 
 *
 * Used for breaking all value transposition symmetries
 *
 */
class GC_ThreadingSymmBreaker : public Gecode::Search::SymmBreaker<MyViewType,int> {
public:
	typedef biu::IPointVec GlobalShiftVec;
private:
	const biu::LatticeFrame *latticeFrame;
	const GlobalShiftVec* shiftVec;
public:

	GC_ThreadingSymmBreaker(Gecode::Space *home,
	const biu::LatticeFrame *latticeFrame_, 
	const GlobalShiftVec* shiftVec_)
	    : Gecode::Search::SymmBreaker<MyViewType,int>(home,latticeFrame_->getDescriptor()->getAutomorphisms().size()),
		  latticeFrame(latticeFrame_), shiftVec(shiftVec_)
		{
		}
	
	GC_ThreadingSymmBreaker(Gecode::Space *home,
	const biu::LatticeFrame *latticeFrame_, 
	const GlobalShiftVec* shiftVec_,
	Gecode::BoolVarArray &init_cp)
	    : Gecode::Search::SymmBreaker<MyViewType,int>(home,latticeFrame_->getDescriptor()->getAutomorphisms().size(), init_cp),
		  latticeFrame(latticeFrame_), shiftVec(shiftVec_)
		{
		}
	
	// copy constructor
	GC_ThreadingSymmBreaker(Gecode::Space *home, bool share, const GC_ThreadingSymmBreaker &sym)
	    : Gecode::Search::SymmBreaker<MyViewType,int>(home,share,sym), 
	      latticeFrame(sym.latticeFrame), shiftVec(sym.shiftVec)
		{}
	
	virtual ~GC_ThreadingSymmBreaker() {}

	virtual GC_ThreadingSymmBreaker* copy(Gecode::Space* home, bool share) {
		return new (home) GC_ThreadingSymmBreaker(home, share, *this);
	}
	
	
protected:
	
		//! Returns the symmetric point index of oldInd using Automorphism symmetry.
	static 
	int 
	getSymmetricIndex(int symmetry, int oldInd, 
		const biu::LatticeFrame *latticeFrame,
		const GlobalShiftVec* shiftVec) 
	{
		assertbiu(latticeFrame != NULL && shiftVec != NULL, "parameter == NULL");
		assertbiu(latticeFrame->getDescriptor()->getAutomorphisms().size() == shiftVec->size(), "different sizes of automorphisms and shift vectors");
		assertbiu(symmetry >= 0 && (biu::AutomorphismVec::size_type)symmetry < latticeFrame->getDescriptor()->getAutomorphisms().size(), "symmetry is no valid index ot automorphism vector");
		
			// calculate symmetric point
		biu::IntPoint symm_p = latticeFrame->getPoint(oldInd);
		symm_p = latticeFrame->getDescriptor()->getAutomorphisms()[symmetry] * symm_p;
		
			// returning index of shifted symmetric point
		return latticeFrame->getIndex(symm_p + shiftVec->at(symmetry));
	}
	
    // symmetry 0 has to be the identity
	//
	Gecode::BoolVar 
	reifySymm(Gecode::Space *home,int symmetry, Gecode::ViewArray<MyViewType> &x, int pos, int val) const {

		int sym_val = GC_ThreadingSymmBreaker::getSymmetricIndex(symmetry, val, 
						latticeFrame, shiftVec);
		
	    Gecode::BoolVar b((Gecode::Space*)home,0,1);
	    rel((Gecode::Space*)home,x[pos],Gecode::IRT_EQ,sym_val,b);
	    return b;
	}
	
public:

		//! Generates the inital Cp variables for the constructor, that represent
		//! if a symmetry is already broken or not.
	static 
	Gecode::BoolVarArray 
	initCpByCore(Gecode::Space* home, HCore* hCore, const biu::LatticeFrame *latticeFrame,
		const GlobalShiftVec* shiftVec ) 
	{
		assertbiu(home != NULL && latticeFrame != NULL && shiftVec != NULL, "parameter == NULL");
		assertbiu(latticeFrame->getDescriptor()->getAutomorphisms().size() == shiftVec->size(), "different sizes of automorphisms and shift vectors");

		biu::AutomorphismVec::size_type size = latticeFrame->getDescriptor()->getAutomorphisms().size(), sym=0; 
		Gecode::BoolVarArray cp = Gecode::BoolVarArray(home, size);
		int sym_val = 0;
			// get shifted core data
		const biu::LatticeFrame::index_set core = hCore->getIndexedHull(latticeFrame, 0);
		biu::LatticeFrame::index_set::const_iterator it;
		
		bool isNotBroken = true;
		for (sym=0; sym<size; sym++) {	// for each symmetry
			isNotBroken = true;
				// for all core points --> test if symm. point is in core again
			for (it = core.begin(); isNotBroken == true && it!= core.end(); it++) {
				sym_val = GC_ThreadingSymmBreaker::getSymmetricIndex(sym, *it, latticeFrame, shiftVec);
				isNotBroken = isNotBroken && (core.find(sym_val) != core.end());  
			}
				// init cp variable with true or false
			cp[sym] = Gecode::BoolVar(home, isNotBroken, isNotBroken);
		}
		return cp;
	}
	
		//! Generates the shift vectors, that are needed calculate the 
		//! symmetries. 
	static
	GlobalShiftVec*
	generateGlobalShiftVec(HCore* hCore, const biu::LatticeFrame *latticeFrame,
		GlobalShiftVec* shiftVec = NULL)
	{
		assertbiu(hCore != NULL && latticeFrame != NULL, "paramter == NULL");
			// init shiftVec
		biu::AutomorphismVec symms = latticeFrame->getDescriptor()->getAutomorphisms();
		biu::AutomorphismVec::size_type i=0, size=symms.size();
		if (shiftVec == NULL)
			shiftVec = new GlobalShiftVec(size, biu::IntPoint(0,0,0));
		else
			shiftVec->resize(size, biu::IntPoint(0,0,0));
		
			// calculate core frame corners
		const biu::IntPoint coreFrameMin =  hCore->getCoreFrame()->getMinCorner() 
			+ latticeFrame->getCenter() - hCore->getCoreFrame()->getCenter();
		const biu::IntPoint coreFrameMax =  hCore->getCoreFrame()->getMaxCorner() 
			+ latticeFrame->getCenter() - hCore->getCoreFrame()->getCenter();
		biu::IntPoint c1, c2;
			// calculate global shift vectors based on the shift of the core frame
		for ( i = 0; i < size; i++)
		{
			c1 = symms[i] * coreFrameMin;
			c2 = symms[i] * coreFrameMax;
												
			shiftVec->at(i) = coreFrameMin 	- biu::IntPoint(	
									std::min(c1.getX(), c2.getX()),
									std::min(c1.getY(), c2.getY()),
									std::min(c1.getZ(), c2.getZ()));
		}
			
		return shiftVec;
	}
	
	

	bool 
	status(const Gecode::Space *home) const {
	    // if all Cps are 0 (except cps[0] which is for identity),
	    // return false
	    
		bool allbroken=true;
	    for(int i=1; i<cps.size(); i++) {
			if (cps[i].min()!=cps[i].max() || cps[i].val()==1) {
				allbroken=false;
				break;
			}
	    }
	    return !allbroken;
	}
};

} } // namespaces
// ------------------------------------------------------------
#endif
