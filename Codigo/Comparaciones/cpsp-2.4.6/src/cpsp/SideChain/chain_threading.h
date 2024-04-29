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

#ifndef THREADING_H_
#define THREADING_H_
#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/minimodel.hh>

#include <gecode/int/branch.hh>

#include <string>
#include <vector>
#include <cpsp/HCore.hh>

#include <cpsp/gecode/GC_LatticeNeighbored2.hh>
#include <cpsp/gecode/GC_ValNearCenter.hh>

#include <cpsp/gecode/gecodeExtensions.hh>

#include <biu/LatticeFrame.hh>

#include "cpsp/gecode/GC_ThreadingSuperSpace.hh"


using namespace Gecode;
namespace cpsp{
//side chain threading
namespace scth{
class SideChainThreading : public cpsp::gecode::SuperSpace {
private:
	
	//! ranks for the view selection ordering (returned by getRank(..))
	std::vector<unsigned int> ranks;

	///The H side chains domain, indexed version
	IntSet HDomain;

	///Thd domain of the H-backbones, also indexed integer values
	IntSet BHDomain;

	///The P side chains and BP-domain (Ps'backbones domains), indexed version
	IntSet OutCoreDomain;
	
	//Here are the constraints
	/** Neighbouring constrains for the backbones, each backbone should be a neighbor to the backbone of the following amino acid in the sequence.
	 */
	void BsNeighbours(const std::string *sequence,const biu::LatticeFrame::index_set* neighVecs);

	/** Neighbourig constrains between the backbones and the side chains of Ps
	 */
	void B_PsNeighbours(const biu::LatticeFrame::index_set* neighVecs);

	/** Neighbourig constrains between the backbones and side chains os Hs
	 */
	void B_HsNeighbours(const biu::LatticeFrame::index_set* neighVecs);

	/** All H's are in the HCore, this constraint will make sure of this
	 */
	void inCore(cpsp::HCore* );

	///All the back bones of H's are 1 distance unit far from the H-Core. This function propagates the necessary contraints for this
	void attachedToCore(cpsp::HCore* );

	///All none H are out of the HCore
	void outCore(const std::string *sequence,cpsp::HCore* );

	///Self avoiding
	void selfAvoid();

	//Here are helping functions
	///Setting domain of the side chains of the Hs
	void setHDomain(const biu::LatticeFrame* latFrame,cpsp::HCore* hCore);

	///Setting the out-core domain. The out core contains all backbones and side chains of Ps
	void setOutCoreDomain(const std::string *sequence,const biu::LatticeFrame* latFrame,cpsp::HCore* hCore);

	///Settiing the BH domain (the backbones of Hs)
	void setBHDomain(const biu::LatticeFrame* latFrame,cpsp::HCore* hCore);

	/** Get the maximum hull. This is a global maximum. The range for each variable is set afterward.
	 *@return the maximum possible value for the out core domain
	 */
	int getMaxHull(const std::string *sequence);

	///Printing the solution points This function print the solutions as coordinates.
	void printSol(const std::string* seq,const biu::LatticeFrame* latFrame);

	/** Printing as moves and return those moves
	 *@param seq The sequence
	 *@param latFrame the lattice frame (Cubic or Lattice)
	 *@return the moves
	 */
	std::string printAsMoves(const std::string* seq,const biu::LatticeFrame* latFrame);

	///Printing for drawing. This function is justs for testing. It is private do not use it!
	void printAsDrowMoves(const std::string* seq,const biu::LatticeFrame* latFrame);
public:
	/** The constructor for threading one sequence to one Core
	 * @param sequence the sequence string 
	 * @param latFrame the lattice (FCC or Cubic)
	 * @param neighVecs the neighboring vector, for effeciency this is a parameter although it could be got from the lattice
	 * @param hCore the HCore
	 * @param shiftVec This is null as default (Will be adjusted automatically in the program)
	 * @param withBreakingSymmetry The default value for breaking the symmetry is true
	 * @param branchHonly if false: a ViewValBranching on all variables is posted. if true: a ViewValBranching is only posted on the H side chain variables, on the remaining variables a FirsSolBranching is posted. 
	 */
	SideChainThreading(const std::string *sequence,
			const biu::LatticeFrame* latFrame
			, const biu::LatticeFrame::index_set* neighVecs
			, cpsp::HCore* hCore
			, const biu::IPointVec* shiftVec=NULL
			, const bool withBreakingSymmetry=true
			, const bool branchHonly=false
			);

	///The copy constructor
	SideChainThreading(bool , SideChainThreading& );

	///Copy constrcutor for the gecode
	virtual Gecode::Space* copy(bool share);

	/** This is the public function for printing the solutions
	 * 1 to print the solution
	 * 2 for printing the directions (U,D,L,R,F,B)
	 */
	virtual std::string print(const std::string* sequence,const biu::LatticeFrame* latFrame,int option=0);


	/**
	 * The destructor	
	 */
	virtual ~SideChainThreading();

	///The positions of the h-side chain sequence. This array to be branched on.
	IntVarArray Hs;

	///The positions of the p-side chain sequence. This array to be branched on
	IntVarArray Ps;

	///The positions of the back bones of H. This array to be branched on
	IntVarArray BHs;

	///The positions of the back bones of P. This array to be branched on
	IntVarArray BPs;

	
	 /*! returns the rank of the variable with the given index
	  * 
	  * @param index of the variable of interest
	  * @return its rank
	  */ 
	virtual
	unsigned int
	getRank(const int index) const;

	//! adds additional constraints if BAB search is used to ensure
	//! a minimal Hcore difference between the next and all other
	virtual void constrain(SideChainThreading* lastSolution)
	{
		
		  // temporary variables to post the new constraint on
		Gecode::IntArgs hPos(Hs.size());
		
		  // collect Hcore assignments of last solution
		for (int i=0; i<Hs.size(); i++) {
			hPos[i]  = lastSolution->Hs[i].val();
		}
		
		  // ensure that at least one Hcore variable is differently assigned
		Gecode::atmost(this, Hs, hPos, (Hs.size()-1), Gecode::ICL_DOM);
		
		  // propagate the constraint
		this->status();
	}

};
}
}
#endif

