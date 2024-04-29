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

#ifndef HPOPTDEG_HH_
#define HPOPTDEG_HH_

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

#include <cpsp/gecode/GC_HPThreading.hh>

//! compressed HP-sequence representation 
typedef std::vector<unsigned char> CSequence;
typedef std::set<std::string> StringSet;
	
typedef std::map<CSequence,unsigned int> CSeqDegHash;


//! local search for low degenerated HP-sequences	
class HPoptdeg {
public:
	
    typedef std::vector< int > ComponentLabel;
		
		
protected:
	
    //! Manages statistical data of HPoptdeg runs.
    class Statistic {
    public:
	//! the runtime
	double time;
	//! number sequences hashed during the run
	CSeqDegHash::size_type hashEntries;
	//! number of sequences in the neutral net
	Statistic() : time(0.0), hashEntries(0)
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
		
    //! maps compressed sequences to their degeneracy 
    CSeqDegHash seq2deg;
		
    unsigned int limitDeg; //!< limit for guessing random sequences
    unsigned int tgtDeg; //!< stop when this deg is found
	    
    biu::Alphabet alphabetHP;
		
    Statistic * stat;
		
    bool verboseOut;
    
    bool hpnx;
    
    std::string::size_type seqLength;

    std::string startSeq;


    float TEMP; // for log metropolis

    int maxSteps; //!< break after such many steps
    int timeLimit; //!< break after this time limit
    
    //! determine exactly only those degeneracies,
    //! where the metropolis criterion probability is
    //! greater than cutoff. If the degeneracy is
    //! larger, the sequence is not
    //! accepted. (Exceptions due to caching.) This
    //! is used to limit the search time, since we
    //! can stop the search each time too many
    //! structures are found.
    double cutoff; 

    //! handler for the cpsp threading object
    cpsp::gecode::GC_HPThreading threadingHandler;

    //! counter for optimization steps
    int steps;
    
    //! a random double uniformly from [0,1[
    double random_double() const;

    void
    initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText ) const;

    //! initialises the threading object for the next run with the 
    //! given sequence (i.e. reinit of H-core database access)
    bool
    initNextThreading(const std::string & seq, int limitDeg);

    CSequence
    compressSeq(const std::string& str) const;

    std::string
    uncompressSeq(const CSequence& cs, size_t seqLength) const;
    
    //! get the degeneracy of a sequence,
    //! either by running predictor or from cache
    int 
    get_degeneracy(std::string sequence, unsigned int limitDeg);
    
    //! run cpsp threading in order to determine degeneracy
    int
    cpsp_deg(std::string sequence, int limitDeg);

    std::string 
    mutate(const std::string &sequence,bool hpnx=false) const;

    std::string 
    point_mutate(const std::string &sequence, int pos) const;
    
    //! calc value for metropolis criiterion
    //! pre next_deg > deg
    double logmetropolis(int deg,int next_deg) const;
    
    //! return a random hp sequence
    std::string random_hp(int len) const;
    
    //! return a random hpnx sequence
    std::string random_hpnx(int len) const;

    //! get true random numbers
    unsigned int 
    read_dev_urandom() const;
public:
		
    bool
    initThreading(int argc, char** argv);
    
      //! construction
    HPoptdeg();
    
      //! destruction
    virtual 
    ~HPoptdeg();
	    
	  //! performs the optimization process to find low degenerated sequences
    void
    optimizeDegeneracy();

	  //! argument parsing and local search
    int
    run(int argc, char** argv);
		
};

#endif /*HPOPTDEG_HH_*/
