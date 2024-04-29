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

#include "HPsurface.hh"
#include "version.hh"

#include <iomanip>
#include <algorithm>
#include <iterator>

#include <cpsp/Exception.hh>

#include <biu/assertbiu.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>

	
	HPsurface::HPsurface()
	  : latDescr(NULL), seq(""), points(), pointSet()
	{
	}
	
	HPsurface::~HPsurface()
	{
	}

	unsigned int 
	HPsurface::getSurface( std::set<size_t> & surfaceIdx,  biu::IPointSet & hull, const char monomerType) 
	{
		assertbiu( latDescr != NULL, "no lattice descriptor present" );
		assertbiu( seq.size() == points.size(), "sequence and structure differ in length" );
		
		  // init data structures
		unsigned int numOfFreeContacts = 0;
		surfaceIdx.clear();
		hull.clear();
		const biu::LatticeNeighborhood& neighs = latDescr->getNeighborhood();
		  // check all neighbors of all H-monomers
		for ( std::string::size_type i = 0; i< seq.size(); i++) {
			if ( seq[i] != monomerType )	// skip the remaining
				continue;
			bool surfaceIdxNotAdded = true;
			for ( biu::LatticeNeighborhood::const_iterator n = neighs.begin();
					n != neighs.end(); n++ )
			{
				if ( pointSet.find( points[i] + (*n) )  == pointSet.end() ) { // neighbor is no structure element
					numOfFreeContacts++;	// count hull contact
					hull.insert( points[i] + (*n) ); // store hull position
					if (surfaceIdxNotAdded) {
						surfaceIdx.insert(i);
						surfaceIdxNotAdded = false;
					}
				}
			}
		}
		
		return numOfFreeContacts;
	}
	
	
	unsigned int 
	HPsurface::getInternalContacts( char type1, char type2, const biu::LatticeModel& lat ) {
		unsigned int ic = 0;
		
		for (std::string::size_type i=0; i<seq.size(); i++) {
			if (seq[i] != type1 && seq[i] != type2)
				continue;
			for (std::string::size_type j=i+2; j<seq.size(); j++) {
				if (seq[j] != type1 && seq[j] != type2)
					continue;
				if ((seq[i] == type1 && seq[j] == type2) ||
					(seq[i] == type2 && seq[j] == type1) )
					ic += ( lat.areNeighbored(points[i],points[j]) ? 1 : 0 );  
			}
		}
		
		return ic;
	}
	
	void
	HPsurface::initAllowedArguments(biu::OptionMap & allowedArgs,
										std::string &infoText ) const
	{
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"the HP-sequence of the protein"));
		allowedArgs.push_back(biu::COption(	
								"abs", true, biu::COption::STRING, 
								"the structure in absolute move string representation"));
		allowedArgs.push_back(biu::COption(	
								"rel", true, biu::COption::STRING, 
								"the structure in relative move string representation"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (SQR, CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string(" HPsurface calculates several surface properties")
			+ std::string(" of a lattice protein.")
			+ std::string("\n\n Available information are:")
			+ std::string("\n  - energy of the structure")
			+ std::string("\n  - internal contacts")
			+ std::string("\n  - free H/P contacts on the proteins surface")
			+ std::string("\n  - H/P monomers on surface of the protein")
			+ std::string("\n  - surrounding positions contacting H/P monomers")
			+ std::string("\n  - surface/interior (+/-) annotation of the protein positions")
			+ std::string("\n");
	}
	

	
	int 
	HPsurface::run( int argc, char** argv ) {

		using namespace biu;
		
		
		try {
			

			// check user parameters 
			biu::OptionMap allowedArgs;
			std::string infoText;
			initAllowedArguments(allowedArgs,infoText);	// init
			
				// parse programm arguments	
			COptionParser opts = COptionParser(	allowedArgs, argc, 
												(char**)argv, infoText);
				// arguments not parseable or not all mandatory arguments given
			if (!opts.noErrors()) {
				throw cpsp::Exception("",-1);
			}
				// help output
			if (opts.argExist("help")) {
				opts.coutUsage();
				throw cpsp::Exception("",0);
			}
			if (opts.argExist("version")) {
				giveVersion();
				throw cpsp::Exception("",0);
			}
				// set lattice
			std::string latStr = opts.getStrVal("lat");
			if (latStr.compare("SQR") == 0)
				latDescr = new LatticeDescriptorSQR();
			else if (latStr.compare("CUB") == 0)
				latDescr = new LatticeDescriptorCUB();
			else if (latStr.compare("FCC") == 0)
				latDescr = new LatticeDescriptorFCC();
			else
				throw cpsp::Exception("Unknown lattice type '"+latStr+"'",-2);
			LatticeModel lattice(latDescr);

				// get and validate sequence
			seq = opts.getStrVal("seq");
			if ( seq.size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no sequence given",-2);
			biu::Alphabet alph("HP",1);
			if (!alph.isAlphabetString(seq))
				throw cpsp::Exception("Error in arguments: sequence '"+seq+"' is no HP-sequence",-2);
			
				// get and validate structure
			points.clear();
			std::string absStructure = "";
			if (opts.argExist("abs") && opts.argExist("rel")) 
				throw cpsp::Exception("Error in arguments: relative and absolute structure given",-2);
			if (!opts.argExist("abs") && !opts.argExist("rel")) 
				throw cpsp::Exception("Error in arguments: no structure given",-2);
			if (opts.argExist("abs")) {
				std::string str = opts.getStrVal("abs");
				if (!latDescr->getAlphabet()->isAlphabetString(str))
					throw  cpsp::Exception("Error in arguments: structure '"+str+"'is no valid move string",-2);
				absStructure = str;
				points = lattice.absMovesToPoints(latDescr->getSequence(str));
			} else {
				std::string str = opts.getStrVal("rel");
				if (!latDescr->getAlphabet()->isAlphabetString(str))
					throw  cpsp::Exception("Error in arguments: structure '"+str+"'is no valid move string",-2);
				absStructure = latDescr->getString(lattice.relMovesToAbsMoves(latDescr->getSequence(str)));
				points = lattice.relMovesToPoints(latDescr->getSequence(str));
			}
			if (points.size() == 0) 
				throw  cpsp::Exception("Error: something wrong in structure conversion",-3);
			if (points.size() != seq.size()) 
				throw  cpsp::Exception("Error in arguments: sequence and structure differ in size",-2);
				// init point set
			pointSet.clear();
			pointSet.insert( points.begin(), points.end() ); // copy data
			
				// check for verbose output
			bool verboseOut = opts.argExist("v");
			
			if (verboseOut) {
				std::cout	<<"\n==================================="
							<<"\n   HPsurface - CPSP-tool-library"
							<<"\n==================================="
							<<"\n"
							<<"\n  sequence   : " << seq
							<<"\n  abs. moves : " << absStructure
							<<"\n  lattice    : " << latDescr->getName()
							<<"\n";
			}
			
			  // do surface investigation
			biu::IPointSet hHull;	// positions in contact to H-monomers
			biu::IPointSet pHull;	// positions in contact to P-monomers
			std::set<size_t> hSurfaceIdx; // index positions of points of H-monomers on the surface
			std::set<size_t> pSurfaceIdx; // index positions of points of H-monomers on the surface
			unsigned int numOfHmonomers = 0;
			for (size_t i = 0; i<seq.size(); i++) {if(seq[i]=='H') numOfHmonomers++;}
			  // get surface information
			unsigned int numOfFreeHcontacts = getSurface( hSurfaceIdx, hHull, 'H');
			unsigned int numOfFreePcontacts = getSurface( pSurfaceIdx, pHull, 'P');
			  // calculate hull positions in contact with H- and P-monomers
			biu::IPointSet hpHull;
			std::set_intersection(	hHull.begin(), hHull.end(),
									pHull.begin(), pHull.end(),
									std::insert_iterator< biu::IPointSet >(hpHull, hpHull.begin()) );
			  // number of structure positions not on the surface
			unsigned int noSurfacePos = seq.size() - hSurfaceIdx.size() - pSurfaceIdx.size();
			  // all surface indices
			std::set<size_t> surfaceIdx(hSurfaceIdx);
			surfaceIdx.insert(pSurfaceIdx.begin(), pSurfaceIdx.end());
			  // generate surface annotation
			std::string surfaceAnnotation(seq.size(),'-');
			for (size_t i=0; i<surfaceAnnotation.size(); i++) {
				if (surfaceIdx.find(i) != surfaceIdx.end())
					surfaceAnnotation[i] = '+';
			}
			
			  // output
			std::cout	<<"\n------------------ sequence information ------------------------\n"
						<<"\n HP-sequence                : " <<seq
						<<"\n absolute move string       : " <<lattice.getString(lattice.pointsToAbsMoves(points))
						<<"\n lattice                    : " <<latDescr->getName()
						<<"\n energy                     : " <<std::setw(7) <<-(int)getInternalContacts('H','H', lattice)
						<<"\n length                     : " <<std::setw(7) <<seq.size();

			std::cout	<<"\n\n------------------- internal contacts --------------------------\n"
						<<"\n # of H-H-contacts          : " <<std::setw(7) <<getInternalContacts('H','H', lattice)
						<<"\n # of H-P-contacts          : " <<std::setw(7) <<getInternalContacts('H','P', lattice)
						<<"\n # of P-P-contacts          : " <<std::setw(7) <<getInternalContacts('P','P', lattice);

			std::cout	<<"\n\n--------------------- free contacts ----------------------------\n"
						<<"\n # free H-contacts          : " <<std::setw(7) <<numOfFreeHcontacts
						<<"\n # free P-contacts          : " <<std::setw(7) <<numOfFreePcontacts;
			if ( numOfFreeHcontacts==0 || numOfFreePcontacts==0 )
				std::cout	<<"\n   ratio free H/P-contacts  : not calculatable";
			else	
				std::cout<<"\n   ratio free H/P-contacts  : " 
							<<std::setprecision(5) <<std::fixed 
							<<(double(numOfFreeHcontacts)/double(numOfFreePcontacts));
					
			std::cout	<<"\n\n-------------- structure positions on surface ------------------\n"
						<<"\n # H-monomers on surface    : " <<std::setw(7) <<hSurfaceIdx.size() 
						<<" of " <<std::setw(4) <<numOfHmonomers
						<<"\n # P-monomers on surface    : " <<std::setw(7) <<pSurfaceIdx.size()
						<<" of " <<std::setw(4) <<(seq.size()-numOfHmonomers)
						<<"\n # no surface monomers      : " <<std::setw(7) <<noSurfacePos
						<<" of " <<std::setw(4) <<seq.size();
			if ( hSurfaceIdx.size()==0 || pSurfaceIdx.size()==0 ) 
				std::cout	<<"\n   ratio surface H/P-mon.   : not calculatable";
			else
				std::cout	<<"\n   ratio surface H/P-mon.   : " 
							<<std::setprecision(5) <<std::fixed 
							<<(double(hSurfaceIdx.size())/double(pSurfaceIdx.size()));
				
			std::cout	<<"\n\n------------------ surrounding positions -----------------------\n"
						<<"\n # overall surface pos.     : " <<std::setw(7) <<(hHull.size()+pHull.size()-hpHull.size())
						<<"\n # H-contacting positions   : " <<std::setw(7) <<hHull.size()
						<<"\n # P-contacting positions   : " <<std::setw(7) <<pHull.size();
			if ( hHull.size()==0 || pHull.size()==0 ) 
				std::cout	<<"\n   ratio H/P-cont. pos.     : not calculatable";
			else
				std::cout	<<"\n   ratio H/P-cont. pos.     : " 
							<<std::setprecision(5) <<std::fixed 
							<<(double(hHull.size())/double(pHull.size()))
							<<"\n";
			std::cout	<<"\n # shared H-/P-cont. pos.   : " <<std::setw(7) <<hpHull.size();

			std::cout	<<"\n\n-------------------- surface annotation ------------------------\n"
						<<(verboseOut?"\n explanation  : '+' = surface, '-' = interior":"")
						<<"\n\n" <<seq
						<<"\n" <<surfaceAnnotation
						<<"\n";

			std::cout	<< std::endl;
			
			if (verboseOut) {
				std::cout	<<"\n====================================\n"
							<<std::endl;
			}
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}
		
		if (latDescr != NULL) delete latDescr;

		return 0;
	}

	int main(int argc, char** argv) {
		HPsurface hps;
			// run surface investigation
		return hps.run(argc, argv);
	}
