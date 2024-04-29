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

#ifndef HPVIEW_HH_
#define HPVIEW_HH_

#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>

	/**
	 * HPview creates an output in CML-file format of a sequence/structure that
	 * can be viewed with molecule viewers like Chime or JMol.
	 */
	class HPview
	{
	 protected:
	 
	 		//! types of viewing programs supported
	 	enum VIEW_TYPE {JMOL, NONE};

			//! the script for jmol viewing	 	
	 	const std::string JMOL_SCRIPT;
	 
			//! Inits the allowed parameters for HPview.
		void initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const;
											
			//! Prints absolute positions to stream in CML-file format.									
		void points2cml(	const std::string & seq,
							const biu::IPointVec & points, 
							const biu::LatticeModel * lattice,
							std::ostream & out = std::cout,
							int distance = 3) const;
						
	public:
	
		HPview();
		
		virtual ~HPview();

	 		//! does conversion and printing according to the given program parameters
	 		//! @return the execution status of parameter parsing and conversion
	 	int convert(int argc, char** argv);
	};

#endif /*HPVIEW_HH_*/
