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

#ifndef SCBB_VIEWER_H
#define SCBB_VIEWER_H

#include <biu/OptionParser.hh>
#include <biu/LatticeModel.hh>
	//Side chain-Backbone lattice model viewer (in order to enable reading the new string with brackets
	class SCBBLatticeModelViewer : public biu::LatticeModel {
	  public:
		SCBBLatticeModelViewer(biu::LatticeDescriptor* lat):biu::LatticeModel(lat){};
		virtual biu::IPointVec absMovesToPoints( const biu::MoveSequence& absMoveSeq ) const;
		virtual ~SCBBLatticeModelViewer();
	};

	/**
	 * HPview creates an output in CML-file format of a sequence/structure that
	 * can be viewed with molecule viewers like Chime or JMol.
	 */
	//Side chain-Backbone model viewer 
	class SCBBViewer
	{
	 private:
	    //Just deleting the brackets of the input string and make it pure moves
	    void deleteBrackets(std::string& absMoves);
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
	
		SCBBViewer();
		
		virtual ~SCBBViewer();

	 		//! does conversion and printing according to the given program parameters
	 		//! @return the execution status of parameter parsing and conversion
	 	int convert(int argc, char** argv);
	};
#endif
