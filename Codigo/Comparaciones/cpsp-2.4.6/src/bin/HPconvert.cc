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


#include "HPconvert.hh"
#include "version.hh"

#include <sstream>
#include <fstream>
#include <limits.h>

#include <biu/util/Util_String.h>
#include <biu/assertbiu.hh>
#include <biu/Alphabet.hh>
#include "cpsp/Exception.hh"
#include "cpsp/PDBsupport.hh"

#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>


	HPconvert::HPconvert()
	{
	}
	
	HPconvert::~HPconvert()
	{
	}
	
	int
	HPconvert::convert(int argc, char** argv)
	{
		using namespace biu;
		using namespace util;
		
			// the lattice descriptor
		LatticeDescriptor* descr = NULL;
		
		try {
			
			std::string* seq = NULL;

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
				descr = new LatticeDescriptorSQR();
			else if (latStr.compare("CUB") == 0)
				descr = new LatticeDescriptorCUB();
			else if (latStr.compare("FCC") == 0)
				descr = new LatticeDescriptorFCC();
			else
				throw cpsp::Exception("Unknown lattice type '"+latStr+"'",-2);
			LatticeModel lattice(descr);
				// conversion directions (mode)
			Mode mode = NONE;
			std::string modeStr = opts.getStrVal("m");
			if (modeStr.compare("a2r") == 0) {
				mode = A2R;
			} else if (modeStr.compare("a2p") == 0) {
				mode = A2P;
			} else if (modeStr.compare("r2a") == 0) {
				mode = R2A;
			} else if (modeStr.compare("r2p") == 0) {
				mode = R2P;
			} else if (modeStr.compare("p2a") == 0) {
				mode = P2A;
			} else if (modeStr.compare("p2n") == 0) {
				mode = P2N;
			} else if (modeStr.compare("p2r") == 0) {
				mode = P2R;
			} else if (modeStr.compare("a2c") == 0) {
				mode = A2C;
			} else if (modeStr.compare("a2n") == 0) {
				mode = A2N;
			} else if (modeStr.compare("r2n") == 0) {
				mode = R2N;
			} else if (modeStr.compare("r2c") == 0) {
				mode = R2C;
			} else if (modeStr.compare("a2x") == 0) {
				mode = A2X;
			} else 
				throw cpsp::Exception("Unsupported mode '"+modeStr+"'",-2);
				// check structure input
			std::string structure = opts.getStrVal("str");
			if ( structure.size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no structure given",-2);
			switch (mode) {
				case P2A : // point file 
				case P2N :
				case P2R : {
					break;
				} 
				case A2R : // move string
				case A2P :
				case R2A :
				case A2N :
				case R2N :
				case A2X :
				case R2P : {
					if (!lattice.getDescriptor()->getAlphabet()->isAlphabetString(structure))
						throw  cpsp::Exception("Structure '"+structure+"'is no valid move string",-2);
					if (opts.argExist("seq")) {
						seq = new std::string(opts.getStrVal("seq"));
						if (seq->size() != lattice.getDescriptor()->getAlphabet()->getSequence(structure).size()+1)
							throw  cpsp::Exception("Sequence '"+(*seq)+"'is too short for given move string",-2);
						biu::Alphabet hpAlph("HP",1);
						if (!hpAlph.isAlphabetString(*seq)) 
							throw  cpsp::Exception("Sequence '"+(*seq)+"'is no HP sequence",-2);
					}
					break;
				} 
				default : assertbiu(false,"should never happen (unsupported conversion identifier found)"); // should never happen
			}
			
				// check for verbose output
			bool verboseOut = opts.argExist("v");
			
				// the outstream to write to
			std::ostream& out = std::cout;
			
			if (verboseOut) {
				out	<<"\n==================================="
					<<"\n   HPconvert - CPSP-tool-library"
					<<"\n==================================="
					<<"\n\n"
					<<" conversion of ";
				switch (mode) {
					case A2R : std::cout <<"absolute to relative move string \n\n"; break;
					case A2P : std::cout <<"absolute move string to xyz position file format\n\n"; break;
					case R2A : std::cout <<"relative to absolute move string \n\n"; break;
					case R2P : std::cout <<"relative move string to xyz position file format\n\n"; break;
					case P2A : std::cout <<"xyz positions from file to absolute move string\n\n"; break;
					case P2N : std::cout <<"xyz positions from file to normalized absolute move string\n\n"; break;
					case P2R : std::cout <<"xyz positions from file to relative move string\n\n"; break;
					case A2C : std::cout <<"absolute move string to CML position file format\n\n"; break;
					case A2N : std::cout <<"absolute move string to normalization\n\n"; break;
					case R2N : std::cout <<"relative move string to normalization\n\n"; break;
					case R2C : std::cout <<"relative move string to CML position file format\n\n"; break;
					case A2X : std::cout <<"absolute move string to PDB position file format\n\n"; break;
					default : assertbiu(false,"should never happen (unsupported conversion identifier found)"); // should never happen
				}
			}
			
				// do conversion
			switch (mode) {
				case A2R : abs2rel(structure, &lattice, out); break;
				case A2P : abs2xyz(structure, &lattice, out); break;
				case R2A : rel2abs(structure, &lattice, out); break;
				case R2P : rel2xyz(structure, &lattice, out); break;
				case P2A : xyz2abs(structure, &lattice, false, out); break;
				case P2N : xyz2abs(structure, &lattice, true, out); break;
				case P2R : xyz2rel(structure, &lattice, out); break;
				case A2C : abs2cml(structure, &lattice, out); break;
				case A2N : mov2norm(structure, &lattice, out); break;
				case R2N : mov2norm(structure, &lattice, out); break;
				case R2C : rel2cml(structure, &lattice, out); break;
				case A2X : abs2pdb(structure, &lattice, out, seq); break;
				default : assertbiu(false,"should never happen (unsupported conversion identifier found)"); // should never happen
			}
			
			if (verboseOut) {
				out	<<"\n====================================\n"
					<<std::endl;
			}
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}
		
		if (descr != NULL) delete descr;

		return 0;
	}
	
	void
	HPconvert::abs2rel(	const std::string & absMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
		out	<<lattice->getString( 
								  lattice->absMovesToRelMoves(
									lattice->parseMoveString(absMoves)
								  )
								) 
			<<std::endl;
	}
		
	void
	HPconvert::rel2abs(	const std::string & relMoves, 
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
		out	<<lattice->getString( 
								  lattice->relMovesToAbsMoves(
									lattice->parseMoveString(relMoves)
								  )
								) 
			<<std::endl;
	}
	
	void
	HPconvert::abs2xyz(	const std::string &absMoves,
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
			// point data
		biu::IPointVec points = lattice->absMovesToPoints( lattice->parseMoveString(absMoves));
			// write output
		writeXYZ(points, lattice, out);
			// clear data
		points.clear();
	}
	
	void
	HPconvert::abs2cml(	const std::string &absMoves,
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
			// point data
		biu::IPointVec points = lattice->absMovesToPoints( lattice->parseMoveString(absMoves));
			// write output
		writeCML(points, lattice, out);
			// clear data
		points.clear();
	}
	
	void
	HPconvert::abs2pdb(	const std::string &absMoves,
						const biu::LatticeModel * lattice,
						std::ostream & out,
						const std::string* const seqStr) const
	{
			// point data
		biu::IPointVec ip = lattice->absMovesToPoints( lattice->parseMoveString(absMoves));
		biu::DPointVec dp;
		const double bondLength = 3.6; // Angstroem
		for (size_t i=0; i< ip.size(); i++) {
			dp.push_back(biu::DblPoint(ip[i].getX() * bondLength, ip[i].getY() * bondLength, ip[i].getZ() * bondLength));
		}
			// write output
		std::string paramInfo = "";
		std::vector<std::string> seq;
		if (seqStr == NULL)
			seq.resize(dp.size(),"LEU");
		else {
			for (size_t i=0; i<seqStr->size(); i++)
				if (seqStr->at(i) == 'H')
					seq.push_back("LEU");
				else
					seq.push_back("HIS");
		}
		cpsp::writePDB(paramInfo, seq, dp, *lattice, out, false, *lattice);
			// clear data
		ip.clear();
		dp.clear();
	}
	
	void
	HPconvert::rel2cml(	const std::string &relMoves,
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
		abs2cml(	lattice->getString(
						lattice->relMovesToAbsMoves(
							lattice->parseMoveString(  relMoves )))
					,lattice
					,out);
	}
	
	void
	HPconvert::rel2xyz(	const std::string &relMoves,
						const biu::LatticeModel * lattice,
						std::ostream & out) const
	{
			// point data
		biu::IPointVec points = lattice->relMovesToPoints( lattice->parseMoveString(relMoves));
			// write output
		writeXYZ(points, lattice, out);
			// clear data
		points.clear();
	}
	
	void
	HPconvert::xyz2abs(	const std::string &xyzFile,
						const biu::LatticeModel * lattice,
						const bool normalize,
						std::ostream & out ) const
	{
			// read point data
		biu::IPointVec points;
		std::ifstream fin(xyzFile.c_str());
		if (fin.good()) {
			readXYZ(fin,points);
		} else {
			std::cout <<"\nERROR: cant read XYZ-file '"+xyzFile+"' !\n" <<std::endl;
			return;
		}
			// convert to absolute move string
		if ( points.size() > 1 ) {
			if (normalize) {
				std::cout	<< lattice->getString(
								lattice->getDescriptor()->normalizeSequence(
									lattice->pointsToAbsMoves(points) ) );
			} else {
				std::cout	<< lattice->getString(
								lattice->pointsToAbsMoves(points) );
			}
			std::cout << std::endl;
		} else {
			std::cout <<"\nERROR: only " <<points.size() <<" points read from '"
					<<xyzFile <<"' !\n --> check line delimiters, only '\\n' supported !\n --> maybe use 'dos2unix INFILE' first !" <<std::endl;
		}
			// clear data
		fin.close();
		points.clear();
	}
	
	void
	HPconvert::xyz2rel(	const std::string &xyzFile,
						const biu::LatticeModel * lattice,
						std::ostream & out) const 
	{
			// read point data
		biu::IPointVec points;
		std::ifstream fin(xyzFile.c_str());
		if (fin.good()) {
			readXYZ(fin,points);
		} else {
			std::cout <<"\nERROR: cant read XYZ-file '"+xyzFile+"' !\n" <<std::endl;
			return;
		}
			// convert to absolute move string
		if ( points.size() > 1 ) {
			std::cout	<< lattice->getString( lattice->pointsToRelMoves(points) )
						<< std::endl;
		} else {
			std::cout <<"\nERROR: only " <<points.size() <<" points read from '"
					<<xyzFile <<"' !\n --> check line delimiters, only '\\n' supported !\n --> maybe use 'dos2unix INFILE' first !" <<std::endl;
		}
			// clear data
		fin.close();
		points.clear();
	}
	
	void 
	HPconvert::writeXYZ(	const biu::IPointVec & points,
							const biu::LatticeModel * lattice,
							std::ostream & out ) const
	{
		out	<<"# lattice protein in absolute positions"
			<<"\n# lattice model  : " <<lattice->getDescriptor()->getName()
			<<"\n# created by HPconvert -- CPSP-tools package"
			<<"\n#"
			<<std::endl;
		for (biu::IPointVec::const_iterator p = points.begin(); p!=points.end(); p++) {
			out <<p->getX() <<" " <<p->getY() <<" " <<p->getZ() <<"\n";
		}
		out	<<"#\n# EOF #" <<std::endl;
		out.flush();
	}
	
	void 
	HPconvert::writeCML(	const biu::IPointVec & points,
							const biu::LatticeModel * lattice,
							std::ostream & out ) const
	{
		int length = 3;
		int shift = points.size()*length;
		out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
			<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
			<<"\n  <list>"
			<<"\n    <molecule id=\""<<lattice->getDescriptor()->getName()<<"-lattice-protein\">"
			<<"\n      <atomArray>"; 
		for (biu::IPointVec::size_type i = 0; i < points.size(); i++) {
			out <<"\n        <atom id=\"a"<<i+1<<"\" elementType=\"C\""
				<<" x3=\""<<(shift + (length * points[i].getX()))<<"\""
				<<" y3=\""<<(shift + (length * points[i].getY()))<<"\""
				<<" z3=\""<<(shift + (length * points[i].getZ()))<<"\""
				<<" />";
		}
		out	<<"\n      </atomArray>" 
			<<"\n      <bondArray>"; 
		for (biu::IPointVec::size_type i = 1; i < points.size(); i++) {
			out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\"a"<<i<<" a"<<i+1<<"\" order=\"S\"/>";
		}
		out	<<"\n      </bondArray>"
			<<"\n    </molecule>"
			<<"\n  </list>"
			<<"\n</list>"
			<<std::endl;
		out.flush();
	}


	void 
	HPconvert::mov2norm(	const std::string & moves, 
							const biu::LatticeModel * lattice,
							std::ostream & out) const
	{
		const biu::LatticeDescriptor* const ld = lattice->getDescriptor();
		out	<<	ld->getString(
					ld->normalizeSequence(
						ld->getSequence(moves)
					)
				)
			<< std::endl;
	}
	
	int
	HPconvert::readXYZ(	std::istream & in, biu::IPointVec & points) const
	{
			// clear data
		points.clear();
			// one line of stream and the allowed alphabet for data lines
		std::string line = "", lineAlph = "-0123456789 \t";
			// positions for reading
		int x=0, y=0, z=0;
		while (std::getline(in, line, '\n')) {
				// check if empty line or comment 
			if (line.empty() || (line[0]=='#'))
				continue;
				// check if right alphabet
			else if (line.find_first_not_of(lineAlph,0) == line.npos) {
				std::istringstream iss(line); // feed string stream for conversion
				z = INT_MIN; // init for failure test
				iss >>x >>y >>z; // read absolute lattice position
				if (z != INT_MIN) { // read was successful
					points.push_back(biu::IntPoint(x,y,z)); // write new position
				}
			}
		}
		return 0;
	}

		
	void
	HPconvert::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const
	{
		allowedArgs.push_back(biu::COption(	
								"m", false, biu::COption::STRING, 
								"conversion mode ( [ar]2[arpc], p2[ar]) - see information below"));
		allowedArgs.push_back(biu::COption(	
								"str", false, biu::COption::STRING, 
								"the structure to convert"));
		allowedArgs.push_back(biu::COption(	
								"seq", true, biu::COption::STRING, 
								"the sequence of the molecule"));
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

		infoText = std::string("HPconvert allows the conversion between different")
			+ std::string(" lattice structure representations and prints them")
			+ std::string(" to STDOUT.")
			+ std::string("\nThe supported modes are (m-parameter):\n")
			+ std::string("\n a2n : normalizing an absolute move string")
			+ std::string("\n a2r : absolute to relative move string")
			+ std::string("\n a2p : absolute move string to xyz-point file")
			+ std::string("\n a2c : absolute move string to CML file for viewing")
			+ std::string("\n a2x : absolute move string to PDB file for viewing")
			+ std::string("\n p2a : xyz-point file to absolute move string")
			+ std::string("\n p2n : xyz-point file to normalized absolute move string")
			+ std::string("\n p2r : xyz-point file to relative move string")
			+ std::string("\n r2a : relative to absolute move string")
			+ std::string("\n r2n : normalizing a relative move string")
			+ std::string("\n r2p : relative move string to xyz-point file")
			+ std::string("\n r2c : relative move string to CML file for viewing")
			+ std::string("\n\nIf the input is NO move string, a file name")
			+ std::string(" has to be given to the '-str=' parameter instead.")
			+ std::string("\n\nThe structure is NOT validated")
			+ std::string(" (if connected and selfavoiding). If it is invalid")
			+ std::string(" normal execution cant be guaranteed.")
			+ std::string("\n\nThe moves are encoded using:")
			+ std::string("\n F/B : +-x\n L/R : +-y\n U/D : +-z")
			+ std::string("\n\nFor an XYZ-file format example run 'HPconvert")
			+ std::string(" -m=a2p -s=FLUBBURD'. A '#' marks a comment line.");
	}
	



	int main(int argc, char** argv) {
		HPconvert hpc;
			// run conversion
		return hpc.convert(argc, argv);
	}
