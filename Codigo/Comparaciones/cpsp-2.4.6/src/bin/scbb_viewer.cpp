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

#include "scbb_viewer.h"
#include "version.hh"

#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>

#include <cpsp/Exception.hh>

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <ctime>

	using namespace biu;

	biu::IPointVec SCBBLatticeModelViewer::absMovesToPoints( const biu::MoveSequence& absMoveSeq ) const{
	    assertbiu(latDescriptor != NULL,  "LatticeModel has no LatticeDescriptor");
	    IPointVec points(absMoveSeq.size()+1);  // rueckgabe vector
	    points[0] = biu::IntPoint(0,0,0); // init mit lattice center
	    for (MoveSequence::size_type i=0; i< absMoveSeq.size(); i++){
		// from a backbone towards the side chain
		if(i%2==0)
		   points[i+1] = applyAbsMove( points[i], absMoveSeq[i] );
		// from a backbone towards another backbone
		else
		   points[i+1] = applyAbsMove( points[i-1], absMoveSeq[i] );
		
       	    }
            return points;
	}
	

	SCBBLatticeModelViewer::~SCBBLatticeModelViewer(){
	}
	
	SCBBViewer::SCBBViewer()
	  : JMOL_SCRIPT("select all; color bonds grey;")
	{
	}
	SCBBViewer::~SCBBViewer()
	{
	}

	void SCBBViewer::deleteBrackets(std::string& absMoves){
	   for(size_t i=0;i<absMoves.size();){
		if(absMoves[i]=='(' || absMoves[i]==')' ){
		   absMoves.erase(i,2);
		}
		i=i+2;
	   }
	}

	void 
	SCBBViewer::points2cml(	const std::string & seq,
						const biu::IPointVec & points, 
						const biu::LatticeModel * lattice,
						std::ostream & out, 
						int distance) const 
	{
		int shift = points.size()*distance;	// shift output to display with positive dimensions only
		
		out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
			<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
			<<"\n<!-- "<<JMOL_SCRIPT<<" -->"
			<<"\n\n<!-- created with CPSP-package : http://www.bioinf.uni-freiburg.de/sw/cpsp/ -->"
			<<"\n"
			<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
			<<"\n  <list>"
			<<"\n    <molecule id=\"lattice-protein\">"
			<<"\n    <!-- model    : HP-model -->"
			<<"\n    <!-- sequence : " <<seq <<" -->"
			<<"\n    <!-- lattice  : " <<lattice->getDescriptor()->getName() <<" -->"
			<<"\n      <atomArray>"; 
		int k=0;
		for (biu::IPointVec::size_type i = 0; i < points.size(); i++) {
			//Here side chain
			if(i%2==1){
			out <<"\n        <atom id=\""<<k+1<<"_"<<seq[k]<<"\""
				<<" elementType=\"C"<<(seq[k]=='H'?"l":"")<<"\""
				<<" x3=\""<<(shift + (distance * points[i].getX()))<<"\""
				<<" y3=\""<<(shift + (distance * points[i].getY()))<<"\""
				<<" z3=\""<<(shift + (distance * points[i].getZ()))<<"\""
				<<" />";
				k++;
			}
			//For backbone
			else{
			   out <<"\n        <atom id=\""<<k+1<<"_b"<<seq[k]<<"\""
				<<" elementType=\"C"<<"2"<<"\""
				<<" x3=\""<<(shift + (distance * points[i].getX()))<<"\""
				<<" y3=\""<<(shift + (distance * points[i].getY()))<<"\""
				<<" z3=\""<<(shift + (distance * points[i].getZ()))<<"\""
				<<" />";
			}
		}
		out	<<"\n      </atomArray>" 
			<<"\n      <bondArray>"; 
		for (biu::IPointVec::size_type i = 1; i < seq.size(); i++) {
			//backbones connections
			out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_b"<<seq[i-1]<<" "<<i+1<<"_b"<<seq[i]<<"\" order=\"S\"/>";
			//backbones-side chains connections
			out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_b"<<seq[i-1]<<" "<<i<<"_"<<seq[i-1]<<"\" order=\"S\"/>";
			
		}
	      //last backbone-side chain connection
	        out     <<"\n        <bond id=\"b"<<seq.size()<<"\" atomRefs2=\""<<seq.size()<<"_b"<<seq[seq.size()-1]<<" "<<seq.size()<<"_"<<seq[seq.size()-1]<<"\" order=\"S\"/>";
		out	<<"\n      </bondArray>"
			<<"\n    </molecule>"
			<<"\n  </list>"
			<<"\n</list>"
			<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
			<<"\n<!-- "<<JMOL_SCRIPT<<" -->"
			<<std::endl;
		out.flush();	
	}
	
	int
	SCBBViewer::convert(int argc, char** argv)
	{
		using namespace biu;
		using namespace biu;
		
			// the lattice descriptor
		LatticeDescriptor* descr = NULL;
		std::ofstream *ofile = NULL;
		std::ofstream *scriptfile = NULL;
		
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
				descr = new LatticeDescriptorSQR();
			else if (latStr.compare("CUB") == 0)
				descr = new LatticeDescriptorCUB();
			else if (latStr.compare("FCC") == 0)
				descr = new LatticeDescriptorFCC();
			else
				throw cpsp::Exception("Unknown lattice type '"+latStr+"'",-2);
			//The lattice model of the sidechain model
			SCBBLatticeModelViewer lattice(descr);

				// get and validate structure
			biu::IPointVec points;
			if (opts.argExist("abs")) {
				std::string str = opts.getStrVal("abs");
				{
					int openBracket = 0;
					std::string tmp;
					for (size_t i=0; i<str.size(); i++) {
						switch(str.at(i)) {
						case '(' : 
							if (openBracket != 0) 
								throw cpsp::Exception("Error in arguments: abs sequence has two successive '('",-2);
							openBracket++;
							break;
						case ')' :
							if (openBracket != 1) 
								throw cpsp::Exception("Error in arguments: abs sequence has two successive ')'",-2);
							openBracket--;
							break;
						default :
							tmp.push_back(str.at(i));
						}
					}
					if (openBracket != 0) 
						throw cpsp::Exception("Error in arguments: in abs sequence the final ')' is missing",-2);
					str = tmp;
				}
				if (!lattice.getDescriptor()->getAlphabet()->isAlphabetString(str))
					throw  cpsp::Exception("Error in arguments: structure '"+str+"'is no valid move string",-2);
				points = lattice.absMovesToPoints(lattice.getDescriptor()->getSequence(str));
			} else {
				throw  cpsp::Exception("Error: no structure given",-2);
			}
			if (points.size() == 0) 
				throw  cpsp::Exception("Error: something wrong in structure conversion",-3);
			
				// get and validate sequence
			std::string seq;
			if (opts.argExist("seq")) {
				seq = opts.getStrVal("seq");
				if ( seq.size() == 0 ) 
					throw cpsp::Exception("Error in arguments: no sequence given",-2);
				biu::Alphabet alph("HP",1);
				if (!alph.isAlphabetString(seq))
					throw cpsp::Exception("Error in arguments: sequence '"+seq+"' is no HP-sequence",-2);
			} else {
				// generate dummy sequence in case none is given 
				seq = std::string(points.size()/2, 'P');
			}
			if (points.size() != 2*seq.size()) 
				throw  cpsp::Exception("Error in arguments: sequence and structure differ in size",-2);
			
			bool printToFile = false;
			std::string outFile = "";
			
				// output file
			if (opts.argExist("o")) {
				outFile = opts.getStrVal("o");
				printToFile = true;
			}
			
				// viewing types 
			VIEW_TYPE viewType = NONE;
			std::string homeDir = "";
			std::string script = "";
			
				// JMol settings
			if (opts.getBoolVal("jmol")) {
				viewType = JMOL;
				script = JMOL_SCRIPT;
				if (opts.argExist("jmolHome")) {
					homeDir = opts.getStrVal("jmolHome");
				} else {
					char * tmp = getenv("JMOL_HOME");
					if (tmp != NULL)
						homeDir = std::string(tmp); 
				}
			}
			
				// temporary data and script file creation
			if (viewType != NONE) {
				if (!printToFile) {
						// temporary output file name
					std::ostringstream  oss("");
					oss <<"HPview.tmp." <<time(NULL) <<".cml";
					outFile = oss.str();
					printToFile = true;
				}
					// temporary script file 
				scriptfile = new std::ofstream(std::string(outFile+".script").c_str());
				if (!(*scriptfile) || !scriptfile->is_open())
					throw  cpsp::Exception("Error: cant open temporary script file '"+outFile+".script'",-3);
				switch (viewType) {
				case JMOL:
				  {
				  	(*scriptfile) <<JMOL_SCRIPT <<std::endl;
				  	break;
				  }
				case NONE:
				default: break;
				}
				scriptfile->close();
			}
			
				// outstream	
			std::ostream * out = &(std::cout);
				// file stream
			if (printToFile) {
				ofile = new std::ofstream(outFile.c_str());
				if (!(*ofile) || !ofile->is_open())
					throw  cpsp::Exception("Error: cant open output file '"+outFile+"'",-3);
				out = ofile;
			}			
				// atom distance in output
			int atomDist = 3;
			if (opts.argExist("d")) {
				atomDist = opts.getIntVal("d");
				if (atomDist < 1)
					throw  cpsp::Exception("Error in arguments: atom distance is smaller than 1",-2);
			}
			
				// check for verbose output
			bool verboseOut = opts.argExist("v");
			
			if (verboseOut) {
				std::cout	<<"\n==================================="
							<<"\n   HPview - CPSP-tool-library"
							<<"\n==================================="
							<<"\n\n";
			}
			
			std::cout <<(verboseOut?"\n - write CML file":""); std::cout.flush();
			points2cml( seq, points, &lattice, *out, atomDist);
			
				// run viewer if selected
			std::cout <<(verboseOut?"\n - run viewer":""); std::cout.flush();
			switch (viewType) {
				case JMOL : 
				  {
					system(std::string("java -Xmx512m -jar "+homeDir+"/Jmol.jar -s "+outFile+".script "+outFile).c_str());
					system(std::string("rm -f "+outFile+".script").c_str());
					break;
				  }
				case NONE :
				default : break;
			}
			
			std::cout <<(verboseOut?"\n - delete temporary files":""); std::cout.flush();
			if (!opts.argExist("o")) {
					// remove temporary output file
				system(std::string("rm -f "+outFile).c_str());
			}
						
			if (verboseOut) {
				std::cout	<<"\n====================================\n"
							<<std::endl;
			}
			std::cout.flush();
			
		// catch errors in user parameters etc.
		} catch (cpsp::Exception& ex) {
			std::cerr <<"\n\t"+ex.toString()+"\n\n" <<std::endl;
			return ex.retValue;
		}
		
		if (scriptfile != NULL) {if (scriptfile->is_open()) scriptfile->close(); delete scriptfile;}
		if (ofile != NULL) {ofile->close(); delete ofile;}
		if (descr != NULL) delete descr;

		return 0;
	}

			
	void
	SCBBViewer::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const
	{
		allowedArgs.push_back(biu::COption(	
								"seq", true, biu::COption::STRING, 
								"the HP sequence"));
		allowedArgs.push_back(biu::COption(	
								"abs", false, biu::COption::STRING, 
								"the structure in absolute moves"));
		allowedArgs.push_back(biu::COption(	
								"o", true, biu::COption::STRING, 
								"name of the output file"));
		allowedArgs.push_back(biu::COption(	
								"lat", true, biu::COption::STRING, 
								"lattice (SQR, CUB, FCC)", "CUB"));
		allowedArgs.push_back(biu::COption(	
								"d", true, biu::COption::INT, 
								"atom distance (>=1)", "3"));
		allowedArgs.push_back(biu::COption(	
								"jmol", true, biu::COption::BOOL, 
								"open and show in Jmol viewer"));
		allowedArgs.push_back(biu::COption(	
								"jmolHome", true, biu::COption::STRING, 
								"home dir of Jmol viewer (or use $JMOL_HOME)"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("SCH_HPview creates an output in CML-file format"
			" of a sequence/structure that can be viewed with"
			" molecule viewers like Chime or Jmol."
			"\n\nThe structure is NOT validated"
			" (if connected and selfavoiding). If it is invalid"
			" normal execution cant be guaranteed."
			"\n\nA side chain move string is encoded like (U)R(U)L(B)"
			" where in brackets the side chain positioning via moves is given"
			" the backbone moves inbetween."
			"\n\nThe moves are encoded using:"
			"\n F/B : +-x\n L/R : +-y\n U/D : +-z"
			"\n\nCurrently supported viewers are:"
			"\n - Jmol (http://jmol.sourceforge.net)"
			"");
	}
	



	int main(int argc, char** argv) {
			// create object
		SCBBViewer hpv;
			// run conversion
		return hpv.convert(argc, argv);
	}

