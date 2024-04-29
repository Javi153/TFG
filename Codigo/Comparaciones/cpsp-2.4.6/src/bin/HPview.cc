#include "HPview.hh"
#include "version.hh"

#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>

#include "cpsp/Exception.hh"

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <ctime>


	HPview::HPview()
	  :	JMOL_SCRIPT("select all; color bonds grey;")
	{
	}
	
	HPview::~HPview()
	{
	}


	void 
	HPview::points2cml(	const std::string & seq,
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
			<<"\n    <!-- encoding : H = Cl-atom, P = C-atom  -->"
			<<"\n    <!-- lattice  : " <<lattice->getDescriptor()->getName() <<" -->"
			<<"\n      <atomArray>"; 
		for (biu::IPointVec::size_type i = 0; i < points.size(); i++) {
			out <<"\n        <atom id=\""<<i+1<<"_"<<seq[i]<<"\""
				<<" elementType=\"C"<<(seq[i]=='H'?"l":"")<<"\""
				<<" x3=\""<<(shift + (distance * points[i].getX()))<<"\""
				<<" y3=\""<<(shift + (distance * points[i].getY()))<<"\""
				<<" z3=\""<<(shift + (distance * points[i].getZ()))<<"\""
				<<" />";
		}
		out	<<"\n      </atomArray>" 
			<<"\n      <bondArray>"; 
		for (biu::IPointVec::size_type i = 1; i < points.size(); i++) {
			out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_"<<seq[i-1]<<" "<<i+1<<"_"<<seq[i]<<"\" order=\"S\"/>";
		}
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
	HPview::convert(int argc, char** argv)
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
			LatticeModel lattice(descr);

				// get and validate sequence
			std::string seq = opts.getStrVal("seq");
			if ( seq.size() == 0 ) 
				throw cpsp::Exception("Error in arguments: no sequence given",-2);
			biu::Alphabet alph("HP",1);
			if (!alph.isAlphabetString(seq))
				throw cpsp::Exception("Error in arguments: sequence '"+seq+"' is no HP-sequence",-2);
			
				// get and validate structure
			biu::IPointVec points;
			if (opts.argExist("abs") && opts.argExist("rel")) 
				throw cpsp::Exception("Error in arguments: relative and absolute structure given",-2);
			if (!opts.argExist("abs") && !opts.argExist("rel")) 
				throw cpsp::Exception("Error in arguments: no structure given",-2);
			if (opts.argExist("abs")) {
				std::string str = opts.getStrVal("abs");
				if (!lattice.getDescriptor()->getAlphabet()->isAlphabetString(str))
					throw  cpsp::Exception("Error in arguments: structure '"+str+"'is no valid move string",-2);
				points = lattice.absMovesToPoints(lattice.getDescriptor()->getSequence(str));
			} else {
				std::string str = opts.getStrVal("rel");
				if (!lattice.getDescriptor()->getAlphabet()->isAlphabetString(str))
					throw  cpsp::Exception("Error in arguments: structure '"+str+"'is no valid move string",-2);
				points = lattice.relMovesToPoints(lattice.getDescriptor()->getSequence(str));
			}
			if (points.size() == 0) 
				throw  cpsp::Exception("Error: something wrong in structure conversion",-3);
			if (points.size() != seq.size()) 
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
	HPview::initAllowedArguments(biu::OptionMap & allowedArgs,
											std::string &infoText ) const
	{
		allowedArgs.push_back(biu::COption(	
								"seq", false, biu::COption::STRING, 
								"the HP sequence"));
		allowedArgs.push_back(biu::COption(	
								"abs", true, biu::COption::STRING, 
								"the structure in absolute moves"));
		allowedArgs.push_back(biu::COption(	
								"rel", true, biu::COption::STRING, 
								"the structure in relative moves"));
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
								"home dir of Jmol viewer (or use $JMOL_HOME)", "$JMOL_HOME"));
		allowedArgs.push_back(biu::COption(	
								"v", true, biu::COption::BOOL, 
								"verbose output"));
		allowedArgs.push_back(biu::COption(	
								"help", true, biu::COption::BOOL, 
								"program parameters and help"));
		allowedArgs.push_back(biu::COption(	
								"version", true, biu::COption::BOOL, 
								"version information of this program"));

		infoText = std::string("HPview creates an output in CML-file format")
			+ std::string(" of a sequence/structure that can be viewed with")	
			+ std::string(" molecule viewers like Chime or Jmol.")
			+ std::string("\n\nThe structure is NOT validated")
			+ std::string(" (if connected and selfavoiding). If it is invalid")
			+ std::string(" normal execution cant be guaranteed.")
			+ std::string("\n\nThe moves are encoded using:")
			+ std::string("\n F/B : +-x\n L/R : +-y\n U/D : +-z")
			+ std::string("\n\nCurrently supported viewers are:")
			+ std::string("\n - Jmol (http://jmol.sourceforge.net)")
			+ std::string("");
	}
	



	int main(int argc, char** argv) {
			// create object
		HPview hpv;
			// run conversion
		return hpv.convert(argc, argv);
	}
