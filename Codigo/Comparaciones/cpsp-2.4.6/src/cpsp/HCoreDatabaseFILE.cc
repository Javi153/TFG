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

#include "cpsp/HCoreDatabaseFILE.hh"
#include "cpsp/Exception.hh"

#include <algorithm>	// fuer <sort> template
#include <functional>	// fuer <greater> template

#ifdef WIN32			// fuer systemspezifische verzeichnis-inhaltsauflistung
#	include <io.h>
#else
#	include <dirent.h>
#endif

#include <biu/util/Util_String.h>

#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>

namespace cpsp
{

	//////////////////////////////////////////////////////////////////////
	// Konstanten
	//////////////////////////////////////////////////////////////////////
	
	const std::string ROOT_CUB = std::string("Cubic/");
	const std::string ROOT_FCC = std::string("Fcc/");
	
	const std::string PREF_CUB = std::string("cub.");
	const std::string PREF_FCC = std::string("fcc.");
	
	const std::string POSTFIX = std::string(".txt");
	

	HCoreDatabaseFILE::HCoreDatabaseFILE(const std::string& dbRootPath) :	
		coreSize(0),  minHH(0), maxHH(UINT_MAX-2), 
			// check for '/' char at the end of the path
		DBROOTPATH(dbRootPath+(dbRootPath[dbRootPath.size()-1] == '/'?"":"/")),
		actFile(NULL), actFileNamePref(""), actFileHHcontInd(-1)
	{
		actHHcontacts.clear();
	}
	
	HCoreDatabaseFILE::HCoreDatabaseFILE(const HCoreDatabaseFILE& toCopy) :
		coreSize(0),  minHH(toCopy.minHH), maxHH(toCopy.maxHH),
		DBROOTPATH(toCopy.DBROOTPATH), actFile(NULL), actFileNamePref(""),
		actFileHHcontInd(-1)
	{
	}
	
	HCoreDatabaseFILE::~HCoreDatabaseFILE()
	{
		if (actFile != NULL) {
			actFile->close();
			delete actFile;
		}
		actFile = NULL;
	}


	//////////////////////////////////////////////////////////////////////
	// inits database access
	//////////////////////////////////////////////////////////////////////
	
	bool HCoreDatabaseFILE::initCoreAccess(const biu::LatticeDescriptor& latDescr, const unsigned int size,
					const unsigned int minHH_, const unsigned int maxHH_)
	{

	// (1) init data
		this->coreSize = size;
		this->minHH = minHH_;
		this->maxHH = maxHH_;
	
		if (actFile != NULL) { delete actFile; actFile = NULL; }
		actFileHHcontInd = -1;
		actFileNamePref = "";
	
	// (2) check available energy files and set <actEnergies>
	
		// create path name
		std::string actDBdirName = DBROOTPATH, prefix = "";
		// lattice root + size dirName
		if (latDescr == biu::LatticeDescriptorCUB()) {
			actDBdirName += ROOT_CUB; prefix += PREF_CUB;
		} else if (latDescr == biu::LatticeDescriptorFCC()) {
			actDBdirName += ROOT_FCC; prefix += PREF_FCC;
		} else {
			std::cerr <<"\n\tError: No core files for lattice type '"+latDescr.getName()+"' available!\n\n";
			return false;
		}
		actDBdirName += (biu::util::Util_String::int2str(size))+"/";	// add folder name for core size
		prefix += (biu::util::Util_String::int2str(size)) + ".";	// prefix of wanted core data file
	
		actHHcontacts.clear();					// clear old data
	
		std::string filename;					// current file name
	
		// open folder and go through entries
	#ifdef WIN32
		long hFile;
		struct _finddata_t c_file;				// folder entry
		hFile = _findfirst ((actDBdirName+prefix+"*"+POSTFIX).c_str(), &c_file);	// open
		if (hFile != -1L) {						// folder exists
	
			do {
				filename = c_file.name;
	#else	// linux
		DIR *actDBdir = opendir(actDBdirName.c_str());	// open folder
		if(actDBdir != NULL)					// folder exists
		{
	
			struct dirent *dirEntry;			// folder entry
	
			while((dirEntry=readdir(actDBdir)))
			{
				filename = dirEntry->d_name;
	#endif
				// get energy out of filename
				if (filename.substr(0,prefix.size()).compare(prefix) == 0) { // correct file
					filename = filename.erase(0,prefix.size());	// cut prefix
					actHHcontacts.push_back(biu::util::Util_String::str2int(filename.erase(filename.find(POSTFIX)))); // cut postfix and convert name
				}
			}
	#ifdef WIN32
			while( _findnext( hFile, &c_file ) == 0 );
	#endif
		} else {
			std::cerr	<<"\n\tError: No core file for size '"<<size <<"' and\n\t       lattice type '" 
						<<latDescr.getName() <<"' available!\n\t==> no directory '"<<actDBdirName<<"'\n\n";
			return false;
			cpsp::Exception::errorCount++;
		}
	
		// close folder
	#ifdef WIN32
		_findclose( hFile );
	#else
		closedir(actDBdir);
	#endif
	
	// (3) sort decreasing
	
		std::sort(actHHcontacts.begin(), actHHcontacts.end()  , std::greater<int>());
	
	// (4) set <actFileEnergyInd> and open file if available
	
			 // go through file until maxenergy found
		unsigned int i=0;
		while ( i<actHHcontacts.size() && actHHcontacts[i]>0 &&((unsigned int)actHHcontacts[i]>maxHH || (unsigned int)actHHcontacts[i]<minHH)) {
//			std::cerr <<" actHH[i] = " <<(int)actHHcontacts[i] <<" size = " <<size <<" minHH = " <<minHH <<" maxHH = " <<maxHH <<"\n";
			i++;
		}
	
		if (i >= actHHcontacts.size()) {// dann KEINE entsprechende datei vorhanden
//			std::cerr <<"no files available\n";
			return false;	
		}
	
		actFileHHcontInd = i;
		actFileNamePref = actDBdirName+prefix;
		if (actFile != NULL) { 
			delete actFile; actFile = NULL;
		}
		actFile = new std::ifstream((actFileNamePref+(biu::util::Util_String::int2str(actHHcontacts[actFileHHcontInd]))+POSTFIX).c_str());
		
		return true;
	} // initCoreAccess()

	bool HCoreDatabaseFILE::getNextCore(HCore& toFill) {
	
		toFill.reset(INT_MAX, INT_MAX);		// invalid dummy core
		if (actFile == NULL) // read only core if db initialised and file already open
			return false;
	
		std::string dataline = "";
		int i=0;
		int p[3];

		do {	////////////////////////////////////
				// READ NEXT CORE
				////////////////////////////////////
			while (std::getline(*actFile, dataline, '#')) { // datablocks seperated by '#'
				if (!dataline.empty() && (dataline.substr(0,8).compare("COREDATA")==0)){	// after keyword 'COREDATA' follows the core data
					try {
						std::istringstream iss(dataline.substr(9)); // remaining string into stringstream
						toFill.reset(actHHcontacts[actFileHHcontInd], coreSize); // new core
						for (i=0; i<(int)coreSize; i++) {
							p[2] = INT_MIN; // init for failure analysis
							iss >>p[0] >>p[1] >>p[2]; // read core positions
							if (p[2] == INT_MIN) // than no complete point coordinates available
								throw ExcEOF();
							toFill.addPoint(HCore::CorePoint(p[0], p[1], p[2])); // add position to core
						}
					} catch (ExcEOF) // for avoiding invalid data in <cores> 
					{
						toFill.reset(INT_MAX,INT_MAX);
						return false;
					}
					return toFill.isValid();
				}
			}
			actFile->close();		// close if everything read
			delete actFile;			// clean up
			actFile = NULL;
	
			// check if other cores with lower energy to read
			if ((++actFileHHcontInd) < (int)actHHcontacts.size()  && actHHcontacts[actFileHHcontInd]>=(int)minHH) {
				// if so open file
				actFile = new std::ifstream((actFileNamePref+(biu::util::Util_String::int2str(actHHcontacts[actFileHHcontInd]))+POSTFIX).c_str());
			}
		} while (actFile != NULL);
	
		return false;
	}


} // namespace cpsp
