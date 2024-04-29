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


#include <iostream>

#include "cpsp/HCoreDatabase.hh"
#include <biu/LatticeDescriptorCUB.hh>
void testClassHCoreDatabase(cpsp::HCoreDatabase& db) {
	using namespace cpsp;
	
	biu::LatticeDescriptorCUB ld;
	std::cout <<"LatticeDescriptor = ld = "<<ld.getName() <<std::endl;
	
	db.initCoreAccess(ld,6,5,7);						// initCoreAccess()
	std::cout	<< "db.initCoreAccess(ld,6,5,7) = ("
				<< db.getActCoreSize() <<","			// getActCoreSize()
				<< db.getActMinHHcontacts() <<","		// getActMinHHcontacts()
				<< db.getActMaxHHcontacts() <<")"		// getActMaxHHcontacts()
				<< std::endl;
				
		HCore core;	
		unsigned int count = 0, invalid = 0;
	
	std::cout << "number of valid cores with <= 6 contacts (21) = ";
	while (db.getNextCore(core)) {
		count++;
		if (!core.isValid())							// kerne muessen gueltig sein
			invalid++;
		if (core.getContacts() == 6)					// setActMinHHcontacts()
			db.setActMinHHcontacts(core.getContacts());
	}
	std::cout <<(count-invalid) <<std::endl;			// anzahl von getNextCore()
	
}

#include "cpsp/HCoreDatabaseFILE.hh"
int main(int argc, char** argv) {
	using namespace cpsp;
	
	HCoreDatabase* db = NULL;
	
	std::string dbPath = "../CoreDB/";

	std::cout	<< "file based database\n" <<std::endl;
	
	db = new HCoreDatabaseFILE(dbPath);

	testClassHCoreDatabase(*db);

	delete db;

	return 0;
}
