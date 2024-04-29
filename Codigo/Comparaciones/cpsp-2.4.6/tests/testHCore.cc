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

#include "cpsp/HCore.hh"
#include <biu/LatticeFrame.hh>
#include <biu/LatticeDescriptorCUB.hh>


int main(int argc, char** argv) {
	
	using namespace cpsp;
	using namespace biu;
	
	HCore core(4,4);
	std::cout	<< "HCore(4,4)           = ("
				<< core.getSize() <<","
				<< core.getContacts() <<")" <<std::endl;
	std::cout	<< "currently invalid    = " <<(core.isValid()?"true":"false")
				<<std::endl;
	
		core.addPoint(IntPoint(0,0,0));
		core.addPoint(IntPoint(1,0,0));
		core.addPoint(IntPoint(1,1,0));
		core.addPoint(IntPoint(1,2,0));
	std::cout	<< "valid after filling  = " <<(core.isValid()?"true":"false")
				<<std::endl;
	
	IPointSet points = core.getPoints();		// getPoints()
	std::cout	<<"points should be     = (0,0,0) (1,0,0) (1,1,0) (1,2,0)\n"
				<<"       and are       =" ;
	for (biu::IPointSet::const_iterator it = points.begin(); it != points.end(); it++) {
		std::cout <<" (" <<*it <<")";
	}
	std::cout <<std::endl;
	
	std::cout	<<"center (0,1,0)       = " <<core.getCenter() <<std::endl;

	LatticeDescriptorCUB ld;
	LatticeFrame *latFrame = new LatticeFrame(&ld, 10);
	std::cout	<<"LatticeFrame(cub,10) = (" 
				<<latFrame->getDescriptor()->getName() <<","
				<<latFrame->getFrameSize() <<")" 
				<<std::endl
				<<"frame center (5,5,5) = (" <<latFrame->getCenter() <<")" <<std::endl;
	
	biu::LatticeFrame::index_set indSet;
	biu::LatticeFrame::index_set::const_iterator iit;

	indSet = core.getIndexedHull(latFrame, 0);	
	std::cout <<"shiftet core positions (to frame center) = \n\t";
	for (iit = indSet.begin(); iit != indSet.end(); iit++) {
		std::cout <<" (" <<latFrame->getPoint(*iit) <<")";
	}
	std::cout <<std::endl;

	indSet = core.getIndexedPSinglets(latFrame);
	std::cout	<<"P-Singlets (5,5,5)   = ";
	for (iit = indSet.begin(); iit != indSet.end(); iit++) {
		std::cout <<" (" <<latFrame->getPoint(*iit) <<")";
	}
	std::cout <<std::endl;

	std::cout	<<"H-Singlet number (4) = "
				<<core.getIndexedHSinglets(latFrame).size()
				<<std::endl
				<<"max Dimension    (3) = "
				<<core.getMaxDimension()
				<<std::endl
				<<"min Even/Odd     (2) = "
				<<core.getMinEvenOdd()
				<<std::endl;
	
	std::cout	<<"reset(0,1) + (1,1,1) = ";
	core.reset(0,1);
	core.addPoint(IntPoint(1,1,1));
	std::cout	<<(core.isValid()?"ok":"failed") <<std::endl;
	
	std::cout	<<"hull 2 elems    (24) = "
				<<core.getIndexedHull(latFrame, 2).size()
				<<std::endl;
	
	
	delete latFrame;
	
	return 0;

}
