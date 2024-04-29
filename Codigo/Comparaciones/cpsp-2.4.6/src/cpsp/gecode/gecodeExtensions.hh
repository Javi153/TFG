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

#ifndef GECODEEXTENSIONS_HH_
#define GECODEEXTENSIONS_HH_


#include <gecode/int/count.hh>

namespace Gecode {

  void
  count(Space* home, const IntVarArgs& xa, const IntArgs& ya,
	IntRelType r, int z, IntConLevel icl);
	
  void
  atmost(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl=ICL_DEF);
  
  void
  atleast(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl=ICL_DEF);
    
  void
  exactly(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl=ICL_DEF); 


  inline
  void
  count(Space* home, const IntVarArgs& xa, const IntArgs& ya,
	IntRelType r, int z, IntConLevel icl) {
    using namespace Int;
    if (home->failed()) return;
    ViewArray<OffsetView> x(home,xa.size());
    for (int i = ya.size(); i--; )
      if ((-ya[i] < Limits::Int::int_min) || (-ya[i] > Limits::Int::int_max))
        throw NumericalOverflow("Int::count");
      else
        x[i].init(xa[i],-ya[i]);
    ConstIntView yv(0);
    switch (r) {
    case IRT_EQ:
      GECODE_ES_FAIL(home,(Count::EqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z)));
      break;
    case IRT_NQ:
      GECODE_ES_FAIL(home,(Count::NqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z)));
      break;
    case IRT_LE:
      GECODE_ES_FAIL(home,(Count::LqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z-1)));
      break;
    case IRT_LQ:
      GECODE_ES_FAIL(home,(Count::LqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z)));
      break;
    case IRT_GR:
      GECODE_ES_FAIL(home,(Count::GqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z+1)));
      break;
    case IRT_GQ:
      GECODE_ES_FAIL(home,(Count::GqInt<OffsetView,ConstIntView>
      			   ::post(home,x,yv,z)));
      break;
    default:
      throw UnknownRelation("Int::count");
    }
  }
  
  inline
  void
  atmost(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl) 
  {
    count(home,x,ya,IRT_LQ,m,icl);
  }
  
  inline
  void
  atleast(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl) 
  {
    count(home,x,ya,IRT_GQ,m,icl);
  }

  inline  
  void
  exactly(Space* home, const IntVarArgs& x, const IntArgs& ya, int m, 
    IntConLevel icl) 
  {
    count(home,x,ya,IRT_EQ,m,icl);
  }
  
} // namespace Gecode




#endif /*GECODEEXTENSIONS_HH_*/
