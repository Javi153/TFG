#ifndef GC_THREADINGSUPERSPACE_HH_
#define GC_THREADINGSUPERSPACE_HH_


#include <gecode/kernel.hh>



namespace cpsp
{
  namespace gecode
  {

	class SuperSpace : public Gecode::Space {
		
	public:
		
		SuperSpace () : Gecode::Space() {}
		virtual ~SuperSpace () {}
		SuperSpace (bool share, Space &s) : Gecode::Space(share,s) {}
		
		 /*! returns the rank of the variable with the given index
		  * 
		  * @param index of the variable of interest
		  * @return its rank
		  */ 
		virtual
		unsigned int
		getRank(const int index) const = 0;
	};
  
  }
}
#endif /*GC_THREADINGSUPERSPACE_HH_*/
