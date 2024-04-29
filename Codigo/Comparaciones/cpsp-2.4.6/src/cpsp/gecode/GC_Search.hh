#ifndef GC_SEARCH_HH_
#define GC_SEARCH_HH_


namespace cpsp {
  namespace gecode {
	
		//! an IntView that provides its index in the original viewarray
	class IdxIntView : public Gecode::Int::IntView {
	  protected:
	    int index;
	  public:
	    IdxIntView() : index(-1) {}
	    IdxIntView(const IdxIntView& toCopy) : IntView(toCopy), index(toCopy.index) {}
	    IdxIntView(const Gecode::IntVar& x, int idx): IntView(x), index(idx) {}
	    void update(Gecode::Space* home, bool share, IdxIntView& x) {
	    	Gecode::Int::IntView::update(home, share, x);
	    	index = x.index;
	    }
		int getIndex(void) const {return index;}
	};
	
	typedef IdxIntView MyViewType;
	
	
  }
}


#endif /*GC_SEARCH_HH_*/
