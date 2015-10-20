#ifndef __DOMAIN_INDEX_H__
#define __DOMAIN_INDEX_H__

#include<string>
#include<RAJA/RangeSegment.hxx>


template<typename IndexTag>
struct IndexTraits {};

template<typename IndexTag>
class Index {

  public:
    explicit Index(int v) : value(v) {}
    
    inline int operator*(void) const {return value;}

    inline static int size(Grid_Data *domain, int sdom_id){
      std::string name = IndexTraits<IndexTag>::getName();
      
      return domain->subdomains[sdom_id].index_size[name];      
    }
    
    inline Index<IndexTag> operator+(int a) const { return Index<IndexTag>(value+a); }

    
    inline static RAJA::RangeSegment range(Grid_Data *domain, int sdom_id){

      int len = Index<IndexTag>::size(domain, sdom_id);
            
      return RAJA::RangeSegment(0, len);
    }

    inline static std::string getName(void){
      return IndexTraits<IndexTag>::getName();
    }

  
  private:
    int value;

};



#define DEF_INDEX(NAME) \
  struct NAME##_TAG {};\
  template<>\
  struct IndexTraits<NAME##_TAG>{\
    static inline std::string getName(void){return #NAME;}\
  };\
  typedef Index<NAME##_TAG> NAME;


template<typename T, typename L, typename IdxI>
struct TView1d {
    inline TView1d(T* ptr, Grid_Data *domain, int sdom_id) :
      view(ptr, 
           IdxI::size(domain, sdom_id))
    {}

    inline T &operator()(IdxI i) const{
      return view(*i);
    }

    View1d<T, L> const view;
};

template<typename T, typename L, typename IdxI, typename IdxJ>
struct TView2d {
    inline TView2d(T* ptr, Grid_Data *domain, int sdom_id) :
      view(ptr, 
           IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id))
    {}

    inline T &operator()(IdxI i, IdxJ j) const{
      return view(*i, *j);
    }

    View2d<T, L> const view;
};


template<typename T, typename L, typename IdxI, typename IdxJ, typename IdxK>
struct TView3d {
    inline TView3d(T* ptr, Grid_Data *domain, int sdom_id) :
      view(ptr, 
           IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id),
           IdxK::size(domain, sdom_id))
    {}

    inline T &operator()(IdxI i, IdxJ j, IdxK k) const{
      return view(*i, *j, *k);
    }

    View3d<T, L> const view;
};

template<typename T, typename L, typename IdxI, typename IdxJ, typename IdxK, typename IdxL>
struct TView4d {
    inline TView4d(T* ptr, Grid_Data *domain, int sdom_id) :
      view(ptr, 
           IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id),
           IdxK::size(domain, sdom_id),
           IdxL::size(domain, sdom_id))
    {}

    inline T &operator()(IdxI i, IdxJ j, IdxK k, IdxL l) const{
      return view(*i, *j, *k, *l);
    }

    View4d<T, L> const view;
};


template<typename T, typename L, typename IdxI>
struct TLayout1d {
    inline TLayout1d(Grid_Data *domain, int sdom_id) :
      layout(IdxI::size(domain, sdom_id))
    {}

    inline T operator()(IdxI i) const{
      return T(layout(*i));
    }

    Layout1d<L> const layout;
};

template<typename T, typename L, typename IdxI, typename IdxJ, typename IdxK>
struct TLayout3d {
    inline TLayout3d(Grid_Data *domain, int sdom_id) :
      layout(IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id),
           IdxK::size(domain, sdom_id))
    {}

    inline T operator()(IdxI i, IdxJ j, IdxK k) const{
      return T(layout(*i, *j, *k));
    }

    Layout3d<L> const layout;
};


template<typename POL, typename IdxI, typename IdxJ, typename BODY>
void forall2T(Grid_Data *domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = IdxI::range(domain, sdom_id);
  RAJA::RangeSegment seg_j = IdxJ::range(domain, sdom_id);

  // Call underlying forall, extracting ranges from domain
  forall2<POL>(seg_i, seg_j,  
    [=](int i, int j){
      // cast indicies to index types
      body(IdxI(i), IdxJ(j));    
    }
  );      
}


template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename TI, typename TJ, typename TK, typename BODY>
void forall3T(TI const &is_i, TJ const &is_j, TK const &is_k, BODY const &body){

  // Call underlying forall, extracting ranges from domain
  forall3<POL, TI, TJ, TK>(is_i, is_j, is_k,  
    [=](int i, int j, int k){
      // cast indicies to index types
      body(IdxI(i), IdxJ(j), IdxK(k));    
    }
  );      
}



template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename IdxL, typename BODY>
void forall4T(Grid_Data *domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = IdxI::range(domain, sdom_id);
  RAJA::RangeSegment seg_j = IdxJ::range(domain, sdom_id);
  RAJA::RangeSegment seg_k = IdxK::range(domain, sdom_id);
  RAJA::RangeSegment seg_l = IdxL::range(domain, sdom_id);

  // Call underlying forall, extracting ranges from domain
  forall4<POL>(seg_i, seg_j, seg_k, seg_l, 
    [=](int i, int j, int k, int l){
      // cast indicies to index types
      body(IdxI(i), IdxJ(j), IdxK(k), IdxL(l));    
    }
  );
}

#endif



