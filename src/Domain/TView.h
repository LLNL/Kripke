#ifndef __DOMAIN_TVIEW_H__
#define __DOMAIN_TVIEW_H__

#include<string>
#include<RAJA/RangeSegment.hxx>
#include<Domain/Index.h>


template<typename T, typename L, typename IdxI>
struct TLayout1d {
  typedef T    IndexLinear;
  typedef IdxI IndexI;

  inline TLayout1d(Grid_Data *domain, int sdom_id) :
    layout(IdxI::size(domain, sdom_id))
  {}

  inline T operator()(IdxI i) const{
    return T(layout(*i));
  }

  Layout1d<L> const layout;
};


template<typename T, typename L, typename IdxI, typename IdxJ>
struct TLayout2d {
  typedef T    IndexLinear;
  typedef IdxI IndexI;
  typedef IdxJ IndexJ;

  inline TLayout2d(Grid_Data *domain, int sdom_id) :
    layout(IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id))
  {}

  inline T operator()(IdxI i, IdxJ j) const{
    return T(layout(*i, *j));
  }

  Layout2d<L> const layout;
};

template<typename T, typename L, typename IdxI, typename IdxJ, typename IdxK>
struct TLayout3d {
  typedef T    IndexLinear;
  typedef IdxI IndexI;
  typedef IdxJ IndexJ;
  typedef IdxK IndexK;

  inline TLayout3d(Grid_Data *domain, int sdom_id) :
    layout(IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id),
           IdxK::size(domain, sdom_id))
  {}

  inline T operator()(IdxI i, IdxJ j, IdxK k) const{
    return T(layout(convertIndex<int>(i), convertIndex<int>(j), convertIndex<int>(k)));
  }

  Layout3d<L> const layout;
};

template<typename T, typename L, typename IdxI, typename IdxJ, typename IdxK, typename IdxL>
struct TLayout4d {
  typedef T    IndexLinear;
  typedef IdxI IndexI;
  typedef IdxJ IndexJ;
  typedef IdxK IndexK;
  typedef IdxL IndexL;

  inline TLayout4d(Grid_Data *domain, int sdom_id) :
    layout(IdxI::size(domain, sdom_id),
           IdxJ::size(domain, sdom_id),
           IdxK::size(domain, sdom_id),
           IdxL::size(domain, sdom_id))
  {}

  inline T operator()(IdxI i, IdxJ j, IdxK k, IdxL l) const{
    return T(layout(*i, *j, *k, *l));
  }

  Layout4d<L> const layout;
};



template<typename DataType, typename TL>
struct TView1d {
  typedef typename TL::IndexI IndexI;

  inline TView1d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    layout(domain, sdom_id),
    data(ptr)
  {}

  inline DataType &operator()(IndexI i) const{
    return data[layout(i)];
  }

  TL const layout;
  DataType *data;
};


template<typename DataType, typename TL>
struct TView2d {
  typedef typename TL::IndexI IndexI;
  typedef typename TL::IndexJ IndexJ;

  inline TView2d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    layout(domain, sdom_id),
    data(ptr)
  {}

  inline DataType &operator()(IndexI i, IndexJ j) const{
    return data[layout(i, j)];
  }

  TL const layout;
  DataType *data;
};


template<typename DataType, typename TL>
struct TView3d {
  typedef typename TL::IndexI IndexI;
  typedef typename TL::IndexJ IndexJ;
  typedef typename TL::IndexK IndexK;

  inline TView3d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    layout(domain, sdom_id),
    data(ptr)
  {}

  inline DataType &operator()(IndexI i, IndexJ j, IndexK k) const{
    return data[layout(i, j, k)];
  }

  TL const layout;
  DataType *data;
};

template<typename DataType, typename TL>
struct TView4d {
  typedef typename TL::IndexI IndexI;
  typedef typename TL::IndexJ IndexJ;
  typedef typename TL::IndexK IndexK;
  typedef typename TL::IndexL IndexL;

  inline TView4d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    layout(domain, sdom_id),
    data(ptr)
  {}

  inline DataType &operator()(IndexI i, IndexJ j, IndexK k, IndexL l) const{
    return data[layout(i, j, k, l)];
  }

  TL const layout;
  DataType *data;
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



