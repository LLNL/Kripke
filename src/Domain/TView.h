#ifndef __DOMAIN_TVIEW_H__
#define __DOMAIN_TVIEW_H__

#include<string>
#include<RAJA/RangeSegment.hxx>
#include<Domain/Index.h>
#include<Domain/View.h>


template<typename DataType, typename Layout>
struct DView1d : public View1d<DataType, Layout> {

  typedef typename Layout::Permutation Permutation;
  typedef typename Layout::IndexI IndexI;

  inline DView1d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    View1d<DataType, Layout>::View1d(
        ptr,
        IndexI::size(domain, sdom_id))
  {}
};

template<typename DataType, typename Layout>
struct DView2d : public View2d<DataType, Layout> {

  typedef typename Layout::Permutation Permutation;
  typedef typename Layout::IndexI IndexI;
  typedef typename Layout::IndexJ IndexJ;

  inline DView2d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    View2d<DataType, Layout>::View2d(
        ptr,
        IndexI::size(domain, sdom_id),
        IndexJ::size(domain, sdom_id))
  {}
};

template<typename DataType, typename Layout>
struct DView3d : public View3d<DataType, Layout> {

  typedef typename Layout::Permutation Permutation;
  typedef typename Layout::IndexI IndexI;
  typedef typename Layout::IndexJ IndexJ;
  typedef typename Layout::IndexK IndexK;

  inline DView3d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    View3d<DataType, Layout>::View3d(
        ptr,
        IndexI::size(domain, sdom_id),
        IndexJ::size(domain, sdom_id),
        IndexK::size(domain, sdom_id))
  {}
};

template<typename DataType, typename Layout>
struct DView4d : public View4d<DataType, Layout> {

  typedef typename Layout::Permutation Permutation;
  typedef typename Layout::IndexI IndexI;
  typedef typename Layout::IndexJ IndexJ;
  typedef typename Layout::IndexK IndexK;
  typedef typename Layout::IndexL IndexL;

  inline DView4d(DataType *ptr, Grid_Data *domain, int sdom_id) :
    View4d<DataType, Layout>::View4d(
        ptr,
        IndexI::size(domain, sdom_id),
        IndexJ::size(domain, sdom_id),
        IndexK::size(domain, sdom_id),
        IndexL::size(domain, sdom_id))
  {}
};


/**
 * Wrapper around Layout1d that provides accessors to Index sizes
 */
template<typename Perm, typename IdxLin, typename IdxI>
struct DLayout1d : public Layout1d<Perm, IdxLin, IdxI>{

  inline DLayout1d(Grid_Data *domain, int sdom_id) :
    Layout1d<Perm, IdxLin, IdxI>::Layout1d(
          IdxI::size(domain, sdom_id))
  {}

  inline DLayout1d(int ni) :
    Layout1d<Perm, IdxLin, IdxI>::Layout1d(ni)
  {}

};


/**
 * Wrapper around Layout2d that provides accessors to Index sizes
 */
template<typename Perm, typename IdxLin, typename IdxI, typename IdxJ>
struct DLayout2d : public Layout2d<Perm, IdxLin, IdxI, IdxJ>{

  inline DLayout2d(Grid_Data *domain, int sdom_id) :
    Layout2d<Perm, IdxLin, IdxI, IdxJ>::Layout2d(
          IdxI::size(domain, sdom_id),
          IdxJ::size(domain, sdom_id))
  {}

  inline DLayout2d(int ni, int nj) :
    Layout2d<Perm, IdxLin, IdxI, IdxJ>::Layout2d(ni, nj)
  {}

};

/**
 * Wrapper around Layout3d that provides accessors to Index sizes
 */
template<typename Perm, typename IdxLin, typename IdxI, typename IdxJ, typename IdxK>
struct DLayout3d : public Layout3d<Perm, IdxLin, IdxI, IdxJ, IdxK>{

  inline DLayout3d(Grid_Data *domain, int sdom_id) :
    Layout3d<Perm, IdxLin, IdxI, IdxJ, IdxK>::Layout3d(
        IdxI::size(domain, sdom_id),
        IdxJ::size(domain, sdom_id),
        IdxK::size(domain, sdom_id))
  {}

  inline DLayout3d(int ni, int nj, int nk) :
    Layout3d<Perm, IdxLin, IdxI, IdxJ, IdxK>::Layout3d(ni, nj, nk)
  {}

};


/**
 * Wrapper around Layout4d that provides accessors to Index sizes
 */
template<typename Perm, typename IdxLin, typename IdxI, typename IdxJ, typename IdxK, typename IdxL>
struct DLayout4d : public Layout4d<Perm, IdxLin, IdxI, IdxJ, IdxK, IdxL>{

  inline DLayout4d(Grid_Data *domain, int sdom_id) :
    Layout4d<Perm, IdxLin, IdxI, IdxJ, IdxK, IdxL>::Layout4d(
        IdxI::size(domain, sdom_id),
        IdxJ::size(domain, sdom_id),
        IdxK::size(domain, sdom_id),
        IdxL::size(domain, sdom_id))
  {}

  inline DLayout4d(int ni, int nj, int nk, int nl) :
    Layout4d<Perm, IdxLin, IdxI, IdxJ, IdxK, IdxL>::Layout4d(ni, nj, nk, nl)
  {}

};


template<typename POL, typename IdxI, typename IdxJ, typename BODY>
void dForall2(Grid_Data *domain, int sdom_id, BODY const &body){

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
void dForall4(Grid_Data *domain, int sdom_id, BODY const &body){

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



