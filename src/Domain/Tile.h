#ifndef DOMAIN_TILE_H__
#define DOMAIN_TILE_H__

#include<RAJA/RAJA.hxx>


// Policy for no tiling
struct tile_none {};


// Policy to tile by given block size
template<int TileSize>
struct tile_fixed {};


template<typename TI, typename BODY>
void forall_tile(tile_none, TI const &is, BODY const &body){
  body(is);
}


template<int TileSize, typename BODY>
void forall_tile(tile_fixed<TileSize>, RAJA::RangeSegment const &is, BODY const &body){
  // tile loop
  int i_begin = is.getBegin();
  int i_end = is.getEnd();
  for(int i0 = i_begin;i0 < i_end;i0 += TileSize){
  
    // Create a new tile
    int i1 = std::min(i0+TileSize, i_end);
    RAJA::RangeSegment is_tile(i0, i1);
      
    // Pass tile index set to body        
    body(is_tile);
  }  
}




#endif

