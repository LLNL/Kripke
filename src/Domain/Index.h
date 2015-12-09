#ifndef __DOMAIN_INDEX_H__
#define __DOMAIN_INDEX_H__

#include<string>
#include<RAJA/RAJA.hxx>



template<typename IndexTag>
struct IndexTraits {
};

/**
 * Strongly typed Index class.
 *
 * Allows integers to be associated with a type, and dissalows automatic
 * conversion.
 *
 * Use the DEF_INDEX(NAME) macro to define new indices.
 */
template<typename IndexTag>
class Index {

  public:
    RAJA_HOST_DEVICE explicit Index(int v) : value(v) {}
    
    RAJA_HOST_DEVICE inline int operator*(void) const {return value;}

    RAJA_HOST_DEVICE inline Index<IndexTag> operator+(int a) const { return Index<IndexTag>(value+a); }

    inline static std::string getName(void){
      return IndexTraits<IndexTag>::getName();
    }

  
  private:
    int value;

};

/**
 * Function provides a way to take either an int or any Index<> type, and
 * convert it to another type, possibly another Index or an int.
 */
template<typename TO>
RAJA_HOST_DEVICE inline TO convertIndex(int val){return TO(val);}


#define DEF_INDEX(NAME) \
  struct NAME##_TAG {};\
  template<>\
  struct IndexTraits<NAME##_TAG>{\
    static inline std::string getName(void){return #NAME;}\
  };\
  typedef Index<NAME##_TAG> NAME;\
  template<typename TO> \
  RAJA_HOST_DEVICE inline TO convertIndex(Index<NAME##_TAG> const &val){return TO(*val);}



#endif



