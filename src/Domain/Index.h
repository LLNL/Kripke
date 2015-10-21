#ifndef __DOMAIN_INDEX_H__
#define __DOMAIN_INDEX_H__

#include<string>
#include<RAJA/RangeSegment.hxx>


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

/**
 * Function provides a way to take either an int or any Index<> type, and
 * convert it to another type, possibly another Index or an int.
 */
template<typename TO>
inline TO convertIndex(int val){return TO(val);}


#define DEF_INDEX(NAME) \
  struct NAME##_TAG {};\
  template<>\
  struct IndexTraits<NAME##_TAG>{\
    static inline std::string getName(void){return #NAME;}\
  };\
  typedef Index<NAME##_TAG> NAME;\
  template<typename TO> \
  inline TO convertIndex(Index<NAME##_TAG> const &val){return TO(*val);}



#endif



