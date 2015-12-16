#ifndef RAJA_INDEX_HXX__
#define RAJA_INDEX_HXX__

#include<RAJA/RangeSegment.hxx>
#include<string>


namespace RAJA {

/**
 * Strongly typed Index class.
 *
 * Allows integers to be associated with a type, and dissalows automatic
 * conversion.
 *
 * Use the DEF_INDEX(NAME) macro to define new indices.
 *
 * Yes, this uses the curiosly-recurring template pattern.
 */
template<typename TYPE>
class IndexValue {

  public:
    inline explicit IndexValue(int v) : value(v) {}

    /**
     * Dereference provides cast-to-integer.
     */
    inline int operator*(void) const {return value;}

    inline TYPE operator+(int a) const { return TYPE(value+a); }

    //static std::string getName(void);
  
  private:
    int value;

};


/**
 * Helper class for convertIndex, since functions cannot be partially
 * specialized
 */
template<typename TO, typename FROM>
struct ConvertIndexHelper {
  static inline TO convert(FROM val){
    return TO(*val);
  }
};

template<typename TO>
struct ConvertIndexHelper<TO, int> {
  static inline TO convert(int val){
    return TO(val);
  }
};

/**
 * Function provides a way to take either an int or any Index<> type, and
 * convert it to another type, possibly another Index or an int.
 */
template<typename TO, typename FROM>
inline TO convertIndex(FROM val){
  return ConvertIndexHelper<TO, FROM>::convert(val);
}




} // namespace RAJA



/**
 * Helper Macro to create new Index types.
 */
#define DEF_INDEX(NAME) \
  class NAME : public RAJA::IndexValue<NAME>{ \
  public: \
    inline explicit NAME(int v) : RAJA::IndexValue<NAME>::IndexValue(v) {} \
    static inline std::string getName(void){return #NAME;} \
  };


#endif



