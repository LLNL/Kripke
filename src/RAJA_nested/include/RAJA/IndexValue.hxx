#ifndef RAJA_INDEXVALUE_HXX__
#define RAJA_INDEXVALUE_HXX__

#include<RAJA/int_datatypes.hxx>
#include<RAJA/RangeSegment.hxx>
#include<string>


namespace RAJA {

/**
 * Strongly typed "integer" class.
 *
 * Allows integers to be associated with a type, and disallows automatic
 * conversion.
 *
 * Use the RAJA_INDEX_VALUE(NAME) macro to define new indices.
 *
 * Yes, this uses the curiously-recurring template pattern.
 */
template<typename TYPE>
class IndexValue {

  public:
    /**
     * Default constructor.
     * Initializes value to 0.
     */
    inline IndexValue() : value(0) {}

    /**
     * Explicit constructor.
     * \param v   Initial value
     */
    inline explicit IndexValue(Index_type v) : value(v) {}

    /**
     * Dereference provides cast-to-integer.
     */
    inline Index_type operator*(void) const {return value;}


    /**
     * Increment by one.
     */
    inline TYPE &operator++(int){
      value++;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Increment by one.
     */
    inline TYPE &operator++(){
      value++;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Decrement by one.
     */
    inline TYPE &operator--(int){
      value--;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Decrement by one.
     */
    inline TYPE &operator--(){
      value--;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Addition.
     */
    inline TYPE operator+(Index_type a) const { return TYPE(value+a); }

    /**
     * Addition.
     */
    inline TYPE operator+(TYPE a) const { return TYPE(value+a.value); }

    /**
     * Subtraction.
     */
    inline TYPE operator-(Index_type a) const { return TYPE(value-a); }

    /**
     * Subtraction.
     */
    inline TYPE operator-(TYPE a) const { return TYPE(value-a.value); }

    /**
     * Multiply.
     */
    inline TYPE operator*(Index_type a) const { return TYPE(value*a); }

    /**
     * Multiply.
     */
    inline TYPE operator*(TYPE a) const { return TYPE(value*a.value); }

    /**
     * Divide.
     */
    inline TYPE operator/(Index_type a) const { return TYPE(value/a); }

    /**
     * Divide.
     */
    inline TYPE operator/(TYPE a) const { return TYPE(value/a.value); }



    /**
     * Add to.
     */
    inline TYPE &operator+=(Index_type x){
      value += x;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Add to.
     */
    inline TYPE &operator+=(TYPE x){
      value += x.value;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Subtract from.
     */
    inline TYPE &operator-=(Index_type x){
      value -= x;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Subtract from.
     */
    inline TYPE &operator-=(TYPE x){
      value -= x.value;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Multiply by.
     */
    inline TYPE &operator*=(Index_type x){
      value *= x;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Multiply by.
     */
    inline TYPE &operator*=(TYPE x){
      value *= x.value;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Divide by.
     */
    inline TYPE &operator/=(Index_type x){
      value /= x;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Divide by.
     */
    inline TYPE &operator/=(TYPE x){
      value /= x.value;
      return *static_cast<TYPE *>(this);
    }

    /**
     * Less-than.
     */
    inline bool operator<(Index_type x) const {
      return( value < x);
    }

    /**
     * Less-than.
     */
    inline bool operator<(TYPE x) const {
      return( value < x.value);
    }

    /**
     * Less-than or equals.
     */
    inline bool operator<=(Index_type x) const {
      return( value <= x);
    }

    /**
     * Less-than or equals.
     */
    inline bool operator<=(TYPE x) const {
      return( value <= x.value);
    }

    /**
     * Greater-than.
     */
    inline bool operator>(Index_type x) const {
      return( value > x);
    }

    /**
     * Greater-than.
     */
    inline bool operator>(TYPE x) const {
      return( value > x.value);
    }

    /**
     * Greater-than or equals.
     */
    inline bool operator>=(Index_type x) const {
      return( value >= x);
    }

    /**
     * Greater-than or equals.
     */
    inline bool operator>=(TYPE x) const {
      return( value >= x.value);
    }

    /**
     * Equal to.
     */
    inline bool operator==(Index_type x) const {
      return( value == x);
    }

    /**
     * Equal to.
     */
    inline bool operator==(TYPE x) const {
      return( value == x.value);
    }

    /**
     * Not equal to.
     */
    inline bool operator!=(Index_type x) const {
      return( value != x);
    }

    /**
     * Not equal to.
     */
    inline bool operator!=(TYPE x) const {
      return( value != x.value);
    }


    // This is not implemented... but should be by the derived type
    // this is done by the macro
    static std::string getName(void);
  
  private:
    Index_type value;

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
struct ConvertIndexHelper<TO, Index_type> {
  static inline TO convert(Index_type val){
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
#define RAJA_INDEX_VALUE(NAME) \
  class NAME : public RAJA::IndexValue<NAME>{ \
  public: \
    inline explicit NAME(RAJA::Index_type v) : RAJA::IndexValue<NAME>::IndexValue(v) {} \
    static inline std::string getName(void){return #NAME;} \
  };


#endif



