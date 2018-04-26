/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#ifndef KRIPKE_CORE_VARLAYOUT_H__
#define KRIPKE_CORE_VARLAYOUT_H__

namespace Kripke {
namespace Core {

template<camp::idx_t N, typename Order, typename T>
struct GetOrderHelper;

template<camp::idx_t N, typename OrderT, typename ... OrderTRest, typename T>
struct GetOrderHelper<N, camp::list<OrderT, OrderTRest...>, T>{
  static const camp::idx_t value = GetOrderHelper<N, camp::list<OrderTRest...>, T>::value;
};

template<camp::idx_t N, typename OrderT, typename ... OrderTRest>
struct GetOrderHelper<N, camp::list<OrderT, OrderTRest...>, OrderT>{
  static const camp::idx_t value = (N-1) - sizeof...(OrderTRest);
};


template<typename Order, typename T>
constexpr camp::idx_t getOrder(){
  return GetOrderHelper<camp::size<Order>::value, Order, T>::value;
}

template<camp::idx_t i, typename Order, typename Types>
struct ExtractStrideOne;

template<camp::idx_t i, typename Order, typename ... Types>
struct ExtractStrideOne<i, Order, camp::list<Types...>>{
  using LTypes = camp::list<Types...>;
  using T = camp::at_v<LTypes, i>;
  using next_t = ExtractStrideOne<i+1, Order, LTypes>;

  static constexpr camp::idx_t our_value = getOrder<Order, T>();

  static constexpr camp::idx_t value =
    next_t::value >= our_value ? next_t::value : our_value;  
  
  static constexpr camp::idx_t arg =
    next_t::value >= our_value ? next_t::arg : i;  
};

template<typename Order, typename ... Types>
struct ExtractStrideOne<sizeof...(Types), Order, camp::list<Types...>>{
  
  static constexpr camp::idx_t value = 0;
  static constexpr camp::idx_t arg = 0;
};

template<typename Order, typename ... T>
struct ArgsToOrder {

  static constexpr camp::idx_t num_types = sizeof...(T);
  using type = camp::idx_seq<getOrder<Order, T>()...>;

  using array_t = std::array<camp::idx_t, sizeof...(T)>;


  static constexpr camp::idx_t stride_one = 
    ExtractStrideOne<0, Order, camp::list<T...> >::arg;    

  template<camp::idx_t ... RangeInts, camp::idx_t ... OrderInts>
  static array_t toArray_expanded(bool debug, camp::idx_seq<RangeInts...>, camp::idx_seq<OrderInts...>){
    using pair_t = std::pair<camp::idx_t, camp::idx_t>;
    using parray_t = std::array<pair_t, sizeof...(T)>;

    parray_t p{{pair_t{RangeInts, OrderInts}...}};

    std::sort(p.begin(), p.end(), 
      [=](pair_t const & a, pair_t const & b){
        return a.second < b.second;
      });

    if(debug){
      array_t a{{(p[RangeInts].second)...}};
      return a;
    }
    else{
      array_t a{{(p[RangeInts].first)...}};
      return a;
    }
  }

  static array_t toArray(bool debug = false){
    return toArray_expanded(debug, camp::make_idx_seq_t<sizeof...(T)>{}, camp::idx_seq<getOrder<Order,T>()...>{});
  }


  static void print(){
    array_t a = toArray(true);

    for(camp::idx_t i = 0;i < (camp::idx_t)sizeof...(T);++i){
      printf("%d ", (int)a[i]);
    }
    printf(" [stride-one arg=%d]\n", (int)stride_one);
  }
  
};



/*
 * Default layout is canonical ordering.
 *
 * This class is specialized for fields that needs data layouts to change.
 */

template<typename Order, typename ... IndexTypes>
struct LayoutInfo {
  
  using args_to_order_t = ArgsToOrder<Order, IndexTypes...>;

  // Default stride-one-index is the right-most index
  constexpr static ptrdiff_t num_dims = sizeof...(IndexTypes);
  constexpr static ptrdiff_t stride_one_dim = args_to_order_t::stride_one;

  using Layout = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IndexTypes...>, stride_one_dim>;

  static std::array<RAJA::Index_type, num_dims> getPermutation(){
    return args_to_order_t::toArray();
  }
};


template<typename Order, typename ... IndexTypes>
using LayoutType = typename LayoutInfo<Order, IndexTypes...>::Layout;

template<typename Order, typename ElementType, typename ... IndexTypes>
using ViewType = RAJA::View<ElementType, LayoutType<Order, IndexTypes...>>;



} // namespace Core
} // namespace Kripke

#endif
