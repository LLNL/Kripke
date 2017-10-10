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

#ifndef KRIPKE_CORE_DATASTORE_H__
#define KRIPKE_CORE_DATASTORE_H__

#include <Kripke.h>
#include <Kripke/Core/BaseVar.h>
#include <map>
#include <stdexcept>
#include <utility>

namespace Kripke {
namespace Core {

/**
 * Container to store variables by name
 */
class DataStore {
  public:
    DataStore();
    ~DataStore();
    DataStore(DataStore const &) = delete;
    DataStore &operator=(DataStore const &) = delete;

    void addVariable(std::string const &name, Kripke::Core::BaseVar *);

    template<typename T, typename ... CTOR_ARGS>
    RAJA_INLINE
    T &newVariable(std::string const &name, CTOR_ARGS &&... ctor_args){
      T *new_var = new T(ctor_args...);
      addVariable(name, new_var);
      return *new_var;
    }

    void deleteVariable(std::string const &name);

    template<typename T>
    RAJA_INLINE
    T &getVariable(std::string const &name){

      // Perform lookup by name
      auto it = m_vars.find(name);
      if(it == m_vars.end()){
        throw std::domain_error("Cannot find '" + name + "' in DataStore");
      }

      // Cast from BaseVar* and check for correctness
      T *var_ptr = dynamic_cast<T*>(it->second);
      KRIPKE_ASSERT(var_ptr != nullptr, "Error casting '%s'", name.c_str());

      return *var_ptr;
    }

    template<typename T>
    RAJA_INLINE
    T const &getVariable(std::string const &name) const{
      return const_cast<DataStore *>(this)-> template getVariable<T>(name);
    }
    

    template<typename T>
    RAJA_INLINE
    bool isVariableType(std::string const &name) const{

      // Perform lookup by name
      auto it = m_vars.find(name);
      if(it == m_vars.end()){
        return false;
      }

      // Cast from BaseVar* to see if it's correct type
      T *var_ptr = dynamic_cast<T*>(it->second);

      return var_ptr != nullptr;
    }

    std::vector<std::string> getVariableList() const;

  private:
    std::map<std::string, Kripke::Core::BaseVar *> m_vars;

};

} } // namespace

#endif
