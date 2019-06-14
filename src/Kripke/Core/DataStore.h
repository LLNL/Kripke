//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

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
    
    std::string getVariableName(BaseVar const &var) const;


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
