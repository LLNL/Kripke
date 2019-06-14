//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>

using namespace Kripke;
using namespace Kripke::Core;

DataStore::DataStore(){}

DataStore::~DataStore(){

  
  while(m_vars.size()){
    auto it = m_vars.begin();
    deleteVariable(it->first);
  }

}

void DataStore::addVariable(std::string const &name,
  Kripke::Core::BaseVar *var)
{
  if(m_vars.find(name) != m_vars.end()){
    throw std::domain_error("Variable '" + name + "' already exists");
  }

  m_vars[name] = var;

  var->setParent(this);
}


void DataStore::deleteVariable(std::string const &name){
  auto it = m_vars.find(name);
  if(it == m_vars.end()){
    throw std::domain_error("Variable '" + name + "' does not exist");
  }

  // destroy object
  //printf("Deleting %s\n", name.c_str());
  delete it->second;

  // remove from map
  m_vars.erase(it);
}


std::vector<std::string> DataStore::getVariableList() const{
  std::vector<std::string> var_list;

  for(auto &i : m_vars){
    var_list.push_back(i.first);
  }

  return var_list;
}


std::string DataStore::getVariableName(BaseVar const &var) const{
  for(auto &i : m_vars){
    if(i.second == &var){
      return i.first;
    }
  }
  return "===";
}
