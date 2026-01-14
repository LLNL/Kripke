//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>

#ifdef KRIPKE_USE_CHAI
#define DEBUG
#include <umpire/Umpire.hpp>
#include <umpire/strategy/QuickPool.hpp>
#include <chai/ManagedArray.hpp>
#undef DEBUG
#endif

using namespace Kripke;
using namespace Kripke::Core;

DataStore::DataStore(int dev_pool_size){
#ifdef KRIPKE_USE_CHAI
  auto &rm = umpire::ResourceManager::getInstance();
  const char * allocator_name = "KRIPKE_DEVICE_POOL";
  size_t umpire_dev_pool_size = ((size_t) dev_pool_size) * 1024 * 1024 * 1024;
  size_t umpire_dev_block_size = 512;
  auto dev_pool_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(allocator_name, rm.getAllocator("DEVICE"), umpire_dev_pool_size, umpire_dev_block_size);
  auto chai_resource_manager = chai::ArrayManager::getInstance();
  chai_resource_manager->setAllocator(chai::GPU, dev_pool_allocator);
  // force allocation of GPU memory pool
  auto tmp = new chai::ManagedArray<int>(100, chai::GPU);
  tmp->free(chai::GPU);
  delete tmp;
#endif // KRIPKE_USE_CHAI
}

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
