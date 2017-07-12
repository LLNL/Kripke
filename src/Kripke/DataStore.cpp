
#include<Kripke.h>
#include<Kripke/DataStore.h>

using namespace Kripke;

DataStore::DataStore(){}

DataStore::~DataStore(){

  
  while(m_vars.size()){
    auto it = m_vars.begin();
    deleteVariable(it->first);
  }

}

void DataStore::addVariable(std::string const &name,
  BaseVar *var)
{
  if(m_vars.find(name) != m_vars.end()){
    throw std::domain_error("Variable '" + name + "' already exists");
  }

  m_vars[name] = var;
}


void DataStore::deleteVariable(std::string const &name){
  auto it = m_vars.find(name);
  if(it == m_vars.end()){
    throw std::domain_error("Variable '" + name + "' does not exist");
  }

  // destroy object
  printf("Deleting %s\n", name.c_str());
  delete it->second;

  // remove from map
  m_vars.erase(it);
}


