
#include<Kripke.h>
#include<Kripke/DataStore.h>


DataStore::DataStore(){}

DataStore::~DataStore(){

  
  while(vars.size()){
    auto it = vars.begin();
    deleteVariable(it->first);
  }

}

void DataStore::addVariable(std::string const &name,
  BaseVar *var)
{
  if(vars.find(name) != vars.end()){
    throw std::domain_error("Variable '" + name + "' already exists");
  }

  vars[name] = var;
}


void DataStore::deleteVariable(std::string const &name){
  auto it = vars.find(name);
  if(it == vars.end()){
    throw std::domain_error("Variable '" + name + "' does not exist");
  }

  // destroy object
  delete it->second;

  // remove from map
  vars.erase(it);
}


