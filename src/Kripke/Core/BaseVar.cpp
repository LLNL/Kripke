//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include<Kripke/Core/BaseVar.h>

#include<Kripke/Core/DataStore.h>

using namespace Kripke::Core;

BaseVar::BaseVar() : m_parent(nullptr){

}

void BaseVar::setParent(DataStore *parent){
  m_parent = parent;
}

std::string BaseVar::getName() const {
  if(m_parent){
    return m_parent->getVariableName(*this);
  }
  else{
    return "---";
  }
}


