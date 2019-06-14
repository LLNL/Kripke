//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_BASE_VAR_H__
#define KRIPKE_CORE_BASE_VAR_H__

#include <string>

namespace Kripke {
namespace Core {

class DataStore;

/**
 * Variable base class for DataStore class
 */
class BaseVar {
  public:
    BaseVar();
    virtual ~BaseVar() = default;

    void setParent(DataStore *parent);

    std::string getName() const;

  private:
    DataStore *m_parent;
};


} } // namespace

#endif
