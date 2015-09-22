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

#ifndef KERNEL_VARIABLE_POLICY_H__
#define KERNEL_VARIABLE_POLICY_H__

#include<Kripke.h>
#include<Domain/Layout.h>
#include<Domain/Forall.h>

struct FixedVariablePolicy {
  typedef LAYOUT_JI LayoutEll;
  typedef LAYOUT_IJ LayoutEllPlus;
};


template<typename T>
struct VariablePolicy{};

template<>
struct VariablePolicy<NEST_DGZ_T> : public FixedVariablePolicy {
  typedef LAYOUT_IJK LayoutPsi;
  typedef LAYOUT_IJK LayoutPhi;
  typedef LAYOUT_IJKL LayoutSigS;
};

template<>
struct VariablePolicy<NEST_DZG_T> : public FixedVariablePolicy {
  typedef LAYOUT_IKJ LayoutPsi;
  typedef LAYOUT_IKJ LayoutPhi;
  typedef LAYOUT_ILJK LayoutSigS;
};

template<>
struct VariablePolicy<NEST_GDZ_T> : public FixedVariablePolicy {
  typedef LAYOUT_JIK LayoutPsi;
  typedef LAYOUT_JIK LayoutPhi;
  typedef LAYOUT_JKIL LayoutSigS;
};

template<>
struct VariablePolicy<NEST_GZD_T> : public FixedVariablePolicy {
  typedef LAYOUT_JKI LayoutPsi;
  typedef LAYOUT_JKI LayoutPhi;
  typedef LAYOUT_JKLI LayoutSigS;
};

template<>
struct VariablePolicy<NEST_ZDG_T> : public FixedVariablePolicy {
  typedef LAYOUT_KIJ LayoutPsi;
  typedef LAYOUT_KIJ LayoutPhi;
  typedef LAYOUT_LIJK LayoutSigS;
};

template<>
struct VariablePolicy<NEST_ZGD_T> : public FixedVariablePolicy {
  typedef LAYOUT_KJI LayoutPsi;
  typedef LAYOUT_KJI LayoutPhi;
  typedef LAYOUT_LJKI LayoutSigS;  
};

template<typename T>
struct VariableView {
  typedef View3d<double, typename T::LayoutPsi> Psi; // D, G, Z
  typedef View3d<double, typename T::LayoutPsi> Phi; // NM, G, Z
  typedef View2d<double, typename T::LayoutEll> Ell; // D, NM
  typedef View2d<double, typename T::LayoutEllPlus> EllPlus; // D, NM  
  typedef View4d<double, typename T::LayoutSigS> SigS; // N, G, Gp, material  
};

#endif
