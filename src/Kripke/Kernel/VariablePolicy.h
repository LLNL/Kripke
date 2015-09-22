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

struct FixedLayoutPolicy {
  typedef LAYOUT_JI Ell;
  typedef LAYOUT_IJ EllPlus;
  
  typedef Layout3d<LAYOUT_KJI> Zone;
  typedef Layout2d<LAYOUT_JI> Face;
};


template<typename T>
struct LayoutPolicy{};

template<>
struct LayoutPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef LAYOUT_IJK Psi;
  typedef LAYOUT_IJK Phi;
  typedef LAYOUT_IJKL SigS;
  typedef LAYOUT_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef LAYOUT_IKJ Psi;
  typedef LAYOUT_IKJ Phi;
  typedef LAYOUT_ILJK SigS;
  typedef LAYOUT_JI SigT;
};

template<>
struct LayoutPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef LAYOUT_JIK Psi;
  typedef LAYOUT_JIK Phi;
  typedef LAYOUT_JKIL SigS;
  typedef LAYOUT_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef LAYOUT_JKI Psi;
  typedef LAYOUT_JKI Phi;
  typedef LAYOUT_JKLI SigS;
  typedef LAYOUT_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef LAYOUT_KIJ Psi;
  typedef LAYOUT_KIJ Phi;
  typedef LAYOUT_LIJK SigS;
  typedef LAYOUT_JI SigT;
};

template<>
struct LayoutPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef LAYOUT_KJI Psi;
  typedef LAYOUT_KJI Phi;
  typedef LAYOUT_LJKI SigS;
  typedef LAYOUT_JI SigT;  
};


template<typename T>
struct ViewPolicy {
  typedef View3d<double, typename T::Psi> Psi; // D, G, Z
  typedef View3d<double, typename T::Psi> Face; // D, G, Z
  typedef View3d<double, typename T::Phi> Phi; // NM, G, Z  
  typedef View2d<double, typename T::Ell> Ell; // D, NM
  typedef View2d<double, typename T::EllPlus> EllPlus; // D, NM  
  typedef View4d<double, typename T::SigS> SigS; // N, G, Gp, material  
  typedef View2d<double, typename T::SigT> SigT; // G, Z
};

#endif
