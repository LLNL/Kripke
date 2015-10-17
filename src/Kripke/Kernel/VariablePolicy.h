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
#include<Domain/Index.h>


struct FixedLayoutPolicy {
  typedef PERM_JI Ell;     // d, nm
  typedef PERM_IJ EllPlus; // d, nm
  
  typedef Layout3d<PERM_KJI> Zone;
  typedef Layout2d<PERM_JI> Face;
};


template<typename T>
struct LayoutPolicy{};

template<>
struct LayoutPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef PERM_IJK Psi;
  typedef PERM_IJK Phi;
  typedef PERM_IJKL SigS;
  typedef PERM_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef PERM_IKJ Psi;
  typedef PERM_IKJ Phi;
  typedef PERM_ILJK SigS;
  typedef PERM_JI SigT;
};

template<>
struct LayoutPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef PERM_JIK Psi;
  typedef PERM_JIK Phi;
  typedef PERM_JKIL SigS;
  typedef PERM_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef PERM_JKI Psi;
  typedef PERM_JKI Phi;
  typedef PERM_JKLI SigS;
  typedef PERM_IJ SigT;
};

template<>
struct LayoutPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef PERM_KIJ Psi;
  typedef PERM_KIJ Phi;
  typedef PERM_LIJK SigS;
  typedef PERM_JI SigT;
};

template<>
struct LayoutPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef PERM_KJI Psi;
  typedef PERM_KJI Phi;
  typedef PERM_LJKI SigS;
  typedef PERM_JI SigT;
};


DEF_INDEX(IMaterial);
DEF_INDEX(ILegendre);
DEF_INDEX(IMoment);
DEF_INDEX(IDirection);
DEF_INDEX(IGlobalGroup);
DEF_INDEX(IGroup);
DEF_INDEX(IZone);
DEF_INDEX(IMix);
DEF_INDEX(IFaceI);
DEF_INDEX(IFaceJ);
DEF_INDEX(IFaceK);

template<typename T>
struct ViewPolicy {
  typedef View3d<double, typename T::Psi> Psi; // D, G, Z
  typedef View3d<double, typename T::Psi> Face; // D, G, Z
  typedef View3d<double, typename T::Phi> Phi; // NM, G, Z  
  typedef View2d<double, typename T::Ell> Ell; // D, NM
  typedef View2d<double, typename T::EllPlus> EllPlus; // D, NM  
  typedef View4d<double, typename T::SigS> SigS; // N, G, Gp, material  
  typedef View2d<double, typename T::SigT> SigT; // G, Z
  
  typedef TView3d<double, typename T::Psi, IDirection, IGroup, IZone> TPsi;
  typedef TView3d<double, typename T::Psi, IDirection, IGroup, IFaceI> TFaceI;
  typedef TView3d<double, typename T::Psi, IDirection, IGroup, IFaceJ> TFaceJ;
  typedef TView3d<double, typename T::Psi, IDirection, IGroup, IFaceK> TFaceK;
  typedef TView3d<double, typename T::Phi, IMoment, IGlobalGroup, IZone> TPhi;  
  typedef TView2d<double, typename T::Ell, IDirection, IMoment> TEll;
  typedef TView2d<double, typename T::EllPlus, IDirection, IMoment> TEllPlus;
  typedef TView4d<double, typename T::SigS, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> TSigS;
  typedef TView2d<double, typename T::EllPlus, IGroup, IZone> TSigT;
};

#endif
