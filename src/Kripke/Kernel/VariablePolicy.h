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
#include<Kripke/Directions.h>
#include<Domain/Layout.h>
#include<Domain/Forall.h>
#include<Domain/Index.h>


DEF_INDEX(IMaterial);
DEF_INDEX(ILegendre);
DEF_INDEX(IMoment);
DEF_INDEX(IDirection);
DEF_INDEX(IGlobalGroup);
DEF_INDEX(IGroup);
DEF_INDEX(IZone);
DEF_INDEX(IZoneIdx);
DEF_INDEX(IMix);
DEF_INDEX(IZoneI);
DEF_INDEX(IZoneJ);
DEF_INDEX(IZoneK);

struct FixedLayoutPolicy {
  typedef PERM_JI Ell;     // d, nm
  typedef PERM_IJ EllPlus; // d, nm
  
  typedef Layout3d<PERM_KJI> Zone;
  typedef TLayout3d<IZone, PERM_KJI, IZoneI, IZoneJ, IZoneK> TZone;      
};


template<typename T>
struct LayoutPolicy{};

template<>
struct LayoutPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef PERM_IJK Psi;
  typedef PERM_IJK Phi;
  typedef PERM_IJKL SigS;
  typedef PERM_IJ SigT;
  
  typedef PERM_IJLK FaceI; // d, g, j, k
  typedef PERM_IJLK FaceJ; // d, g, i, k
  typedef PERM_IJLK FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef PERM_IKJ Psi;
  typedef PERM_IKJ Phi;
  typedef PERM_ILJK SigS;
  typedef PERM_JI SigT;
  
  typedef PERM_ILKJ FaceI; // d, g, j, k
  typedef PERM_ILKJ FaceJ; // d, g, i, k
  typedef PERM_ILKJ FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef PERM_JIK Psi;
  typedef PERM_JIK Phi;
  typedef PERM_JKIL SigS;
  typedef PERM_IJ SigT;
  
  typedef PERM_JILK FaceI; // d, g, j, k
  typedef PERM_JILK FaceJ; // d, g, i, k
  typedef PERM_JILK FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef PERM_JKI Psi;
  typedef PERM_JKI Phi;
  typedef PERM_JKLI SigS;
  typedef PERM_IJ SigT;
  
  typedef PERM_JLKI FaceI; // d, g, j, k
  typedef PERM_JLKI FaceJ; // d, g, i, k
  typedef PERM_JLKI FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef PERM_KIJ Psi;
  typedef PERM_KIJ Phi;
  typedef PERM_LIJK SigS;
  typedef PERM_JI SigT;
  
  typedef PERM_LKIJ FaceI; // d, g, j, k
  typedef PERM_LKIJ FaceJ; // d, g, i, k
  typedef PERM_LKIJ FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef PERM_KJI Psi;
  typedef PERM_KJI Phi;
  typedef PERM_LJKI SigS;
  typedef PERM_JI SigT;
  
  typedef PERM_LKJI FaceI; // d, g, j, k
  typedef PERM_LKJI FaceJ; // d, g, i, k
  typedef PERM_LKJI FaceK; // d, g, i, j
};

struct FixedViewPolicy {
  typedef TView1d<double, PERM_I, IZoneI> Tdx;  
  typedef TView1d<double, PERM_I, IZoneJ> Tdy;  
  typedef TView1d<double, PERM_I, IZoneK> Tdz;
  typedef TView1d<Directions, PERM_I, IDirection> TDirections;
  
  typedef TView1d<IZoneI, PERM_I, IZoneIdx> TIdxToI;
  typedef TView1d<IZoneJ, PERM_I, IZoneIdx> TIdxToJ;
  typedef TView1d<IZoneK, PERM_I, IZoneIdx> TIdxToK;
};

template<typename T>
struct ViewPolicy : public FixedViewPolicy {
  typedef TView3d<double, typename T::Psi, IDirection, IGroup, IZone> TPsi;
  typedef TView4d<double, typename T::FaceI, IDirection, IGroup, IZoneJ, IZoneK> TFaceI;
  typedef TView4d<double, typename T::FaceJ, IDirection, IGroup, IZoneI, IZoneK> TFaceJ;
  typedef TView4d<double, typename T::FaceK, IDirection, IGroup, IZoneI, IZoneJ> TFaceK;
  typedef TView3d<double, typename T::Phi, IMoment, IGlobalGroup, IZone> TPhi;  
  typedef TView2d<double, typename T::Ell, IDirection, IMoment> TEll;
  typedef TView2d<double, typename T::EllPlus, IDirection, IMoment> TEllPlus;
  typedef TView4d<double, typename T::SigS, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> TSigS;
  typedef TView2d<double, typename T::EllPlus, IGroup, IZone> TSigT;
};

#endif
