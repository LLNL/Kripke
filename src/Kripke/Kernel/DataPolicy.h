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
  typedef PERM_JI Layout_Ell;     // d, nm
  typedef PERM_IJ Layout_EllPlus; // d, nm
  
  typedef Layout3d<PERM_KJI> Layout_Zone;
  typedef TLayout3d<IZone, PERM_KJI, IZoneI, IZoneJ, IZoneK> TLayout_Zone;
};


template<typename T>
struct LayoutPolicy{};

template<>
struct LayoutPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef PERM_IJK Layout_Psi;
  typedef PERM_IJK Layout_Phi;
  typedef PERM_IJKL Layout_SigS;
  typedef PERM_IJ Layout_SigT;
  
  typedef PERM_IJLK Layout_FaceI; // d, g, j, k
  typedef PERM_IJLK Layout_FaceJ; // d, g, i, k
  typedef PERM_IJLK Layout_FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef PERM_IKJ Layout_Psi;
  typedef PERM_IKJ Layout_Phi;
  typedef PERM_ILJK Layout_SigS;
  typedef PERM_JI Layout_SigT;
  
  typedef PERM_ILKJ Layout_FaceI; // d, g, j, k
  typedef PERM_ILKJ Layout_FaceJ; // d, g, i, k
  typedef PERM_ILKJ Layout_FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef PERM_JIK Layout_Psi;
  typedef PERM_JIK Layout_Phi;
  typedef PERM_JKIL Layout_SigS;
  typedef PERM_IJ Layout_SigT;
  
  typedef PERM_JILK Layout_FaceI; // d, g, j, k
  typedef PERM_JILK Layout_FaceJ; // d, g, i, k
  typedef PERM_JILK Layout_FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef PERM_JKI Layout_Psi;
  typedef PERM_JKI Layout_Phi;
  typedef PERM_JKLI Layout_SigS;
  typedef PERM_IJ Layout_SigT;
  
  typedef PERM_JLKI Layout_FaceI; // d, g, j, k
  typedef PERM_JLKI Layout_FaceJ; // d, g, i, k
  typedef PERM_JLKI Layout_FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef PERM_KIJ Layout_Psi;
  typedef PERM_KIJ Layout_Phi;
  typedef PERM_LIJK Layout_SigS;
  typedef PERM_JI Layout_SigT;
  
  typedef PERM_LKIJ Layout_FaceI; // d, g, j, k
  typedef PERM_LKIJ Layout_FaceJ; // d, g, i, k
  typedef PERM_LKIJ Layout_FaceK; // d, g, i, j
};

template<>
struct LayoutPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef PERM_KJI Layout_Psi;
  typedef PERM_KJI Layout_Phi;
  typedef PERM_LJKI Layout_SigS;
  typedef PERM_JI Layout_SigT;
  
  typedef PERM_LKJI Layout_FaceI; // d, g, j, k
  typedef PERM_LKJI Layout_FaceJ; // d, g, i, k
  typedef PERM_LKJI Layout_FaceK; // d, g, i, j
};

struct FixedViewPolicy {
  typedef TView1d<double, PERM_I, IZoneI> View_dx;
  typedef TView1d<double, PERM_I, IZoneJ> View_dy;
  typedef TView1d<double, PERM_I, IZoneK> View_dz;
  typedef TView1d<Directions, PERM_I, IDirection> View_Directions;
  
  typedef TView1d<IZoneI, PERM_I, IZoneIdx> View_IdxToI;
  typedef TView1d<IZoneJ, PERM_I, IZoneIdx> View_IdxToJ;
  typedef TView1d<IZoneK, PERM_I, IZoneIdx> View_IdxToK;
};

template<typename T>
struct ViewPolicy : public FixedViewPolicy {
  typedef TView3d<double, typename T::Layout_Psi, IDirection, IGroup, IZone> View_Psi;
  typedef TView4d<double, typename T::Layout_FaceI, IDirection, IGroup, IZoneJ, IZoneK> View_FaceI;
  typedef TView4d<double, typename T::Layout_FaceJ, IDirection, IGroup, IZoneI, IZoneK> View_FaceJ;
  typedef TView4d<double, typename T::Layout_FaceK, IDirection, IGroup, IZoneI, IZoneJ> View_FaceK;
  typedef TView3d<double, typename T::Layout_Phi, IMoment, IGlobalGroup, IZone> View_Phi;
  typedef TView2d<double, typename T::Layout_Ell, IDirection, IMoment> View_Ell;
  typedef TView2d<double, typename T::Layout_EllPlus, IDirection, IMoment> View_EllPlus;
  typedef TView4d<double, typename T::Layout_SigS, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> View_SigS;
  typedef TView2d<double, typename T::Layout_EllPlus, IGroup, IZone> View_SigT;
};

template<typename T>
struct DataPolicy : public LayoutPolicy<T>, public ViewPolicy<LayoutPolicy<T> > {

};

#endif
