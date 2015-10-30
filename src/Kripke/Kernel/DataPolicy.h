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
#include<Kripke/DView.h>
#include<Domain/Layout.h>
#include<Domain/Forall.h>
#include<Domain/Index.h>


/*
 * Define strongly-typed indices used in Kripke
 */
DEF_INDEX(IMaterial);     // Material ID
DEF_INDEX(ILegendre);     // Legendre expansion coefficient
DEF_INDEX(IMoment);       // Spherical harmonic moment
DEF_INDEX(IDirection);    // Local direction
DEF_INDEX(IGlobalGroup);  // Global energy group
DEF_INDEX(IGroup);        // Local energy group
DEF_INDEX(IZone);         // Cannonical zone number
DEF_INDEX(IZoneIdx);      // Mapped zone index (sequential in hyperplane)
DEF_INDEX(IMix);          // Mixed element slot
DEF_INDEX(IZoneI);        // zone on the I boundary face
DEF_INDEX(IZoneJ);        // zone on the K boundary face
DEF_INDEX(IZoneK);        // zone on the K boundary face



/**
 * Layout policies that don't change with nesting.
 */
struct FixedLayoutPolicy {
  typedef DLayout2d<PERM_JI, int, IDirection, IMoment> Layout_Ell;
  typedef DLayout2d<PERM_IJ, int, IDirection, IMoment> Layout_EllPlus;

  typedef DLayout3d<PERM_KJI, IZone, IZoneI, IZoneJ, IZoneK> TLayout_Zone;
};


/**
 * Layout policies tied directly to nesting.
 */
template<typename T>
struct NestingPolicy{};

template<>
struct NestingPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_IJK, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_IJK, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_IJKL, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_IJ, int, IGroup, IZone> Layout_SigT;
  
  typedef DLayout4d<PERM_IJLK, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_IJLK, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_IJLK, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_IKJ, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_IKJ, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_ILJK, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_JI, int, IGroup, IZone> Layout_SigT;

  typedef DLayout4d<PERM_ILKJ, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_ILKJ, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_ILKJ, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_JIK, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_JIK, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_JKIL, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_IJ, int, IGroup, IZone> Layout_SigT;

  typedef DLayout4d<PERM_JILK, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_JILK, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_JILK, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_JKI, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_JKI, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_JKLI, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_IJ, int, IGroup, IZone> Layout_SigT;

  typedef DLayout4d<PERM_JLKI, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_JLKI, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_JLKI, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_KIJ, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_KIJ, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_LIJK, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_JI, int, IGroup, IZone> Layout_SigT;

  typedef DLayout4d<PERM_LKIJ, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_LKIJ, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_LKIJ, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef DLayout3d<PERM_KJI, int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout3d<PERM_KJI, int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout4d<PERM_LJKI, int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout2d<PERM_JI, int, IGroup, IZone> Layout_SigT;

  typedef DLayout4d<PERM_LKJI, int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout4d<PERM_LKJI, int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout4d<PERM_LKJI, int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};


/**
 * Views that have fixed policies
 */
struct FixedViewPolicy {
  typedef DView1d<double, DLayout1d<PERM_I, int, IZoneI> >View_dx;
  typedef DView1d<double, DLayout1d<PERM_I, int, IZoneJ> > View_dy;
  typedef DView1d<double, DLayout1d<PERM_I, int, IZoneK> > View_dz;
  typedef DView1d<Directions, DLayout1d<PERM_I, int, IDirection> > View_Directions;
  
  typedef DView1d<IZoneI, DLayout1d<PERM_I, int, IZoneIdx> > View_IdxToI;
  typedef DView1d<IZoneJ, DLayout1d<PERM_I, int, IZoneIdx> > View_IdxToJ;
  typedef DView1d<IZoneK, DLayout1d<PERM_I, int, IZoneIdx> > View_IdxToK;

  typedef DView1d<IZone, DLayout1d<PERM_I, int, IMix> > View_MixedToZones;
  typedef DView1d<IMaterial, DLayout1d<PERM_I, int, IMix> > View_MixedToMaterial;
  typedef DView1d<double, DLayout1d<PERM_I, int, IMix> > View_MixedToFraction;
  typedef DView1d<ILegendre, DLayout1d<PERM_I, int, IMoment> > View_MomentToCoeff;
};

/**
 * Views with policies that vary between nestings.
 */
template<typename T>
struct ViewPolicy : public FixedViewPolicy {
  // Discrete and Moment Unknowns
  typedef DView3d<double, typename T::Layout_Psi> View_Psi;
  typedef DView3d<double, typename T::Layout_Phi> View_Phi;

  // Spatial domain face indices
  typedef DView4d<double, typename T::Layout_FaceI> View_FaceI;
  typedef DView4d<double, typename T::Layout_FaceJ> View_FaceJ;
  typedef DView4d<double, typename T::Layout_FaceK> View_FaceK;

  // L and L+ matrices
  typedef DView2d<double, typename T::Layout_Ell> View_Ell;
  typedef DView2d<double, typename T::Layout_EllPlus> View_EllPlus;

  // Data tables
  typedef DView4d<double, typename T::Layout_SigS> View_SigS;
  typedef DView2d<double, typename T::Layout_SigT> View_SigT;
};



/**
 * Combined Policies for Layouts, Views.
 *
 * A convenience class: makes it easier to include in application.
 */
struct FixedDataPolicy {
  static const int memory_alignment = 64;
};

template<typename T>
struct DataPolicy : public FixedDataPolicy, public NestingPolicy<T>, public ViewPolicy<NestingPolicy<T> >
{
};

#endif
