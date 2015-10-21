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
#include<Domain/TView.h>

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
  typedef TLayout2d<int, PERM_JI, IDirection, IMoment> Layout_Ell;
  typedef TLayout2d<int, PERM_IJ, IDirection, IMoment> Layout_EllPlus;

  typedef TLayout3d<IZone, PERM_KJI, IZoneI, IZoneJ, IZoneK> TLayout_Zone;
};


/**
 * Layout policies tied directly to nesting.
 */
template<typename T>
struct NestingPolicy{};

template<>
struct NestingPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_IJK, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_IJK, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_IJKL, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_IJ, IGroup, IZone> Layout_SigT;
  
  typedef TLayout4d<int, PERM_IJLK, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_IJLK, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_IJLK, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_IKJ, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_IKJ, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_ILJK, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_JI, IGroup, IZone> Layout_SigT;

  typedef TLayout4d<int, PERM_ILKJ, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_ILKJ, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_ILKJ, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_JIK, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_JIK, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_JKIL, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_IJ, IGroup, IZone> Layout_SigT;

  typedef TLayout4d<int, PERM_JILK, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_JILK, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_JILK, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_JKI, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_JKI, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_JKLI, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_IJ, IGroup, IZone> Layout_SigT;

  typedef TLayout4d<int, PERM_JLKI, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_JLKI, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_JLKI, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_KIJ, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_KIJ, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_LIJK, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_JI, IGroup, IZone> Layout_SigT;

  typedef TLayout4d<int, PERM_LKIJ, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_LKIJ, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_LKIJ, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  typedef TLayout3d<int, PERM_KJI, IDirection, IGroup, IZone>    Layout_Psi;
  typedef TLayout3d<int, PERM_KJI, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef TLayout4d<int, PERM_LJKI, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef TLayout2d<int, PERM_JI, IGroup, IZone> Layout_SigT;

  typedef TLayout4d<int, PERM_LKJI, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef TLayout4d<int, PERM_LKJI, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef TLayout4d<int, PERM_LKJI, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};


/**
 * Views that have fixed policies
 */
struct FixedViewPolicy {
  typedef TView1d<double, TLayout1d<int, PERM_I, IZoneI> >View_dx;
  typedef TView1d<double, TLayout1d<int, PERM_I, IZoneJ> > View_dy;
  typedef TView1d<double, TLayout1d<int, PERM_I, IZoneK> > View_dz;
  typedef TView1d<Directions, TLayout1d<int, PERM_I, IDirection> > View_Directions;
  
  typedef TView1d<IZoneI, TLayout1d<int, PERM_I, IZoneIdx> > View_IdxToI;
  typedef TView1d<IZoneJ, TLayout1d<int, PERM_I, IZoneIdx> > View_IdxToJ;
  typedef TView1d<IZoneK, TLayout1d<int, PERM_I, IZoneIdx> > View_IdxToK;

  typedef TView1d<IZone, TLayout1d<int, PERM_I, IMix> > View_MixedToZones;
  typedef TView1d<IMaterial, TLayout1d<int, PERM_I, IMix> > View_MixedToMaterial;
  typedef TView1d<double, TLayout1d<int, PERM_I, IMix> > View_MixedToFraction;
  typedef TView1d<ILegendre, TLayout1d<int, PERM_I, IMoment> > View_MomentToCoeff;
};

/**
 * Views with policies that vary between nestings.
 */
template<typename T>
struct ViewPolicy : public FixedViewPolicy {
  // Discrete and Moment Unknowns
  typedef TView3d<double, typename T::Layout_Psi> View_Psi;
  typedef TView3d<double, typename T::Layout_Phi> View_Phi;

  // Spatial domain face indices
  typedef TView4d<double, typename T::Layout_FaceI> View_FaceI;
  typedef TView4d<double, typename T::Layout_FaceJ> View_FaceJ;
  typedef TView4d<double, typename T::Layout_FaceK> View_FaceK;

  // L and L+ matrices
  typedef TView2d<double, typename T::Layout_Ell> View_Ell;
  typedef TView2d<double, typename T::Layout_EllPlus> View_EllPlus;

  // Data tables
  typedef TView4d<double, typename T::Layout_SigS> View_SigS;
  typedef TView2d<double, typename T::Layout_SigT> View_SigT;
};



/**
 * Combined Policies for Layouts, Views.
 *
 * A convenience class: makes it easier to include in application.
 */
template<typename T>
struct DataPolicy : public NestingPolicy<T>, public ViewPolicy<NestingPolicy<T> >
{};

#endif
