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

#ifndef KRIPKE_VARTYPES_H__
#define KRIPKE_VARTYPES_H__

#include <Kripke.h>
#include <Kripke/Set.h>
#include <Kripke/Field.h>

namespace Kripke {

  RAJA_INDEX_VALUE(Dimension, "Dimension");
  RAJA_INDEX_VALUE(Direction, "Direction");
  RAJA_INDEX_VALUE(GlobalGroup, "GlobalGroup");
  RAJA_INDEX_VALUE(Group, "Group");
  RAJA_INDEX_VALUE(Legendre, "Legendre");
  RAJA_INDEX_VALUE(Material, "Material");
  RAJA_INDEX_VALUE(MixElem, "MixElem");
  RAJA_INDEX_VALUE(Moment, "Moment");
  RAJA_INDEX_VALUE(Zone, "Zone");
  RAJA_INDEX_VALUE(ZoneI, "ZoneI");
  RAJA_INDEX_VALUE(ZoneJ, "ZoneJ");
  RAJA_INDEX_VALUE(ZoneK, "ZoneK");

  using Field_Flux = Kripke::Field<double, Direction, Group, Zone>;
  using Field_Moments = Kripke::Field<double, Moment, Group, Zone>;

  using Field_IPlane = Kripke::Field<double, Direction, Group, ZoneJ, ZoneK>;
  using Field_JPlane = Kripke::Field<double, Direction, Group, ZoneI, ZoneK>;
  using Field_KPlane = Kripke::Field<double, Direction, Group, ZoneI, ZoneJ>;

  using Field_Ell = Kripke::Field<double, Moment, Direction>;

  using Field_Speed  = Kripke::Field<double, Material, GlobalGroup>;
  using Field_SigmaT = Kripke::Field<double, Material, GlobalGroup>;
  using Field_SigmaS = Kripke::Field<double, Material, Legendre, GlobalGroup, GlobalGroup>;

  using Field_Direction2Double = Kripke::Field<double, Direction>;
  using Field_Direction2Int    = Kripke::Field<int, Direction>;

  using Field_Adjacency        = Kripke::Field<GlobalSdomId, Dimension>;

  using Field_Moment2Legendre  = Kripke::Field<Legendre, Moment>;

  using Field_ZoneI2Double  = Kripke::Field<double, ZoneI>;
  using Field_ZoneJ2Double  = Kripke::Field<double, ZoneJ>;
  using Field_ZoneK2Double  = Kripke::Field<double, ZoneK>;
  using Field_Zone2Double   = Kripke::Field<double, Zone>;
  using Field_Zone2Int      = Kripke::Field<int, Zone>;
  using Field_Zone2MixElem  = Kripke::Field<MixElem, Zone>;

  using Field_MixElem2Double   = Kripke::Field<double, MixElem>;
  using Field_MixElem2Material = Kripke::Field<Material, MixElem>;
  using Field_MixElem2Zone     = Kripke::Field<Zone, MixElem>;

  using Field_SigmaTZonal = Kripke::Field<double, Group, Zone>;
}

#endif
