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

#ifndef KERNEL_SWEEP_POLICY_H__
#define KERNEL_SWEEP_POLICY_H__

#include<Kripke.h>
#include<Domain/Layout.h>
#include<Domain/Forall.h>


template<typename T>
struct SweepPolicy{}; // d, g, z

template<>
struct SweepPolicy<NEST_DGZ_T> : ForallPolicy3<PERM_IJK, seq_pol, seq_pol, sweep_seq_pol>
{};

template<>
struct SweepPolicy<NEST_DZG_T> : ForallPolicy3<PERM_IKJ, seq_pol, seq_pol, sweep_seq_pol>
{};

template<>
struct SweepPolicy<NEST_GDZ_T> : ForallPolicy3<PERM_JIK, seq_pol, seq_pol, sweep_seq_pol>
{};

template<>
struct SweepPolicy<NEST_GZD_T> : ForallPolicy3<PERM_JKI, seq_pol, seq_pol, sweep_seq_pol>
{};

template<>
struct SweepPolicy<NEST_ZDG_T> : ForallPolicy3<PERM_KIJ, seq_pol, seq_pol, sweep_seq_pol>
{};

template<>
struct SweepPolicy<NEST_ZGD_T> : ForallPolicy3<PERM_KJI, seq_pol, seq_pol, sweep_seq_pol>
{};


#endif
