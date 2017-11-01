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

#ifndef KRIPKE_CORE_COMM_H__
#define KRIPKE_CORE_COMM_H__

#include <Kripke.h>
#include <Kripke/Core/BaseVar.h>

#ifdef KRIPKE_USE_MPI
#include <mpi.h>
#endif

namespace Kripke {
namespace Core {
/**
 * An interprocess communicator.
 *
 * Used as an abstraction layer around MPI... mainly to allow compilation w/o
 * MPI
 */
class Comm : public Kripke::Core::BaseVar {
  public:


#ifdef KRIPKE_USE_MPI
    RAJA_INLINE
    static void init(int *argc, char ***argv){
      MPI_Init(argc, argv);
    }
#else
    RAJA_INLINE
    static void init(int *, char ***){
    }
#endif

    RAJA_INLINE
    static void finalize() {
#ifdef KRIPKE_USE_MPI
      MPI_Finalize();
#endif
    }


    RAJA_INLINE
    static Comm getSelf() {
#ifdef KRIPKE_USE_MPI
      return Comm(MPI_COMM_SELF);
#else
      return Comm();
#endif
    }

#ifdef KRIPKE_USE_MPI
    RAJA_INLINE
    Comm() :
      m_comm(MPI_COMM_WORLD),
      m_rank(0),
      m_size(0)
    {
      int r, s;
      MPI_Comm_rank(m_comm, &r);
      MPI_Comm_size(m_comm, &s);
      m_rank = r;
      m_size = s;
    }

    RAJA_INLINE
    Comm(MPI_Comm c) :
      m_comm(c),
      m_rank(0),
      m_size(0)
    {
      int r, s;
      MPI_Comm_rank(m_comm, &r);
      MPI_Comm_size(m_comm, &s);
      m_rank = r;
      m_size = s;
    }
#else
    RAJA_INLINE
    Comm() :
      m_rank(0),
      m_size(1)
    {}
#endif


    virtual ~Comm() = default;

    RAJA_INLINE
    size_t size() const {
      return m_size; 
    }
    
    RAJA_INLINE
    size_t rank() const {
      return m_rank;
    }
    
    RAJA_INLINE
#ifdef KRIPKE_USE_MPI
    Comm split(int color, int key) const {
      MPI_Comm split_comm;
      MPI_Comm_split(m_comm, color, key, &split_comm);
      return Comm(split_comm);
#else
    Comm split(int , int ) const {
      return Comm();
#endif
    }

    /**
     * Allreduce SUM a single value.
     * Without MPI, this is a NOP
     */
    RAJA_INLINE
    long allReduceSumLong(long value) const {
#ifdef KRIPKE_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_LONG, MPI_SUM, m_comm);
#endif
      return value;
    }

    /**
     * Allreduce SUM an array, in-place
     * Without MPI, this is a NOP
     */
    RAJA_INLINE
#ifdef KRIPKE_USE_MPI
    void allReduceSumLong(long *value, size_t len) const {
      MPI_Allreduce(MPI_IN_PLACE, value, len, MPI_LONG, MPI_SUM, m_comm);
    }
#else
    void allReduceSumLong(long *, size_t ) const {}
#endif


    /**
     * Allreduce SUM an array, in-place
     * Without MPI, this is a NOP
     */
    RAJA_INLINE
#ifdef KRIPKE_USE_MPI
    void allReduceSumInt(int *value, size_t len) const {
      MPI_Allreduce(MPI_IN_PLACE, value, len, MPI_INT, MPI_SUM, m_comm);
    }
#else
    void allReduceSumInt(int *, size_t) const {}
#endif

    /**
     * Allreduce SUM a single value.
     * Without MPI, this is a NOP
     */
    RAJA_INLINE
    double allReduceSumDouble(double value) const {
#ifdef KRIPKE_USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, m_comm);
#endif
      return value;
    }

    /**
     * Allreduce SUM an array, in-place.
     * Without MPI, this is a NOP
     */
#ifdef KRIPKE_USE_MPI
    RAJA_INLINE
    void allReduceSumDouble(double *value, size_t len) const {
      MPI_Allreduce(MPI_IN_PLACE, value, len, MPI_DOUBLE, MPI_SUM, m_comm);
    }
#else
    RAJA_INLINE
    void allReduceSumDouble(double *, size_t) const {
    }
#endif

    /**
     * Prefix scan SUM a single value.
     * Without MPI, this is a NOP
     */
    RAJA_INLINE
    long scanSumLong(long value) const {
#ifdef KRIPKE_USE_MPI
      MPI_Scan(MPI_IN_PLACE, &value, 1, MPI_LONG, MPI_SUM, m_comm);
#endif
      return value;
    }

  private:
#ifdef KRIPKE_USE_MPI
    MPI_Comm m_comm;
#endif
    size_t m_rank;
    size_t m_size;
};



} } // namespace

#endif
