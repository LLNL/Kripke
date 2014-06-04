/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

#include <Kripke.h>

struct User_Data;
struct Grid_Data;
struct SubTVec;
struct Group_Dir_Set;

/**
 * This is the Kernel base-class and interface definition.
 * This abstracts the storage of Psi, Phi, L, L+ from the rest of the code,
 * providing data-layout specific routines.
 */
class Kernel {
  public:
    Kernel();
    virtual ~Kernel();

    virtual Nesting_Order nestingPsi(void) const = 0;
    virtual Nesting_Order nestingPhi(void) const = 0;

    // Variable Creation
    void allocateStorage(User_Data *user_data);

    // Computational Kernels
    virtual void LTimes(Grid_Data *grid_data) = 0;
    virtual void LPlusTimes(Grid_Data *grid_data) = 0;
    virtual void sweep(Grid_Data *grid_data, Group_Dir_Set *ga_set, double *i_plane_ptr, double *j_plane_ptr, double *k_plane_ptr) = 0;
};


// Factory to create correct kernel object
Kernel *createKernel(Nesting_Order, int num_dims);

#endif
