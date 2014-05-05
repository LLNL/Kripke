/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

enum Nesting_Order {
  NEST_GDZ,
  NEST_DGZ,
  NEST_ZDG,
  NEST_DZG,
  NEST_ZGD,
  NEST_GZD
};

struct User_Data;
struct Grid_Data;
struct SubTVec;
struct Group_Dir_Set;

class Kernel {
  public:
    Kernel();
    virtual ~Kernel();

    virtual Nesting_Order nesting(void) const = 0;

    // Variable Creation
    void allocateStorage(Grid_Data *grid_data);

    // Computational Kernels
    virtual void evalSigmaTot(User_Data *grid_data, Group_Dir_Set *ga_set) = 0;
    virtual void evalSigmaS(Grid_Data *grid_data) = 0;
    virtual void LTimes(Grid_Data *grid_data) = 0;
    virtual void LPlusTimes(Grid_Data *grid_data) = 0;
    virtual void sweep(Grid_Data *grid_data, Group_Dir_Set *ga_set, double *i_plane_ptr, double *j_plane_ptr, double *k_plane_ptr) = 0;
};


// Factory to create correct kernel object
Kernel *createKernel(Nesting_Order, int num_dims);

#endif
