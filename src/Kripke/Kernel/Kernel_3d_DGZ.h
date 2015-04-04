#ifndef KRIPKE_KERNEL_3D_DGZ_H__
#define KRIPKE_KERNEL_3D_DGZ_H__

#include<Kripke/Kernel.h>

class Kernel_3d_DGZ : public Kernel {
  public:
    Kernel_3d_DGZ();
    virtual ~Kernel_3d_DGZ();

    virtual Nesting_Order nestingPsi(void) const;
    virtual Nesting_Order nestingPhi(void) const;
    virtual Nesting_Order nestingSigt(void) const;
    virtual Nesting_Order nestingEll(void) const;
    virtual Nesting_Order nestingEllPlus(void) const;

    virtual void LTimes(Grid_Data *grid_data);
    virtual void LPlusTimes(Grid_Data *grid_data);
    virtual void sweep(Subdomain *ga_set);
};

#endif
