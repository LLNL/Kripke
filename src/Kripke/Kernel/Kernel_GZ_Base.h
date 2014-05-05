#ifndef KRIPKE_KERNEL_GZ_BASE_H__
#define KRIPKE_KERNEL_GZ_BASE_H__

#include<Kripke/Kernel.h>

class Kernel_GZ_Base : public Kernel {
  public:
    Kernel_GZ_Base();
    virtual ~Kernel_GZ_Base();

    // Computational Kernels
    virtual void evalSigmaTot(User_Data *grid_data, Group_Dir_Set *ga_set);
};

#endif
