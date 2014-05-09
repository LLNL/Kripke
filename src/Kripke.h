#ifndef KRIPKE_H__
#define KRIPKE_H__
/* Prototypes */
struct User_Data;

#define KRIPKE_USE_PAPI 1
#define KRIPKE_USE_PERFTOOLS 1

/* directions.c */
void InitDirections(User_Data *grid_data, int num_directions_per_octant);

/* driver.c */
void Driver(User_Data *user_data);

/* sweep_solver.c */
int SweepSolver(User_Data *user_data);
void CreateBufferInfo(User_Data *user_data);


enum Nesting_Order {
  // Nestings for Psi and Phi
  // D referes to directions OR moments, depending on context
  NEST_GDZ,
  NEST_DGZ,
  NEST_ZDG,
  NEST_DZG,
  NEST_ZGD,
  NEST_GZD,

  // Nestings for L and L+ matrices
  NEST_DNM,
  NEST_NMD
};


#endif

