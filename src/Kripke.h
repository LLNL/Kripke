#ifndef KRIPKE_H__
#define KRIPKE_H__
/* Prototypes */
struct User_Data;

/* directions.c */
void InitDirections(User_Data *grid_data, int num_directions_per_octant);

/* driver.c */
void Driver(User_Data *user_data);

/* sweep_solver.c */
int SweepSolver(User_Data *user_data);
void CreateBufferInfo(User_Data *user_data);


#endif

