/*BHEADER**********************************************************************
 * (c) 1998   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision: 1.2 $
 *********************************************************************EHEADER*/

/*--------------------------------------------------------------------------
 * Header file for the BC_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef included_bc_data
#define included_bc_data

/*--------------------------------------------------------------------------
 * Define the BC_Data structure.
 *
 * bc_types  : The boundary condition type values for the six faces
 * bc_values : The boundary condition values
 *--------------------------------------------------------------------------*/

typedef struct {
  int bc_types[6];
  double bc_values[6];

} BC_Data;

#endif
