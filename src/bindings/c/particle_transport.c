
#include "particle_transport.h"

void solve_particles_decoupled_c(double *initial_concentration, double *inlet_concentration, double *particle_size);

void solve_particles_decoupled(double initial_concentration, double inlet_concentration, double particle_size)
{
  solve_particles_decoupled_c(&initial_concentration, &inlet_concentration, &particle_size);
}


//void write_airway_c(int *ne_field, const char *EXELEMFILE, int *EXELEMFILE_LEN,
//                            const char *group_name, int *group_name_len, const char *field_name, int *field_name_len );
//void write_airway(int ne_field, const char *EXELEMFILE, const char *group_name, const char *field_name )
//{
//  int filename_len = strlen(EXELEMFILE);
//  int group_name_len = strlen(group_name);
//  int field_name_len = strlen(field_name);
//
//  write_airway_c(&ne_field, EXELEMFILE, &filename_len, group_name, &group_name_len, field_name, &field_name_len);
//}
//
//void write_terminal_c(const char *EXNODEFILE, int *EXNODEFILE_LEN, const char *name, int *name_len);
//void write_terminal(const char *EXNODEFILE, const char *name)
//{
//  int filename_len = strlen(EXNODEFILE);
//  int name_len = strlen(name);
//
//  write_terminal_c(EXNODEFILE, &filename_len, name, &name_len);
//}