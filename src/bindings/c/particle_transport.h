#ifndef AETHER_PARTICLE_TRANSPORT_H
#define AETHER_PARTICLE_TRANSPORT_H

#include "symbol_export.h"

SHO_PUBLIC void solve_particles_decoupled(double initial_concentration, double inlet_concentration, double particle_size);
//SHO_PUBLIC void write_airway(int ne_field, const char *EXELEMFILE, const char *group_name, const char *field_name );
//SHO_PUBLIC void write_terminal(const char *EXNODEFILE, const char *name);

#endif /* AETHER_PARTICLE_TRANSPORT_H */
