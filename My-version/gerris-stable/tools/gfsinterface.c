/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2009 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <stdlib.h>
#include "config.h"
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include "init.h"
#include "simulation.h"
#include "vof.h"

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi_boundary.h"
#endif /* HAVE_MPI */


static void interface_pos(FttCell *cell, GfsVariableTracerVOF *t)
{ 
  FttDirection d;
  int count = 0;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    double voidf = 0.;
    if(neighbor) voidf = GFS_VALUE(neighbor,t->alpha);
    if(voidf*(1.-voidf) > 1.E-8 ) count++;
  }
 
  if(count>1){
    FttVector p[2], m;
    if (gfs_vof_facet (cell, t, p, &m) == 2) {
      FttVector p1 = p[0];
      FttVector p2 = p[1];
      gfs_simulation_map_inverse (gfs_object_simulation(t), &p1);
      gfs_simulation_map_inverse (gfs_object_simulation(t), &p2);
      printf("%lf %lf\n",(p1.x + p2.x)/2.0,(p1.y + p2.y)/2.);
    }
  }

}

static void gfs_interface (GfsSimulation * sim)
{
  g_return_if_fail (sim != NULL);
  FttVector p[2];
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsVariableTracerVOF *t = NULL;

  if (domain->variables) {
    GfsVariable * v = NULL;
    GSList * i = domain->variables;
    
    while (i && !v) {
      if (GFS_IS_VARIABLE_TRACER_VOF (i->data))
	v = i->data;
      i = i->next;
    }
    t = GFS_VARIABLE_TRACER_VOF(v);
 
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) interface_pos, t);

  }
}

int main (int argc, char * argv[])
{

  int c = 0;
  gdouble time;
  gint np;
  GtsFile * fp;
  FILE * fptr;
  gint maxlevel = -2;

  gfs_init (&argc, &argv);

    /*ParentFile*/
  int optind = 1;

    if (optind >= argc) { /* missing FILE */
        fprintf (stderr,
                "gerris: missing FILE\n"
                "Try `gerris --help' for more information.\n");
        return 1; /* failure */
    }
    
    fptr = fopen (argv[optind], "r");

    if (fptr == NULL) {
        fprintf (stderr, "gerris: unable to open file `%s'\n", argv[optind]);
        return 1;
    }

    fp = gts_file_new (fptr);
    GfsSimulation * simulation;

    if (!(simulation = gfs_simulation_read (fp))) {
        fprintf (stderr,
                "gerris: file `%s' is not a valid simulation file\n"
                "%s:%d:%d: %s\n",
                argv[optind], argv[optind],
                fp->line, fp->pos, fp->error);
        return 1;
    } 


    gfs_interface (simulation);
    gts_file_destroy (fp);
    fclose(fptr);
}

