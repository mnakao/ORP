#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "orp.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static bool accept(const int hosts, const int switches, const int current_diameter, const int diameter,
                   const double current_ASPL, const double ASPL, double *max_diff_energy)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }

  //  diameter == current_diameter
  double diff = (current_ASPL-ASPL)*switches*hosts;
  *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);
    
  return (ASPL <= current_ASPL);
}

double calc_min_temp()
{
  return -2.0 / log(0.0001);
}

double calc_max_temp(const int hosts, const int switches, const int radix, const int seed)
{
  int lines, diameter, current_diameter, ncalcs = 100;
  long sum;
  double ASPL, current_ASPL, max_diff_energy = 0;
  ORP_Restore r;

  ORP_Srand(seed);
  int *h_degree  = malloc(sizeof(int) * switches);
  int *s_degree  = malloc(sizeof(int) * switches);
  int (*edge)[2] = ORP_Generate_random(hosts, switches, radix, false, &lines, h_degree, s_degree);

  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);
  
  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  current_diameter = diameter;
  current_ASPL     = ASPL;

  for(int i=0;i<ncalcs;i++){
    if(uniform_rand() > 0.5)
      ORP_Swap_adjacency(switches, radix, s_degree, &r, adjacency);
    else
      ORP_Swing_adjacency(switches, radix, h_degree, s_degree, &r, adjacency);
    
    ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
    
    if(accept(hosts, switches, current_diameter, diameter, current_ASPL, ASPL, &max_diff_energy)){
      current_diameter = diameter;
      current_ASPL     = ASPL;
    }
    else{
      ORP_Restore_adjacency(r, radix, h_degree, s_degree, adjacency);
    }
  }
  ORP_Finalize_aspl();

  return (-1.0 * max_diff_energy) / log(0.5);
}
