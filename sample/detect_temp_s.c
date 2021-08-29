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

static bool accept_s(const int hosts, const int switches, const int current_diameter, const int diameter,
                     const double current_ASPL, const double ASPL, const int symmetries, double *max_diff_energy)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }

  //  diameter == current_diameter
  double diff = (current_ASPL-ASPL)*switches*hosts/symmetries;
  *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);
    
  return (ASPL <= current_ASPL);
}

double calc_min_temp_s()
{
  return -2.0 / log(0.0001);
}

double calc_max_temp_s(const int hosts, const int switches, const int radix, const int seed, const int symmetries, const bool assign_evenly)
{
  int lines, diameter, current_diameter, ncalcs = 100, based_switches = switches/symmetries;
  long sum;
  double ASPL, current_ASPL, max_diff_energy = 0;
  ORP_Restore r;

  ORP_Srand(seed);
  int *h_degree  = malloc(sizeof(int) * based_switches);
  int *s_degree  = malloc(sizeof(int) * based_switches);
  int (*edge)[2] = ORP_Generate_random_s(hosts, switches, radix, true, symmetries, &lines, h_degree, s_degree);

  int (*adjacency)[radix] = malloc(sizeof(int) * based_switches * radix);
  ORP_Conv_edge2adjacency_s(hosts, switches, radix, lines, edge, symmetries, adjacency);

  // Save and unset value of ORP_PROFILE
  char *val = getenv("ORP_PROFILE");
  if(val) unsetenv("ORP_PROFILE");
  
  ORP_Init_aspl_s(hosts, switches, radix, symmetries);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  current_diameter = diameter;
  current_ASPL     = ASPL;

  for(int i=0;i<ncalcs;i++){
    if(assign_evenly || uniform_rand() > 0.5)
      ORP_Swap_adjacency_s(switches, radix, s_degree, symmetries, &r, adjacency);
    else
      ORP_Swing_adjacency_s(switches, radix, symmetries, h_degree, s_degree, &r, adjacency);
    
    ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
    
    if(accept_s(hosts, switches, current_diameter, diameter, current_ASPL, ASPL, symmetries, &max_diff_energy)){
      current_diameter = diameter;
      current_ASPL     = ASPL;
    }
    else{
      ORP_Restore_adjacency(r, radix, h_degree, s_degree, adjacency);
    }
  }
  ORP_Finalize_aspl();

  // Undo value of ORP_ASPL
  if(val) setenv("ORP_PROFILE", val, 1);
  
  return (-1.0 * max_diff_energy) / log(0.5);
}
