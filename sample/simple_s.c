#include <stdio.h>
#include <stdlib.h>
#include "orp.h"

int main(int argc, char **argv)
{
  int hosts, switches, radix, lines, diameter = -1, low_diameter = -1;
  long sum = -1;
  double ASPL = -1, low_ASPL = -1;

  if(argc != 3){
    fprintf(stderr, "./%s xxx.edges symmetries\n", argv[0]);
    exit(1);
  }

  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  int symmetries = atoi(argv[2]);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ORP_Read_edge(argv[1], edge);
  if(hosts % symmetries != 0){
    fprintf(stderr, "hosts(%d) must be divisible by symmetries(%d)\n", hosts, symmetries);
    exit(1);
  }
  else if(switches % symmetries != 0){
    fprintf(stderr, "switches(%d) must be divisible by symmetries(%d)\n", switches, symmetries);
    exit(1);
  }

  int based_switches = switches / symmetries;
  int *s_degree = malloc(sizeof(int) * based_switches);
  int *h_degree = malloc(sizeof(int) * based_switches);
  ORP_Set_host_degree_s  (hosts, switches, lines, edge, symmetries, h_degree);
  ORP_Set_switch_degree_s(hosts, switches, lines, edge, symmetries, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * based_switches * radix); // int adjacency[based_switches][radix];
  ORP_Conv_edge2adjacency_s(hosts, switches, radix, lines, edge, symmetries, adjacency);

  ORP_Init_aspl_s(hosts, switches, radix, symmetries);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Finalize_aspl();

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  printf("Hosts = %d, Switches %d, Radix = %d\n", hosts, switches, radix);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
 
  free(adjacency);
  free(h_degree);
  free(s_degree);
  free(edge);

  return 0;
}
