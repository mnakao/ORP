#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "orp.h"

int main(int argc, char *argv[])
{
  int hosts, switches, radix, lines, diameter, low_diameter;
  long sum;
  double ASPL, low_ASPL;
  
  if(argc != 2){
    fprintf(stderr, "%s <edgefile>\n", argv[0]);
    exit(1);
  }
  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  
  int *h_degree  = malloc(sizeof(int) * switches);
  int *s_degree  = malloc(sizeof(int) * switches);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  
  ORP_Read_edge(argv[1], edge);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);
  
  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Finalize_aspl();

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  printf("Hosts = %d, Switches %d, Radix = %d\n", hosts, switches, radix);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  printf("--\n");
  
  int *hist = malloc((radix+1) * sizeof(int));
  for(int i=0;i<radix+1;i++)  hist[i] = 0;
  for(int i=0;i<switches;i++) hist[s_degree[i]]++;

  int s = 0;
  printf("#Num\tconnected to switches\n");
  for(int i=1;i<radix+1;i++){
    s += hist[i];
    printf("%d\t%d\n", i, hist[i]);
  }
  printf("---\n");
  printf("SUM : %d\n", s);
  
  if(s != switches || hist[0] != 0){
    printf("Error\n");
    exit(1);
  }
  
  return 0;
}
