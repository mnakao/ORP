#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "orp.h"

int main(int argc, char **argv)
{
  int hosts, switches, radix, lines, diameter = -1, low_diameter = -1, new_diameter = -1;
  long sum = -1, new_sum = -1;
  double ASPL = -1, low_ASPL = -1, new_ASPL = -1;

  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
    
  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ORP_Read_edge(argv[1], edge);
  
  int *s_degree = malloc(sizeof(int) * switches);
  int *h_degree = malloc(sizeof(int) * switches);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix); // int adjacency[switches][radix];
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);

  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  //  printf("%f\n", ASPL - low_ASPL);

  int (*new_edge)[2] = malloc(sizeof(int) * lines * 2);
  int *new_s_degree  = malloc(sizeof(int) * switches);
  int *new_h_degree  = malloc(sizeof(int) * switches);
  int new_lines      = lines - 1;
  char fname[100];
  for(int i=hosts;i<lines;i++){
    memcpy(new_edge, edge, sizeof(int) * lines * 2);
    if(new_edge[i][0] < hosts || new_edge[i][1] < hosts) exit(1);
    new_edge[i][0] = new_edge[lines-1][0];
    new_edge[i][1] = new_edge[lines-1][1];
    ORP_Set_degrees(hosts, switches, new_lines, new_edge, new_h_degree, new_s_degree);
    ORP_Conv_edge2adjacency(hosts, switches, radix, new_lines, new_edge, adjacency);
    ORP_Set_aspl(new_h_degree, new_s_degree, adjacency, &new_diameter, &new_sum, &new_ASPL);
    if(new_ASPL <= ASPL){
      printf("Delete %d -- %d\n", edge[i][0], edge[i][1]);
      sprintf(fname, "%d.edges", i);
      ORP_Write_edge(hosts, switches, radix, new_lines, new_edge, fname);
      //      printf("Found\n");
      //      exit(0);
    }
  }
  
  ORP_Finalize_aspl();
 
  free(adjacency);
  free(h_degree);
  free(s_degree);
  free(edge);
  return 0;
}
