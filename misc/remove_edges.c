#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "orp.h"

int main(int argc, char **argv)
{
  int hosts = -1, switches = -1, radix = -1, lines = -1, diameter = -1, low_diameter = -1, new_diameter = -1, num = 0;
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

  // Remove self loop
  for(int i=0;i<lines;i++){
    if(edge[i][0] >= hosts && edge[i][1] >= hosts){
      if(edge[i][0] == edge[i][1]){
        int u = edge[i][0], v = edge[i][1];
        edge[i][0] = edge[lines-1][0];
        edge[i][1] = edge[lines-1][1];
        printf("Remove self loop        %d -- %d\n", u, v);
        i--;
        lines--;
      }
    }
  }

  // Remove multiple edges
  for(int i=0;i<lines;i++){
    if(edge[i][0] >= hosts && edge[i][1] >= hosts){
      for(int j=i+1;j<lines;j++){
        if(edge[j][0] >= hosts && edge[j][1] >= hosts){
          if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) || (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0])){
            int u = edge[i][0], v = edge[i][1];
            edge[j][0] = edge[lines-1][0];
            edge[j][1] = edge[lines-1][1];
            printf("Remove multiple edges   %d -- %d\n", u, v);
            j--;
            lines--;
          }
        }
      }
    }
  }
  
  // Remove unnecessary edge
  for(int i=0;i<lines;i++){
    if(edge[i][0] >= hosts && edge[i][1] >= hosts){
      int u = edge[i][0], v = edge[i][1];
      edge[i][0] = edge[lines-1][0];
      edge[i][1] = edge[lines-1][1];
      ORP_Set_degrees(hosts, switches, lines-1, edge, h_degree, s_degree);
      ORP_Conv_edge2adjacency(hosts, switches, radix, lines-1, edge, adjacency);
      ORP_Set_aspl(h_degree, s_degree, adjacency, &new_diameter, &new_sum, &new_ASPL);
      if(new_ASPL <= ASPL && new_diameter == diameter){
        char fname[100];
        printf("Remove unnecessary edge %d -- %d : Created %d.edges\n", u, v, num);
        sprintf(fname, "%d.edges", num++);
        ORP_Write_edge(hosts, switches, radix, lines-1, (void *)edge, fname);
      }

      edge[i][0] = u;
      edge[i][1] = v;
    }
  }

  ORP_Finalize_aspl();
 
  return 0;
}
