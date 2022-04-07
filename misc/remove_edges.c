#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "orp.h"
#define NOT_DEFINED -1

static bool is_edge_in_switches(const int i, const int hosts, const int (*edge)[2])
{
  return (edge[i][0] >= hosts && edge[i][1] >= hosts);
}

int main(int argc, char **argv)
{
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, lines = NOT_DEFINED;
  int diameter = NOT_DEFINED, low_diameter = NOT_DEFINED, new_diameter = NOT_DEFINED;
  long sum = NOT_DEFINED, new_sum = NOT_DEFINED;
  double ASPL = NOT_DEFINED, low_ASPL = NOT_DEFINED, new_ASPL = NOT_DEFINED;
  char fname[100];

  if(argc != 2){
    fprintf(stderr, "Please specify an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
    
  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2); // int edge[lines][2];
  ORP_Read_edge(argv[1], edge);
  
  int *s_degree = malloc(sizeof(int) * switches);
  int *h_degree = malloc(sizeof(int) * switches);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix); // int adjacency[switches][radix];
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);

  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);

  // Remove self-loop
  bool flag = false;
  for(int i=0;i<lines;i++){
    if(is_edge_in_switches(i, hosts, edge)){
      if(edge[i][0] == edge[i][1]){
        flag = true;
        printf("Remove self-loop        %d -- %d\n", edge[i][0], edge[i][1]);
        edge[i][0] = edge[lines-1][0];
        edge[i][1] = edge[lines-1][1];
        lines--;
        i--;  // Since edge[lines-1][0:2] may have a self-loop, return i by one.
      }
    }
  }

  // Remove multiple edges
  for(int i=0;i<lines;i++){
    if(is_edge_in_switches(i, hosts, edge)){
      for(int j=i+1;j<lines;j++){
        if(is_edge_in_switches(j, hosts, edge)){
          if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) || (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0])){
            flag = true;
            printf("Remove multiple edges   %d -- %d\n", edge[j][0], edge[j][1]);
            edge[j][0] = edge[lines-1][0];
            edge[j][1] = edge[lines-1][1];
            lines--;
            j--; // Since edge[i][0:2] and edge[lines-1][0:2] may be multiple edges, return j by one.
          }
        }
      }
    }
  }

  int num = 0;
  if(flag){
    printf("Created %d.edges\n", num);
    sprintf(fname, "%d.edges", num++);
    ORP_Write_edge(hosts, switches, radix, lines, (void *)edge, fname);
  }
  
  // Remove unnecessary edge
  for(int i=0;i<lines;i++){
    if(is_edge_in_switches(i, hosts, edge)){
      int u = edge[i][0], v = edge[i][1]; // backup
      edge[i][0] = edge[lines-1][0];
      edge[i][1] = edge[lines-1][1];
      ORP_Set_degrees(hosts, switches, lines-1, edge, h_degree, s_degree);
      ORP_Conv_edge2adjacency(hosts, switches, radix, lines-1, edge, adjacency);
      ORP_Set_aspl(h_degree, s_degree, adjacency, &new_diameter, &new_sum, &new_ASPL);
      if(new_ASPL == ASPL && new_diameter == diameter){
        printf("Remove unnecessary edge %d -- %d : Created %d.edges\n", u, v, num);
        sprintf(fname, "%d.edges", num++);
        ORP_Write_edge(hosts, switches, radix, lines-1, (void *)edge, fname);
      }
      edge[i][0] = u; edge[i][1] = v; // restore
    }
  }

  ORP_Finalize_aspl();

  free(edge);
  free(s_degree);
  free(h_degree);
  
  return 0;
}
