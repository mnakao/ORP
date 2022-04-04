#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -f input -o output [-G symmetries]\n", argv);
  fprintf(stderr, "  -f : Input file\n");
  fprintf(stderr, "  -o : Output file\n");
  fprintf(stderr, "  -G : Numer of symmetries (Default: 1)\n");
  exit(1);
}

static void set_args(const int argc, char **argv, char **infname, char **outfname, int *symmetries)
{
  int result;
  while((result = getopt(argc,argv,"f:o:G:"))!=-1){
    switch(result){
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      *outfname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*outfname, optarg);
      break;
    case 'G':
      *symmetries = atoi(optarg);
      if(*symmetries < 1)
        ERROR("-G value >= 1\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static int tmp[2][2];
static void backup(const int i, const int j, const int radix, const int (*adjacency)[radix])
{
  tmp[0][0] = i; tmp[0][1] = tmp[1][0] = j;
  
}

static int search_index(const int v, const int target, const int exclusion, const int *s_degree,
                        const int radix, const int (*adjacency)[radix])
{
  if(v == target){
    for(int i=0;i<s_degree[v];i++){
      if(adjacency[v][i] == target && i != exclusion){
        return i;
      }
    }
  }
  else{
    for(int i=0;i<s_degree[v];i++){
      if(adjacency[v][i] == target){
        return i;
      }
    }
  }

  ERROR("Something Wrong (id=2)\n");
  return -1; // dummy
}

int main(int argc, char **argv)
{
  char *infname = NULL, *outfname = NULL;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, lines = NOT_DEFINED;
  int diameter = NOT_DEFINED, low_diameter = NOT_DEFINED, new_diameter = NOT_DEFINED, symmetries = 1;
  long sum = NOT_DEFINED, new_sum = NOT_DEFINED;
  double ASPL = NOT_DEFINED, low_ASPL = NOT_DEFINED, new_ASPL = NOT_DEFINED;
  bool flag = false;

  set_args(argc, argv, &infname, &outfname, &symmetries);
  if(infname == NULL || outfname == NULL)
    ERROR("./%s -f input -o output\n", argv[0]);
  
  ORP_Read_property(infname, &hosts, &switches, &radix, &lines);
  printf("Hosts = %d, Switches = %d, Radix = %d\n", hosts, switches, radix);

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ORP_Read_edge(infname, edge);
  
  int *s_degree = malloc(sizeof(int) * switches);
  int *h_degree = malloc(sizeof(int) * switches);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix); // int adjacency[switches][radix];
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);

  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  printf("--\n");
  for(int u=0;u<switches/symmetries;u++){
    for(int u_d=0;u_d<s_degree[u];u_d++){
      int v   = adjacency[u][u_d];
      int v_d = search_index(v, u, u_d, s_degree, radix, adjacency);
      adjacency[u][u_d] = adjacency[u][s_degree[u]-1];
      s_degree[u]--;
      adjacency[v][v_d] = adjacency[v][s_degree[v]-1];
      s_degree[v]--;
      ORP_Set_aspl(h_degree, s_degree, adjacency, &new_diameter, &new_sum, &new_ASPL);
      if((new_ASPL <= ASPL && new_diameter <= diameter) || new_diameter < diameter){
        flag = true;
        printf("Delete : %d -- %d\n", u+hosts, v+hosts);
        printf(" Diameter Gap = %d (%d - %d)\n", new_diameter - low_diameter, new_diameter, low_diameter);
        printf(" ASPL Gap     = %.10f (%.10f - %.10f)\n", new_ASPL - low_ASPL, new_ASPL, low_ASPL);
      }
      
      adjacency[u][u_d] = v; s_degree[u]++;
      adjacency[v][v_d] = u; s_degree[v]++;
    }
  }
  ORP_Finalize_aspl();
  
  if(!flag)
    printf("Not found.\n");
  
  free(adjacency); free(h_degree);
  free(s_degree);  free(edge);
  return 0;
}
