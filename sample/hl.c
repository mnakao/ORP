#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept(const int current_diameter, const int diameter, const double current_ASPL, const double ASPL)
{
  if(diameter < current_diameter){
    return true;
  }
  else if(diameter > current_diameter){
    return false;
  }
  else
    return (ASPL <= current_ASPL);
}

static double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

static void print_help(char *argv)
{
  ERROR("%s [-H hosts] [-S switches] [-R radix] [-f input] [-o output] [-s seed] [-n calcs] [-E]\n", argv);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, char **infname, char **outfname,
		     int *seed, long *ncalcs, bool *assign_evenly)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:f:o:s:n:E"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts <= 3)
        ERROR("-H value > 3\n");
      break;
    case 'S':
      *switches = atoi(optarg);
      if(*switches <= 2)
        ERROR("-S value > 2\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix <= 3)
        ERROR("-R value > 3\n");
      break;
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
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atol(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
      break;
    case 'E':
      *assign_evenly = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool assign_evenly = false;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, seed = 0;
  int lines, diameter, current_diameter, best_diameter, low_diameter;
  int (*edge)[2], *h_degree, *s_degree;
  long sum, best_sum, ncalcs = 10000;
  double ASPL, current_ASPL, best_ASPL, low_ASPL;
  ORP_Restore r;

  set_args(argc, argv, &hosts, &switches, &radix, &infname, &outfname, &seed, &ncalcs, &assign_evenly);
  
  ORP_Srand(seed);
  if(infname){
    ORP_Read_property(infname, &hosts, &switches, &radix, &lines);
    h_degree = malloc(sizeof(int) * switches);
    s_degree = malloc(sizeof(int) * switches);
    edge     = malloc(sizeof(int) * lines * 2);
    ORP_Read_edge(infname, edge);
    ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  }
  else{
    if(hosts == NOT_DEFINED || switches == NOT_DEFINED || radix == NOT_DEFINED)
      print_help(argv[0]);
    h_degree = malloc(sizeof(int) * switches);
    s_degree = malloc(sizeof(int) * switches);
    edge     = ORP_Generate_random(hosts, switches, radix, assign_evenly, &lines, h_degree, s_degree);
  }
  
  printf("Hosts = %d, Switches = %d, Degrees = %d\n", hosts, switches, radix);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);

  int (*adjacency)[radix]      = malloc(sizeof(int) * switches * radix);
  int (*best_adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  int *best_h_degree           = malloc(sizeof(int) * switches);
  int *best_s_degree           = malloc(sizeof(int) * switches);

  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);
  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int) * switches * radix);
  memcpy(best_h_degree,  h_degree,  sizeof(int) * switches);
  memcpy(best_s_degree,  s_degree,  sizeof(int) * switches);

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  double hl_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    int	interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Diameter Gap : ASPL Gap\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0)
	printf("%ld\t%d\t%f\n", i, best_diameter-low_diameter, best_ASPL-low_ASPL);

      if(assign_evenly || uniform_rand() > 0.5)
        ORP_Swap_adjacency(switches, radix, s_degree, &r, adjacency);
      else
        ORP_Swing_adjacency(switches, radix, h_degree, s_degree, &r, adjacency);
      
      ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
      
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int) * switches * radix);
        memcpy(best_h_degree,  h_degree,  sizeof(int) * switches);
        memcpy(best_s_degree,  s_degree,  sizeof(int) * switches);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept(current_diameter, diameter, current_ASPL, ASPL)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
        ORP_Restore_adjacency(r, radix, h_degree, s_degree, adjacency);
      }
    }
  }
  hl_time = get_time() - hl_time;
  ORP_Finalize_aspl();
  ORP_Conv_adjacency2edge(hosts, switches, radix, best_h_degree, best_s_degree, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", hl_time);
  printf("Assign evenly?  = %s\n", (assign_evenly)? "Yes" : "No");

  if(outfname)
    ORP_Write_edge(hosts, switches, radix, lines, edge, outfname);

  return 0;
}
