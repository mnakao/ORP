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
#define DEFAULT_SEED 0
#define DEFAULT_NCALCS 10000
extern double calc_max_temp_s(const int hosts, const int switches, const int radix, const int seed, const int symmetries, const bool assign_evenly);
extern double calc_min_temp_s();

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

bool accept_s(const int hosts, const int switches, const int current_diameter, const int diameter, const double current_ASPL,
              const double ASPL, const double temp, const bool ASPL_priority, const int symmetries)
{
  if(diameter < current_diameter && !ASPL_priority){
    return true;
  }
  else if(diameter > current_diameter && !ASPL_priority){
    return false;
  }
  else{ //  diameter == current_diameter
    if(ASPL <= current_ASPL){
      return true;
    }
    else{
      double diff = ((current_ASPL-ASPL)*switches*(hosts-1)/symmetries);
      return exp(diff/temp) > uniform_rand();
    }
  }
}

static double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

static void print_help(char *argv)
{
  fprintf(stderr, "%s -H hosts -S switches -R radix -G symmetries [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-A] [-E]\n", argv);
  fprintf(stderr, "  -H : Number of hosts (Required when -f option is not specified)\n");
  fprintf(stderr, "  -S : Number of switches (Set automatically when -f and -S are not specified).\n");
  fprintf(stderr, "  -R : Radix (Required when -f is not specified)\n");
  fprintf(stderr, "  -G : Numer of symmetries (Default: 1)\n");
  fprintf(stderr, "  -f : Input file\n");
  fprintf(stderr, "  -o : Output file\n");
  fprintf(stderr, "  -s : Random seed (Default: %d)\n", DEFAULT_SEED);
  fprintf(stderr, "  -n : Number of calculations (Default: %d)\n", DEFAULT_NCALCS);
  fprintf(stderr, "  -w : Max temperature\n");
  fprintf(stderr, "  -c : Min temperature\n");
  fprintf(stderr, "  -A : ASPL takes precedence over Diameter\n");
  fprintf(stderr, "  -E : Assign hosts evenly to switches\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, int *symmetries, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *assign_evenly, bool *ASPL_priority)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:G:f:o:s:n:w:c:AE"))!=-1){
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
    case 'G':
      *symmetries = atoi(optarg);
      if(*symmetries <= 0)
        ERROR("-G value > 0\n");
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
    case 'w':
      *max_temp = atof(optarg);
      if(*max_temp <= 0)
        ERROR("-w value > 0\n");
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("MIN value > 0\n");
      break;
    case 'A':
      *ASPL_priority = true;
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
  bool ASPL_priority = false, assign_evenly = false;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, seed = DEFAULT_SEED;
  int symmetries = 1, based_switches = NOT_DEFINED;
  int lines, diameter, current_diameter, best_diameter, low_diameter;
  int (*edge)[2], *h_degree, *s_degree;
  long sum, best_sum, ncalcs = DEFAULT_NCALCS;
  double max_temp = NOT_DEFINED, min_temp = NOT_DEFINED, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ORP_Restore r;

  set_args(argc, argv, &hosts, &switches, &radix, &symmetries, &infname, &outfname, &seed,
	   &ncalcs, &max_temp, &min_temp, &assign_evenly, &ASPL_priority);

  ORP_Srand(seed);
  if(infname){
    if(hosts != NOT_DEFINED || radix != NOT_DEFINED || switches != NOT_DEFINED)
      ERROR("When using -f option, you cannot use -H, -R, and -S.\n");
    ORP_Read_property(infname, &hosts, &switches, &radix, &lines);
    if(hosts%symmetries != 0 || switches%symmetries != 0)
      ERROR("hosts and switches must be even numbers\n ");
    
    based_switches = switches/symmetries;
    h_degree = malloc(sizeof(int) * based_switches);
    s_degree = malloc(sizeof(int) * based_switches);
    edge     = malloc(sizeof(int) * lines * 2);
    ORP_Read_edge(infname, edge);
    ORP_Set_degrees_s(hosts, switches, lines, edge, symmetries, h_degree, s_degree);
  }
  else{
    if(hosts == NOT_DEFINED || radix == NOT_DEFINED)
      print_help(argv[0]);
    else if(switches == NOT_DEFINED)
      switches = ORP_Optimize_switches(hosts, radix);

    if(hosts%symmetries != 0 || switches%symmetries != 0)
      ERROR("hosts and switches must be divisible by symmetries\n ");
    based_switches = switches/symmetries;
    h_degree = malloc(sizeof(int) * based_switches);
    s_degree = malloc(sizeof(int) * based_switches);
    edge     = ORP_Generate_random_s(hosts, switches, radix, true, symmetries, &lines, h_degree, s_degree);
  }

  if(max_temp == NOT_DEFINED)
    max_temp = calc_max_temp_s(hosts, switches, radix, seed, symmetries, assign_evenly);
  
  if(min_temp == NOT_DEFINED)
    min_temp = calc_min_temp_s();
  
  printf("Hosts = %d, Switches = %d, Radix = %d, Symmetries = %d\n", hosts, switches, radix, symmetries);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
  if(infname)
    printf("Input file name = %s\n", infname);
  if(outfname)
    printf("Output file name = %s\n", outfname);

  int (*adjacency)[radix]      = malloc(sizeof(int) * based_switches * radix);
  int (*best_adjacency)[radix] = malloc(sizeof(int) * based_switches * radix);
  int *best_h_degree           = malloc(sizeof(int) * based_switches);
  int *best_s_degree           = malloc(sizeof(int) * based_switches);

  ORP_Conv_edge2adjacency_s(hosts, switches, radix, lines, edge, symmetries, adjacency);
  ORP_Init_aspl_s(hosts, switches, radix, symmetries);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);

  best_diameter = current_diameter = diameter;
  best_sum      = sum;
  best_ASPL     = current_ASPL = ASPL;
  memcpy(best_adjacency, adjacency, sizeof(int) * based_switches * radix);
  memcpy(best_h_degree,  h_degree,  sizeof(int) * based_switches);
  memcpy(best_s_degree,  s_degree,  sizeof(int) * based_switches);

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    long interval = (ncalcs < 100)? 1 : ncalcs/100;
    printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    for(long i=0;i<ncalcs;i++){
      if(i%interval == 0)
        printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
               current_ASPL-low_ASPL, current_diameter-low_diameter,
               best_ASPL-low_ASPL, best_diameter-low_diameter);

      if(assign_evenly || uniform_rand() > 0.5)
        ORP_Swap_adjacency_s(switches, radix, s_degree, symmetries, &r, adjacency);
      else
        ORP_Swing_adjacency_s(switches, radix, symmetries, h_degree, s_degree, &r, adjacency);
      
      ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
      
      if(diameter < best_diameter || (diameter == best_diameter && ASPL < best_ASPL)){
	best_diameter = diameter;
	best_sum      = sum;
	best_ASPL     = ASPL;
	memcpy(best_adjacency, adjacency, sizeof(int) * based_switches * radix);
        memcpy(best_h_degree,  h_degree,  sizeof(int) * based_switches);
        memcpy(best_s_degree,  s_degree,  sizeof(int) * based_switches);
	if(diameter == low_diameter && ASPL == low_ASPL){
	  printf("Find optimum solution\n");
	  break;
	}
      }
      
      if(accept_s(hosts, switches, current_diameter, diameter, current_ASPL, ASPL, temp, ASPL_priority, symmetries)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
        ORP_Restore_adjacency(r, radix, h_degree, s_degree, adjacency);
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ORP_Finalize_aspl();
  ORP_Conv_adjacency2edge_s(hosts, switches, radix, best_h_degree, best_s_degree, best_adjacency, symmetries, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", sa_time);
  printf("Assign evenly?  = %s\n", (assign_evenly)? "Yes" : "No");
  printf("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  printf("BIAS of host?   = %s\n", (ORP_is_bias())? "Yes" : "No");
  printf("Verify?         = %s\n", (ORP_Verify_edge(hosts, switches, radix, lines, edge))? "Yes" : "No");

  if(outfname)
    ORP_Write_edge(hosts, switches, radix, lines, edge, outfname);

  return 0;
}
