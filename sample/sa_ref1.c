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
#define DEFAULT_MAX_TEMP 100.00
#define DEFAULT_MIN_TEMP 0.22

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static int search_index(const int v, const int target, const int exclusion,
                        const int *s_degree, const int radix, const int (*adjacency)[radix])
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

  ERROR("Something Wrong (id=0)\n");
  return -1; // dummy
}

static void set_random_edge(const int total_edges, const int switches, const int s_degree[switches], const int radix,
                            const int (*adjacency)[radix], int *u, int *v, int *u_d, int *v_d)
{
  int r = get_random(total_edges);
  *u = NOT_DEFINED;
  for(int i=0;i<switches;i++){
    r -= s_degree[i];
    if(r < 0){
      *u = i;
      break;
    }
  }
  if(*u == NOT_DEFINED) ERROR("Something Wrong (id=1)\n");
  
  *u_d = get_random(s_degree[*u]);
  *v   = adjacency[*u][*u_d];
  *v_d = search_index(*v, *u, *u_d, s_degree, radix, adjacency);
}

static void set_random_edge_bias(const int hosts, const int switches, const int h_degree[switches], const int s_degree[switches],
                                 const int radix, const int (*adjacency)[radix], int *u, int *v, int *u_d, int *v_d)
{
  int r = get_random(hosts+switches);
  *u = NOT_DEFINED;
  for(int i=0;i<switches;i++){
    r -= (h_degree[i] + 1);
    if(r < 0){
      *u = i;
      break;
    }
  }
  if(*u == NOT_DEFINED) ERROR("Something Wrong (id=2)\n");
  
  *u_d = get_random(s_degree[*u]);
  *v   = adjacency[*u][*u_d];
  *v_d = search_index(*v, *u, *u_d, s_degree, radix, adjacency);
}

bool accept(const int switches, const int current_diameter, const int diameter, const double current_ASPL,
            const double ASPL, const double temp, const bool ASPL_priority)
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
      double diff = ((current_ASPL-ASPL)*switches*switches);
      return (exp(diff/temp) > uniform_rand());
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
  fprintf(stderr, "%s [-H hosts] [-S switches] [-R radix] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-A] [-B]\n", argv);
  fprintf(stderr, "  -H : Number of hosts (Required when -f option is not specified)\n");
  fprintf(stderr, "  -S : Number of switches (Set automatically when -f and -S are not specified).\n");
  fprintf(stderr, "  -R : Radix (Required when -f is not specified)\n");
  fprintf(stderr, "  -f : Input file\n");
  fprintf(stderr, "  -o : Output file\n");
  fprintf(stderr, "  -s : Random seed (Default: %d)\n", DEFAULT_SEED);
  fprintf(stderr, "  -n : Number of calculations (Default: %d)\n", DEFAULT_NCALCS);
  fprintf(stderr, "  -w : Max temperature (Default: %.2f)\n", (double)DEFAULT_MAX_TEMP);
  fprintf(stderr, "  -c : Min temperature (Default: %.2f)\n", (double)DEFAULT_MIN_TEMP);
  fprintf(stderr, "  -A : ASPL takes precedence over Diameter\n");
  fprintf(stderr, "  -B : Increase the bias in the number of hosts\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, char **infname, char **outfname,
		     int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *ASPL_priority, bool *bias_of_host)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:f:o:s:n:w:c:AB"))!=-1){
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
    case 'B':
      *bias_of_host = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool ASPL_priority = false, bias_of_host = false;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, seed = DEFAULT_SEED;
  int lines, diameter, current_diameter, best_diameter, low_diameter;
  int (*edge)[2], *h_degree, *s_degree;
  long sum, best_sum, ncalcs = DEFAULT_NCALCS;
  double max_temp = DEFAULT_MAX_TEMP, min_temp = DEFAULT_MIN_TEMP, ASPL, current_ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &hosts, &switches, &radix, &infname, &outfname, &seed,
           &ncalcs, &max_temp, &min_temp, &ASPL_priority, &bias_of_host);
  
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
    if(hosts == NOT_DEFINED || radix == NOT_DEFINED)
      print_help(argv[0]);

    if(switches == NOT_DEFINED)
      switches = ORP_Optimize_switches(hosts, radix);

    h_degree = malloc(sizeof(int) * switches);
    s_degree = malloc(sizeof(int) * switches);
    edge     = ORP_Generate_random(hosts, switches, radix, false, &lines, h_degree, s_degree);
  }
    
  printf("Hosts = %d, Switches = %d, Radix = %d\n", hosts, switches, radix);
  printf("Random seed = %d\n", seed);
  printf("Number of calculations = %ld\n", ncalcs);
  printf("Max, Min temperature = %f, %f\n", max_temp, min_temp);
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
  double sa_time = get_time();
  if(diameter == low_diameter && ASPL == low_ASPL){
    printf("Find optimum solution\n");
  }
  else{
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs);
    double temp = max_temp;
    int total_edges = 0;
    for(int i=0;i<switches;i++) total_edges += s_degree[i];
    long interval = (ncalcs < 100)? 1 : ncalcs/100;
    long j = 0;
    printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    for(long i=0;i<ncalcs;i++){
      if(i/interval == j){
        printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
               current_ASPL-low_ASPL, current_diameter-low_diameter,
               best_ASPL-low_ASPL, best_diameter-low_diameter);
        j++;
      }

      bool enable_swing = true;
      int u[2], v[2], u_d[2], v_d[2];
      while(1){
        if(bias_of_host)
          set_random_edge_bias(hosts, switches, h_degree, s_degree, radix, adjacency, &u[0], &v[0], &u_d[0], &v_d[0]);
        else
          set_random_edge(total_edges, switches, s_degree, radix, adjacency, &u[0], &v[0], &u_d[0], &v_d[0]);
        set_random_edge(total_edges, switches, s_degree, radix, adjacency, &u[1], &v[1], &u_d[1], &v_d[1]);
        if(u[0] == u[1] || s_degree[u[0]] == 1) continue;
        //        if(v[0] == u[1]) continue;
        if(/*v[1] == u[0] || */v[0] == v[1]) continue;
        else if(v[0] == u[1] && v[1] == u[0]) continue;
        else if(h_degree[u[0]] == 0 && h_degree[u[1]] == 0 && h_degree[v[0]] == 0 && h_degree[v[1]] == 0)
          enable_swing = false;
	else if(h_degree[u[1]] == 0)     continue;
        
        break;
      }
      
      // SWING
      if(enable_swing){
        adjacency[v[0]][v_d[0]]           = u[1];
        adjacency[u[0]][u_d[0]]           = adjacency[u[0]][s_degree[u[0]]-1];
        adjacency[u[0]][s_degree[u[0]]-1] = NOT_DEFINED;
        adjacency[u[1]][s_degree[u[1]]]   = v[0];
        h_degree[u[0]]++; s_degree[u[0]]--; h_degree[u[1]]--; s_degree[u[1]]++;

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
      }
      
      if(enable_swing && accept(switches, current_diameter, diameter, current_ASPL, ASPL, temp, ASPL_priority)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
      }
      else{
        if(enable_swing){ // UNDO
          h_degree[u[0]]--; s_degree[u[0]]++; h_degree[u[1]]++; s_degree[u[1]]--;
          adjacency[u[0]][s_degree[u[0]]-1] = adjacency[u[0]][u_d[0]];
          adjacency[v[0]][v_d[0]]           = u[0];
          adjacency[u[0]][u_d[0]]           = v[0];
          adjacency[u[1]][s_degree[u[1]]]   = NOT_DEFINED;
          temp *= cooling_rate;
          i++;
        }

        // SWAP
        if(get_random(2)){
          adjacency[u[0]][u_d[0]] = v[1];
          adjacency[u[1]][u_d[1]] = v[0];
          adjacency[v[0]][v_d[0]] = u[1];
          adjacency[v[1]][v_d[1]] = u[0];
        }
        else{
          adjacency[u[0]][u_d[0]] = u[1];
          adjacency[u[1]][u_d[1]] = u[0];
          adjacency[v[0]][v_d[0]] = v[1];
          adjacency[v[1]][v_d[1]] = v[0];
        }
        
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
        
        if(accept(switches, current_diameter, diameter, current_ASPL, ASPL, temp, ASPL_priority)){
          current_diameter = diameter;
          current_ASPL     = ASPL;
        }
        else{  // UNDO
          adjacency[u[0]][u_d[0]] = v[0];
          adjacency[u[1]][u_d[1]] = v[1];
          adjacency[v[0]][v_d[0]] = u[0];
          adjacency[v[1]][v_d[1]] = u[1];
        }
      }
      temp *= cooling_rate;
    }
  }
  sa_time = get_time() - sa_time;
  ORP_Finalize_aspl();
  ORP_Conv_adjacency2edge(hosts, switches, radix, best_h_degree, best_s_degree, best_adjacency, edge);
  
  printf("---\n");
  printf("Diameter        = %d\n", best_diameter);
  printf("Diameter Gap    = %d (%d - %d)\n", best_diameter - low_diameter, best_diameter, low_diameter);
  printf("ASPL            = %.10f (%ld/%.0f)\n", best_ASPL, best_sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);
  printf("Time            = %f sec.\n", sa_time);
  printf("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  //  printf("Verify?         = %s\n", (ORP_Verify_edge(hosts, switches, radix, lines, edge))? "Yes" : "No");

  if(outfname)
    ORP_Write_edge(hosts, switches, radix, lines, edge, outfname);

  return 0;
}
