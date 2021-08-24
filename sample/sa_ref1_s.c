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
extern double calc_max_temp_s(const int hosts, const int switches, const int radix, const int seed, const int symmetries);
extern double calc_min_temp_s();

static void SWAP(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

static bool IS_DIAMETER(const int u, const int v, const int switches, const int symmetries)
{
  return (symmetries%2 == 0 && abs(u-v) == switches/2);
}

static int NORM(int x, const int switches)
{
  while(x < 0 || x >= switches)
    x = (x < 0)? x + switches : x - switches;
  return x;
}

static int LOCAL_INDEX(const int v, const int position, const int switches, const int symmetries)
{
  int based_switches = switches/symmetries;
  int n = v - (position/based_switches)*based_switches;
  return NORM(n, switches);
}

// return adjacency[v][d];
static int GLOBAL_ADJ(const int switches, const int radix, const int symmetries,
                      const int adjacency[switches/symmetries][radix], const int v, const int d)
{
  int based_swtiches = switches/symmetries;
  int n = adjacency[v%based_swtiches][d] + (v/based_swtiches)*based_swtiches;
  return NORM(n, switches);
}

static int search_index_s(const int v, const int target, const int exclusion, const int *s_degree,
                          const int radix, const int switches, const int (*adjacency)[radix], const int symmetries)
{
  int based_switches = switches/symmetries;
  if(v == target){
    for(int i=0;i<s_degree[v%based_switches];i++){
      if(GLOBAL_ADJ(switches, radix, symmetries, adjacency, v, i) == target && i != exclusion){
	return i;
      }
    }
  }
  else if(symmetries%2 == 0 && abs(target-v) == switches/2){
    return exclusion;
  }
  else{
    for(int i=0;i<s_degree[v%based_switches];i++){
      if(GLOBAL_ADJ(switches, radix, symmetries, adjacency, v, i) == target){
        return i;
      }
    }
  }

  ERROR("Something Wrong (id=3)\n");
  return -1; // dummy
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static void set_random_edge_s(const int based_total_edges, const int switches, const int *s_degree, const int radix,
                              const int (*adjacency)[radix], const int symmetries, int *u, int *v, int *u_d, int *v_d)
{
  int based_switches = switches/symmetries;
  int r = get_random(based_total_edges);

  *u = NOT_DEFINED;
  for(int i=0;i<based_switches;i++){
    r -= s_degree[i];
    if(r < 0){
      *u = i;
      break;
    }
  }
  if(*u == NOT_DEFINED) ERROR("Something Wrong (id=1)\n");
  
  *u   = get_random(symmetries) * based_switches + (*u);
  *u_d = get_random(s_degree[(*u)%based_switches]);
  *v   = GLOBAL_ADJ(switches, radix, symmetries, adjacency, *u, *u_d); // v = adjacency[u][u_d];
  *v_d = search_index_s(*v, *u, *u_d, s_degree, radix, switches, adjacency, symmetries);
}

#if 0
static void set_random_edge_bias_s(const int based_hosts, const int switches, const int h_degree[switches], const int s_degree[switches],
                                   const int radix, const int (*adjacency)[radix], const int symmetries, int *u, int *v, int *u_d, int *v_d)
{
  int based_switches = switches/symmetries;
  int r = get_random(based_hosts+based_switches);

  *u = NOT_DEFINED;
  for(int i=0;i<based_switches;i++){
    r -= (h_degree[i] + 1);
    if(r < 0){
      *u = i;
      break;
    }
  }
  if(*u == NOT_DEFINED) ERROR("Something Wrong (id=2)\n");
  
  *u   = get_random(symmetries) * based_switches + (*u);
  *u_d = get_random(s_degree[(*u)%based_switches]);
  *v   = GLOBAL_ADJ(switches, radix, symmetries, adjacency, *u, *u_d); // v = adjacency[u][u_d];
  *v_d = search_index_s(*v, *u, *u_d, s_degree, radix, switches, adjacency, symmetries);
}
#endif

static bool check_rotated_edges_overlap(const int u0, const int v0, const int u1, const int v1,
                                        const int switches, const int symmetries)
{
  int diff0 = (u0 > v0)? v0 - u0 + switches : v0 - u0;
  int diff1 = (u1 > v1)? v1 - u1 + switches : v1 - u1;
  int diff2 = (u1 < v1)? u1 - v1 + switches : u1 - v1;
  int based_switches = switches/symmetries;

  if(diff0 == diff1 && u0%based_switches == u1%based_switches)
    return true;
  else if(diff0 == diff2 && u0%based_switches == v1%based_switches)
    return true;
  else
    return false;
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
      double diff = ((current_ASPL-ASPL)*switches*(hosts-1))/symmetries;
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
  fprintf(stderr, "%s [-H hosts] [-S switches] [-R radix] [-G symmetries] [-f input] [-o output] [-s seed] [-n calcs] [-w max_temp] [-c min_temp] [-A] [-B]\n", argv);
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
  fprintf(stderr, "  -B : Increase the bias in the number of hosts\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, int *symmetries, char **infname,
                     char **outfname, int *seed, long *ncalcs, double *max_temp, double *min_temp, bool *ASPL_priority, bool *bias_of_host)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:G:f:o:s:n:w:c:AB"))!=-1){
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
    case 'B':
      *bias_of_host = true;
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static bool swap_adjacency_1opt_s(const int u, const int u_d, const int switches, const int radix,
                                  const int *s_degree, const int symmetries, int adjacency[switches/symmetries][radix])
{
  int based_switches = switches/symmetries;
  int v = GLOBAL_ADJ(switches, radix, symmetries, adjacency, u, u_d);
  if(IS_DIAMETER(u, v, switches, symmetries)) return false;
  int v_d = search_index_s(v, u, u_d, s_degree, radix, switches, adjacency, symmetries);
  int new_v;
  if((u-v)%based_switches == 0){
    if(symmetries <= 4) return false;
    while(1){
      new_v = v + based_switches * get_random(symmetries);
      if(u == new_v || v == new_v) continue;
      else if(NORM(u-v, switches) == NORM(new_v-u, switches)) continue;
      else if(switches%2 == 1) break;
      else if(/* nodes%2 == 0 && */ (u-new_v)%(switches/2) != 0) break;
    }
  }
  else{
    int rnd   = (symmetries%2 == 1)? get_random(symmetries-1) : get_random(symmetries);
    new_v = (rnd != symmetries-1)? v + based_switches*(rnd+1) : u + based_switches*(symmetries/2);
    //  int rnd   = get_random(symmetries-1);
    //  new_v = v + based_nodes*(rnd+1);
  }
  int new_u = v - (new_v - u);
  int tmp[2] = {LOCAL_INDEX(new_v, u, switches, symmetries), LOCAL_INDEX(new_u, v, switches, symmetries)};

  adjacency[u%based_switches][u_d] = tmp[0];
  adjacency[v%based_switches][v_d] = tmp[1];
  return true;
}

int main(int argc, char *argv[])
{
  char *infname = NULL, *outfname = NULL;
  bool ASPL_priority = false, bias_of_host = false;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, seed = DEFAULT_SEED;
  int symmetries = 1, based_switches = NOT_DEFINED;;
  int lines, diameter, current_diameter, best_diameter, low_diameter;
  int (*edge)[2], *h_degree, *s_degree;
  long sum, best_sum, ncalcs = DEFAULT_NCALCS;
  double max_temp = NOT_DEFINED, min_temp = NOT_DEFINED, ASPL, current_ASPL, best_ASPL, low_ASPL;
  ORP_Restore r;

  set_args(argc, argv, &hosts, &switches, &radix, &symmetries, &infname, &outfname, &seed,
           &ncalcs, &max_temp, &min_temp, &ASPL_priority, &bias_of_host);

  if(bias_of_host)
    setenv("ORP_BIAS", "1", 0);
  
  ORP_Srand(seed);
  if(infname){
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
    edge     = ORP_Generate_random_s(hosts, switches, radix, false, symmetries, &lines, h_degree, s_degree);
  }

  if(max_temp == NOT_DEFINED)
    max_temp = calc_max_temp_s(hosts, switches, radix, seed, symmetries);
  
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
    double cooling_rate = pow(min_temp/max_temp, (double)1.0/ncalcs), temp = max_temp;
    int u[2], v[2], u_d[2], v_d[2], tmp[4];
    int based_total_edges = 0;
    for(int i=0;i<based_switches;i++)
      based_total_edges += s_degree[i];

    long interval = (ncalcs < 100)? 1 : ncalcs/100;
    long i = 0, j = 0;
    printf("Ncalcs : Temp : current ASPL Gap ( Dia. ) : Best ASPL Gap ( Dia. )\n");
    while(i < ncalcs){
      if(i/interval == j){
        printf("%ld\t%f\t%f ( %d )\t%f ( %d )\n", i, temp,
               current_ASPL-low_ASPL, current_diameter-low_diameter,
               best_ASPL-low_ASPL, best_diameter-low_diameter);
        j++;
      }

      bool enable_swing = true, enable_swap = true;
      while(1){
        if(bias_of_host){
          u[0]   = get_random(switches);
          u_d[0] = get_random(s_degree[u[0]%based_switches]);
          v[0]   = GLOBAL_ADJ(switches, radix, symmetries, adjacency, u[0], u_d[0]); // v[0] = adjacency[u[0]][u_d[0]];
          v_d[0] = search_index_s(v[0], u[0], u_d[0], s_degree, radix, switches, adjacency, symmetries);
          //
          v[1]   = get_random(switches);
          v_d[1] = get_random(s_degree[v[1]%based_switches]);
          u[1]   = GLOBAL_ADJ(switches, radix, symmetries, adjacency, v[1], v_d[1]); // u[1] = adjacency[v[1]][v_d[1]];
          u_d[1] = search_index_s(u[1], v[1], v_d[1], s_degree, radix, switches, adjacency, symmetries);
        }
        else{
          set_random_edge_s(based_total_edges, switches, s_degree, radix, adjacency, symmetries, &u[0], &v[0], &u_d[0], &v_d[0]);
          set_random_edge_s(based_total_edges, switches, s_degree, radix, adjacency, symmetries, &u[1], &v[1], &u_d[1], &v_d[1]);
        }
        if(u[0] == u[1] || s_degree[u[0]%based_switches] == 1) continue;
        // if(v[0] == u[1]) continue;
        if(/*v[1] == u[0] || */v[0] == v[1]) continue;
        else if(v[0] == u[1] && v[1] == u[0]) continue;
        else if(h_degree[u[0]%based_switches] == 0 && h_degree[u[1]%based_switches] == 0 &&
                h_degree[v[0]%based_switches] == 0 && h_degree[v[1]%based_switches] == 0)
          enable_swing = false;
	else if(h_degree[u[1]%based_switches] == 0) continue;
        break;
      }
      
      // SWING
      if(enable_swing){
        if(symmetries%2 == 0 && abs(u[0]-v[0]) == switches/2){
          int tmp = GLOBAL_ADJ(switches, radix, symmetries, adjacency, u[0], s_degree[u[0]%based_switches]-1);
          adjacency[u[0]%based_switches][u_d[0]] = LOCAL_INDEX(tmp, u[0], switches, symmetries);
          adjacency[u[0]%based_switches][s_degree[u[0]%based_switches]-1] = NOT_DEFINED;
          s_degree[u[0]%based_switches]--;
          adjacency[u[1]%based_switches][s_degree[u[1]%based_switches]] = LOCAL_INDEX(u[1]+switches/2, u[1], switches, symmetries);
          s_degree[u[1]%based_switches]++;
          h_degree[u[0]%based_switches]++; h_degree[u[1]%based_switches]--;
          // adjacency[u[0]][u_d[0]]           = adjacency[u[0]][s_degree[u[0]]-1];
          // adjacency[u[0]][s_degree[u[0]]-1] = NOT_DEFINED;
          // adjacency[u[1]][s_degree[u[1]]]   = u[1] + switches/2;
        }
        else{
          adjacency[v[0]%based_switches][v_d[0]] = LOCAL_INDEX(u[1], v[0], switches, symmetries);
          int tmp = GLOBAL_ADJ(switches, radix, symmetries, adjacency, u[0], s_degree[u[0]%based_switches]-1);
          adjacency[u[0]%based_switches][u_d[0]] = LOCAL_INDEX(tmp, u[0], switches, symmetries);
          adjacency[u[0]%based_switches][s_degree[u[0]%based_switches]-1] = NOT_DEFINED;
          s_degree[u[0]%based_switches]--;
          adjacency[u[1]%based_switches][s_degree[u[1]%based_switches]] = LOCAL_INDEX(v[0], u[1], switches, symmetries);
          s_degree[u[1]%based_switches]++;
          h_degree[u[0]%based_switches]++; h_degree[u[1]%based_switches]--;
        }
            
        ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
        i++;

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
      } // end if(enable_swing)
      
      if(enable_swing && accept_s(hosts, switches, current_diameter, diameter, current_ASPL, ASPL, temp, ASPL_priority, symmetries)){
	current_diameter = diameter;
	current_ASPL     = ASPL;
        temp *= cooling_rate;
      }
      else{
        if(enable_swing){ // UNDO
          h_degree[u[0]%based_switches]--; s_degree[u[0]%based_switches]++;
          h_degree[u[1]%based_switches]++; s_degree[u[1]%based_switches]--;
          int tmp = GLOBAL_ADJ(switches, radix, symmetries, adjacency, u[0], u_d[0]);
          adjacency[u[0]%based_switches][s_degree[u[0]%based_switches]-1] = LOCAL_INDEX(tmp, u[0], switches, symmetries);
          adjacency[u[0]%based_switches][u_d[0]] = LOCAL_INDEX(v[0], u[0], switches, symmetries);
          adjacency[v[0]%based_switches][v_d[0]] = LOCAL_INDEX(u[0], v[0], switches, symmetries);
          adjacency[u[1]%based_switches][s_degree[u[1]%based_switches]] = NOT_DEFINED;      
          temp *= cooling_rate;
        }

        // SWAP
        if(IS_DIAMETER(u[0], v[0], switches, symmetries) && IS_DIAMETER(u[1], v[1], switches, symmetries)){
          if((u[0] - u[1])%based_switches == 0){
            enable_swap = false;
          }
          else{
            if(get_random(2)){
              tmp[0] = LOCAL_INDEX(u[1], u[0], switches, symmetries);
              tmp[1] = LOCAL_INDEX(u[0], u[1], switches, symmetries);
            }
            else{
              tmp[0] = LOCAL_INDEX(v[1], u[0], switches, symmetries);
              tmp[1] = LOCAL_INDEX(v[0], u[1], switches, symmetries);
            }
            adjacency[u[0]%based_switches][u_d[0]] = tmp[0];
            adjacency[u[1]%based_switches][u_d[1]] = tmp[1];
          }
        }
        else if(IS_DIAMETER(u[0], v[0], switches, symmetries) || IS_DIAMETER(u[1], v[1], switches, symmetries)){
          if(IS_DIAMETER(u[1], v[1], switches, symmetries)){
            SWAP(&u[0], &u[1]); SWAP(&u_d[0], &u_d[1]);
            SWAP(&v[0], &v[1]); SWAP(&v_d[0], &v_d[1]);
          }
          
          int opposite = switches/2, rnd = get_random(4);
          if(rnd == 0){ // u[0]--v[1], u[1]--u[1]', v[0]--v[1]'
            tmp[0] = LOCAL_INDEX(v[1],          u[0], switches, symmetries);
            tmp[1] = LOCAL_INDEX(u[1]+opposite, u[1], switches, symmetries);
            tmp[2] = LOCAL_INDEX(v[1]+opposite, v[0], switches, symmetries);
            tmp[3] = LOCAL_INDEX(u[0],          v[1], switches, symmetries);
          }
          else if(rnd == 1){ // u[0]--v[1]', v[0]--v[1], u[1]--u[1]'
            tmp[0] = LOCAL_INDEX(v[1]+opposite, u[0], switches, symmetries);
            tmp[1] = LOCAL_INDEX(u[1]+opposite, u[1], switches, symmetries);
            tmp[2] = LOCAL_INDEX(v[1],          v[0], switches, symmetries);
            tmp[3] = LOCAL_INDEX(v[0],          v[1], switches, symmetries);
          }
          else if(rnd == 2){ // u[0]--u[1], v[0]--u[1]', v[1]--v[1]'
            tmp[0] = LOCAL_INDEX(u[1],          u[0], switches, symmetries);
            tmp[1] = LOCAL_INDEX(u[0],          u[1], switches, symmetries);
            tmp[2] = LOCAL_INDEX(u[1]+opposite, v[0], switches, symmetries);
            tmp[3] = LOCAL_INDEX(v[1]+opposite, v[1], switches, symmetries);
          }
          else if(rnd == 3){ // u[0]--u[1]', u[1]--v[0], v[1]--v[1]'
            tmp[0] = LOCAL_INDEX(u[1]+opposite, u[0], switches, symmetries);
            tmp[1] = LOCAL_INDEX(v[0],          u[1], switches, symmetries);
            tmp[2] = LOCAL_INDEX(u[1],          v[0], switches, symmetries);
            tmp[3] = LOCAL_INDEX(v[1]+opposite, v[1], switches, symmetries);
          }
          
          adjacency[u[0]%based_switches][u_d[0]] = tmp[0];
          adjacency[u[1]%based_switches][u_d[1]] = tmp[1];
          adjacency[v[0]%based_switches][v_d[0]] = tmp[2];
          adjacency[v[1]%based_switches][v_d[1]] = tmp[3];
        }
        else if(check_rotated_edges_overlap(u[0], v[0], u[1], v[1], switches, symmetries)){
          // Two selected edges are symmetrical
          enable_swap = swap_adjacency_1opt_s(u[0], u_d[0], switches, radix, s_degree, symmetries, adjacency);
        }
        else{
          if(get_random(2)){ // u[0]--v[1], v[0]--u[1]
            if(IS_DIAMETER(u[0], v[1], switches, symmetries) || IS_DIAMETER(v[0], u[1], switches, symmetries)){
              enable_swap = false;
            }
            else if(check_rotated_edges_overlap(u[0], v[1], v[0], u[1], switches, symmetries)){
              enable_swap = false;
            }
            else{
              tmp[0] = LOCAL_INDEX(v[1], u[0], switches, symmetries);
              tmp[1] = LOCAL_INDEX(v[0], u[1], switches, symmetries);
              tmp[2] = LOCAL_INDEX(u[1], v[0], switches, symmetries);
              tmp[3] = LOCAL_INDEX(u[0], v[1], switches, symmetries);
            }
          }
          else{ // u[0]--u[1], v[0]--v[1]
            if(IS_DIAMETER(u[0], u[1], switches, symmetries) || IS_DIAMETER(v[0], v[1], switches, symmetries)){
              enable_swap = false;
            }
            else if(check_rotated_edges_overlap(u[0], u[1], v[0], v[1], switches, symmetries)){
              enable_swap = false;
            }
            else{
              tmp[0] = LOCAL_INDEX(u[1], u[0], switches, symmetries);
              tmp[1] = LOCAL_INDEX(u[0], u[1], switches, symmetries);
              tmp[2] = LOCAL_INDEX(v[1], v[0], switches, symmetries);
              tmp[3] = LOCAL_INDEX(v[0], v[1], switches, symmetries);
            }
          }
          if(enable_swap){
            adjacency[u[0]%based_switches][u_d[0]] = tmp[0];
            adjacency[u[1]%based_switches][u_d[1]] = tmp[1];
            adjacency[v[0]%based_switches][v_d[0]] = tmp[2];
            adjacency[v[1]%based_switches][v_d[1]] = tmp[3];
          }
        }

        if(enable_swap){
          ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
          i++;
          
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
          else{  // UNDO
            adjacency[u[0]%based_switches][u_d[0]] = LOCAL_INDEX(v[0], u[0], switches, symmetries);
            adjacency[u[1]%based_switches][u_d[1]] = LOCAL_INDEX(v[1], u[1], switches, symmetries);
            adjacency[v[0]%based_switches][v_d[0]] = LOCAL_INDEX(u[0], v[0], switches, symmetries);
            adjacency[v[1]%based_switches][v_d[1]] = LOCAL_INDEX(u[1], v[1], switches, symmetries);
          }
          temp *= cooling_rate;
        }
      }
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
  printf("ASPL priority?  = %s\n", (ASPL_priority)? "Yes" : "No");
  //  printf("Verify?         = %s\n", (ORP_Verify_edge(hosts, switches, radix, lines, edge))? "Yes" : "No");

  if(outfname)
    ORP_Write_edge(hosts, switches, radix, lines, edge, outfname);

  return 0;
}
