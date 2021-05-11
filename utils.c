#include "common.h"
#ifdef _OPENMP
static int *_local_frontier;
#pragma omp threadprivate(_local_frontier)

void ORP_declare_local_frontier(const int switches)
{
#pragma omp parallel
  {
    _local_frontier = malloc(sizeof(int) * switches);
  }
}

void ORP_free_local_frontier()
{
#pragma omp parallel
  {
    free(_local_frontier);
  }
}
#endif

void ORP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2)
{
#if defined(__ARM_NEON) || defined(__FUJITSU)
  posix_memalign((void **)a, ALIGN_VALUE, s);
#else
  if(enable_avx2)
    *a = _mm_malloc(s, ALIGN_VALUE);
  else
    posix_memalign((void **)a, ALIGN_VALUE, s);
#endif
}

void ORP_Free(uint64_t *a, const bool enable_avx2)
{
#if defined(__ARM_NEON) || defined(__FUJITSU)
  free(a);
#else
  if(enable_avx2)
    _mm_free(a);
  else
    free(a);
#endif
}

double moore_bound(const double nodes, const double degree)
{
  if(degree + 1 >= nodes)
    return 1;
  
  double diam = -1, n = 1, r = 1, aspl = 0.0, prev_tmp;
  while(1){
    double tmp = n + degree * pow(degree-1, r-1);
    if(tmp >= nodes || (r > 1 && prev_tmp == tmp))
      break;
    
    n = tmp;
    aspl += r * degree * pow(degree-1, r-1);
    diam = r++;
    prev_tmp = tmp;
  }

  diam++;
  aspl += diam * (nodes - n);
  aspl /= (nodes - 1);
  
  return aspl;
}

void CHECK_PARAMETERS(const int hosts, const int switches, const int radix)
{
  if(hosts < 3)
    ERROR("Hosts (%d) >= 3\n", hosts);
  else if(switches < 4)
    ERROR("Switches (%d) >= 4\n", switches); // It is troublesome to make the case of switches = 1,2,3.
  else if(radix < 3)
    ERROR("Radix (%d) >= 3\n", radix);
  else if(switches*radix-2*(switches-1) < hosts)
    ERROR("Switches (%d) * Radix (%d) - 2 * (Switches (%d) - 1) >= Host (%d))\n",
          switches, radix, switches, hosts);
}

// This function is based on the following paper.
//   Ryota Yasudo, Michihiro Koibuchi, Koji Nakano.
//   Designing High-Performance Interconnection Networks with Host-Switch Graphs.
//   IEEE Transactions on Parallel and Distributed Systems, 30(2), 315-330. 2019.
//   https://doi.org/10.1109/TPDS.2018.2864286
static double continuous_moore_bound(const int hosts, const int switches, const int radix)
{
  double h = hosts;
  double s = switches;
  double r = radix;
  return moore_bound(s,r-h/s)*(s*h-h)/(s*h-s)+2;
}

int ORP_Optimize_switches(const int hosts, const int radix)
{
  if(hosts < 3)
    ERROR("Hosts (%d) >= 3\n", hosts);
  else if(radix < 3)
    ERROR("Radix (%d) >= 3\n", radix);
  
  int s = 3;
  double prev = DBL_MAX;
  while(1){
    if(s*radix-2*(s-1) >= hosts){
      double tmp = continuous_moore_bound(hosts, s, radix);
      if(prev <= tmp)
        break;
      prev = tmp;
    }
    s++;
  }

  return s - 1;
}

static bool IS_HOST(const int v, const int hosts)
{
  return (v < hosts);
}

static bool IS_SWITCH(const int v, const int hosts)
{
  return (v >= hosts);
}

void ORP_Set_degrees(const int hosts, const int switches, const int lines, const int (*edge)[2],
                     int h_degree[switches], int s_degree[switches])
{
  for(int i=0;i<switches;i++)
    h_degree[i] = s_degree[i] = 0;

  for(int i=0;i<lines;i++){
    if(IS_HOST(edge[i][0], hosts)){
      h_degree[edge[i][1]-hosts]++;
    }
    else if(IS_HOST(edge[i][1], hosts)){
      h_degree[edge[i][0]-hosts]++;
    }
    else{ // Both vertices are switches
      s_degree[edge[i][0]-hosts]++;
      s_degree[edge[i][1]-hosts]++;
    }
  }
}

bool ORP_Verify_edge(const int hosts, const int switches, const int radix, const int lines, const int (*edge)[2])
{
  int local_hosts[hosts];
  for(int i=0;i<hosts;i++)
    local_hosts[i] = 0;
  
  for(int i=0;i<lines;i++){
    if(IS_HOST(edge[i][0], hosts))
      local_hosts[edge[i][0]]++;
    else if(IS_HOST(edge[i][1], hosts))
      local_hosts[edge[i][1]]++;
    else if(IS_HOST(edge[i][0], hosts) && IS_HOST(edge[i][1], hosts)){
      fprintf(stderr, "Both vertices are hosts\n");
      return false;
    }
  }
  for(int i=0;i<hosts;i++){
    if(local_hosts[i] != 1){
       fprintf(stderr, "hosts is wrong.\n");
       return false;
    }
  }
  
  int h_degree[switches], s_degree[switches];
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  for(int i=0;i<switches;i++){
    if(radix < h_degree[i] + s_degree[i]){
      fprintf(stderr, "radix is too large.\n");
      return false;
    }
  }

  return true;
}

void ORP_Srand(const unsigned int seed)
{
  srand(seed);
}

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

void ORP_Set_host_degree(const int hosts, const int switches, const int lines, const int (*edge)[2],
                         int h_degree[switches])
{
  for(int i=0;i<switches;i++)
    h_degree[i] = 0;

  for(int i=0;i<lines;i++){
    if(IS_HOST(edge[i][0], hosts))
      h_degree[edge[i][1]-hosts]++;
    else if(IS_HOST(edge[i][1], hosts))
      h_degree[edge[i][0]-hosts]++;
  }
}

void ORP_Set_switch_degree(const int hosts, const int switches, const int lines, const int (*edge)[2],
                           int s_degree[switches])
{
  for(int i=0;i<switches;i++)
    s_degree[i] = 0;

  for(int i=0;i<lines;i++){
    if(IS_SWITCH(edge[i][0],hosts) && IS_SWITCH(edge[i][1],hosts)){
      s_degree[edge[i][0]-hosts]++;
      s_degree[edge[i][1]-hosts]++;
    }
  }
}

void ORP_Print_degree(const int hosts, const int switches, const int degree[switches])
{
  for(int i=0;i<switches;i++)
    printf("%d : %d\n", i + hosts, degree[i]);
}

void ORP_Read_property(const char* fname, int* hosts, int* switches, int* radix, int *lines)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  fscanf(fp, "%d %d %d", hosts, switches, radix);

  *lines = 0;
  int n1, n2;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF)
    (*lines)++;

  fclose(fp);
}

void ORP_Read_edge(const char* fname, int (*edge)[2])
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int hosts, switches, radix;
  fscanf(fp, "%d %d %d", &hosts, &switches, &radix); // Skip the first line
  
  int n1, n2, i = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    edge[i][0] = n1;
    edge[i][1] = n2;
    i++;
  }

  fclose(fp);
}

void ORP_Write_edge(const int hosts, const int switches, const int radix, const int lines,
                    const int edge[lines][2], char *fname)
{
  FILE *fp = NULL;

  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);

  fprintf(fp, "%d %d %d\n", hosts, switches, radix);
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);
  
  fclose(fp);
}

void ORP_Conv_edge2adjacency(const int hosts, const int switches, const int radix,
                             const int lines, const int (*edge)[2], int (*adjacency)[radix])
{
  int s_degree[switches];
  for(int i=0;i<switches;i++)
    s_degree[i] = 0;

  for(int i=0;i<switches;i++)
    for(int j=0;j<radix;j++)
      adjacency[i][j] = NOT_DEFINED;

  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    if(IS_SWITCH(n1,hosts) && IS_SWITCH(n2,hosts)){
      int s1 = n1 - hosts;
      int s2 = n2 - hosts;
      adjacency[s1][s_degree[s1]++] = s2;
      adjacency[s2][s_degree[s2]++] = s1;
    }
  }
}

void ORP_Conv_adjacency2edge(const int hosts, const int switches, const int radix, const int *h_degree,
                             const int *s_degree, const int (*adjacency)[radix], int (*edge)[2])
{
  int lines = 0;
  for(int i=0;i<switches;i++){
    for(int j=0;j<h_degree[i];j++){
      edge[lines][0] = lines;
      edge[lines][1] = hosts + i;
      lines++;
    }
  }

  for(int i=0;i<switches;i++){
    int loop_count = 0;
    for(int j=0;j<s_degree[i];j++){
      int v = adjacency[i][j];
      if(i < v){
        edge[lines][0] = i + hosts;
        edge[lines][1] = v + hosts;
        lines++;
      }
      else if(i == v){
        loop_count++;
        if(loop_count%2 == 0){
          edge[lines][0] = i + hosts;
          edge[lines][1] = v + hosts;
          lines++;
        }
      }
    }
    if(loop_count%2 == 1) ERROR("Something Wrong. (id=3)\n");
  }
}

void ORP_Print_adjacency(const int hosts, const int switches, const int radix, const int s_degree[switches],
                         const int (*adjacency)[radix])
{
  for(int i=0;i<switches;i++){
    printf("%d :", i);
    for(int j=0;j<s_degree[i];j++)
      printf("%3d ", adjacency[i][j]);
    printf("\n");
  }
}

void ORP_Print_info(const int hosts, const int switches, const int radix, const int h_degree[switches],
                    const int s_degree[switches], const int (*adjacency)[radix])
{
  printf("hosts = %d, swithces = %d, radix = %d\n", hosts, switches, radix);
  printf("   h_degree s_degree h_degree + s_degree\n");
  for(int i=0;i<switches;i++)
    printf("%2d : %7d %7d %7d\n", i, h_degree[i], s_degree[i], h_degree[i]+s_degree[i]);
  printf("   adjacency\n");
  for(int i=0;i<switches;i++){
    printf("%2d : ", i);
    for(int j=0;j<s_degree[i];j++)
      printf("%3d ", adjacency[i][j]);
    printf("\n");
  }
}

void ORP_Print_edge(const int lines, const int (*edge)[2])
{
  for(int i=0;i<lines;i++)
    printf("%d %d\n", edge[i][0], edge[i][1]);
}

void ORP_Print_switch(const int hosts, const int lines, const int (*edge)[2])
{
  for(int i=0;i<lines;i++)
    if(IS_SWITCH(edge[i][0],hosts) && IS_SWITCH(edge[i][1],hosts))
      printf("%d %d\n", edge[i][0]-hosts, edge[i][1]-hosts);
}

void ORP_Write_switch(const int hosts, const int lines, const int edge[lines][2], char *fname)
{
  FILE *fp = NULL;
  
  if((fp = fopen(fname, "w")) == NULL)
    ERROR("Cannot open %s\n", fname);
  
  for(int i=0;i<lines;i++)
    if(IS_SWITCH(edge[i][0],hosts) && IS_SWITCH(edge[i][1],hosts))
      fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);

  fclose(fp);
}

double ORP_Get_mem_usage(const int kind, const int switches, const int symmetries)
{
  double mem;
  if(kind == ASPL_MATRIX){
    mem  = ((double)switches * switches * sizeof(uint64_t) * 2)/UINT64_BITS;
    mem += (double)switches * switches * sizeof(char);
  }
  else{ // kind == ASPL_BFS
    mem  = (double)switches * 3 * sizeof(int);
#ifdef _OPENMP
    mem += (double)switches * sizeof(int) * omp_get_max_threads();
#endif
  }

  return mem/(1024*1024); // to Mbytes
}

int ORP_Get_kind()
{
  int kind  = ASPL_MATRIX;
  char *val = getenv("ORP_ASPL");
  if(val)
    if(strcmp(val, "BFS") == 0)
      kind = ASPL_BFS;
    else
      ERROR("Unknown ORP_ASPL value (%s)\n", val);

  return kind;
}

bool ORP_Check_profile()
{
  static bool first = true, enable_profile = false;
  
  if(first){
    first = false;
    char *val = getenv("ORP_PROFILE");
    if(val)
      if(atoi(val) == 1)
        enable_profile = true;
      else if(atoi(val) == 0)
        enable_profile = false;
      else
        ERROR("Unknown ORP_PROFILE value (%d)\n", atoi(val));
  }
  
  return enable_profile;
}

void ORP_Profile(const char* name, const int kind, const int switches, const int symmetries, 
                 const double elapsed_time, const unsigned int times)
{
  char kind_name[7], hostname[MAX_HOSTNAME_LENGTH];
  if(kind == ASPL_MATRIX) strcpy(kind_name, "MATRIX");
  else /* ASPL_BFS */     strcpy(kind_name, "BFS");
  gethostname(hostname, sizeof(hostname));
  time_t t = time(NULL);

  printf("------ Profile for SET_ASPL ------\n");
  printf("Date            = %s", ctime(&t));
  printf("Hostname        = %s\n", hostname);
  printf("Number of Times = %d\n", times);
  printf("Total Time      = %f sec.\n", elapsed_time);
  printf("Average Time    = %f sec.\n", elapsed_time/times);
  printf("Algorithm       = %s (%s)\n", kind_name, name);
  printf("Memory Usage    = %.3f MB\n", ORP_Get_mem_usage(kind, switches, symmetries));
  if(symmetries != 1)
    printf("Symmetries      = %d\n", symmetries);
#ifdef _OPENMP
  printf("Num of Threads  = %d\n", omp_get_max_threads());
#else
  printf("Num of Threads  = %d\n", 1);
#endif
  printf("--------- End of Profile ---------\n");
}

double ORP_Get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

#ifdef __AVX2__
static void matmul_avx2(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                        const int *restrict s_degree, const int *restrict adjacency, const int elements,
                        const int symmetries)
{
  int quarter_elements = elements/4;
  if(symmetries == 1){
#pragma omp parallel for
    for(int i=0;i<switches;i++){
      __m256i *b = (__m256i *)(B + i*elements);
      for(int j=0;j<s_degree[i];j++){
        int n = *(adjacency + i * radix + j);  // int n = adjacency[i][j];
        __m256i *a = (__m256i *)(A + n*elements);
        for(int k=0;k<quarter_elements;k++){
          __m256i aa = _mm256_load_si256(a+k);
          __m256i bb = _mm256_load_si256(b+k);
          _mm256_store_si256(b+k, _mm256_or_si256(aa, bb));
        }
      }
    }
  }
  else{ /* NOT IMPLEMENTED */
  }
}
#endif

static void matmul(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                   const int *restrict s_degree, const int *restrict adjacency, const int elements,
                   const int symmetries)
{
  if(symmetries == 1){
#pragma omp parallel for
    for(int i=0;i<switches;i++){
      for(int j=0;j<s_degree[i];j++){
        int n = *(adjacency + i * radix + j);  // int n = adjacency[i][j];
        for(int k=0;k<elements;k++)
          B[i*elements+k] |= A[n*elements+k];
      }
    }
  }
  else{	/* NOT IMPLEMENTED */
  }
}

void ORP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                const int *restrict s_degree, const int *restrict adjacency, const int elements,
                const int symmetries, const bool enable_avx2)
{
#ifdef __AVX2__
  if(enable_avx2) matmul_avx2(A, B, switches, radix, s_degree, adjacency, elements, symmetries);
  else            matmul     (A, B, switches, radix, s_degree, adjacency, elements, symmetries);
#else
  matmul(A, B, switches, radix, s_degree, adjacency, elements, symmetries);
#endif
}

static int simple_top_down_step(const int switches, const int num_frontier, const int radix,
                                int *s_degree, int (*adjacency)[radix],
                                const int* restrict frontier, int* restrict next, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<s_degree[v];j++){
      int n = adjacency[v][j];
      if(bitmap[n] == NOT_VISITED){
        bitmap[n] = VISITED;
        next[count++] = n;
      }
    }
  }
  return count;
}

static bool simple_bfs(const int switches, const int radix, int* s_degree, int (*adjacency)[radix])
{
  char *bitmap  = malloc(sizeof(char) * switches);
  int *frontier = malloc(sizeof(int)  * switches);
  int *next     = malloc(sizeof(int)  * switches);
  int num_frontier = 1, root = 0;

  for(int i=0;i<switches;i++)
    bitmap[i] = NOT_VISITED;

  frontier[0]  = root;
  bitmap[root] = VISITED;

  while(1){
    num_frontier = simple_top_down_step(switches, num_frontier, radix, s_degree,
                                        adjacency, frontier, next, bitmap);
    if(num_frontier == 0) break;

    int *tmp = frontier;
    frontier = next;
    next     = tmp;
  }

  bool flag = false;
  for(int i=0;i<switches;i++)
    if(bitmap[i] == NOT_VISITED)
      flag = true;

  free(frontier);
  free(next);
  free(bitmap);

  return flag;
}

#ifdef _OPENMP
int ORP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency, 
                      const int switches, const int radix, const int* restrict s_degree, const int symmetries,
                      int* restrict frontier, int* restrict next, int* restrict distance)
{
  int count = 0;
  if(symmetries == 1){
#pragma omp parallel
    {
      int local_count = 0;
#pragma omp for nowait
      for(int i=0;i<num_frontier;i++){
        int v = frontier[i];
        for(int j=0;j<s_degree[v];j++){
          int n = *(adjacency + v * radix + j);
          if(distance[n] == NOT_USED){
            distance[n] = level;
            _local_frontier[local_count++] = n;
          }
      	}
      }  // end for i
#pragma omp critical
      {
      	memcpy(&next[count], _local_frontier, local_count*sizeof(int));
        count += local_count;
      }
    }
  }
  else{ /* NOT_IMPLEMENTED */
  }
  return count;
}
#else
int ORP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
                      const int switches, const int radix, const int* restrict s_degree, const int symmetries,
                      int* restrict frontier, int* restrict next, int* restrict distance, char* restrict bitmap)
{
  int count = 0;
  if(symmetries == 1){
    for(int i=0;i<num_frontier;i++){
      int v = frontier[i];
      for(int j=0;j<s_degree[v];j++){
        int n = *(adjacency + v * radix + j);
        if(distance[n] == NOT_USED){
          distance[n] = level;
          next[count++] = n;
	}
      }
    }
  }
  else{	/* NOT IMPLEMENTED */
  }
  return count;
}
#endif

static void backup_restore(const int u[2], const int u_d[2], const int v[2], const int v_d[2],
                           const int op, ORP_Restore *r)
{
  if(r == NULL) return;

  r->op   = op;
  for(int i=0;i<2;i++){
    r->u[i]   = u[i];
    r->u_d[i] = u_d[i];
    r->v[i]   = v[i];
    r->v_d[i] = v_d[i];
  }
}

void ORP_Restore_adjacency(const ORP_Restore r, const int radix, int *h_degree, int *s_degree,
                           int (*adjacency)[radix])
{
  if(r.op == OP_SWAP){
    for(int i=0;i<2;i++){
      adjacency[r.u[i]][r.u_d[i]] = r.v[i];
      adjacency[r.v[i]][r.v_d[i]] = r.u[i];
    }
  }
  else{ // OP_SWING
    adjacency[r.u[0]][s_degree[r.u[0]]]   = r.v[0];
    adjacency[r.v[0]][r.v_d[0]]           = r.u[0];
    adjacency[r.u[1]][s_degree[r.u[1]]-1] = NOT_DEFINED;
    h_degree[r.u[0]]--; s_degree[r.u[0]]++; h_degree[r.u[1]]++; s_degree[r.u[1]]--;
  }
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
}

bool swing_adjacency(const int switches, const int radix, int h_degree[switches], int s_degree[switches],
                     ORP_Restore *r, int adjacency[switches][radix])
{
  int u[2], v[2], u_d[2], v_d[2]; // u_d[1], v[1], and v_d[1] are not used because it is a host.

  while(1){
    u[0] = get_random(switches);
    u[1] = get_random(switches);
    if(u[0] == u[1] || h_degree[u[1]] == 0) continue;

    u_d[0] = get_random(s_degree[u[0]]);
    v[0]   = adjacency[u[0]][u_d[0]];
    if(v[0] == u[1]) continue;
    break;
  }
  
  // search index
  v_d[0] = search_index(v[0], u[0], u_d[0], s_degree, radix, adjacency);

  // u[0]--v[0], u[1]--h[0] -> u[0]--h[0], u[1]--v[0]
  if(s_degree[u[0]] == 1) return false;
  backup_restore(u, u_d, v, v_d, OP_SWING, r);
  adjacency[v[0]][v_d[0]]           = u[1];
  adjacency[u[0]][u_d[0]]           = adjacency[u[0]][s_degree[u[0]]-1];
  adjacency[u[0]][s_degree[u[0]]-1] = NOT_DEFINED;
  adjacency[u[1]][s_degree[u[1]]]   = v[0];
  h_degree[u[0]]++; s_degree[u[0]]--; h_degree[u[1]]--; s_degree[u[1]]++;

  return true;
}

void ORP_Swing_adjacency(const int switches, const int radix, int h_degree[switches], int s_degree[switches],
                         ORP_Restore *r, int adjacency[switches][radix])
{
  while(1){
    if(swing_adjacency(switches, radix, h_degree, s_degree, r, adjacency))
      break;
  }
}

void ORP_Swap_adjacency(const int switches, const int radix, const int s_degree[switches],
                        ORP_Restore *r, int adjacency[switches][radix])
{
  int u[2], v[2], u_d[2], v_d[2];

  while(1){
    u[0] = get_random(switches);
    u[1] = get_random(switches);
    if(u[0] == u[1]) continue;
    
    u_d[0] = get_random(s_degree[u[0]]);
    v[0]   = adjacency[u[0]][u_d[0]];
    if(v[0] == u[1]) continue;
    
    u_d[1] = get_random(s_degree[u[1]]);
    v[1]   = adjacency[u[1]][u_d[1]];
    if(v[1] == u[0] || v[0] == v[1]) continue;
    break;
  }

  // search index
  v_d[0] = search_index(v[0], u[0], u_d[0], s_degree, radix, adjacency);
  v_d[1] = search_index(v[1], u[1], u_d[1], s_degree, radix, adjacency);

  // backup for restore
  backup_restore(u, u_d, v, v_d, OP_SWAP, r);

  // u[0]--v[0], u[1]--v[1] -> u[0]--v[1], u[1]--v[0]
  adjacency[u[0]][u_d[0]] = v[1];
  adjacency[u[1]][u_d[1]] = v[0];
  adjacency[v[0]][v_d[0]] = u[1];
  adjacency[v[1]][v_d[1]] = u[0];
}

void* ORP_Generate_random(const int hosts, const int switches, const int radix, const bool assign_evenly,
                          int *lines, int h_degree[switches], int s_degree[switches])
{
  CHECK_PARAMETERS(hosts, switches, radix);

  // malloc edge
  *lines = (switches * radix - hosts)/2 + hosts;
  int (*edge)[2] = malloc(sizeof(int) * (*lines) * 2); // int edge[*lines][2];

  for(int i=0;i<switches;i++)
    h_degree[i] = hosts/switches;

  for(int i=0;i<hosts%switches;i++)
    h_degree[i]++;

  // connect switch-switch
  int degree[switches];
  for(int i=0;i<switches;i++){
    degree[i] 	= h_degree[i];
    s_degree[i] = 0;
  }
  
  int tmp_lines = 0;
  for(int i=0;i<switches-1;i++){
    int d = radix - degree[i];
    for(int j=0;j<d;j++){
      edge[tmp_lines][0] = i+hosts;
      edge[tmp_lines][1] = (i+1)+hosts;
      s_degree[i]++;
      s_degree[i+1]++;
      degree[i]++;
      degree[i+1]++;
      tmp_lines++;
    }
  }

  // Add loop to the last switch
  int last = switches - 1;
  if(degree[last] < radix-1){
    int loops = (radix - degree[last])/2;
    for(int i=0;i<loops;i++){
      edge[tmp_lines][0] = edge[tmp_lines][1] = last + hosts;
      s_degree[last] += 2;
      tmp_lines++;
    }
  }
  if(tmp_lines + hosts != *lines) ERROR("Something Wrong (id=1)\n");
                                    
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  ORP_Conv_edge2adjacency(hosts, switches, radix, tmp_lines, edge, adjacency);

  // Give randomness
  for(int i=0;i<tmp_lines*GEN_GRAPH_ITERS;i++){
    ORP_Swap_adjacency(switches, radix, s_degree, NULL, adjacency);
    if(!assign_evenly) ORP_Swing_adjacency(switches, radix, h_degree, s_degree, NULL, adjacency);
  }

  // Repeat until there are no unreachable vertices
  while(simple_bfs(switches, radix, s_degree, adjacency)){
    ORP_Swap_adjacency(switches, radix, s_degree, NULL, adjacency);
    if(!assign_evenly) ORP_Swing_adjacency(switches, radix, h_degree, s_degree, NULL, adjacency);
  }

  ORP_Conv_adjacency2edge(hosts, switches, radix, h_degree, s_degree, adjacency, edge);
  free(adjacency);
  
  return edge;
}

// https://github.com/r-ricdeau/host-switch-aspl
void ORP_Set_lbounds(const int hosts, const int radix, int *low_diameter, double *low_ASPL)
{
  int i, leaf_n = hosts - 1, total_n, total_m = 1, current_layer, next_layer, wiener, d_lb;
  double aspl_lb;

  for(i=1;;i++){
    if ((int)pow((radix - 1), i) > leaf_n) break;
    total_m += (int)pow((radix - 1), i);
  }

  i--;
  total_n  = (int)pow((radix - 1), i);
  total_m -= (int)pow((radix - 1), i);
  current_layer = total_n;

  next_layer = ceil((double)(leaf_n - total_n) / (double)(radix - 2));
  current_layer -= next_layer;
  total_m += next_layer;
  next_layer += leaf_n - total_n;
  wiener = current_layer * (i + 1) + next_layer * (i + 2);

  aspl_lb = (double)wiener / (double)(hosts - 1);
  d_lb = (next_layer > 0)? i + 2 : i + 1;

  *low_ASPL     = aspl_lb;
  *low_diameter = d_lb;
}
