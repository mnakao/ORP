#include "common.h"
static bool _enable_avx2 = false, _is_profile = false, _enable_disconnected = false;
static char* _bitmap = NULL;
static int _hosts, _switches, _based_switches, _radix, _symmetries, _kind, _elements, _times;
static int *_frontier = NULL, *_distance = NULL, *_next = NULL;
static uint64_t *_A, *_B;
static double _elapsed_time;

static void printb(unsigned int v) {
  unsigned int mask = (int)1 << (sizeof(v) * CHAR_BIT - 1);
  do putchar(mask & v ? '1' : '0');
  while (mask >>= 1);
  printf("\n");
}

void ORP_Profile(const char* name, const int kind, const int nodes, const int symmetries,
                 const double elapsed_time, const unsigned int times);
void ORP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                const int *restrict ports, const int *restrict adjacency, const int elements,
                const bool enable_avx2);
void ORP_Matmul_s(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                  const int *restrict ports, const int *restrict adjacency, const int elements,
                  const bool enable_avx2, const int symmetries);
void ORP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2);
void ORP_Free(uint64_t *a, const bool enable_avx2);
void ORP_declare_local_frontier(int swithces);
void ORP_free_local_frontier();
bool ORP_Check_profile();
int ORP_Get_kind();
int ORP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
                      const int switches, const int radix, const int* restrict ports,
                      int* restrict frontier, int* restrict next, int* restrict distance);
int ORP_top_down_step_s(const int level, const int num_frontier, const int* restrict adjacency,
                        const int switches, const int radix, const int* restrict ports,
                        int* restrict frontier, int* restrict next, int* restrict distance, const int symmetries);
extern double ORP_Get_time();

static bool CHECK_DISCONNECTED()
{
  char *val = getenv("ORP_DISCONNECTED");
  if(!val){
    return false;
  }
  else{
    if(atoi(val) == 1)
      return true;
    else if(atoi(val) == 0)
      return false;
    else
      ERROR("Unknown ORP_DISCONNECTED value (%d)\n", atoi(val));
  }

  return false; // dummy
}

static void aspl_mat(const int* restrict h_degree, const int* restrict s_degree, const int* restrict adjacency,
                     int *diameter, long *sum, double *ASPL)
{
#pragma omp parallel for
  for(int i=0;i<_switches*_switches;i++)
    _bitmap[i] = NOT_VISITED;

#pragma omp parallel for
  for(int i=0;i<_switches*_elements;i++)
    _A[i] = _B[i] = 0;
  
  long k = 0, stop_k = ((long)_switches*_switches-_switches)/2, local_sum = 0;
  int local_diameter = 0;
#pragma omp parallel for
  for(int i=0;i<_switches;i++){
    unsigned int offset = i*_elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *diameter = 1;
  long level = 2;
  for(int kk=0;kk<_switches;kk++){
    ORP_Matmul(_A, _B, _switches, _radix, s_degree, adjacency, _elements, _enable_avx2);

    level++;
#pragma omp parallel for reduction(+:k,local_sum) reduction(max:local_diameter)
    for(int i=0;i<_switches;i++){
      for(int j=i+1;j<_switches;j++){
        int ii = i*_switches+j;
        if(_bitmap[ii] == NOT_VISITED){
          if(_B[i*_elements+(j/UINT64_BITS)] & (0x1ULL<<(j%UINT64_BITS))){
            _bitmap[ii] = VISITED;
            k++;
            if(h_degree[i] != 0 && h_degree[j] != 0){
              local_diameter = MAX(local_diameter, level-2);
              local_sum += level * h_degree[i] * h_degree[j];
            }
          }
          else if(_enable_disconnected && (h_degree[i] == 0 || h_degree[j] == 0)){
            _bitmap[ii] = VISITED;
            k++;
          }
        }
      }
    }
    *diameter = local_diameter;
    if(k == stop_k) break;

    // swap A <-> B
    uint64_t* tmp = _A;
    _A = _B;
    _B = tmp;

    if(kk == _switches-1)
      (*diameter) = _switches;
  }

#pragma omp parallel for reduction(+:local_sum)
  for(int i=0;i<_switches;i++)
    local_sum += (long)h_degree[i] * (h_degree[i] - 1);

  *ASPL = local_sum / (double)(((long)_hosts*(_hosts-1))/2);
  *sum  = local_sum;
  *diameter += 2;
}

static void aspl_mat_s(const int* restrict h_degree, const int* restrict s_degree, const int* restrict adjacency,
                       int *diameter, long *sum, double *ASPL)
{
#pragma omp parallel for
  for(int i=0;i<_switches*_based_switches;i++)
    _bitmap[i] = NOT_VISITED;

#pragma omp parallel for
  for(int i=0;i<_switches*_elements;i++)
    _A[i] = _B[i] = 0;

  long k = 0, stop_k = ((long)_switches-_based_switches)*_based_switches+(_based_switches * (_based_switches-1)/2), local_sum = 0;
  int local_diameter = 0;
#pragma omp parallel for
  for(int i=0;i<_based_switches;i++){
    unsigned int offset = i*_elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *diameter = 1;
  long level = 2;
  for(int kk=0;kk<_switches;kk++){
    ORP_Matmul_s(_A, _B, _switches, _radix, s_degree, adjacency, _elements, _enable_avx2, _symmetries);

    level++;
#pragma omp parallel for reduction(+:k,local_sum) reduction(max:local_diameter)
    for(int i=0;i<_switches;i++){
      int ib = i%_based_switches;
      int end = (i < _based_switches)? i : _based_switches;
      int ss  = (i < _based_switches)? _symmetries * 2 : _symmetries;
      for(int j=0;j<end;j++){
        int ii = i*_based_switches+j;
        if(_bitmap[ii] == NOT_VISITED){
          if(_B[i*_elements+(j/UINT64_BITS)] & (0x1ULL<<(j%UINT64_BITS))){
            _bitmap[ii] = VISITED;
            k++;
            if(h_degree[ib] != 0 && h_degree[j] != 0){
              local_diameter = MAX(local_diameter, level-2);
              local_sum += level * h_degree[ib] * h_degree[j] * ss;
            }
          }
          else if(_enable_disconnected && (h_degree[ib] == 0 || h_degree[j] == 0)){
            _bitmap[ii] = VISITED;
            k++;
          }
        }
      }
    }
    *diameter = local_diameter;
    if(k == stop_k) break;

    // swap A <-> B
    uint64_t* tmp = _A;
    _A = _B;
    _B = tmp;

    if(kk == _switches-1)
      (*diameter) = _switches;
  }
  
  local_sum = local_sum / 2;
#pragma omp parallel for reduction(+:local_sum)
  for(int i=0;i<_based_switches;i++)
    local_sum += (long)h_degree[i] * (h_degree[i] - 1) * _symmetries;

  *ASPL = local_sum / (double)(((long)_hosts*(_hosts-1))/2);
  *sum  = local_sum;
  *diameter += 2;
}

void ORP_Init_aspl_s(const int hosts, const int switches, const int radix, const int symmetries)
{
  if(hosts % symmetries != 0)
    ERROR("hosts(%d) must be divisible by symmetries(%d)\n", hosts, symmetries);
  else if(switches % symmetries != 0)
    ERROR("switches(%d) must be divisible by symmetries(%d)\n", switches, symmetries);

  _enable_disconnected = CHECK_DISCONNECTED();
    _kind              = ORP_Get_kind(switches, symmetries);
  _based_switches      = switches/symmetries;
  _elements            = (_based_switches+(UINT64_BITS-1))/UINT64_BITS;
#ifdef __AVX2__
  if(_elements >= 4){ // For performance
    _enable_avx2 = true;
    _elements = ((_elements+3)/4)*4;  // _elements must be multiple of 4
  }
#endif
  
  if(_kind == ASPL_MATRIX){
    ORP_Malloc(&_A, switches*_elements*sizeof(uint64_t), _enable_avx2); // uint64_t A[switches][_elements];
    ORP_Malloc(&_B, switches*_elements*sizeof(uint64_t), _enable_avx2); // uint64_t B[switches][_elements];
    _bitmap = malloc(sizeof(char) * switches * _based_switches);        // char _bitmap[switches][_based_switches];
  }
  else{ // _kind == ASPL_BFS
    _frontier = malloc(sizeof(int)  * switches);
    _distance = malloc(sizeof(int)  * switches);
    _next     = malloc(sizeof(int)  * switches);
#ifdef _OPENMP
    ORP_declare_local_frontier(switches);
#endif
  }
  
  _hosts        = hosts;
  _switches     = switches;
  _radix        = radix;
  _symmetries   = symmetries;
  _is_profile   = ORP_Check_profile();
  _elapsed_time = 0;
  _times        = 0;
}

void ORP_Init_aspl(const int hosts, const int switches, const int radix)
{
  ORP_Init_aspl_s(hosts, switches, radix, 1);
}

static void aspl_bfs(const int* restrict h_degree, const int* restrict s_degree, const int* restrict adjacency,
                     int* diameter, long *sum, double* ASPL)
{
  *diameter = 0;
  *sum      = 0;
  for(int s=0;s<_switches;s++){
    bool flag = true;
    if(h_degree[s] == 0) continue;
    
    int num_frontier = 1, level = 1;
    for(int i=0;i<_switches;i++)
      _distance[i] = NOT_USED;
    
    _frontier[0] = s;
    _distance[s] = level;
    
    while(1){
      num_frontier = ORP_top_down_step(level++, num_frontier, adjacency, _switches, _radix, s_degree, 
				       _frontier, _next, _distance);
      if(num_frontier == 0) break;

      int *tmp = _frontier;
      _frontier = _next;
      _next     = tmp;
    }

    if(flag){
      flag = false;
      if(_enable_disconnected){
        for(int i=0;i<_switches;i++){
          if(_distance[i] == NOT_USED && h_degree[i] != 0){
            *diameter = INT_MAX;
            return;
          }
        }
      }
      else{
        for(int i=0;i<_switches;i++){
          if(_distance[i] == NOT_USED){
            *diameter = INT_MAX;
            return;
          }
        }
      }
    }

    for(int i=s+1;i<_switches;i++){
      if(h_degree[i] != 0){
        *sum += (long)(_distance[i] + 2) * h_degree[i] * h_degree[s];
        *diameter = MAX(*diameter, _distance[i]);
      }
    }
  }

  for(int s=0;s<_switches;s++)
    *sum += (long)h_degree[s] * (h_degree[s] - 1);
  
  *ASPL = *sum / (double)(((long)_hosts*(_hosts-1))/2);
  *diameter += 2;
}

static void aspl_bfs_s(const int* restrict h_degree, const int* restrict s_degree, const int* restrict adjacency,
                       int* diameter, long *sum, double* ASPL)
{
  *diameter = 0;
  *sum      = 0;
  for(int s=0;s<_based_switches;s++){
    bool flag =	true;
    if(h_degree[s] == 0) continue;

    int num_frontier = 1, level = 1;
    for(int i=0;i<_switches;i++)
      _distance[i] = NOT_USED;

    _frontier[0] = s;
    _distance[s] = level;

    while(1){
      num_frontier = ORP_top_down_step_s(level++, num_frontier, adjacency, _switches, _radix, s_degree,
                                         _frontier, _next, _distance, _symmetries);
      if(num_frontier == 0) break;

      int *tmp = _frontier;
      _frontier = _next;
      _next     = tmp;
    }

    for(int i=s+1;i<_switches;i++)
      if(h_degree[i%_based_switches] != 0)
        *diameter = MAX(*diameter, _distance[i]);

    if(flag){
      flag = false;
      if(_enable_disconnected){
        for(int i=0;i<_switches;i++){
          if(_distance[i] == NOT_USED && h_degree[i%_based_switches] != 0){
            *diameter = INT_MAX;
            return;
          }
        }
      }
      else{
        for(int i=0;i<_switches;i++){
          if(_distance[i] == NOT_USED){
            *diameter = INT_MAX;
            return;
          }
        }
      }
    }

    for(int i=0;i<_switches;i++)
      if(i!=s)
        *sum += (long)(_distance[i] + 2) * h_degree[i%_based_switches] * h_degree[s];
  }
  
  *sum = *sum * _symmetries / 2;
  for(int s=0;s<_based_switches;s++)
    *sum += (long)h_degree[s] * (h_degree[s] - 1) * _symmetries;
  
  *ASPL = *sum / (double)(((long)_hosts*(_hosts-1))/2);
  *diameter += 2;
}

void ORP_Finalize_aspl()
{
  if(_kind == ASPL_MATRIX){
    ORP_Free(_A, _enable_avx2);
    ORP_Free(_B, _enable_avx2);
    free(_bitmap);
  }
  else{ // _kind == ASPL_BFS
    free(_frontier);
    free(_distance);
    free(_next);
#ifdef _OPENMP
    ORP_free_local_frontier();
#endif
  }
  
  if(_is_profile){
#ifdef _OPENMP
    ORP_Profile("THREADS", _kind, _switches, _symmetries, _elapsed_time, _times);
#else
    ORP_Profile("SERIAL",  _kind, _switches, _symmetries, _elapsed_time, _times);
#endif
  }
}

void ORP_Set_aspl(const int* restrict h_degree, const int* restrict s_degree, const int* restrict adjacency,
                  int *diameter, long *sum, double *ASPL)
{
  double t = ORP_Get_time();

  if(_symmetries == 1){
    if(_kind == ASPL_MATRIX)
      aspl_mat(h_degree, s_degree, adjacency, diameter, sum, ASPL);
    else // _kind == ASPL_MATRIX_BFS
      aspl_bfs(h_degree, s_degree, adjacency, diameter, sum, ASPL);
  }
  else{
    if(_kind == ASPL_MATRIX)
      aspl_mat_s(h_degree, s_degree, adjacency, diameter, sum, ASPL);
    else // _kind == ASPL_MATRIX_BFS
      aspl_bfs_s(h_degree, s_degree, adjacency, diameter, sum, ASPL);
  }
  _elapsed_time += ORP_Get_time() - t;

  if(*diameter >= _switches+2){
    *diameter = INT_MAX;
    *sum      = LONG_MAX;
    *ASPL     = DBL_MAX;
  }
  
  _times++;
}
