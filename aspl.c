 #include "common.h"
static bool _enable_avx2 = false, _is_profile = false;
static char* _bitmap = NULL;
static int _hosts, _switches, _radix, _symmetries, _kind, _elements, _times;
static int *_frontier = NULL, *_distance = NULL, *_next = NULL;
static uint64_t *_A, *_B;
static double _elapsed_time;

void ORP_Profile(const char* name, const int kind, const int nodes, const int symmetries,
                 const double elapsed_time, const unsigned int times);
void ORP_Matmul(const uint64_t *restrict A, uint64_t *restrict B, const int switches, const int radix,
                const int *restrict ports, const int *restrict adjacency, const int elements,
                const int symmetries, const bool enable_avx2);
void ORP_Malloc(uint64_t **a, const size_t s, const bool enable_avx2);
void ORP_Free(uint64_t *a, const bool enable_avx2);
void ORP_declare_local_frontier(int swithces);
void ORP_free_local_frontier();
bool ORP_Check_profile();
int ORP_Get_kind();
int ORP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
                      const int switches, const int radix, const int* restrict ports, const int symmetries,
                      int* restrict frontier, int* restrict next, int* restrict distance);
extern double ORP_Get_time();

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
#pragma omp parallel for
  for(int i=0;i<_switches/_symmetries;i++){
    unsigned int offset = i*_elements+i/UINT64_BITS;
    _A[offset] = _B[offset] = (0x1ULL << (i%UINT64_BITS));
  }

  *diameter = 1;
  long level = 2;
  for(int kk=0;kk<_switches;kk++){
    ORP_Matmul(_A, _B, _switches, _radix, s_degree, adjacency, _elements, _symmetries, _enable_avx2);

    level++;
#pragma omp parallel for reduction(+:k,local_sum)
    for(int i=0;i<_switches;i++){
      for(int j=i+1;j<_switches;j++){
        int ii = i*_switches+j;
        if(_bitmap[ii] == NOT_VISITED && (_B[i*_elements+(j/UINT64_BITS)] & (0x1ULL<<(j%UINT64_BITS)))){
          _bitmap[ii] = VISITED;
          local_sum += level * h_degree[i] * h_degree[j];
          k++;
        }
      }
    }

    if(k == stop_k) break;

    // swap A <-> B
    uint64_t* tmp = _A;
    _A = _B;
    _B = tmp;

    (*diameter) += 1;
  }

#pragma omp parallel for reduction(+:local_sum)
  for(int i=0;i<_switches;i++)
    local_sum += (long)h_degree[i] * (h_degree[i] - 1);

  *ASPL = local_sum / (double)(((long)_hosts*(_hosts-1))/2);
  *sum  = local_sum;
  *diameter += 2;
}

void ORP_Init_aspl_s(const int hosts, const int switches, const int radix, const int symmetries)
{
  _kind     = ORP_Get_kind();
  _elements = (switches/symmetries+(UINT64_BITS-1))/UINT64_BITS;
#ifdef __AVX2__
  if(_elements >= 4){ // For performance
    _enable_avx2 = true;
    _elements = ((_elements+3)/4)*4;  // _elements must be multiple of 4
  }
#endif
  
  if(_kind == ASPL_MATRIX){
    ORP_Malloc(&_A, switches*_elements*sizeof(uint64_t), _enable_avx2); // uint64_t A[switches][_elements];
    ORP_Malloc(&_B, switches*_elements*sizeof(uint64_t), _enable_avx2); // uint64_t B[switches][_elements];
    _bitmap = malloc(sizeof(char) * switches * switches);               // char _bitmap[switches][switches];
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
  int based_nodes = _switches/_symmetries;
  bool reached = true;
  *diameter = 0;
  *sum      = 0;
  for(int s=0;s<based_nodes;s++){
    if(h_degree[s] == 0) continue;
    
    int num_frontier = 1, level = 0;
    for(int i=0;i<_switches;i++)
      _distance[i] = NOT_USED;
    
    _frontier[0] = s;
    _distance[s] = level;
    
    while(1){
      num_frontier = ORP_top_down_step(level++, num_frontier, adjacency, _switches, _radix, s_degree, 
				       _symmetries, _frontier, _next, _distance);
      if(num_frontier == 0) break;

      int *tmp = _frontier;
      _frontier = _next;
      _next     = tmp;
    }

    *diameter = MAX(*diameter, level-1);

    if(s == 0){
      for(int i=1;i<_switches;i++)
        if(_distance[i] == NOT_USED)
          reached = false;
      
      if(!reached){
        *diameter = INT_MAX;
        return;
      }
    }

    *sum += (long)h_degree[s] * (h_degree[s] - 1);
    for(int i=s+1;i<_switches;i++)
      *sum += (long)(_distance[i] + 3) * h_degree[i] * h_degree[s] * _symmetries;
  }

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

  if(_kind == ASPL_MATRIX)
    aspl_mat(h_degree, s_degree, adjacency, diameter, sum, ASPL);
  else // _kind == ASPL_MATRIX_BFS
    aspl_bfs(h_degree, s_degree, adjacency, diameter, sum, ASPL);

  _elapsed_time += ORP_Get_time() - t;

  if(*diameter >= _switches+2){
    *diameter = INT_MAX;
    *sum      = LONG_MAX;
    *ASPL     = DBL_MAX;
  }
  
  _times++;
}
