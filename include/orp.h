#include <stdbool.h>
typedef struct {
  int u[2], v[2], u_d[2], v_d[2], rand, op;
} ORP_Restore;

void ORP_Conv_adjacency2edge(int hosts, int switches, int radix, int *h_degree, int *s_degree, void *adjacency, void *edge);
void ORP_Conv_edge2adjacency(int hosts, int switches, int radix, int lines, void *edge, void *adjacency);
void ORP_Finalize_aspl();
void ORP_Init_aspl(int hosts, int switches, int radix);
void ORP_Print_adjacency(int hosts, int switches, int radix, int *s_degree, void *adjacency);
void ORP_Print_degree(int hosts, int switches, int degree[switches]);
void ORP_Print_edge(int lines, void *edge);
void ORP_Print_switch(int hosts, int lines, void *edge);
void ORP_Read_edge(char* fname, void *edge);
void ORP_Read_property(char* fname, int* host, int* switches, int* radix, int *lines);
void ORP_Restore_adjacency(ORP_Restore r, int radix, int *h_degree, int *s_degree, void *adjacency);
void ORP_Set_degrees(int hosts, int switches, int lines, void *edge, int h_degree[switches], int s_degree[switches]);
void ORP_Set_host_degree(int hosts, int switches, int lines, void *edge, int h_degree[switches]);
void ORP_Set_switch_degree(int hosts, int switches, int lines, void *edge, int s_degree[switches]);
void ORP_Set_degrees_s(int hosts, int switches, int lines, void *edge, int symmetries, int h_degree[switches/symmetries], int s_degree[switches/symmetries]);
void ORP_Set_host_degree_s(int hosts, int switches, int lines, void *edge, int symmetries, int h_degree[switches/symmetries]);
void ORP_Set_switch_degree_s(int hosts, int switches, int lines, void *edge, int symmetries, int s_degree[switches/symmetries]);
void ORP_Set_lbounds(int hosts, int radix, int *low_diameter, double *low_ASPL);
void ORP_Set_aspl(int* h_degree, int* s_degree, void* adjacency, int *diameter, long *sum, double *ASPL);
void ORP_Srand(unsigned int seed);
void ORP_Swap_adjacency(int switches, int radix, int *s_degree, ORP_Restore *r, void *adjacency);
void ORP_Swing_adjacency(int switches, int radix, int *h_degree, int *s_degree, ORP_Restore *r, void *adjacency);
void ORP_Write_edge(int hosts, int switches, int radix, int lines, void *edge, char *fname);
void ORP_Write_switch(int hosts, int lines, void *edge, char *fname);
void* ORP_Generate_random(int hosts, int switches, int radix, bool assign_evenly, int *lines, int *h_degree, int *s_degree);
bool ORP_Verify_edge(int hosts, int switches, int radix, int lines, void *edge);
int ORP_Optimize_switches(int hosts, int radix);

