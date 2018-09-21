// Read in prepared list files regarding parts of the cco2 SDF files
// Extract the corresponding thermal trajectory data to new files

// Last modified 9/21/18 by Greg Vance

#include <stdio.h>
#include <stdlib.h>

// Fetch lists prepared ahead-of-time for this program
#define PID_FILE "pid_list"
#define SDF_FILE "sdf_list"

// New files this program creates as output
#define TIME_FILE "time_fetched.dat"
#define RHO_FILE "rho_fetched.dat"
#define TEMP_FILE "temp_fetched.dat"

// Particle struct taken from the cco2 SDF file headers
typedef struct {
    double x, y, z;             /* position of body */
    float mass;                 /* mass of body */
    float vx, vy, vz;           /* velocity of body */
    float u;                    /* internal energy */
    float h;                    /* smoothing length */
    float rho;                  /* density */
    float drho_dt;              /* time derivative of rho */
    float udot;                 /* time derivative of u */
    float ax, ay, az;           /* acceleration */
    float lax, lay, laz;        /* acceleration at tpos-dt */
    float phi;                  /* potential */
    float idt;                  /* timestep */
    float pr;           /* pressure */
    unsigned int nbrs;          /* number of neighbors */
    unsigned int ident;         /* unique identifier */
    unsigned int windid;        /* wind id */
    float temp;                 /* temperature */
    float Y_el;                  /* for alignment */
    float mfp;                  /* mean free path */
    float f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20;
} particle;

int * read_pid_file(int * n_id);

int main()
{
	int n_id, n_sdf, sdf, iter;
	int * id;
	FILE * listp, * sdfp;
	char tpos[100], name[500];
	
	id = read_pid_file(&n_id);
	listp = fopen(SDF_FILE, "r");
	
	fscanf(listp, "n_sdf %d\niter tpos name", &n_sdf);
	
	for (sdf = 0; sdf < n_sdf; ++sdf)
	{
		fscanf(listp, "%d %s %s", iter, tpos, name);
		
		
	}
	
	fclose(listp);
	free(id);
	
	return 0;
}

