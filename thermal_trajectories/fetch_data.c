// Read in prepared list files regarding parts of the cco2 SDF files
// Extract the corresponding thermal trajectory data to new files
// This data includes times, densities, temperatures, and velocities

// Last modified 10/8/18 by Greg Vance

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Fetch lists prepared ahead-of-time for this program
#define PID_FILE "pid_list"
#define SDF_FILE "sdf_list"

// New files this program creates as output
#define TIME_FILE "time_fetched.dat"
#define RHO_FILE "rho_fetched.dat"
#define TEMP_FILE "temp_fetched.dat"
#define VRAD_FILE "vrad_fetched.dat"

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

// Function prototypes
int * read_pid_file(int * n_id);
int get_offset(char * sdf_name);
int get_nobj(char * sdf_name, int offset);
int find_index(int x, int * set, int n_set);
float hypot3d(float x, float y, float z);
float ** alloc_float_2d(int rows, int cols);
void free_float_2d(float ** ptr, int rows, int cols);
void write_outfile(char * name, int * id, float ** dat, int n_id, int n_sdf);

int main()
{
	int n_id, n_sdf, sdf, iter, offset, nobj, k, i;
	int * id;
	FILE * listfp, * timefp, * sdfp;
	char tpos[100], name[500];
	particle part;
	float ** rho, ** temp, ** vrad;
	
	printf("reading input selection files\n");
	
	id = read_pid_file(&n_id);
	
	listfp = fopen(SDF_FILE, "r");
	fscanf(listfp, "n_sdf %d\niter tpos name", &n_sdf);
	
	printf("setting up for reading SDF data\n");
	
	timefp = fopen(TIME_FILE, "w");
	
	rho = alloc_float_2d(n_id, n_sdf);
	temp = alloc_float_2d(n_id, n_sdf);
	vrad = alloc_float_2d(n_id, n_sdf);
	
	printf("reading data from SDF files\n");
	
	for (sdf = 0; sdf < n_sdf; ++sdf)
	{
		fscanf(listfp, "%d %s %s", &iter, tpos, name);
		
		printf("reading SDF file %d/%d\n", sdf + 1, n_sdf);
		
		offset = get_offset(name);
		nobj = get_nobj(name, offset);
		
		fprintf(timefp, "%d %s\n", iter, tpos);
		
		sdfp = fopen(name, "r");
		fseek(sdfp, offset, SEEK_SET);
		
		for (k = 0; k < nobj; ++k)
		{
			fread(&part, sizeof(particle), 1, sdfp);
			
			if ((i = find_index(part.ident, id, n_id)) != -1)
			{
				rho[i][sdf] = part.rho;
				temp[i][sdf] = part.temp;
				vrad[i][sdf] = hypot3d(part.vx, part.vy, part.vz);
			}
		}
		
		fclose(sdfp);
	}
	
	fclose(listfp);
	fclose(timefp);
	
	printf("writing output data files\n");
	
	write_outfile(RHO_FILE, id, rho, n_id, n_sdf);
	write_outfile(TEMP_FILE, id, temp, n_id, n_sdf);
	write_outfile(VRAD_FILE, id, vrad, n_id, n_sdf);
	
	printf("cleaning up\n");
	
	free(id);
	free_float_2d(rho, n_id, n_sdf);
	free_float_2d(temp, n_id, n_sdf);
	free_float_2d(vrad, n_id, n_sdf);
	
	return 0;
}

// Read the list of particle IDs to fetch from file
int * read_pid_file(int * n_id)
{
	int * id;
	FILE * idfp;
	int l;
	
	idfp = fopen(PID_FILE, "r");
	fscanf(idfp, "n_id %d", n_id);
	
	id = malloc(*n_id * sizeof(int));
	if (id == NULL)
	{
		printf("couldn't allocate id array\n");
		exit(1);
	}
	
	for (l = 0; l < *n_id; ++l)
		fscanf(idfp, "%d", &(id[l]));
	
	fclose(idfp);
	return id;
}

// Find the header length of an SDF file given its name
int get_offset(char * sdf_name)
{
	int offset, matched, length;
	char eoh[] = "\n# SDF-EOH";
	FILE * sdfp;
	
	length = strlen(eoh);
	sdfp = fopen(sdf_name, "r");
	
	// Read characters one at a time until the header end is  matched
	for (offset = matched = 0; matched < length; ++offset)
	{
		if (fgetc(sdfp) == eoh[matched])
			++matched;
		else
			matched = 0;
	}
	
	// Account for additional newline characters
	for (; fgetc(sdfp) != '\n'; ++offset);
	
	fclose(sdfp);
	return offset + 1;
}

// Find the number of particles in an SDF file
int get_nobj(char * sdf_name, int offset)
{
	FILE * sdfp;
	int sz;
	
	sdfp = fopen(sdf_name, "r");
	fseek(sdfp, 0, SEEK_END);
	sz = ftell(sdfp);
	
	fclose(sdfp);
	return (sz - offset) / sizeof(particle);
}

// Use binary search to return index of x in set or -1 on failure
int find_index(int x, int * set, int n_set)
{
	int lo, hi, mid;
	
	lo = 0;
	hi = n_set - 1;
	
	while (lo <= hi)
	{
		mid = (lo + hi) / 2;
		
		if (x < set[mid])
			hi = mid - 1;
		else if (x > set[mid])
			lo = mid + 1;
		else
			return mid;
	}
	
	return -1;
}

// Return the sum of three floats added in quadrature
float hypot3d(float x, float y, float z)
{
	return hypotf(hypotf(x, y), z);
}

// Allocate a zeroed 2d array of floats with given size
float ** alloc_float_2d(int rows, int cols)
{
	float ** ptr;
	int r, c;
	
	ptr = malloc(rows * sizeof(float *));
	if (ptr == NULL)
	{
		printf("couldn't allocate rows for 2d float array\n");
		exit(1);
	}
	
	for (r = 0; r < rows; ++r)
	{
		ptr[r] = malloc(cols * sizeof(float));
		if (ptr[r] == NULL)
		{
			printf("couldn't allocate column for 2d float array\n");
			exit(1);
		}
	}
	
	for (r = 0; r < rows; ++r)
		for(c = 0; c < cols; ++c)
			ptr[r][c] = 0.;
	
	return ptr;
}

// Free a previously allocated 2d array of floats
void free_float_2d(float ** ptr, int rows, int cols)
{
	int r;
	
	for (r = 0; r < rows; ++r)
		free(ptr[r]);
	
	free(ptr);
}

// Write rho, temp, or vrad data to an output file
void write_outfile(char * name, int * id, float ** dat, int n_id, int n_sdf)
{
	FILE * ofp;
	int i, t;
	
	ofp = fopen(name, "w");
	
	for (i = 0; i < n_id; ++i)
	{
		fprintf(ofp, "%d", id[i]);
		
		for (t = 0; t < n_sdf; ++t)
			fprintf(ofp, " %.10e", dat[i][t]);
		
		fprintf(ofp, "\n");
	}
	
	fclose(ofp);
}

