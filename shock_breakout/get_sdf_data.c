// get_sdf_data.c

// Start by reading a prepared list of all the cco2 SDF files
// Extract the values r, vr, and rho for each particle
// Values for r and vr are calcuated from x, y, z, vx, vy, vz
// Write out one file with this data for each SDF processed

// Last modified 2/25/19 by Greg Vance

#include <stdio.h>
#include <string.h>
#include <math.h>

// List of SDF files prepared ahead-of-time for this program
#define INPUT_LIST "sdf_list"

// Directory for this program to store its output files
#define OUTPUT_DIR "data_from_sdfs/"

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
int get_offset(char * sdf_name);
int get_nobj(char * sdf_name, int offset);
float hypot3d(float x, float y, float z);

int main()
{
	int n_sdf, sdf, iter, offset, nobj, k;
	FILE * listfp, * sdfp, * outfp;
	char tpos[100], sdf_name[500], out_name[500];
	particle part;
	float rad, vrad, rho, mass;
	
	// Open list file and read its header to prepare for iteration
	listfp = fopen(INPUT_LIST, "r");
	fscanf(listfp, "n_sdf %d\niter tpos name", &n_sdf);
	
	for (sdf = 0; sdf < n_sdf; ++sdf)
	{
		printf("now processsing SDF %d/%d...\n", sdf + 1, n_sdf);
		fscanf(listfp, "%d %s %s", &iter, tpos, sdf_name);
		
		// Set up the next output file to recieve the data
		sprintf(out_name, "%ssdf_iter%06d.dat", OUTPUT_DIR, iter);
		outfp = fopen(out_name, "w");
		fprintf(outfp, "iter %d\ntpos %s\n", iter, tpos);
		fprintf(outfp, "rad vrad rho mass\n");
		
		// Prepare to read data from the next SDF file
		offset = get_offset(sdf_name);
		nobj = get_nobj(sdf_name, offset);
		sdfp = fopen(sdf_name, "r");
		fseek(sdfp, offset, SEEK_SET);
		
		for (k = 0; k < nobj; ++k)
		{
			fread(&part, sizeof(particle), 1, sdfp);
			
			rad = hypot3d(part.x, part.y, part.z);
			vrad = hypot3d(part.vx, part.vy, part.vz);
			rho = part.rho;
			mass = part.mass;
			
			fprintf(outfp, "%.10e %.10e %.10e %.10e\n", rad, vrad, rho, mass);
		}
		
		fclose(sdfp);
		fclose(outfp);
	}
	
	fclose(listfp);
	
	return 0;
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

// Return the sum of three floats added in quadrature
float hypot3d(float x, float y, float z)
{
	return hypotf(hypotf(x, y), z);
}

