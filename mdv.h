/***********************************************************************
  particle-collision.h: an include file for particle-collision.c
***********************************************************************/
#define Ratom 1.0         /* Default RGB color of an atom */
#define Gatom 0.0
#define Batom 0.0

#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define RCUT 0.4
#define PI 3.141592653589793
/* Constants for the random number generator */
#define D2P31M 2147483647.0
#define DMUL 16807.0

#include <stdbool.h>

int nlon=18, nlat=9;      /* Number of polygons for a sphere in the longitudinal & lateral directions */
float atom_radius = 0.3;  /* Atomic radius  */
int winx=640, winy=640;   /* Window size */
float min_ext[3], max_ext[3];  /* Range of atomic coordinates:
                             (left,lower,back), (right,top,front) */
float eye[3];             /* position of eye point */
float center[3];          /* position of look reference point */
float up[3];              /* up direction for camera */

/* Starting parameters  *******************************************************/
double DeltaT = 0.005;      /* Size of a time step (in reduced unit) */
int StepLimit = 1000000;      /* Number of time steps to be simulated */

/* Constants ******************************************************************/

double Region[3] = {10,10,10};  /* MD (Simulated usable area) box lengths */
double DeltaTH;                 /* Half the time step */

/* Variables ******************************************************************/

int nAtom = 50;           /* Number of atoms initially set to 50*/
double r[NMAX][3];        /* r[i][0|1|2] is the x|y|z coordinate of atom i */
double rv[NMAX][3];       /* Atomic velocities */
double ra[NMAX][3];       /* Acceleration on atoms */
int stepCount;            /* Current time step */
double gravity = -9.81;   /* Acceleration of gravity */
float camera_dist;        /* Distance from camera to center */
double angle_delta;       /* Amount to rotate camera viewing angle */

/* Interactive Parameters *****************************************************/
int clickedAtom = -1;             /* The atom the user has clicked */
double clickedAtomNewVelocity[3]; /* To be the new velocity of atom once released */
float collisionDamping = 0.8;     /* Damping coefficient for collisions */
bool mapVelToColor = true;        /* Boolean to toggle on/off mapping velocities to 3D color cube */

/******************************************************************************/
