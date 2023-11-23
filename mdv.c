/***********************************************************************
  Program mdv.c--ball representation of atoms.
  Required files
    mdv.h:   Include file
    md.conf:   MD configuration file containing atomic coordinates
***********************************************************************/
#include "mdv.h"
#include <stdio.h>
#include <math.h>
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>    /* Header file for the OpenGL library */
#include <OpenGL/glu.h>   /* Header file for the GLu library */
#include <GLUT/glut.h>    /* Header file for the GLut library */
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>                // for gettimeofday()
// #include <stdbool.h>

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

GLuint sphereid;          /* display-list id of atom sphere geom */
GLuint atomsid;           /* display-list id of all atoms */
GLuint boundsid;          /* display-list id for boundary lines */
GLuint buttonsid;
GLdouble fovy, aspect, near_clip, far_clip;  
                          /* parameters for gluPerspective() */
FILE *fp;                 /* pointer to open an MD-configuration file */

/* Function prototypes ************************************************/
void reshape(int, int);
void makeFastNiceSphere(GLuint, double);
void makeBoundaryLines(GLuint);
void makeAtoms(void);
void makeCurframeGeom(void);
void drawScene(void);
void display(void);
void initView(float *, float *);
void readConf(void);
void rotateCamera (double delta);

double SignR(double v,double x) {if (x > 0) return v; else return -v;}
double Dmod(double a, double b) {
	int n;
	n = (int) (a/b);
	return (a - b*n);
}
double RandR(double *seed) {
	*seed = Dmod(*seed*DMUL,D2P31M);
	return (*seed/D2P31M);
}
void RandVec3(double *p, double *seed) {
	double x,y,s = 2.0;
	while (s > 1.0) {
		x = 2.0*RandR(seed) - 1.0; y = 2.0*RandR(seed) - 1.0; s = x*x + y*y;
	}
	p[2] = 1.0 - 2.0*s; s = 2.0*sqrt(1.0 - s); p[0] = s*x; p[1] = s*y;
}

double dot_product(double *v, double *u, int n)
{
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += (double)v[i] * (double)u[i];
    return result;
}

double distance(double *p1, double *p2, int n) {
  double dist = 0.0;
  for (int i = 0; i < n; i++){
    dist = sqrt(pow(p1[0]-p2[0],2.0)+pow(p1[1]-p2[1],2.0)+pow(p1[2]-p2[2],2.0));
  }
  return dist;
}

// void InitParams();
void InitConf();
// void ComputeAccel();
void ComputeAcceleration();
void SingleStep();
void HalfKick();
void ApplyBoundaryCond();
// void EvalProps();

void animate(void);
// void moveAtom(void);

float * mapVelocityToColor();

/*----------------------------------------------------------------------------*/
// void InitParams() {
// /*------------------------------------------------------------------------------
// 	Initializes parameters.
// ------------------------------------------------------------------------------*/
// 	int k;
// 	double rr,ri2,ri6,r1;

// 	/* Reads control parameters */
// 	// scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
// 	// scanf("%le",&Density);
// 	// scanf("%le",&InitTemp);
// 	// scanf("%le",&DeltaT);
// 	// scanf("%d",&StepLimit);
// 	// scanf("%d",&StepAvg);

// 	/* Computes basic parameters */
// 	DeltaTH = 0.5*DeltaT;
// 	for (k=0; k<3; k++) {
// 		Region[k] = InitUcell[k]/pow(Density/4.0,1.0/3.0);
// 		RegionH[k] = 0.5*Region[k];

//     // printf("REGION: %f\n", Region[k]);
// 	}



// 	/* Constants for potential truncation */
// 	rr = RCUT*RCUT; ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1=sqrt(rr);
// 	Uc = 4.0*ri6*(ri6 - 1.0);
// 	Duc = -48.0*ri6*(ri6 - 0.5)/r1;
// }

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to face-centered cubic (fcc) lattice positions.  
	rv are initialized with a random velocity corresponding to Temperature.  
------------------------------------------------------------------------------*/
	// double c[3],gap[3],e[3],vSum[3],vMag;

  // double rLoc; 

	int j,n,k;
	double seed, e[3];
	/* FCC atoms in the original unit cell */
	// double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
	//                          {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

	/* Sets up a face-centered cubic (fcc) lattice */
	// for (k=0; k<3; k++) gap[k] = Region[k]/InitUcell[k];

  // printf("The gap is: %f, %f, %f \n", gap[0], gap[1], gap[2]);

  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);


  // for (k=0; k<3; k++) {
  //   rLoc[k] = ((double)rand()/(double)(RAND_MAX)) * Region[k];
  //   printf("rLoc[%f] = %f", k, rLoc[k]);
  // }


  int nAtomTotal = 50;
  nAtom= 0;
  for (int i = 0; i<nAtomTotal;i++){
    // printf("In first loop\n");
    for (k=0; k<3; k++) {
      // rLoc = ((double)rand()/(double)(RAND_MAX)) * Region[k];
      // printf("rLoc[%d] = %f\n", k, rLoc);
      // printf("REGION: %f\n", Region[k]);
      r[nAtom][k] = ((double)rand()/(double)(RAND_MAX)) * Region[k];
      // r[nAtom][k] = rLoc;
      // printf("r[%d][%d] = %f \n", i, k, r[nAtom][k]);

    }
    nAtom++;
  }





	// nAtom = 0;
	// for (nZ=0; nZ<InitUcell[2]; nZ++) {
	// 	c[2] = nZ*gap[2];
	// 	for (nY=0; nY<InitUcell[1]; nY++) {
	// 		c[1] = nY*gap[1];
	// 		for (nX=0; nX<InitUcell[0]; nX++) {
	// 			c[0] = nX*gap[0];
	// 			for (j=0; j<4; j++) {
	// 				for (k=0; k<3; k++) {
	// 					// r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
  //           r[nAtom][k] = ((double)rand()/(double)(RAND_MAX)) * Region[k];
  //         }
	// 				++nAtom;
	// 			}
	// 		}
	// 	}
	// }

	/* Generates random velocities */
	seed = 13597.0;
	// vMag = sqrt(3*InitTemp);
	// for(k=0; k<3; k++) vSum[k] = 0.0;
	for(n=0; n<nAtom; n++) {
		RandVec3(e,&seed);
		for (k=0; k<3; k++) {
			rv[n][k] = 20*e[k];
      // rv[n][k] = ((float)rand()/(float)(RAND_MAX)) * 10.0;
      // printf("rv[n][k] = %f\n", rv[n][k]);
			// vSum[k] = vSum[k] + rv[n][k];
		}
	}
	/* Makes the total momentum zero */
	// for (k=0; k<3; k++) vSum[k] = vSum[k]/nAtom;
	// for (n=0; n<nAtom; n++) for(k=0; k<3; k++) rv[n][k] = rv[n][k] - vSum[k];
}

void ComputeAcceleration() {
  int n,k;

  for (n=0; n<nAtom; n++) {

    if (n != clickedAtom) {

      ra[n][0] = 0;
      ra[n][1] = gravity;
      // ra[n][1] = 0;
      ra[n][2] = 0;
    }

    // printf("In acceleration, position is: %f, %f, %f\n", r[n][0], r[n][1], r[n][2]);
  }

}

// /*----------------------------------------------------------------------------*/
// void ComputeAccel() {
// /*------------------------------------------------------------------------------
// 	Acceleration, ra, are computed as a function of atomic coordinates, r,
// 	using the Lennard-Jones potential.  The sum of atomic potential energies,
// 	potEnergy, is also computed.   
// ------------------------------------------------------------------------------*/
// 	double dr[3],f,fcVal,rrCut,rr,ri2,ri6,r1;
// 	int j1,j2,n,k;

// 	rrCut = RCUT*RCUT;
// 	for (n=0; n<nAtom; n++) for (k=0; k<3; k++) ra[n][k] = 0.0;
// 	potEnergy = 0.0;

// 	/* Doubly-nested loop over atomic pairs */
// 	for (j1=0; j1<nAtom-1; j1++) {
// 		for (j2=j1+1; j2<nAtom; j2++) {
// 			/* Computes the squared atomic distance */
// 			for (rr=0.0, k=0; k<3; k++) {
// 				dr[k] = r[j1][k] - r[j2][k];
// 				/* Chooses the nearest image */
// 				dr[k] = dr[k] - SignR(RegionH[k],dr[k]-RegionH[k])
// 				              - SignR(RegionH[k],dr[k]+RegionH[k]);
// 				rr = rr + dr[k]*dr[k];
// 			}
// 			/* Computes acceleration & potential within the cut-off distance */
// 			if (rr < rrCut) {
// 				ri2 = 1.0/rr; 
//         ri6 = ri2*ri2*ri2; 
//         r1 = sqrt(rr);
// 				fcVal = 48.0*ri2*ri6*(ri6-0.5) + Duc/r1;
// 				for (k=0; k<3; k++) {
// 					f = fcVal*dr[k];
// 					ra[j1][k] = ra[j1][k] + f;
// 					ra[j2][k] = ra[j2][k] - f;
// 				}
// 				potEnergy = potEnergy + 4.0*ri6*(ri6-1.0) - Uc - Duc*(r1-RCUT);
// 			} 
// 		} 
//     // ra[j1][1]+= -9.81;  // Add gravity
// 	}
// }

/*----------------------------------------------------------------------------*/
void SingleStep() {
/*------------------------------------------------------------------------------
	r & rv are propagated by DeltaT in time using the velocity-Verlet method.
------------------------------------------------------------------------------*/
	int n,k;

	HalfKick(); /* First half kick to obtain v(t+Dt/2) */
	for (n=0; n<nAtom; n++) /* Update atomic coordinates to r(t+Dt) */
		for (k=0; k<3; k++) r[n][k] = r[n][k] + DeltaT*rv[n][k];
	ApplyBoundaryCond();
	// ComputeAccel(); /* Computes new accelerations, a(t+Dt) */
  ComputeAcceleration();
	HalfKick(); /* Second half kick to obtain v(t+Dt) */
}

/*----------------------------------------------------------------------------*/
void HalfKick() {
/*------------------------------------------------------------------------------
	Accelerates atomic velocities, rv, by half the time step.
------------------------------------------------------------------------------*/
	// int n,k;
	// for (n=0; n<nAtom; n++)
	// 	for (k=0; k<3; k++) rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];

  DeltaTH = 0.5*DeltaT;

  int n,k;
	for (n=0; n<nAtom; n++) {
		for (k=0; k<3; k++) {
       rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];

      //  if (n == clickedAtom)
      //     printf("For atom %d, DeltaTH: %f, accel: %f, velocity:%f\n", n, DeltaTH, ra[n][k], rv[n][k]);
    }
  }

  
}

/*----------------------------------------------------------------------------*/
void ApplyBoundaryCond() {
/*------------------------------------------------------------------------------
	Applies periodic boundary conditions to atomic coordinates.
------------------------------------------------------------------------------*/
	// double rv_init[NMAX][3];
  // double totSystemSpeedStart[3];
  // // for (n=0; n<nAtom; n++) {
  //   for (int k=0; k<3; k++) {
  //     totSystemSpeedStart[k] = 0;
  //   }
  // // }

  // for (int n=0; n<nAtom; n++) {
  //   for (int k=0; k<3; k++) {
  //     totSystemSpeedStart[k] += fabs(rv[n][k]);
  //     rv_init[n][k] = rv[n][k];
  //   }
  // }

  // printf("\nTotal system velocity start: %f,%f,%f = %f \n\n", totSystemVel[0], totSystemVel[1], totSystemVel[2], (totSystemVel[0] + totSystemVel[1] + totSystemVel[2]));

  
  
  
  
  
  
  int n,k;
  int toPrint = 0;
  // float collisionDamping = 1;
  // float collisionDamping = 0.5;
  // float collisionDamping = 0.8;
  // float collisionDamping = 0.95;

	for (n=0; n<nAtom; n++) {
    // if (r[n][0] > 5.129928 || r[n][1] > 5.129928 || r[n][2] > 5.129928 ||
    //   r[n][0] < 0 || r[n][1] < 0 || r[n][2] < 0){
      // printf("Position Before Boundary condition: %f, %f, %f\n", r[n][0], r[n][1], r[n][2]);
    //   toPrint = 1;
    // }

    // Check box collision
    // printf("Particle %d position is: %f,%f,%f \n\n", n, r[n][0], r[n][1], r[n][2]);
		for (k=0; k<3; k++)  {
      // r[n][k] = r[n][k] - SignR(RegionH[k],r[n][k])
			//                   - SignR(RegionH[k],r[n][k]-Region[k]);

      if (r[n][k] + atom_radius > Region[k] || r[n][k] - atom_radius < 0) {
        // printf("Particle %d wall collision\n", n);

        if (r[n][k] - atom_radius < 0)
          r[n][k] = 0.0 + atom_radius;
        else
          r[n][k] = Region[k] - atom_radius;

        // float rv_old = rv[n][k];
        rv[n][k] = rv[n][k] * -1.0 * collisionDamping;

        // ra[n][k] = (rv[n][k] - rv_old)/DeltaT;
        // ra[n][1] = ra[n][1] - 9.81;
      }

    }

    // double m1 = 10.0;
    // double m2 = 10.0;
    // double force;
    // long double norm[3];
    // long double relative_v[3];
    // long double vel_norm_dot;
    // long double dot_scalar_norm[3];
    // long double norm_sum = 0;

    // Check particle collision
    // for (int j1=0; j1<nAtom-1; j1++) {
		  for (int n2=n+1; n2<nAtom; n2++) {

        int p1Stationary = 1;
        int p2Stationary = 1;
        for (k=0; k<3; k++)  {
          if (rv[n][k] != 0)
            p1Stationary = 0;
          if (rv[n2][k] != 0)
            p2Stationary = 0;
        }
        if (p1Stationary && p2Stationary)
          continue;


        long double norm[3];
        long double relative_v[3];
        long double vel_norm_dot = 0;
        long double dot_scalar_norm[3];
        long double norm_sum = 0;

        // Check wall collision for second particle first
        for (k=0; k<3; k++)  {

          if (r[n2][k] + atom_radius > Region[k] || r[n2][k] - atom_radius < 0) {
            // printf("Particle %d wall collision\n", n2);

            if (r[n2][k] - atom_radius < 0)
              r[n2][k] = 0.0 + atom_radius;
            else
              r[n2][k] = Region[k] - atom_radius;

            rv[n2][k] = rv[n2][k] * -1.0 * collisionDamping;
          }

        }



        long double dst = sqrt(pow(r[n][0]-r[n2][0],2.0)+pow(r[n][1]-r[n2][1],2.0)+pow(r[n][2]-r[n2][2],2.0));
        // printf("Particles have distance %Lf \n", dst);
        if (dst <= 2 * atom_radius) {
          // printf("Particles %d and %d with distance %Lf are colliding \n", n, n2, dst);
          // printf("Particle %d has position %f,%f,%f and velocity %f,%f,%f  \n", n, r[n][0], r[n][1], r[n][2], rv[n][0], rv[n][1], rv[n][2]);
          // printf("Particle %d has position %f,%f,%f and velocity %f,%f,%f  \n", n2, r[n2][0], r[n2][1], r[n2][2], rv[n2][0], rv[n2][1], rv[n2][2]);


          // Get norm vector between two particles
          for (k=0; k<3; k++)  {
            norm[k] = r[n][k] - r[n2][k];
            norm_sum += pow(norm[k],2);

            if (n == clickedAtom)
              relative_v[k] = clickedAtomNewVelocity[k] - rv[n2][k];
            else if (n2 == clickedAtom)
              relative_v[k] = rv[n][k] - clickedAtomNewVelocity[k];
            else
              relative_v[k] = rv[n][k] - rv[n2][k];
          }
          for (k=0; k<3; k++)  {
            norm[k] = norm[k] / sqrt(norm_sum);
          }

          



          double rv_mag = sqrt(pow(rv[n][0],2.0)+pow(rv[n][1],2.0)+pow(rv[n][2],2.0));
          // printf("rv_mag = %f\n", rv_mag);

          double rv_mag2 = sqrt(pow(rv[n2][0],2.0)+pow(rv[n2][1],2.0)+pow(rv[n2][2],2.0));
          // printf("rv_mag2 = %f\n", rv_mag2);

          double distToMove = (2*atom_radius - dst);
          // printf("distToMove = %f\n", distToMove);

          // double new_rv[3];
          // double new_rv2[3];

          // new_rv[0] = (2*atom_radius - dst)*(rv[n][0]/rv_mag);
          // new_rv[1] = (2*atom_radius - dst)*(rv[n][1]/rv_mag);
          // new_rv[2] = (2*atom_radius - dst)*(rv[n][2]/rv_mag);


          // new_rv2[0] = (2*atom_radius - dst)*(rv[n2][0]/rv_mag2);
          // new_rv2[1] = (2*atom_radius - dst)*(rv[n2][1]/rv_mag2);
          // new_rv2[2] = (2*atom_radius - dst)*(rv[n2][2]/rv_mag2);


          // double rv_mag_new = sqrt(pow(new_rv[0],2.0)+pow(new_rv[1],2.0)+pow(new_rv[2],2.0));
          // printf("rv_mag_new = %f\n", rv_mag_new);

          // double rv_mag2_new = sqrt(pow(new_rv2[0],2.0)+pow(new_rv2[1],2.0)+pow(new_rv2[2],2.0));
          // printf("rv_mag2_new = %f\n", rv_mag2_new);


          double shift_dist = 2.0*atom_radius - dst;
          long double p1_move_vec[3];
          long double p2_move_vec[3];

          for (k=0; k<3; k++)  {
            p1_move_vec[k] = norm[k] * shift_dist;
            p2_move_vec[k] = norm[k] * shift_dist * -1;
          }

          int moveP1 = 1;

          for (k=0; k<3; k++)  {
            if (r[n][k] + p1_move_vec[k] > Region[k] || r[n][k] + p1_move_vec[k]) {
              moveP1 = 1;
            }
          }

          for (k=0; k<3; k++)  {
            if (moveP1){
              r[n][k] += p1_move_vec[k];
            }
            else {
              r[n2][k] -= p2_move_vec[k];
            }
          }


          long double newdst = sqrt(pow(r[n][0]-r[n2][0],2.0)+pow(r[n][1]-r[n2][1],2.0)+pow(r[n][2]-r[n2][2],2.0));
          // printf("Particles %d and %d new distance %Lf \n\n", n, n2, newdst);
          




          // printf("norm = %f,%f,%f  \n", norm[0], norm[1], norm[2]);
          // printf("relative_v = %f,%f,%f  \n", relative_v[0], relative_v[1], relative_v[2]);

          for (k=0; k<3; k++)  {
            // printf("vel_norm_dot = %Lf\n", vel_norm_dot);
            vel_norm_dot += relative_v[k]*norm[k];
          }
          for (k=0; k<3; k++)  {
            dot_scalar_norm[k] = vel_norm_dot*norm[k];
          }
          // printf("norm = %Lf,%Lf,%Lf \n", norm[0], norm[1], norm[2]);
          // printf("relative_v = %Lf,%Lf,%Lf \n", relative_v[0], relative_v[1], relative_v[2]);
          // printf("vel_norm_dot = %Lf  \n", vel_norm_dot);
          // printf("dot_scalar_norm = %Lf,%Lf,%Lf  \n", dot_scalar_norm[0], dot_scalar_norm[1], dot_scalar_norm[2]);


          // TODO: Add in logic so that if one of the collision atoms is being clicked, only the other
          // atom get the velocity transfer

          // double velocity_error[3];

          if (n == clickedAtom) {
            // printf("Atom n (%d) is the clicked atom and velocity is %f,%f,%f\n", n, clickedAtomNewVelocity[0], clickedAtomNewVelocity[1], clickedAtomNewVelocity[2]);
            // printf("Atom n2 (%d) is the other atom and velocity is %f,%f,%f\n", n2, rv[n2][0], rv[n2][1], rv[n2][2]);
            // printf("dot_scalar_norm = %Lf,%Lf,%Lf \n", dot_scalar_norm[0], dot_scalar_norm[1], dot_scalar_norm[2]);
            for (k=0; k<3; k++)  {
              rv[n2][k] = (rv[n2][k] + dot_scalar_norm[k]) * collisionDamping;
            }
            // printf("Atom n2 (%d) new velocity is %f,%f,%f\n", n2, rv[n2][0], rv[n2][1], rv[n2][2]);
            
          }
          else if (n2 == clickedAtom) {
            // printf("Atom n2 (%d) is the clicked atom and velocity is %f,%f,%f\n", n2, clickedAtomNewVelocity[0], clickedAtomNewVelocity[1], clickedAtomNewVelocity[2]);
            // printf("Atom n (%d) is the other atom and velocity is %f,%f,%f\n", n, rv[n][0], rv[n][1], rv[n][2]);
            for (k=0; k<3; k++)  {
              
              rv[n][k] = (rv[n][k] - dot_scalar_norm[k]) * collisionDamping;
            }
            // printf("Atom n (%d) new velocity is %f,%f,%f\n", n, rv[n][0], rv[n][1], rv[n][2]);
          }
                  

          else {
            for (k=0; k<3; k++)  {
              long double rvn_old = rv[n][k];
              long double rvn2_old = rv[n2][k];

              rv[n][k] = (rv[n][k] - dot_scalar_norm[k]) * collisionDamping;
              rv[n2][k] = (rv[n2][k] + dot_scalar_norm[k]) * collisionDamping;


              // velocity_error[k] = (rvn_old + rvn2_old) - (rv[n][k] + rv[n2][k]);
              
            }
          }



          rv_mag = sqrt(pow(rv[n][0],2.0)+pow(rv[n][1],2.0)+pow(rv[n][2],2.0));
          // printf("rv_mag = %f\n", rv_mag);

          rv_mag2 = sqrt(pow(rv[n2][0],2.0)+pow(rv[n2][1],2.0)+pow(rv[n2][2],2.0));
          // printf("rv_mag2 = %f\n", rv_mag2);

          // If velocity gets low enough, set to zero
          if (rv_mag < 0.5) {
            rv[n][0] = 0;
            rv[n][1] = 0;
            rv[n][2] = 0;
          }

          if (rv_mag2 < 0.5) {
            rv[n2][0] = 0;
            rv[n2][1] = 0;
            rv[n2][2] = 0;
          }

        }
      
      }

  }
}


// /*----------------------------------------------------------------------------*/
// void EvalProps() {
// /*------------------------------------------------------------------------------
// 	Evaluates physical properties: kinetic, potential & total energies.
// ------------------------------------------------------------------------------*/
// 	double vv;
// 	int n,k;

// 	kinEnergy = 0.0;
// 	for (n=0; n<nAtom; n++) {
// 		vv = 0.0;
// 		for (k=0; k<3; k++)
// 			vv = vv + rv[n][k]*rv[n][k];
// 		kinEnergy = kinEnergy + vv;
// 	}
// 	kinEnergy *= (0.5/nAtom);
// 	potEnergy /= nAtom;
// 	totEnergy = kinEnergy + potEnergy;
// 	temperature = kinEnergy*2.0/3.0;

// 	/* Print the computed properties */
// 	// printf("%9.6f %9.6f %9.6f %9.6f\n",
// 	// stepCount*DeltaT,temperature,potEnergy,totEnergy);
// }

void animate() { /* Callback function for idle events */
    /* Keep updating the scene until the last MD step is reached */
    if (stepCount <= StepLimit) {
        SingleStep(); /* One MD-step integration */
        // if (stepCount%StepAvg == 0) EvalProps(); 
        makeCurframeGeom(); /* Redraw the scene (make a display list) */
        glutPostRedisplay(); 
        ++stepCount;
    }
}

void mouseClick(int button, int state, int x, int y){
  // printf("Mouse was clicked\n");
  // printf("button = %d, state = %d, x = %d, y = %d\n", button, state, x, y);

  // Mouse unclicked so turn gravity back on for particle that was held and set velocity
  if (state) {
    if (clickedAtom > -1) {
      
      

      // printf("Resetting velocity to be: %f,%f,%f \n", clickedAtomNewVelocity[0], clickedAtomNewVelocity[1], clickedAtomNewVelocity[2]);
      for (int k =0; k<3;k++){
        
        rv[clickedAtom][k] = clickedAtomNewVelocity[k];
      }

      ra[clickedAtom][1] = gravity;

      // Reset clickedAtom
      clickedAtom = -1;
    }
  }

  else if (state == 0) {
    double modelMatrix[16];
    double projMatrix[16];
    GLint viewport[4];
    GLdouble mouseVec_start[3];
    GLdouble mouseVec_end[3];

    clickedAtom = -1;

    glGetIntegerv( GL_VIEWPORT, viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
    glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
    GLint realy = viewport[3] - (GLint) y;
    gluUnProject( (GLdouble) x, (GLdouble) realy, 0.0, modelMatrix,
        projMatrix, viewport, &mouseVec_start[0], &mouseVec_start[1], &mouseVec_start[2] );
    gluUnProject( (GLdouble) x, (GLdouble) realy, 1.0, modelMatrix,
        projMatrix, viewport, &mouseVec_end[0], &mouseVec_end[1], &mouseVec_end[2] );

    // printf("mouseVec_start = %f, %f,%f \n", mouseVec_start[0], mouseVec_start[1], mouseVec_start[2]);
    // printf("mouseVec_end = %f,%f,%f \n", mouseVec_end[0], mouseVec_end[1], mouseVec_end[2]);

    double mouseVec[3];
    for (int k=0; k <3;k++){
      mouseVec[k] = mouseVec_end[k] - mouseVec_start[k];
    }  

    double mouseVec_square = dot_product(mouseVec, mouseVec, 3); 

    // printf("mouseVec = %f,%f,%f \n", mouseVec[0], mouseVec[1], mouseVec[2]);

    double mouseStart_toAtom[3];
    double mouseVect_mouseAtomVec_dot[3];
    double mousePoint_mouseVect_proj[3];
    double closestPoint[3];
    int mouseRayAtomCollision = 0;

    for (int n=0; n<nAtom; n++) {
      for (int k=0; k <3;k++){
        mouseStart_toAtom[k] = r[n][k] - mouseVec_start[k];
      }

      double mouseVect_mouseAtomVec_dot = dot_product(mouseStart_toAtom, mouseVec, 3); 

      for (int k=0; k <3;k++){
        mousePoint_mouseVect_proj[k] = mouseVect_mouseAtomVec_dot / mouseVec_square;
      }

      for (int k=0; k <3;k++){
        closestPoint[k] = mouseVec_start[k] + mouseVec[k] * mousePoint_mouseVect_proj[k];
      }

      // printf("closestPoint = %f,%f,%f\n", closestPoint[0], closestPoint[1], closestPoint[2]);

      double len = distance(closestPoint, r[n], 3); 

      // printf("len = %f\n", len);

      mouseRayAtomCollision = 0;
      double epsilon = 0.01;
      if (len < atom_radius + epsilon) {
        mouseRayAtomCollision = 1;
        clickedAtom = n;
        break;
        
      }

    }

    if (mouseRayAtomCollision) {
      // printf("The mouse is grabbing particle %d\n", clickedAtom);
      for(int k=0; k<3;k++) {
        rv[clickedAtom][k] = 0;
        ra[clickedAtom][k] = 0;
      }
    }
    // else
      // printf("The mouse is not grabbing any particle\n");
  }



    // Vec3d AP = P - A; 
    // double ap_dot_ab = DotProduct(AP, AB); 
    // // t is a projection param when we project vector AP onto AB 
    // *t = ap_dot_ab / ab_square; 
    // // calculate the closest point 
    // Vec3d Q = A + AB * (*t); 



  // int n;
  // int k;

  // double minX = 100;
  // double minY = 100;
  // double minDist = 100;
  // int closestAtom;

  // double new_x = x / (winx/2.0) - 1.0 ;
  // double new_y = -1 * (y / (winy/2.0) - 1.0);

  // printf("new_x = %f, new_y = %f \n", new_x, new_y);


  // for (n=0; n<nAtom; n++) {
  //   double xdiff = r[n][0] - x;
  //   double ydiff = r[n][1] - y;
  //   double dist = sqrt(pow(xdiff, 2) + pow(ydiff, 2));
  //    if (dist < minDist) {
  //     closestAtom = n;
  //    }
  // }

  // printf("closes atom is number %d with position: %f,%f,%f\n", closestAtom, r[closestAtom][0], r[closestAtom][1], r[closestAtom][2]);
}


void moveAtom(int x, int y) {

  double atomStartPos[3];


  if (clickedAtom > -1) {

    // TIMER  //
    struct timeval t1, t2;
    double elapsedTime;

    // start timer
    gettimeofday(&t1, NULL);

    double modelMatrix[16];
    double projMatrix[16];
    GLint viewport[4];
    GLdouble mouseVec_start[3];
    GLdouble mouseVec_end[3];

    glGetIntegerv( GL_VIEWPORT, viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
    glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
    GLint realy = viewport[3] - (GLint) y;
    gluUnProject( (GLdouble) x, (GLdouble) realy, 0.0, modelMatrix,
        projMatrix, viewport, &mouseVec_start[0], &mouseVec_start[1], &mouseVec_start[2] );
    gluUnProject( (GLdouble) x, (GLdouble) realy, 1.0, modelMatrix,
        projMatrix, viewport, &mouseVec_end[0], &mouseVec_end[1], &mouseVec_end[2] );

    // printf("mouseVec_start = %f, %f,%f \n", mouseVec_start[0], mouseVec_start[1], mouseVec_start[2]);
    // printf("mouseVec_end = %f,%f,%f \n", mouseVec_end[0], mouseVec_end[1], mouseVec_end[2]);

    double mouseVec[3];
    double newLoc[3];

    for (int k=0; k <3;k++){
      mouseVec[k] = mouseVec_end[k] - mouseVec_start[k];
    }  

    // printf("Currently grabbed particle position: %f,%f,%f \n", r[clickedAtom][0], r[clickedAtom][1], r[clickedAtom][2]);

    GLdouble obj_win_coord[3];
    gluProject(r[clickedAtom][0], r[clickedAtom][1], r[clickedAtom][2], modelMatrix,
        projMatrix, viewport, &obj_win_coord[0], &obj_win_coord[1], &obj_win_coord[2]);

    // printf("Obj pos in win space is %f,%f,%f \n", obj_win_coord[0], obj_win_coord[1], obj_win_coord[2]);


    gluUnProject( (GLdouble) x, (GLdouble) realy, obj_win_coord[2], modelMatrix,
        projMatrix, viewport, &newLoc[0], &newLoc[1], &newLoc[2] );

    
    for(int k =0; k <3;k++){
      atomStartPos[k] = r[clickedAtom][k];
      r[clickedAtom][k]= newLoc[k];
    }



    // stop timer
    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    // printf("%f ms.\n", elapsedTime);

    // Get velocity
    double vel[3];
    for(int k =0; k <3;k++){
        clickedAtomNewVelocity[k] = (r[clickedAtom][k] - atomStartPos[k]) / (elapsedTime);

    }

  }

}

void MyKeyboardFunc(unsigned char Key, int x, int y)
{
  // printf("Key pressed: %c \n", Key);

  switch(Key)
  {
    case 'g': // Add to gravity value
      gravity += 10;
      printf("gravity value is now %f \n", gravity);
      break;
    case 'h': // Subtract from gravity value
      gravity -= 10;
      printf("gravity value is now %f \n", gravity);
      break;
    case 'z': // Increase the size of the particles
      atom_radius += 0.1; 
      sphereid = glGenLists(1);
      makeFastNiceSphere(sphereid,atom_radius);
      makeAtoms();
      break;
    case 'x': // Reduce the size of the particles
      atom_radius -= 0.1;
      sphereid = glGenLists(1);
      makeFastNiceSphere(sphereid,atom_radius);
      makeAtoms();
      break;
    case 'c': // Add a particle to the system from top (or botton if gravity is flipped) with random downward velocity
      nAtom++; 
      r[nAtom-1][0]= ((double)rand()/(double)(RAND_MAX)) * Region[0];
      r[nAtom-1][1]= Region[1];
      r[nAtom-1][2]= ((double)rand()/(double)(RAND_MAX)) * Region[2];

      if (gravity > 0)
        r[nAtom-1][1]=r[nAtom-1][1]*-1;

      rv[nAtom-1][1]= ((double)rand()/(double)(RAND_MAX)) * Region[1];

      makeAtoms();
      break;
    case 'v': // Remove a particle from system
      nAtom--; 
      makeAtoms();
      break;
    case 'q': // Quit program
      exit(1);
      break;
    case 'l': // Pan camera view clockwise
      angle_delta = 0.1;
      rotateCamera(angle_delta);
    break;
    case 'k': // Pan camera view counter-clockwise
      angle_delta = -0.1;
      rotateCamera(angle_delta);
    break;
    case 'd': // Decrease damping coefficient (<1 will make particles slow down with each collision)
      if (collisionDamping - 0.1 > 0)
        collisionDamping -= 0.1;
        printf("New collisionDamping: %f\n", collisionDamping);
    break;
    case 'f': // Increase damping coefficient (>1 will make particles speed up with each collision)
      if (collisionDamping + 0.1 <= 1.5)
        collisionDamping += 0.1;
        printf("New collisionDamping: %f\n", collisionDamping);
    break;
    case 'm': // Toggle on/off mapping the particle velocity to 3D color cube
      mapVelToColor = !mapVelToColor;
    break;
    
  };
  glutPostRedisplay();
}

void rotateCamera (double delta) {
  double thetaX = acos((eye[0] - center[0]) / camera_dist);
  double thetaZ = asin((eye[2] - center[2]) / camera_dist);

  if (thetaZ < 0)
    thetaX *= -1;

  if (thetaX + delta > PI && thetaX <= PI) {
    thetaX = thetaX + delta;
    thetaX = -1 * (2*PI - thetaX);
  }
  else {
    thetaX = thetaX + delta;
  }
  
  double newX = center[0] + camera_dist * cos(thetaX);
  double newZ = center[2] + camera_dist * sin(thetaX);

  eye[0] = newX;
  eye[2] = newZ;

  gluLookAt(
  (GLdouble)eye[0],(GLdouble)eye[1],(GLdouble)eye[2],
  (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
  (GLdouble)up[0],(GLdouble)up[1],(GLdouble)up[2]);
}

/**********************************************************************/
void reshape (int w, int h) {
/***********************************************************************
  Callback for glutReshapeFunc()
***********************************************************************/
  /* set the GL viewport to match the full size of the window */
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  aspect = w/(float)h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy,aspect,near_clip,far_clip);
  glMatrixMode(GL_MODELVIEW);
}

/**********************************************************************/
void makeFastNiceSphere(GLuint listid, double radius) {
/***********************************************************************
Called once to generate and compile sphere geometry into the given
display list id.
***********************************************************************/
  int i,j;
  float lon,lat;
  float loninc,latinc;
  float x,y,z;

  loninc = 2*M_PI/nlon;
  latinc = M_PI/nlat;

  glNewList(listid,GL_COMPILE);

  /* South-pole triangular fan */
  glBegin(GL_TRIANGLE_FAN);
  glNormal3f(0,-1,0);
  glVertex3f(0,-radius,0);
  lon = 0;
  lat = -M_PI/2 + latinc;
  y = sin(lat);
  for (i=0; i<=nlon; i++) {
    x = cos(lon)*cos(lat);
    z = -sin(lon)*cos(lat);
    glNormal3f(x,y,z);
    glVertex3f(x*radius,y*radius,z*radius);
    lon += loninc;
  }
  glEnd();

  /* Quadrilateral stripes to cover the sphere */
  for (j=1; j<nlat-1; j++) {
    lon = 0;
    glBegin(GL_QUAD_STRIP);
      for (i=0; i<=nlon; i++) {
        x = cos(lon)*cos(lat);
        y = sin(lat);
        z = -sin(lon)*cos(lat);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        x = cos(lon)*cos(lat+latinc);
        y = sin(lat+latinc);
        z = -sin(lon)*cos(lat+latinc);
        glNormal3f(x,y,z);
        glVertex3f(x*radius,y*radius,z*radius);
        lon += loninc;
      }
    glEnd();
    lat += latinc;
  }

  /* North-pole triangular fan */
  glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0,1,0);
    glVertex3f(0,radius,0);
    y = sin(lat);
    lon = 0;
    for (i=0; i<=nlon; i++) {
      x = cos(lon)*cos(lat);
      z = -sin(lon)*cos(lat);
      glNormal3f(x,y,z);
      glVertex3f(x*radius,y*radius,z*radius);
      lon += loninc;
    }
  glEnd();


  glEndList();
}


void makeBoundaryLines(GLuint listid){

    glNewList(listid, GL_COMPILE);
    /* Boundary lines? */

    // Region[0], Region[1], Region[2]
    // Region[0], Region[1], 0
    // Region[0], 0,         Region[2]
    // 0,         Region[1], Region[2]
    // Region[0], 0,         0
    // 0,         0,         Region[2]
    // 0,         Region[1], 0
    // 0,         0,         0

    

    glBegin(GL_LINES);
      glLineWidth(50);

      glColor3f(0, 0, 1);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 0, Region[2]);

      glColor3f(0, 1, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(Region[0], 0, 0);

      glColor3f(0, 0.5, 0.5);
      glVertex3f(0, 0, 0);
      glVertex3f(0, Region[1], 0);


      glColor3f(1, 0, 0);

      glVertex3f(Region[0], 0, 0);
      glVertex3f(Region[0], Region[1], 0);

      glVertex3f(Region[0], 0, 0);
      glVertex3f(Region[0], 0,  Region[2]);

      glVertex3f(0, Region[1], 0);
      glVertex3f(0, Region[1], Region[2]);

      glVertex3f(0, Region[1], 0);
      glVertex3f(Region[0], Region[1], 0);

      glVertex3f(0, 0, Region[2]);
      glVertex3f(0, Region[1], Region[2]);

      glVertex3f(0, 0, Region[2]);
      glVertex3f(Region[0], 0, Region[2]);

      glVertex3f(Region[0], Region[1], 0);
      glVertex3f(Region[0], Region[1], Region[2]);

      glVertex3f(Region[0], Region[1], Region[2]);
      glVertex3f(Region[0], 0, Region[2]);

      glVertex3f(Region[0], Region[1], Region[2]);
      glVertex3f(0, Region[1], Region[2]);

      

    glEnd();

    
    glEndList();


}

// void makeButtons(GLuint listid){

//   double modelMatrix[16];
//   double projMatrix[16];
//   GLint viewport[4];
//   GLdouble world_coord_v1[3];
//   GLdouble world_coord_v2[3];
//   GLdouble world_coord_v3[3];

//   glGetIntegerv( GL_VIEWPORT, viewport );
//   glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
//   glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );

//   int x = 50; int y = 50;
//   GLint realy = viewport[3] - (GLint) y;

//   gluUnProject( (GLdouble) x, (GLdouble) realy, 0.0, modelMatrix,
//         projMatrix, viewport, &world_coord_v1[0], &world_coord_v1[1], &world_coord_v1[2] );

//   gluUnProject( (GLdouble) (x - 25), (GLdouble) (realy + 25), 0.0, modelMatrix,
//         projMatrix, viewport, &world_coord_v2[0], &world_coord_v2[1], &world_coord_v2[2] );

//   gluUnProject( (GLdouble) (x + 25) , (GLdouble) (realy + 25), 0.0, modelMatrix,
//         projMatrix, viewport, &world_coord_v3[0], &world_coord_v3[1], &world_coord_v3[2] );


//   printf("world_coord_v1: %f,%f,%f \n", world_coord_v1[0], world_coord_v1[1], world_coord_v1[2]);
//   printf("world_coord_v2: %f,%f,%f \n", world_coord_v2[0], world_coord_v2[1], world_coord_v2[2]);
//   printf("world_coord_v3: %f,%f,%f \n", world_coord_v3[0], world_coord_v3[1], world_coord_v3[2]);

  
//   glNewList(listid, GL_COMPILE);

//   glBegin(GL_TRIANGLES);
//     glColor3f(0,0,1);
//     // glVertex3f( world_coord_v1[0], world_coord_v1[1], world_coord_v1[2] );
//     // glVertex3f( world_coord_v2[0], world_coord_v2[1], world_coord_v2[2] );
//     // glVertex3f( world_coord_v3[0], world_coord_v3[1], world_coord_v3[2] );
//     glVertex3f( 3, 10.7, 22 );
//     glVertex3f( 3.5, 10.7, 22 );
//     glVertex3f( 3.25, 10.45, 22 );

//   glEnd();


    

//   glEndList();

//   // printf("world_coord_v1: %f,%f,%f \n", world_coord_v1[0], world_coord_v1[1], world_coord_v1[2]);
//   // printf("world_coord_v2: %f,%f,%f \n", world_coord_v2[0], world_coord_v2[1], world_coord_v2[2]);
//   // printf("world_coord_v3: %f,%f,%f \n", world_coord_v3[0], world_coord_v3[1], world_coord_v3[2]);
// }



/**********************************************************************/
void makeAtoms() {
/***********************************************************************
  Makes display-list of all atoms in the current frame using spheres.
***********************************************************************/
  int i;
  float *colorVals;
  float rval,gval,bval;


  glNewList(atomsid, GL_COMPILE);
  rval = Ratom; gval = Gatom; bval = Batom;  /* RGB color of an atom */
  // rval = 0; gval = 0; bval = 1;
  // printf("Number of atoms: %d\n", nAtom);
  for (i=0; i < nAtom; i++) {
    glPushMatrix();
    glTranslatef(r[i][0],r[i][1],r[i][2]);
    // if (i == 0)
      // printf(" Curent atom position: %f %f %f\n", r[i][0],r[i][1],r[i][2]);

    // colorVals = mapVelocityToColor(r[i][0],r[i][1],r[i][2]);

    if (mapVelToColor) {
      colorVals = mapVelocityToColor(rv[i][0],rv[i][1],rv[i][2]);
      rval = colorVals[0]; gval = colorVals[1]; bval = colorVals[2];
    }

    glColor3f(rval,gval,bval);
    glCallList(sphereid);
    glPopMatrix();
  }


  // Draw boundary lines
  glCallList(boundsid);

  // glPopMatrix();
  // glPopMatrix();
  // glPopMatrix();
  // glPushMatrix();

  // glTranslatef(0.5,0.5,0.5);
  glCallList(buttonsid);
  // glPushMatrix();
  // glPushMatrix();
  // glPushMatrix();
  // glPopMatrix();



  glEndList();


  // Try to do triangle
  // double projMatrix[16];
  // glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
  // glPushMatrix();
  // glPushMatrix();
  // glTranslatef(0,0,0);
// glMatrixMode(GL_PROJECTION);
// gluOrtho2D(0,400,0,500);

  // glColor3ub(255,255,0);  // yellow
  // glBegin(GL_TRIANGLES);
  // glColor3f(0.5,0,0);
  // glVertex2f( -0.5, -0.5 );
  // glVertex2f( 0.5, -0.5 );
  // glVertex2f( 0, 0.5 );
  // glPopMatrix();
   

    // glVertex2f(300.0,210.0);
    // glVertex2f(340.0,215.0);
    // glVertex2f(320.0,250.0);
  // glEnd();

  //  glMatrixMode(GL_MODELVIEW);

  // glNewList(boundsid, GL_COMPILE);
  // glCallList(boundsid);
  // glEndList();



  // glBegin(GL_LINES);
  //     glVertex3f(0, 0, 0);
  //     glVertex3f(Region[0], 0, 0);
  // glEnd();

  // glLineWidth(10);
  // glBegin(GL_LINES);
  //   glColor3f(rval, gval, bval);
  //   glVertex3f(0, 0, 0);
  //   glVertex3f(5, 5, 5);
  // glEnd();

}

float * mapVelocityToColor(double x, double y, double z ) {

   static float colors[3];
   int i;

  //  /* set the seed */
  //  srand( (unsigned)time( NULL ) );
  
  colors[0] = (x / Region[0]);
  colors[1] = (y / Region[1]);
  colors[2] = (z / Region[2]);

  // printf( "Red: %f, Blue: %f, Green: %f \n", colors[0], colors[1], colors[2]);

  //  for ( i = 0; i < 3; ++i) {
  //     colors[i] = (x / Region[i]);
  //     printf( "colors[%d] = %f\n", i, colors[i]);
  //  }

   return colors;
}

/**********************************************************************/
void makeCurframeGeom() {
/***********************************************************************
  Reads the atoms information for the current time frame and makes the
  display-list of all the atoms' geometry.
***********************************************************************/
  makeAtoms();
}

/**********************************************************************/
void drawScene() {
/***********************************************************************
  Called by display() to draw the view of the current scene.
***********************************************************************/
  /* Define viewing transformation */
  gluLookAt(
    (GLdouble)eye[0],(GLdouble)eye[1],(GLdouble)eye[2],
    (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
    (GLdouble)up[0],(GLdouble)up[1],(GLdouble)up[2]);
  glCallList(atomsid);
}

/**********************************************************************/
void display() {
/***********************************************************************
  Callback for glutDisplayFunc().  It clears the frame and depth 
  buffers and draws the atoms in the current frame.
***********************************************************************/
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glLoadIdentity();
  drawScene();
  glutSwapBuffers();
}

/**********************************************************************/
void initView (float *min_ext, float *max_ext) {
/***********************************************************************
  Initializes global viewing, lighting, and projection values.
***********************************************************************/
  GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position1[] = {0.5, 0.5, 1.0, 0.0};
  float dif_ext[3],dis;
  int i;

  /* Define normal light */
  glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
  glLightfv(GL_LIGHT0,GL_POSITION,light_position1);

  /* Enable a single OpenGL light */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  /* Use depth buffering for hidden surface elimination */
  glEnable(GL_DEPTH_TEST);

  /* get diagonal and average distance of extent */
  for (i=0; i<3; i++){
    min_ext[i] = 0.0;
    max_ext[i] = Region[i];
  }
  for (i=0; i<3; i++) dif_ext[i] = max_ext[i] - min_ext[i];
  dis = 0.0;
  for (i=0; i<3; i++) dis += dif_ext[i]*dif_ext[i];
  dis = (float)sqrt((double)dis);

  /* set center in world space */
  for (i=0; i<3; i++) center[i] = min_ext[i] + dif_ext[i]/2.0;

  /* set initial eye & look at location in world space */
  eye[0] = center[0];
  eye[1] = center[1];
  eye[2] = center[2] + dis;

  camera_dist = dis;

  // printf("eye = %f,%f,%f\n dis = %f\n", eye[0], eye[1], eye[2], dis);
  // printf("center = %f,%f,%f\n", center[0], center[1], center[2]);

  up[0] = 0.0;
  up[1] = 1.0;
  up[2] = 0.0;

  /* set parameters for gluPerspective() */
  /* Near- & far clip-plane distances */
  near_clip = (GLdouble)( 0.5*(dis-0.5*dif_ext[2]) );
  far_clip  = (GLdouble)( 2.0*(dis+0.5*dif_ext[2]) );
  /* Field of view */
  fovy = (GLdouble)( 0.5*dif_ext[1]/(dis-0.5*dif_ext[2]) );
  fovy = (GLdouble)( 2*atan((double)fovy)/M_PI*180.0 );
  fovy = (GLdouble)(1.2*fovy);

  /* Enable the color material mode */
  glEnable(GL_COLOR_MATERIAL);
}

/**********************************************************************/
int main(int argc, char **argv) {
/**********************************************************************/

  glutInit(&argc, argv);

  /* Read atomic coordinates from an MD-configuration file */
    // InitParams();
    InitConf();
    // ComputeAccel();
    ComputeAcceleration();
    stepCount=1;


  /* Set up an window */
  /* Initialize display mode */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  /* Specify window size */
  glutInitWindowSize(winx, winy);
  /* Open window */
  glutCreateWindow("Colliding Particles");

  /* Initialize view */
  initView(min_ext, max_ext);

  /* Set a glut callback functions */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(animate);

  glutMotionFunc(moveAtom);
  glutMouseFunc(mouseClick);
  glutKeyboardFunc(MyKeyboardFunc);

  

  /* generate an OpenGL display list for single sphere */
  sphereid = glGenLists(1);
  makeFastNiceSphere(sphereid,atom_radius);

  /* Draw boundary lines */
  boundsid = glGenLists(1);
  makeBoundaryLines(boundsid);

  // /* Draw buttons */
  // buttonsid = glGenLists(1);
  // makeButtons(buttonsid);
  
  /* generate an OpenGL display list for the atoms' geometry */
  atomsid = glGenLists(1);


  /* make the geometry of the current frame's atoms */
  makeCurframeGeom();

  /* Start main display loop */
  glutMainLoop();
  
  return 0;
}

