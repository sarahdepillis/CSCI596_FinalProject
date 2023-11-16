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

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

GLuint sphereid;          /* display-list id of atom sphere geom */
GLuint atomsid;           /* display-list id of all atoms */
GLuint boundsid;          /* display-list id for boundary lines */
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

void InitParams();
void InitConf();
void ComputeAccel();
void ComputeAcceleration();
void SingleStep();
void HalfKick();
void ApplyBoundaryCond();
void EvalProps();

void animate(void);
// void changeView(void);

float * mapVelocityToColor();

/*----------------------------------------------------------------------------*/
void InitParams() {
/*------------------------------------------------------------------------------
	Initializes parameters.
------------------------------------------------------------------------------*/
	int k;
	double rr,ri2,ri6,r1;

	/* Reads control parameters */
	scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
	scanf("%le",&Density);
	scanf("%le",&InitTemp);
	scanf("%le",&DeltaT);
	scanf("%d",&StepLimit);
	scanf("%d",&StepAvg);

	/* Computes basic parameters */
	DeltaTH = 0.5*DeltaT;
	for (k=0; k<3; k++) {
		Region[k] = InitUcell[k]/pow(Density/4.0,1.0/3.0);
		RegionH[k] = 0.5*Region[k];

    // printf("REGION: %f\n", Region[k]);
	}



	/* Constants for potential truncation */
	rr = RCUT*RCUT; ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1=sqrt(rr);
	Uc = 4.0*ri6*(ri6 - 1.0);
	Duc = -48.0*ri6*(ri6 - 0.5)/r1;
}

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to face-centered cubic (fcc) lattice positions.  
	rv are initialized with a random velocity corresponding to Temperature.  
------------------------------------------------------------------------------*/
	double c[3],gap[3],e[3],vSum[3],vMag;

  double rLoc; 

	int j,n,k,nX,nY,nZ;
	double seed;
	/* FCC atoms in the original unit cell */
	double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
	                         {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

	/* Sets up a face-centered cubic (fcc) lattice */
	for (k=0; k<3; k++) gap[k] = Region[k]/InitUcell[k];

  // printf("The gap is: %f, %f, %f \n", gap[0], gap[1], gap[2]);

  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);
  // printf("Printing random number: %f \n",((double)rand()/(double)(RAND_MAX)) * Region[k]);


  // for (k=0; k<3; k++) {
  //   rLoc[k] = ((double)rand()/(double)(RAND_MAX)) * Region[k];
  //   printf("rLoc[%f] = %f", k, rLoc[k]);
  // }


  int nAtomTotal = 300;
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
	vMag = sqrt(3*InitTemp);
	for(k=0; k<3; k++) vSum[k] = 0.0;
	for(n=0; n<nAtom; n++) {
		RandVec3(e,&seed);
		for (k=0; k<3; k++) {
			rv[n][k] = vMag*e[k];
      // rv[n][k] = ((float)rand()/(float)(RAND_MAX)) * 10.0;
      // printf("rv[n][k] = %f\n", rv[n][k]);
			vSum[k] = vSum[k] + rv[n][k];
		}
	}
	/* Makes the total momentum zero */
	for (k=0; k<3; k++) vSum[k] = vSum[k]/nAtom;
	for (n=0; n<nAtom; n++) for(k=0; k<3; k++) rv[n][k] = rv[n][k] - vSum[k];
}

void ComputeAcceleration() {
  int n,k;
  for (n=0; n<nAtom; n++) {



    ra[n][0] = 0;
    ra[n][1] = -19.81;
    // ra[n][1] = 0;
    ra[n][2] = 0;

    // printf("In acceleration, position is: %f, %f, %f\n", r[n][0], r[n][1], r[n][2]);
  }

}

/*----------------------------------------------------------------------------*/
void ComputeAccel() {
/*------------------------------------------------------------------------------
	Acceleration, ra, are computed as a function of atomic coordinates, r,
	using the Lennard-Jones potential.  The sum of atomic potential energies,
	potEnergy, is also computed.   
------------------------------------------------------------------------------*/
	double dr[3],f,fcVal,rrCut,rr,ri2,ri6,r1;
	int j1,j2,n,k;

	rrCut = RCUT*RCUT;
	for (n=0; n<nAtom; n++) for (k=0; k<3; k++) ra[n][k] = 0.0;
	potEnergy = 0.0;

	/* Doubly-nested loop over atomic pairs */
	for (j1=0; j1<nAtom-1; j1++) {
		for (j2=j1+1; j2<nAtom; j2++) {
			/* Computes the squared atomic distance */
			for (rr=0.0, k=0; k<3; k++) {
				dr[k] = r[j1][k] - r[j2][k];
				/* Chooses the nearest image */
				dr[k] = dr[k] - SignR(RegionH[k],dr[k]-RegionH[k])
				              - SignR(RegionH[k],dr[k]+RegionH[k]);
				rr = rr + dr[k]*dr[k];
			}
			/* Computes acceleration & potential within the cut-off distance */
			if (rr < rrCut) {
				ri2 = 1.0/rr; 
        ri6 = ri2*ri2*ri2; 
        r1 = sqrt(rr);
				fcVal = 48.0*ri2*ri6*(ri6-0.5) + Duc/r1;
				for (k=0; k<3; k++) {
					f = fcVal*dr[k];
					ra[j1][k] = ra[j1][k] + f;
					ra[j2][k] = ra[j2][k] - f;
				}
				potEnergy = potEnergy + 4.0*ri6*(ri6-1.0) - Uc - Duc*(r1-RCUT);
			} 
		} 
    // ra[j1][1]+= -9.81;  // Add gravity
	}
}

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
  // ComputeAcceleration();
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

  int n,k;
	for (n=0; n<nAtom; n++) {
		for (k=0; k<3; k++) {
       rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];
      //  printf("DeltaTH: %f, accel: %f, velocity:%f\n", DeltaTH, ra[n][k], rv[n][k]);
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
  float collisionDamping = 0.8;
  // float collisionDamping = 0.95;

	for (n=0; n<nAtom; n++) {
    // if (r[n][0] > 5.129928 || r[n][1] > 5.129928 || r[n][2] > 5.129928 ||
    //   r[n][0] < 0 || r[n][1] < 0 || r[n][2] < 0){
    //   printf("Position Before Boundary condition: %f, %f, %f\n", r[n][0], r[n][1], r[n][2]);
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


          double velocity_error[3];
          for (k=0; k<3; k++)  {
            long double rvn_old = rv[n][k];
            long double rvn2_old = rv[n2][k];

            rv[n][k] = (rv[n][k] - dot_scalar_norm[k]) * collisionDamping;
            rv[n2][k] = (rv[n2][k] + dot_scalar_norm[k]) * collisionDamping;

            // ra[n][k] = (rv[n][k] - rvn_old)/DeltaT;
            // ra[n][1] = ra[n][1] - 9.81;
            // ra[n2][k] = (rv[n2][k] - rvn2_old)/DeltaT;
            // ra[n][1] = ra[n2][1] - 9.81;


            velocity_error[k] = (rvn_old + rvn2_old) - (rv[n][k] + rv[n2][k]);
            
          }

          // printf("Particle %d has NEW velocity %f,%f,%f  \n", n, rv[n][0], rv[n][1], rv[n][2]);
          
          // rv[n][0] += velocity_error[0];
          // rv[n][1] += velocity_error[1];
          // rv[n][2] += velocity_error[2];

          // printf("velocity_error = %f,%f,%f\n", velocity_error[0], velocity_error[1], velocity_error[2]);

           

          // printf("rv[n][k] = %f,%f,%f  \n", rv[n][0], rv[n][1], rv[n][2]);
          // printf("rv[n2][k] = %f,%f,%f  \n", rv[n2][0], rv[n2][1], rv[n2][2]);

          // double rv_sum = 0;
          // for (k=0; k<3; k++)  {
          //   norm[k] = r[n][k] - r[n2][k];
          //   rv_sum += pow(norm[k],2);
          //   relative_v[k] = rv[n][k] - rv[n2][k];
          // }
          // for (k=0; k<3; k++)  {
          //   norm[k] = norm[k] / sqrt(norm_sum);
          // }

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

          // double rv_mag2 = sqrt(pow(rv[n2][0],2.0)+pow(rv[n2][1],2.0)+pow(rv[n2][2],2.0));
          // printf("rv_mag2 = %f\n", rv_mag2);


          // Add some logic here that checks if the new position will shove the ball outside the bounds

          // double new_rv[3];
          // double new_rv2[3];
          // double newPos1[3];
          // double newPos2[3];

          // new_rv[0] = (2*atom_radius - dst)*(rv[n][0]/rv_mag);
          // new_rv[1] = (2*atom_radius - dst)*(rv[n][1]/rv_mag);
          // new_rv[2] = (2*atom_radius - dst)*(rv[n][2]/rv_mag);

          // newPos1[0] = r[n][0] + new_rv[0];
          // newPos1[1] = r[n][1] + new_rv[1];
          // newPos1[2] = r[n][2] + new_rv[2];


          // new_rv2[0] = (2*atom_radius - dst)*(rv[n2][0]/rv_mag2);
          // new_rv2[1] = (2*atom_radius - dst)*(rv[n2][1]/rv_mag2);
          // new_rv2[2] = (2*atom_radius - dst)*(rv[n2][2]/rv_mag2);

          // newPos2[0] = r[n2][0] + new_rv2[0];
          // newPos2[1] = r[n2][1] + new_rv2[1];
          // newPos2[2] = r[n2][2] + new_rv2[2];

          // printf("REGION: %f,%f,%f\n", Region[0], Region[1], Region[2]);
          // printf("New position 1 would be = %f,%f,%f\n", newPos1[0], newPos1[1], newPos1[2]);
          // printf("New position 2 would be = %f,%f,%f\n", newPos2[0], newPos2[1], newPos2[2]);

          // for (k=0; k<3; k++)  {
          //   if (newPos1[k] > Region[k] || newPos1[k] < 0) {
          //     new_rv2[k] = (2*atom_radius - dst)*(rv[n2][k]/rv_mag2);

          //     r[n2][k] += new_rv2[k];

          //   }
          //   else {
          //     r[n][k] += new_rv[k];
          //   }
          // }

          // r[n][0] += (2*atom_radius - dst)*(rv[n][0]/rv_mag);
          // r[n][1] += (2*atom_radius - dst)*(rv[n][1]/rv_mag);
          // r[n][2] += (2*atom_radius - dst)*(rv[n][2]/rv_mag);

          // double new_rv_mag = sqrt(pow(new_rv[0],2.0)+pow(new_rv[1],2.0)+pow(new_rv[2],2.0));
          // printf("new_rv_mag = %f\n", new_rv_mag);

          // new_rv2[0] = (2*atom_radius - dst/2.0)*(rv[n2][0]/rv_mag2);
          // new_rv2[1] = (2*atom_radius - dst/2.0)*(rv[n2][1]/rv_mag2);
          // new_rv2[2] = (2*atom_radius - dst/2.0)*(rv[n2][2]/rv_mag2);

          // r[n2][0] += new_rv2[0];
          // r[n2][1] += new_rv2[1];
          // r[n2][2] += new_rv2[2];

          // r[n][0] += (2*atom_radius - dst)*(rv[n][0]/rv_mag);
          // r[n][1] += (2*atom_radius - dst)*(rv[n][1]/rv_mag);
          // r[n][2] += (2*atom_radius - dst)*(rv[n][2]/rv_mag);

          // double new_rv_mag2 = sqrt(pow(new_rv2[0],2.0)+pow(new_rv2[1],2.0)+pow(new_rv2[2],2.0));
          // printf("new_rv_mag2 = %f\n", new_rv_mag2);


          // printf("Particle %d has NEW position %f,%f,%f and NEW velocity %f,%f,%f  \n", n, r[n][0], r[n][1], r[n][2], rv[n][0], rv[n][1], rv[n][2]);
          // printf("Particle %d has NEW position %f,%f,%f and NEW velocity %f,%f,%f  \n\n", n2, r[n2][0], r[n2][1], r[n2][2], rv[n2][0], rv[n2][1], rv[n2][2]);




          // long double newdst = sqrt(pow(r[n][0]-r[n2][0],2.0)+pow(r[n][1]-r[n2][1],2.0)+pow(r[n][2]-r[n2][2],2.0));
          // printf("Particles %d and %d new distance %Lf \n\n", n, n2, newdst);
          

        }
      
      }

      

    // }





    // if (toPrint == 1){
    //   printf("Position After Boundary condition:   %f, %f, %f\n", r[n][0], r[n][1], r[n][2]);
    //   printf("Velocity After Boundary condition:   %f, %f, %f\n", rv[n][0], rv[n][1], rv[n][2]);
      
    //   toPrint = 0;
    // }
  }

  // double totSystemSpeed[3];
  // // for (n=0; n<nAtom; n++) {
  //   for (k=0; k<3; k++) {
  //     totSystemSpeed[k] = 0;
  //   }
  // // }

  // for (n=0; n<nAtom; n++) {
  //   for (k=0; k<3; k++) {
  //     totSystemSpeed[k] += fabs(rv[n][k]);

  //     if (fabs(rv_init[n][k]) != fabs(rv[n][k]))
  //       printf("Particle %d speed changed: \n rv_init[%d][%d] = %f\n rv[%d][%d] = %f\n", n, n, k, rv_init[n][k],n,k,rv[n][k]);
  //   }
  // }

  // printf("\nTotal system speed start: %f,%f,%f = %f \n\n", totSystemSpeedStart[0], totSystemSpeedStart[1], totSystemSpeedStart[2], (totSystemSpeedStart[0] + totSystemSpeedStart[1] + totSystemSpeedStart[2]));
  // printf("\nTotal system speed end: %f,%f,%f = %f \n\n", totSystemSpeed[0], totSystemSpeed[1], totSystemSpeed[2], (totSystemSpeed[0] + totSystemSpeed[1] + totSystemSpeed[2]));

}

/*----------------------------------------------------------------------------*/
void EvalProps() {
/*------------------------------------------------------------------------------
	Evaluates physical properties: kinetic, potential & total energies.
------------------------------------------------------------------------------*/
	double vv;
	int n,k;

	kinEnergy = 0.0;
	for (n=0; n<nAtom; n++) {
		vv = 0.0;
		for (k=0; k<3; k++)
			vv = vv + rv[n][k]*rv[n][k];
		kinEnergy = kinEnergy + vv;
	}
	kinEnergy *= (0.5/nAtom);
	potEnergy /= nAtom;
	totEnergy = kinEnergy + potEnergy;
	temperature = kinEnergy*2.0/3.0;

	/* Print the computed properties */
	// printf("%9.6f %9.6f %9.6f %9.6f\n",
	// stepCount*DeltaT,temperature,potEnergy,totEnergy);
}

void animate() { /* Callback function for idle events */
    /* Keep updating the scene until the last MD step is reached */
    if (stepCount <= StepLimit) {
        SingleStep(); /* One MD-step integration */
        if (stepCount%StepAvg == 0) EvalProps(); 
        makeCurframeGeom(); /* Redraw the scene (make a display list) */
        glutPostRedisplay(); 
        ++stepCount;
    }
}

void changeView(int x, int y) {
  printf("changeView function being triggered: x = %d, y = %d \n", x, y);
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
      glLineWidth(20);

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
    // colorVals = mapVelocityToColor(rv[i][0],rv[i][1],rv[i][2]);
    // rval = colorVals[0]; gval = colorVals[1]; bval = colorVals[2];

    glColor3f(rval,gval,bval);
    glCallList(sphereid);
    glPopMatrix();
  }


  // Draw boundary lines
  glCallList(boundsid);


  glEndList();



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

  printf("eye = %f,%f,%f\n dis = %f\n", eye[0], eye[1], eye[2], dis);

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
    InitParams();
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
  glutCreateWindow("Lennard-Jones Atoms");

  /* Initialize view */
  initView(min_ext, max_ext);

  /* Set a glut callback functions */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(animate);
  glutMotionFunc(changeView);

  

  /* generate an OpenGL display list for single sphere */
  sphereid = glGenLists(1);
  makeFastNiceSphere(sphereid,atom_radius);

  boundsid = glGenLists(1);
  makeBoundaryLines(boundsid);
  
  /* generate an OpenGL display list for the atoms' geometry */
  atomsid = glGenLists(1);


  /* make the geometry of the current frame's atoms */
  makeCurframeGeom();

  /* Start main display loop */
  glutMainLoop();
  
  return 0;
}

