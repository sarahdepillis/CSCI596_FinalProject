/***********************************************************************
  Program particle-collision.c -- ball collision simulation.
  Required files
    particle-collision.h:   Include file
***********************************************************************/
#include "particle-collision.h"
#include <stdio.h>
#include <math.h>
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>    /* Header file for the OpenGL library */
#include <OpenGL/glu.h>   /* Header file for the GLu library */
#include <GLUT/glut.h>    /* Header file for the GLut library */

#include <time.h>
#include <stdlib.h>
#include <sys/time.h>     /* for gettimeofday() */ 

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

GLuint sphereid;          /* display-list id of atom sphere geom */
GLuint atomsid;           /* display-list id of all atoms */
GLuint boundsid;          /* display-list id for boundary lines */
GLdouble fovy, aspect, near_clip, far_clip;  /* parameters for gluPerspective() */

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
    dist = dist + pow(p1[i]-p2[i],2.0);
  }
  dist = sqrt(dist);
  return dist;
}
double magnitude(double *v, int n) {
  return sqrt(pow(v[0],2.0)+pow(v[1],2.0)+pow(v[2],2.0));
}


void InitConf();
void ComputeAcceleration();
void SingleStep();
void HalfKick();
void ApplyBoundaryCond();
void animate(void);
float * mapVelocityToColor(double x, double y, double z);

/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
	r are initialized to random positions withing region.  
	rv are initialized with a random velocities.  
------------------------------------------------------------------------------*/

	int j,n,k;
	double seed, e[3];

  /* Generates random positions */
  for (int i = 0; i<nAtom;i++){
    for (k=0; k<3; k++) {
      r[i][k] = ((double)rand()/(double)(RAND_MAX)) * Region[k];
    }
  }

	/* Generates random velocities */
	seed = 13597.0;
	for(n=0; n<nAtom; n++) {
		RandVec3(e,&seed);
		for (k=0; k<3; k++) {
			rv[n][k] = 20*e[k];
		}
	}
}

void ComputeAcceleration() {
  int n,k;

  for (n=0; n<nAtom; n++) {

    if (n != clickedAtom) {
      ra[n][0] = 0;
      ra[n][1] = gravity;
      ra[n][2] = 0;
    }
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
  ComputeAcceleration();
	HalfKick(); /* Second half kick to obtain v(t+Dt) */
}

/*----------------------------------------------------------------------------*/
void HalfKick() {
/*------------------------------------------------------------------------------
	Accelerates velocities, rv, by half the time step.
------------------------------------------------------------------------------*/

  DeltaTH = 0.5*DeltaT;

  int n,k;
	for (n=0; n<nAtom; n++) {
		for (k=0; k<3; k++) {
       rv[n][k] = rv[n][k] + DeltaTH*ra[n][k];
    }
  }

  
}

/*----------------------------------------------------------------------------*/
void ApplyBoundaryCond() {
/*------------------------------------------------------------------------------
	Applies periodic boundary conditions to atomic coordinates.
------------------------------------------------------------------------------*/
  
  int n,k;

	for (n=0; n<nAtom; n++) {

    // Check wall collision
		for (k=0; k<3; k++)  {

      if (r[n][k] + atom_radius > Region[k] || r[n][k] - atom_radius < 0) {
        // Check if outside bounds
        if (r[n][k] - atom_radius < 0)
          r[n][k] = 0.0 + atom_radius;
        else
          r[n][k] = Region[k] - atom_radius;

        rv[n][k] = rv[n][k] * -1.0 * collisionDamping;
      }

    }

    // Check particle collision
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


        double norm[3];
        double relative_v[3];
        long double vel_norm_dot = 0;
        long double dot_scalar_norm[3];
        long double norm_sum = 0;

        // Check wall collision for second particle first
        for (k=0; k<3; k++)  {

          if (r[n2][k] + atom_radius > Region[k] || r[n2][k] - atom_radius < 0) {

            if (r[n2][k] - atom_radius < 0)
              r[n2][k] = 0.0 + atom_radius;
            else
              r[n2][k] = Region[k] - atom_radius;

            rv[n2][k] = rv[n2][k] * -1.0 * collisionDamping;
          }

        }

        long double dst = distance(r[n], r[n2], 3);
        
        if (dst <= 2 * atom_radius) {
          // Get norm vector between two particles


          for (k=0; k<3; k++)  {
            norm[k] = r[n][k] - r[n2][k];
            norm_sum += pow(norm[k],2);
            norm[k] = norm[k] / sqrt(norm_sum);

            if (n == clickedAtom)
              relative_v[k] = clickedAtomNewVelocity[k] - rv[n2][k];
            else if (n2 == clickedAtom)
              relative_v[k] = rv[n][k] - clickedAtomNewVelocity[k];
            else
              relative_v[k] = rv[n][k] - rv[n2][k];
          }

          // Shift particles so they arent overlapping
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
          // End shift particles


          vel_norm_dot = dot_product(relative_v, norm, 3);

          for (k=0; k<3; k++)  {
            dot_scalar_norm[k] = vel_norm_dot*norm[k];
          }

          if (n == clickedAtom) {
            for (k=0; k<3; k++)  {
              rv[n2][k] = (rv[n2][k] + dot_scalar_norm[k]) * collisionDamping;
            }
          }
          else if (n2 == clickedAtom) {
            for (k=0; k<3; k++)  {
              rv[n][k] = (rv[n][k] - dot_scalar_norm[k]) * collisionDamping;
            }
          }
                  

          else {
            for (k=0; k<3; k++)  {
              rv[n][k] = (rv[n][k] - dot_scalar_norm[k]) * collisionDamping;
              rv[n2][k] = (rv[n2][k] + dot_scalar_norm[k]) * collisionDamping;
            }
          }

          // If velocity gets low enough, set to zero
          if (magnitude(rv[n], 3) < 0.5) {
            rv[n][0] = 0;
            rv[n][1] = 0;
            rv[n][2] = 0;
          }

          if (magnitude(rv[n2], 3) < 0.5) {
            rv[n2][0] = 0;
            rv[n2][1] = 0;
            rv[n2][2] = 0;
          }

        }
      
      }

  }
}

void animate() { /* Callback function for idle events */
    /* Keep updating the scene until the last step is reached */
    if (stepCount <= StepLimit) {
        SingleStep(); /* One step integration */
        makeCurframeGeom(); /* Redraw the scene (make a display list) */
        glutPostRedisplay(); 
        ++stepCount;
    }
}

void mouseClick(int button, int state, int x, int y){
  if (state && clickedAtom > -1) { // Button up
      // Update atom velocity and reset gravity on click release
      for (int k =0; k<3;k++){
        rv[clickedAtom][k] = clickedAtomNewVelocity[k];
      }
      ra[clickedAtom][1] = gravity;

      // Reset clickedAtom
      clickedAtom = -1;
  }
  else if (state == 0) {   // Button down
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

    double mouseVec[3];
    for (int k=0; k <3;k++){
      mouseVec[k] = mouseVec_end[k] - mouseVec_start[k];
    }  

    double mouseVec_square = dot_product(mouseVec, mouseVec, 3); 
    double mouseStart_toAtom[3];
    double mouseVect_mouseAtomVec_dot[3];
    double closestPoint[3];
    int mouseRayAtomCollision = 0;
    double epsilon = 0.01;

    for (int n=0; n<nAtom; n++) {
      for (int k=0; k <3;k++){
        mouseStart_toAtom[k] = r[n][k] - mouseVec_start[k];
      }

      double mouseVect_mouseAtomVec_dot = dot_product(mouseStart_toAtom, mouseVec, 3); 
      double mousePoint_mouseVect_proj = mouseVect_mouseAtomVec_dot / mouseVec_square;

      for (int k=0; k <3;k++){
        closestPoint[k] = mouseVec_start[k] + mouseVec[k] * mousePoint_mouseVect_proj;
      }

      double len = distance(closestPoint, r[n], 3); 
      
      if (len < atom_radius + epsilon) {  // Check for collision between mouse ray and atom
        mouseRayAtomCollision = 1;
        clickedAtom = n;
        break;
      }

    }

    // If an atom has been clicked, set velocities and acceleration to 0
    if (mouseRayAtomCollision) {
      for(int k=0; k<3;k++) {
        rv[clickedAtom][k] = 0;
        ra[clickedAtom][k] = 0;
      }
    }
  }
}


void moveAtom(int x, int y) {

  double atomStartPos[3];


  if (clickedAtom > -1) {

    // TIMER setup //
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

    double mouseVec[3];
    double newLoc[3];

    for (int k=0; k <3;k++){
      mouseVec[k] = mouseVec_end[k] - mouseVec_start[k];
    }  

    GLdouble obj_win_coord[3];
    gluProject(r[clickedAtom][0], r[clickedAtom][1], r[clickedAtom][2], modelMatrix,
        projMatrix, viewport, &obj_win_coord[0], &obj_win_coord[1], &obj_win_coord[2]);

    gluUnProject( (GLdouble) x, (GLdouble) realy, obj_win_coord[2], modelMatrix,
        projMatrix, viewport, &newLoc[0], &newLoc[1], &newLoc[2] );

    
    for(int k =0; k <3;k++){
      atomStartPos[k] = r[clickedAtom][k];
      r[clickedAtom][k]= newLoc[k];
    }

    // stop timer
    gettimeofday(&t2, NULL);

    // compute the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

    // Get velocity
    double vel[3];
    for(int k =0; k <3;k++){
        clickedAtomNewVelocity[k] = (r[clickedAtom][k] - atomStartPos[k]) / (elapsedTime);

    }

  }

}

void MyKeyboardFunc(unsigned char Key, int x, int y)
{
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
      if (nAtom > 0)
        nAtom--; 
      makeAtoms();
      break;
    case 'q': // Quit program
      exit(1);
      break;
    case 'l': // Pan camera view clockwise
      angle_delta = 0.05;
      rotateCamera(angle_delta);
    break;
    case 'k': // Pan camera view counter-clockwise
      angle_delta = -0.05;
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

  drawScene();
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
    /* Boundary lines */

    glBegin(GL_LINES);
      // glLineWidth(500);

      glColor3f(0, 0, 1);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 0, Region[2]);

      // glLineWidth(50);
      glColor3f(0, 1, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(Region[0], 0, 0);

      // glLineWidth(50);
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
  for (i=0; i < nAtom; i++) {
    glPushMatrix();
    glTranslatef(r[i][0],r[i][1],r[i][2]);

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

  glEndList();

}

float * mapVelocityToColor(double x, double y, double z ) {

  static float colors[3];
  int i;
  
  colors[0] = (x / Region[0]);
  colors[1] = (y / Region[1]);
  colors[2] = (z / Region[2]);

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
  
  /* generate an OpenGL display list for the atoms' geometry */
  atomsid = glGenLists(1);


  /* make the geometry of the current frame's atoms */
  makeCurframeGeom();

  /* Start main display loop */
  glutMainLoop();
  
  return 0;
}

