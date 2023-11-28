# CSCI596_FinalProject
Colliding Particle System Simulator

Compiling
---------------
To compile mdv.c, type:
cc -o mdv mdv.c -L/System/Library/Frameworks -framework GLUT -framework OpenGL -lm

To run mdv.c, type:
./md

Usage
---------------
Interactive options:
- You can grab and throw a particle with the mouse
- 'g' increases upward gravity
- 'h' increases downward gravity
- 'z' increases size of particles
- 'x' decreases size of particles
- 'c' add particle to the system
- 'v' remove particle from system
- 'q' quit program
- 'l' pan camera view clockwise
- 'k' pan camera view counter-clockwise
- 'd' decrease damping coefficient
- 'f' increase damping coefficient
- 'm' toggle on/off mapping the particle velocity to 3D color cube

Assumptions
---------------
- All particles are equal size and equal mass



![rotating](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/42521c94-61dd-457f-82c8-56d5b8fe2b6f)


