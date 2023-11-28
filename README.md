# CSCI596_FinalProject
Colliding Particle System Simulator

Compiling
---------------
To compile mdv.c, type:
cc -o particle-collision particle-collision.c -L/System/Library/Frameworks -framework GLUT -framework OpenGL -lm

To run particle-collision.c, type:
./particle-collision

Assumptions
---------------
- All particles are equal size and equal mass

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



Rotate View:
![rotating](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/ff2ef9e5-25e0-4979-8125-919f99506631)
Click and Drag Particle:
![drag_particle](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/f36ca340-3294-4c05-ac6b-576a76e2ee94)
Increase Gravity:
![gravity_down](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/c6dd8a7b-ccfe-401f-b37a-ca2fd2ac151d)
Decrease Gravity:
![gravity_up](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/2753f6d9-15ea-4a01-9949-152ede65a6ef)
Damping Coefficient Greater Than One:
![damping_coeeff_gt1](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/51d31f7d-d32d-44e2-b194-de3a5ab852f2)



