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

![drag_particle](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/70e36e10-b472-47da-8da0-371c0dcead61)

Increase Gravity:

![gravity_down](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/9f99571b-fb04-4d85-ae7d-b3e58484eede)

Decrease Gravity:

![gravity_up](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/38cec86c-778a-4c75-938a-6bcf6491476a)

Damping Coefficient Greater Than One:

![damping_coeeff_gt1](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/51d31f7d-d32d-44e2-b194-de3a5ab852f2)



