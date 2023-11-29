# CSCI596_FinalProject
Interactive Colliding Particle System Simulator

Compiling
---------------
To compile mdv.c, type:
cc -o particle-collision particle-collision.c -L/System/Library/Frameworks -framework GLUT -framework OpenGL -lm

To run particle-collision.c, type:
./particle-collision

Assumptions
---------------
- All particles are equal size and equal mass


Collisions
---------------

### Wall Collisions
A particle that lives in a region , $`R = \begin{bmatrix} R_x \\ R_y \\ R_z \end{bmatrix}, `$ has a position, $`s = \begin{bmatrix} x \\ y \\ z\end{bmatrix}`$

When 
```
$`0 \leq x \leq R_x `$ OR
$`0 \leq y \leq R_y `$ OR
$`0 \leq z \leq R_z `$
```
a collsion has occured.


![image](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/3bcc8f04-973a-4c7d-baca-5d40c8f4329a)


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

![drag_particle](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/54a1fcf6-5979-4c70-8408-810cf626ed9e)

Increase/Decrease Gravity:

![gravity_down](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/0b3c781c-9d6f-422a-af30-94ab87c3f31b)
![gravity_up](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/dee8c3f7-c64e-4e7b-8bcc-2ad39d542282)

Add Particles:

![add_particles](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/f79b7a16-97e9-4da6-bd13-b67b2645a11f)

Remove/Add Particles: 

![rm_particles](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/7d19c1e5-b239-4e4f-a7e1-2c0645d11180)
![rm_particles (1)](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/d5c40519-644e-41a2-b33d-be8a8453a292)

Size Change:

![size_change](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/7611ea71-b0e8-4f5e-b214-85f1d4b5e5f3)
![size_change (1)](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/e27302ca-7575-423c-a3b0-c88572b73d95)
![size_change (2)](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/53698618-e518-443c-8909-f6ba40ef018c)



Damping Coefficient Greater Than One:

![damping_coeeff_gt1](https://github.com/sarahdepillis/CSCI596_FinalProject/assets/28903687/51d31f7d-d32d-44e2-b194-de3a5ab852f2)



