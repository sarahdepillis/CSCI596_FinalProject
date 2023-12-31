# CSCI596_FinalProject
Interactive Colliding Particle System Simulator

Compiling
---------------
To compile particle-collision.c, type:
cc -o particle-collision particle-collision.c -L/System/Library/Frameworks -framework GLUT -framework OpenGL -lm

To run particle-collision.c, type:
./particle-collision

Assumptions
---------------
- All particles are equal in size and mass


Collisions
---------------

### Wall Collisions
A particle that lives in a region, $`R = \begin{bmatrix} R_x \\ R_y \\ R_z \end{bmatrix}, `$ has a position, $`s = \begin{bmatrix} x \\ y \\ z\end{bmatrix}`$ and a velocity of $`v = \begin{bmatrix} v_x \\ v_y \\ v_z\end{bmatrix}`$

When  
$`0 \geq x \geq R_x`$ or  
$`0 \geq y \geq R_y`$ or  
$`0 \geq z \geq R_z`$  
a wall collision has occurred. 

If $`0 \geq x \geq R_x`$ then $`v_x = -1 * D * v_x`$  
If $`0 \geq y \geq R_y`$ then $`v_y = -1 * D * v_y`$  
If $`0 \geq z \geq R_z`$ then $`v_z = -1 * D * v_z`$  
Where D is the damping coefficient if one exists


### Particle Collisions
Particle A has a position, $`s_A = \begin{bmatrix} x_1 \\ y_1 \\ z_1\end{bmatrix}`$ and velocity $`\vec{v_A} = \begin{bmatrix} v_x1 \\ v_y1 \\ v_z1\end{bmatrix}`$  
Particle B has a position, $`s_B = \begin{bmatrix} x_2 \\ y_2 \\ z_2\end{bmatrix}`$ and velocity $`\vec{v_B} = \begin{bmatrix} v_x2 \\ v_y2 \\ v_z2\end{bmatrix}`$  
Both particles have radius r.

Find the distance between the two particles: $$d = \sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2}$$

If $`d < 2*r`$, a collision between two particles has occurred.

Now new velocity vectors need to be calculated.

1. Find the normal vector between the two particles: $`\vec{N} = s_A - s_B `$
2. Get the magnitude of the normal vector by taking the square root of the dot product of the normal vector with itself: $`||N|| = \sqrt {\vec{N} \cdot \vec{N}}`$
3. Divide the normal vector by its magnitude to get the unit norm vector: $`\vec{n} = \frac{\vec{N}}{||N||}`$
4. Get the relative velocity between the two particles: $`\vec{V} = v_A - v_B`$
5. Find the dot product of the unit norm vector and the relative velocity vector: $`vn = \vec{V} \cdot  \vec{n}`$
6. Multiply the resulting dot product with the unit norm vector: $`\vec{VN} = vn * \vec{n}`$
7. Find the resulting velocities for each particle:
   - $`v_A = (\vec{v_A} - \vec{VN}) * collisionDamping `$
   - $`v_B = (\vec{v_B} + \vec{VN}) * collisionDamping `$


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



