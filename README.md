# Modeling and Simulation of Robotic 3D Printer

## Introduction

Robotic 3D Printers allow engineers and designers to rapidly prototype new mechanisms, products, or art. There are two related goals of this project. The first is to simulate the kinematics, and dynamics of a robotic 3D printer, and the filament particles it extrudes. The second is to simplify the simulation by removing the electric and drag forces acting on the particles, and optimize this simplified simulation using a Genetic Algorithm. The Genetic Algorithm was tasked with minimizing the difference between a desired path, and the path generated by the simulation. To do this, design strings of 4 parameters, 3 angular velocities, and the initial velocity of extrusion, were randomly generated, tested, ranked, and new strings were bred from the top performing design strings.

## Background and Theory

The motion of the dispenser at the end of the robotic arm is described by:

![eq1](https://user-images.githubusercontent.com/52175303/109432441-1355bd80-79c0-11eb-9204-49f7a8c4483b.png)

![eq2](https://user-images.githubusercontent.com/52175303/109432448-2072ac80-79c0-11eb-90d2-511974e3877e.png)

![eq3](https://user-images.githubusercontent.com/52175303/109432471-36806d00-79c0-11eb-8755-19c2db1d2bff.png)

![eq4](https://user-images.githubusercontent.com/52175303/109432485-46984c80-79c0-11eb-9331-fc1abbc5e72c.png)

![eq5](https://user-images.githubusercontent.com/52175303/109432506-56179580-79c0-11eb-8e97-69e024f44bd9.png)

![eq6](https://user-images.githubusercontent.com/52175303/109432513-6760a200-79c0-11eb-9d31-bf0c335ba8a0.png)

![eq7](https://user-images.githubusercontent.com/52175303/109432524-78111800-79c0-11eb-87ee-b54bef8f5702.png)

![eq8](https://user-images.githubusercontent.com/52175303/109432538-86f7ca80-79c0-11eb-81a6-b158b8b46e43.png)

When a droplet is extruded from the dispenser, it exits with an initial position and velocity:

![eq9](https://user-images.githubusercontent.com/52175303/109432559-9ecf4e80-79c0-11eb-9d27-86d927c72f43.png)

Where ![eq10](https://user-images.githubusercontent.com/52175303/109432579-b575a580-79c0-11eb-981a-2b3bc53716ca.png) is the extrusion velocity, one of the design parameters.

The dynamics of each droplet is described by Newton's second law:

![eq11](https://user-images.githubusercontent.com/52175303/109432615-df2ecc80-79c0-11eb-872a-dfafbd7ab483.png)

The force of gravity acting on each particle is:

![eq12](https://user-images.githubusercontent.com/52175303/109432645-fa014100-79c0-11eb-901d-cd134a459275.png)

The force of the electric field is:

![eq13](https://user-images.githubusercontent.com/52175303/109432662-0f766b00-79c1-11eb-9b95-a0b50ed77ff4.png)

Where ![eq14](https://user-images.githubusercontent.com/52175303/109432673-274def00-79c1-11eb-9470-5149817e4ca2.png) is the number of point charges, ![eq15](https://user-images.githubusercontent.com/52175303/109432681-3896fb80-79c1-11eb-89b3-58326893f337.png) is their charge, and ![eq16](https://user-images.githubusercontent.com/52175303/109432693-4e0c2580-79c1-11eb-9f8e-faefa0f1b63a.png) are the positions of the point charges.

The force of drag is:

![eq17](https://user-images.githubusercontent.com/52175303/109432709-654b1300-79c1-11eb-930d-ad3474b21055.png)

The coefficient of drag depends on the Reynolds number of the particle:

![eq18](https://user-images.githubusercontent.com/52175303/109432726-785de300-79c1-11eb-9164-5aaf57bc190f.png)

![image](https://user-images.githubusercontent.com/52175303/109432758-a04d4680-79c1-11eb-8fe4-4ad6f9b553ea.png)

Once the total force on a member is known, its position and velocity are updated using the Forward Euler method.

![eq19](https://user-images.githubusercontent.com/52175303/109432780-b6f39d80-79c1-11eb-9460-c5acf5b3bd00.png)

![eq20](https://user-images.githubusercontent.com/52175303/109432792-c7a41380-79c1-11eb-8a62-c6cf46014017.png)

The filament used in this simulation is a mixture of two materials, therefore some of its material properties are a mixture of those of two disparate materials. Its effective density, and charge capacity are:

![eq21](https://user-images.githubusercontent.com/52175303/109432815-e86c6900-79c1-11eb-9669-e0bccc7afc68.png)

![eq22](https://user-images.githubusercontent.com/52175303/109432833-f8844880-79c1-11eb-91de-3fa0fa8dd166.png)

Where ![eq23](https://user-images.githubusercontent.com/52175303/109432850-0b971880-79c2-11eb-9e0d-0843c09316bc.png) and ![eq24](https://user-images.githubusercontent.com/52175303/109432863-1ea9e880-79c2-11eb-9f0b-c249a0b1b879.png) are the densities, and ![eq25](https://user-images.githubusercontent.com/52175303/109432888-3d0fe400-79c2-11eb-80d1-8c61982b53e2.png) and ![eq26](https://user-images.githubusercontent.com/52175303/109432899-4b5e0000-79c2-11eb-9789-1554ff37b16a.png) are the charge capacities of the materials in the mixture, and ![eq27](https://user-images.githubusercontent.com/52175303/109432909-59138580-79c2-11eb-8af6-0a169734eea9.png) is the volume fraction of the mixture.

## Results and Discussion

### Full Model with Electrical Forces



