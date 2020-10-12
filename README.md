A repository of past Partial Differential Equation projects. Intro to PDES projects 1-3 were completed in MATLAB in 2018 for the Imperial introduction to computational PDE's course.

### Intro to PDE's Project 1 2018
*  Implemented Centered Finite difference method on a linked diffusion system in polar coordinates.
*  Example contour graphs shown below for some fixed times.

![example_graph](https://user-images.githubusercontent.com/58078485/95724705-94f5ba80-0c6e-11eb-8428-b497f0791ad0.png)

### Intro to PDE's Project 2 2018
Project was  concerned with numerically solving equations of the form on a uniform rectangular grid:
![image](https://user-images.githubusercontent.com/58078485/95731668-b7d89c80-0c77-11eb-9522-632b230ee670.png)
![image](https://user-images.githubusercontent.com/58078485/95730833-b5297780-0c76-11eb-9672-72e0fbbb6525.png)

*  Numerical Solutions for equation (1) were achieved by using a Centred Finite Difference Method along with Gauss Seidel and Multi-Grid method to improve convergence. 
*  Numerical Solutions for equation (10) were achieved by using a Centred Finite Difference and Crank Nicolson Method along with an altered Gauss Seidel and Multi-Grid methods to improve convergence.
*  Below is an example simulation with a known solution for error comparison.

![image](https://user-images.githubusercontent.com/58078485/95731442-71833d80-0c77-11eb-88da-47db4a045187.png)

### Intro to PDE's Project 3 2018
Project focused on the solving of Turing Systems shown below using similar methods to Project 2 with a focus on stabilty analysis and non-regular domains:

![image](https://user-images.githubusercontent.com/58078485/95735471-c07fa180-0c7c-11eb-8c12-818ed77de28f.png)

For example the Gauss Seidel step now solves the following:

![image](https://user-images.githubusercontent.com/58078485/95737121-19503980-0c7f-11eb-8828-b586b1cc9757.png)


* The graphs below show how the numerical simulation can generate a spots pattern which can typically be seen on animals. This was done via domain decomposition. 

![image](https://user-images.githubusercontent.com/58078485/95735141-5a931a00-0c7c-11eb-8092-597a0629c696.png) 

* The graph below shows a stripes pattern and spots pattern with their relative stabilities.

![image](https://user-images.githubusercontent.com/58078485/95736416-20c31300-0c7e-11eb-9000-34ca47935136.png) 


