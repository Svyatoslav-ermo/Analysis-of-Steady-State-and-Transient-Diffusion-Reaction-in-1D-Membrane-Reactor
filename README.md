![image](https://github.com/user-attachments/assets/cc3caab5-f10e-4195-8fa2-1536b267e52a)# Analysis-of-Steady-State-and-Transient-Diffusion-Reaction-in-1D-Membrane-Reactor
The code simulates 4 chemical reactions to analyze 1D concentration profiles of 4 chemical species (A, B, U, F) diffusing in 1 direction inside a membrane reactor. Steady-state and transient-state solutions were obtained using multivariable Newton and implicit Euler methods, respectively.

Run main.m with all files in the same folder.

Uncomment either of the (a), (b), or (c) in main.m to test different conditions.

Graphical results include:
1. 2D plot of chemical species concentrations (Y) in steady-state along membrane reactor length (X)
2. 3D plot of how chemical species concentrations (Y) change over time (t) along membrane reactor length (X)

# Introduction

The following 4 chemical reactions were studied:
1. A + 2B -> D
2. A -> P
3. B -> U
4. U -> F

The reactions took place in a 1D membrane chemical reactor as follows:

![image](https://github.com/user-attachments/assets/8335ee75-317f-4154-8f42-790e81c53603)

where C_AF and C_BF are fixed concentrations of reactants A and B at r = 0 diffusing into the membrane reactor in direction r. It is assumed that the molecular diffusion happens in a one-dimensional horizontal r direction and follows Fick's law. The species deposit at r = L with different deposition rates.

The following partial differential equations describe the change in species concentrations with position r and time t:

![image](https://github.com/user-attachments/assets/a2b5816d-291b-49d1-afa0-06e5766403b2)

where k are chemical reaction rates and D is the diffusion coefficient (the same for all species)

# Dimensionless model

The following dimensionless variables and parameters were defined:

![image](https://github.com/user-attachments/assets/1bf384f2-5926-4cb5-9457-fc53c6e2af82)

The above partial differential equations in Introduction section were written in dimensionless form as follows:

![image](https://github.com/user-attachments/assets/694ed04b-4109-4257-b221-64acaa20cfd9)

The boundary conditions are given as follows:

![image](https://github.com/user-attachments/assets/7399937d-3c64-4fea-8c17-7a0e22c2d36d)

# Methods to Obtain Steady-State Solution

A steady-state solution was obtained by applying the centered finite difference method and dividing the membrane reactor into N nodes. The 4 ordinary differential equations with time derivative set to 0 were discretized and implemented in the Jacobian matrix to solve for steady-state concentrations using multivariable Newton's method. At each Newton's method iteration, the LU decomposition with partial pivoting was performed to solve the multivariable system describing species concentrations. The LU decomposition was chosen as the method to solve the system because it is the most efficient in solving band matrices, and the Jacobian is a band matrix. Newton's method was iterated until the square of vector 2-norm was less than 0.01.

The least number of nodes required to obtain an accurate steady-state solution is 41. To determine the least number of nodes required for accurate steady-state solution, concentrations of species A at x = 0.3 were compared for results with different numbers of nodes. The testing algorithm was developed that started Newton’s method computation of results with N = 11 nodes. The number of nodes was increased by 10 at each iteration of the testing algorithm and the results of the subsequent and previous iterations were compared by calculating the percentage difference. The testing algorithm stopped once the percentage difference became less than 0.01%. 41 is the number of nodes that give a percentage difference of less than 0.01%.

![image](https://github.com/user-attachments/assets/46fd7db6-5ff7-4315-9871-795744b36717)

**Steady-state Concentration Profile for Conditions (a)**

![image](https://github.com/user-attachments/assets/4d0b1440-9ab5-469e-a7f6-288f6205dfc9)

**Steady-state Concentration Profile for Conditions (b)**

![image](https://github.com/user-attachments/assets/9bc5ca70-37fd-445a-a61b-abab0b0f4799)

**Steady-state Concentration Profile for Conditions (c)**

# Methods to Obtain Transient Solution

At the first attempt, explicit Euler was attempted to integrate partial differential equations with respect to time. However, because explicit Euler is unstable for stiff functions, results diverged and an extremely small time step on the order of 10^(-8) was required to compute the transient solution, which required unrealistic computation time. Therefore, implicit Euler was implemented as a method to integrate differential equations with respect to time because implicit Euler is unconditionally stable and works with any time step. The assembled implicit Euler system with the Jacobian was solved as an ODE using multivariable Newton's method and LU decomposition.

The implicit Euler integrated chemical species concentration in the time domain until 99% of the steady-state solutions, determined in the previous section, was reached. To calculate how close transient state is to the steady-state solution, the difference between the current time-step vector and the steady state was taken, and the square of the 2-norm of the difference was calculated. When the square of the 2-norm became less than 0.01, the transient state reached 99% of the steady state.

![image](https://github.com/user-attachments/assets/b23d0ecc-04c0-4c81-9eae-942e79daaa99)

**Change in Chemical Species Concentrations (Y) in Time (t) Along Reactor Length (X) for Conditions (c)**
