# Computational Fluid Dynamics Assignments

3 computational homeworks for Boğaziçi University's **Fluid Mechanics** (ME 353) class:

1. Dimension Analysis
2. Effect of Mass Flow Rate on Major Head Loss
3. Finding Boundary Layer Parameters

## HW1-Dimension Analysis
The purpose of this assignment is to determine the necessary area of a paraglider to land a spacecraft. 

A dimensional analysis is carried out using existing data and a relationship is found between the variables. To show the relationship, curve fitting is applied using three different functions: second degree polynomial, power series, and exponential.

Variables:
* $\mu$: dynamic viscosity of the air
* $\rho$: density of the air
* *m*: mass of the spacecraft
* *A*: area of the paraglider

Power series fit:

$$f(x) = ax^{b}$$ 

where $a=4.21\times 10^6, b=-1.5$

The formula that gives the area:

$$area = (\frac{m^{b+1}\nu}{a\mu d^{b}})^{\frac{2}{3b+2}}$$

Necessary area with respect to different velocities and masses are shown below:

<p align="center">
  <img src="https://github.com/edizferit/CFD_Assignments/blob/main/figures/hw1.jpg?raw=true" width="50%">
</p>

## HW2-Effect of Mass Flow Rate on Major Head Loss

Friction factor can be calculated by either using Colebrook formula which is a transcendental equation, or with Haaland formula which gives an approximate result. In this HW, accuracy of Haaland formula with respect to Colebrook formula is calculated for a certain mass flow rate interval.

* Secant method is used to solve Colebrook formula for a certain mass flow rate interval.
* Major head losses regarding two methods are found by Darcy-Weisbach equation.

Percent error of the results of Haaland formula with respect to the Colebrook formula is plotted on mass flow rate below:

<p align="center">
  <img src="https://github.com/edizferit/CFD_Assignments/blob/main/figures/hw2.jpg?raw=true" width="50%">
</p>

## HW3-Laminar Boundary Layer on a Flat Plate

### Blasius Solution:

The following second order differential equation (Blasius Equation) is transformed to a set of first order differential equations and solved with ode45 solver.

$$2f''' + f\times f'' = 0$$

Result:

<p align="center">
  <img src="https://github.com/edizferit/CFD_Assignments/blob/main/figures/hw31.jpg?raw=true" width="50%">
</p>

### Boundary Layer Thickness:

Followings are calculated by numerical integration with trapezoidal rule and integration errors are estimated:
* Boundary layer thickness $\delta$
* Boundary layer displacement thickness $\delta^{*}$ 
* Boundary layer momentum thickness $\theta$

### Velocity and Shear Stress Distribution in the Boundary Layer:

<p align="center">
  <img src="https://github.com/edizferit/CFD_Assignments/blob/main/figures/hw32.jpg?raw=true" width="80%">
</p>


