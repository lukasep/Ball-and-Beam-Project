# UC Berkeley EE222/ME237 Nonlinear Systems Ball and Beam Project


Progress:
Try 2 has Feedback Linearization with LQR 
Average Tracking Error: 0.0005 
Average Energy Consumption: 0.0152 
Safety Contraint Violation: 0 
Tracking Cost: 0.84 
Energy Cost: 0.08 
Safety Cost: 0.00 
Total Score: 0.91 

Methods: We implement a controller for the ball-and-beam system by combining feedback linearization with an LQR design. The controller is structured as a subclass of matlab.System with a primary method, stepImpl, that is called at every simulation time step, and a constructor that sets up the necessary symbolic derivations and computes the LQR gain. The system dynamics are defined by a set of nonlinear equations where the state vector consists of the ball position, ball velocity, beam angle, and beam angular velocity. Using the ball position as the controlled output, the controller employs feedback linearization by computing successive Lie derivatives (captured in functions such as lie1_func, lie2_func, and lie3_func) to transform the nonlinear system into a linear form. A symbolic control law is derived by solving for the servo voltage input that forces the fourth derivative of the output to equal a new virtual input, v, which is then computed using MATLAB'S LQR function, based on the error between the transformed system outputs and a reference trajectory provided by the function get_ref_traj(t). The LQR gain is calculated from a linearized state-space representation of the system with matrices A and B using weighting matrices Q and R that can be tuned. The virtual input is converted into the actual servo voltage using the derived feedback linearization function, and safety constraints are applied to ensure that the beam angle remains within predefined limits.


Original doc:

EE222/ME237 Nonlinear Systems, Spring 2025 Starter code and instructions for the course project.

## Project Overview

This project involves designing and testing nonlinear controllers for a ball and beam system. The objective is to develop controllers that stabilize the ball at a desired position on the beam. You will first implement your controllers in MATLAB simulations and later test them on physical hardware.

## Understanding the Problem

To gain a full understanding of the problem and project expectations, please refer to the following documents in this repository:

[`EE_222_Course_Project.pdf`](EE_222_Course_Project.pdf) – Overview of the project and system model. (Disregard the due dates and GitHub link in this older document)

[`EE222 Lab Feedback and FAQ.pdf`](EE222_Lab_Feedback_and_FAQ.pdf) – Common issues and recommendations.

[`EE222_Lab_Part_1_Simulation.pdf`](EE222_Lab_Part_1_Simulation.pdf) – Instructions for running simulation.

[`EE222_Lab_Part_2_Hardware_Testing.pdf`](EE222_Lab_Part_2_Hardware_Testing.pdf) – Instructions for hardware testing. (To be updated)

## Code Instructions

### Prerequisites

Install MATLAB and Simulink using the Berkeley academic license.

### Getting Started

Clone or fork this repository.

Run `setup.m` or manually add the repository and its subfolders to the MATLAB path.

Modify only studentControllerInterface.m to implement your controller.

To test your controller:

Run `run_matlab_ball_and_beam.m` for a MATLAB-based simulation.

Run `run_simulink_ball_and_beam.m` for a Simulink-based simulation.
