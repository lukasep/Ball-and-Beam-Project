# UC Berkeley EE222/ME237 Nonlinear Systems Ball and Beam Project

**Lab Project Phase I: Simulations**  
Soomi Lee, Arvind Kruthiventy, Emily Lukas  
4 April 2025  

**Final Controllers:**
- Feedback Linearization with LQR
- Linearization around Trajectory

**Observer:**
- Luenberger Observer

**Results:**  
Controller 1:  
Sine wave animation: https://github.com/lukasep/Ball-and-Beam-Project/blob/main/controller1video.mp4  
Square wave animation: https://github.com/lukasep/Ball-and-Beam-Project/blob/main/controlelr1_sqaure.mp4
-Average Tracking Error: 0.0005  
-Average Energy Consumption: 0.0091  
-Safety Constraint Violation: 0  
-Tracking Cost: 0.94  
-Energy Cost: 0.05  
-Safety Cost: 0.00  
**Total Score: 0.98**


Controller 2: https://github.com/lukasep/Ball-and-Beam-Project/blob/main/controller2video.mp4  
-Average Tracking Error: 0.0004  
-Average Energy Consumption: 0.0309  
-Safety Constraint Violation: 0  
-Tracking Cost: 0.67  
-Energy Cost: 0.15  
-Safety Cost: 0.00  
**Total Score: 0.82**


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
