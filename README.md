# SutureModel
This repository hosts the code for:
the ICRA 2025 paper "Suture Thread Modeling Using Control Barrier Functions for Autonomous Surgery"

This repository contains the MATLAB implementation for real-time control of needle position and suture thread simulation in a hernia repair task. The control system allows interactive movement of the needle using keyboard inputs while simulating the thread’s interaction with a tissue model.

If you have access to a joystick download Sim Playground folder and run main_joystick.

# Main Script
# main.m
Run this script to start the real-time control of the needle position and observe the suture thread simulation.
Use the keyboard controls:
W → Move Up
A → Move Left
S → Move Down
D → Move Right
The thread color changes dynamically based on its interaction with the tissue model.

# Barrier Certificate Functions
These functions implement control barrier certificates to ensure safe and stable suture behavior while accounting for material stiffness.

# create_si_connectivity_barrier_certificate_with_obstacles_stiff.m
Use this function when incorporating suture stiffness properties.
Adjust k2 to modify the material stiffness:
Increase k2 → More stiff thread.
Recommended value: k2 = 1 (empirically tuned for Polyamide suture).
# create_si_connectivity_barrier_certificate_with_obstacles.m
Use this function when suture stiffness is negligible.
Works best for Silk suture, where flexibility dominates.
