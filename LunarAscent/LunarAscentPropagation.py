"""
Copyright (c) 2010-2022, Delft University of Technology
All rights reserved

This file is part of the Tudat. Redistribution and use in source and
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.

AE4866 Propagation and Optimization in Astrodynamics
Lunar Ascent
First name: ***COMPLETE HERE***
Last name: ***COMPLETE HERE***
Student number: ***COMPLETE HERE***

This module computes the dynamics of a Lunar ascent vehicle, according to a simple thrust guidance law.  This file propagates the dynamics
using a variety of integrator and propagator settings. For each run, the differences w.r.t. a benchmark propagation are
computed, providing a proxy for setting quality. The benchmark settings are currently defined semi-randomly, and are to be
analyzed/modified.

The propagtion starts with a small velocity close to the surface of the Moon, and an initial flight path angle of 90
degrees. Making (small) adjustments to this initial state is permitted if properly argued in the report.

The propagation is terminated as soon as one of the following conditions is met:

- Altitude > 100 km
- Altitude < 0 km
- Propagation time > 3600 s
- Vehicle mass < 2250 kg

This propagation assumes only point mass gravity by the Moon and thrust acceleration of the vehicle. Both the
translational dynamics and mass of the vehicle are propagated, using a fixed specific impulse.

The thrust is computed based on a constant thrust magnitude, and a variable thrust direction. The trust direction is defined
on a set of 5 nodes, spread evenly in time. At each node, a thrust angle theta is defined, which gives the angle between
the -z and y angles in the ascent vehicle's vertical frame (see Mooij, 1994, "The motion of a vehicle in a planetary
atmosphere" ). Between the nodes, the thrust is linearly interpolated. If the propagation goes beyond the bounds of
the nodes, the boundary value is used. The thrust profile is parameterized by the values of the vector thrust_parameters.
The thrust guidance is implemented in the LunarAscentThrustGuidance class in the LunarAscentUtilities.py file.

The entries of the vector 'thrust_parameters' contains the following:
- Entry 0: Constant thrust magnitude
- Entry 1: Constant spacing in time between nodes
- Entry 2-6: Thrust angle theta, at nodes 1-5 (in order)

Details on the outputs written by this file can be found:
- benchmark data: comments for 'generateBenchmarks' function
- results for integrator/propagator variations: comments under "RUN SIMULATION FOR VARIOUS SETTINGS"

Frequent warnings and/or errors that might pop up:
* One frequent warning could be the following (mock values):
    "Warning in interpolator, requesting data point outside of boundaries, requested data at 7008 but limit values are
    0 and 7002, applying extrapolation instead."
It can happen that the benchmark ends earlier than the regular simulation, due to the smaller step size. Therefore,
the code will be forced to extrapolate the benchmark states (or dependent variables) to compare them to the
simulation output, producing a warning. This warning can be deactivated by forcing the interpolator to use the boundary
value instead of extrapolating (extrapolation is the default behavior). This can be done by setting:

    interpolator_settings = interpolators.lagrange_interpolation(
        8, boundary_interpolation = interpolators.extrapolate_at_boundary)

* One frequent error could be the following:
    "Error, propagation terminated at t=4454.723896, returning propagation data up to current time."
    This means that an error occurred with the given settings. Typically, this implies that the integrator/propagator
    combination is not feasible. It is part of the assignment to figure out why this happens.

* One frequent error can be one of:
    "Error in RKF integrator, step size is NaN"
    "Error in ABM integrator, step size is NaN"
    "Error in BS integrator, step size is NaN"

This means that a variable time-step integrator wanting to take a NaN time step. In such cases, the selected
integrator settings are unsuitable for the problem you are considering.

NOTE: When any of the above errors occur, the propagation results up to the point of the crash can still be extracted
as normal. It can be checked whether any issues have occured by using the function

dynamics_simulator.integration_completed_successfully

which returns a boolean (false if any issues have occured)

* A frequent issue can be that a simulation with certain settings runs for too long (for instance if the time steo
becomes excessively small). To prevent this, you can add an additional termination setting (on top of the existing ones!)

    cpu_tim_termination_settings = propagation_setup.propagator.cpu_time_termination(
        maximum_cpu_time )

where maximum_cpu_time is a varaiable (float) denoting the maximum time in seconds that your simulation is allowed to
run. If the simulation runs longer, it will terminate, and return the propagation results up to that point.

* Finally, if the following error occurs, you can NOT extract the results up to the point of the crash. Instead,
the program will immediately terminate

    SPICE(DAFNEGADDR) --

    Negative value for BEGIN address: -214731446

This means that a state is extracted from Spice at a time equal to NaN. Typically, this is indicative of a
variable time-step integrator wanting to take a NaN time step, and the issue not being caught by Tudat.
In such cases, the selected integrator settings are unsuitable for the problem you are considering.
"""

###########################################################################
# IMPORT STATEMENTS #######################################################
###########################################################################

import sys
sys.path.insert(0, '/home/dominic/Software/tudat-bundle/build-tudat-bundle-Desktop-Default/tudatpy/')

# General imports
import numpy as np
import os

# Tudatpy imports
from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.math import interpolators

# Problem-specific imports
import LunarAscentUtilities as Util

###########################################################################
# DEFINE GLOBAL SETTINGS ##################################################
###########################################################################

# Load spice kernels
spice_interface.load_standard_kernels()
# Define problem independent variables
thrust_parameters_list = [
[13232.2025345638,47.5319798593,-0.0128010195,0.0507979044,0.6538417737,-0.5560136572,1.0214363814],
[14041.4505556691,10.0102942972,0.0099324952,-0.2090952564,0.6617580954,-0.5865616919,0.8351957248],
[17869.1842977423,21.5312002995,0.0895461222,-0.3786714207,0.4978693228,-0.2725262092,-1.132938021],
[16672.3513533361,22.6348243258,0.0693122983,-0.2407475549,-0.638175925,0.2575758141,0.9782958947],
[12104.00634096,70.3475237079,-0.0664353371,-0.479385498,0.5952144227,0.9919778286,0.8465391759],
[16613.5052975733,47.8996858629,-0.0596513546,0.1723836786,0.6126002189,-0.5908889086,1.2849393354],
[10797.3346707877,47.1284960583,-0.0761030402,0.2135657477,-0.2574372935,0.632335363,0.8361722032],
[8155.7384249754,91.7735954025,-0.0126505475,0.0862529017,0.2613339143,0.9629374784,0.5450343889],
[10372.2825529985,18.6555033061,-0.0292649253,0.3710720506,-0.1977270775,-0.8709862055,-0.6660235576],
[10455.6616011541,21.7025712226,0.0907484445,-0.0453779178,0.17490319,0.0929129929,0.7042905128],
[9875.7084133104,57.2319444711,0.0596400778,-0.4247998863,-0.6336323069,0.2174051679,0.9263492577],
[11994.6615921799,63.0374983558,0.0987704021,0.1910929221,0.4622633364,0.6333203944,0.8949938918],
[19651.4163620304,22.5348709058,-0.0480510497,0.2239391366,0.0206275638,-0.2690444668,-1.147023045]]

for count in range(len(thrust_parameters_list)):

    thrust_parameters=thrust_parameters_list[count]
    print(thrust_parameters)
    # Choose whether benchmark is run
    use_benchmark = True
    run_integrator_analysis = True

    # Choose whether output of the propagation is written to files
    write_results_to_file = True
    # Get path of current directory
    current_dir = os.path.dirname(__file__)

    ###########################################################################
    # DEFINE SIMULATION SETTINGS ##############################################
    ###########################################################################

    # Set simulation start epoch
    simulation_start_epoch = 0.0  # s
    # Vehicle settings
    vehicle_mass = 4.7E3  # kg
    vehicle_dry_mass = 2.25E3  # kg
    constant_specific_impulse = 311.0  # s
    # Fixed simulation termination settings
    maximum_duration = constants.JULIAN_DAY  # s
    termination_altitude = 100.0E3  # m

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Define settings for celestial bodies
    bodies_to_create = ['Moon']
    # Define coordinate system
    global_frame_origin = 'Moon'
    global_frame_orientation = 'ECLIPJ2000'

    # Create body settings
    # N.B.: all the bodies added after this function is called will automatically
    # be placed in the same reference frame, which is the same for the full
    # system of bodies
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create,
        global_frame_origin,
        global_frame_orientation)
    # Create bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Create vehicle object and add it to the existing system of bodies
    bodies.create_empty_body('Vehicle')
    # Set mass of vehicle
    bodies.get_body('Vehicle').mass = vehicle_mass

    ###########################################################################
    # CREATE PROPAGATOR SETTINGS ##############################################
    ###########################################################################

    # Retrieve termination settings
    termination_settings = Util.get_termination_settings(simulation_start_epoch,
                                                         maximum_duration,
                                                         termination_altitude,
                                                         vehicle_dry_mass)
    # Retrieve dependent variables to save
    dependent_variables_to_save = Util.get_dependent_variable_save_settings()
    # Check whether there is any
    are_dependent_variables_to_save = False if not dependent_variables_to_save else True

    ###########################################################################
    # IF DESIRED, GENERATE BENCHMARK ##########################################
    ###########################################################################

    if use_benchmark:

        # Define benchmark interpolator settings to make a comparison between the two benchmarks
        benchmark_interpolator_settings = interpolators.lagrange_interpolation(
            8,boundary_interpolation = interpolators.extrapolate_at_boundary)

        # Create propagator settings for benchmark (Cowell)
        propagator_settings = Util.get_propagator_settings(
            thrust_parameters,
            bodies,
            simulation_start_epoch,
            constant_specific_impulse,
            vehicle_mass,
            termination_settings,
            dependent_variables_to_save)

        benchmark_output_path = current_dir + '/SimulationOutput/benchmarks_' + str(count) + '/' if write_results_to_file else None

        # Generate benchmarks
        benchmark_step_size = 1.0
        benchmark_list = Util.generate_benchmarks(benchmark_step_size,
                                                  simulation_start_epoch,
                                                  bodies,
                                                  propagator_settings,
                                                  are_dependent_variables_to_save,
                                                  benchmark_output_path)

        # Extract benchmark states
        first_benchmark_state_history = benchmark_list[0]
        second_benchmark_state_history = benchmark_list[1]

        # Create state interpolator for first benchmark
        benchmark_state_interpolator = interpolators.create_one_dimensional_vector_interpolator(
            first_benchmark_state_history,
            benchmark_interpolator_settings)

        # Compare benchmark states, returning interpolator of the first benchmark, and writing difference to file if
        # write_results_to_file is set to True
        benchmark_state_difference = Util.compare_benchmarks(first_benchmark_state_history,
                                                             second_benchmark_state_history,
                                                             benchmark_output_path,
                                                             'benchmarks_state_difference.dat')

        # Extract benchmark dependent variables, if present
        if are_dependent_variables_to_save:
            first_benchmark_dependent_variable_history = benchmark_list[2]
            second_benchmark_dependent_variable_history = benchmark_list[3]

            # Create dependent variable interpolator for first benchmark
            benchmark_dependent_variable_interpolator = interpolators.create_one_dimensional_vector_interpolator(
                first_benchmark_dependent_variable_history,
                benchmark_interpolator_settings)

            # Compare benchmark dependent variables, returning interpolator of the first benchmark, and writing difference
            # to file if write_results_to_file is set to True
            benchmark_dependent_difference = Util.compare_benchmarks(first_benchmark_dependent_variable_history,
                                                                     second_benchmark_dependent_variable_history,
                                                                     benchmark_output_path,
                                                                     'benchmarks_dependent_variable_difference.dat')

