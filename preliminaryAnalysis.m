clc
clear all
close all

low_thrust_states = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LowThrust/SimulationOutput/benchmarks/benchmark_1_states.dat');
low_thrust_dependent = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LowThrust/SimulationOutput/benchmarks/benchmark_1_dependent_variables.dat');

figure
plot3(low_thrust_states(:,2),low_thrust_states(:,3),low_thrust_states(:,4))
axis equal

figure
for i=1:3
    plot(low_thrust_dependent(:,1),low_thrust_dependent(:,i+1))
    hold on
end

entry_states = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/ShapeOptimization/SimulationOutput/benchmarks/benchmark_1_states.dat');
entry_dependent = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/ShapeOptimization/SimulationOutput/benchmarks/benchmark_1_dependent_variables.dat');

figure
plot3(entry_states(:,2),entry_states(:,3),entry_states(:,4))
axis equal

figure
for i=1:2
    subplot(2,2,i)
    plot(entry_dependent(:,1),entry_dependent(:,i+1))
    hold on
end

ascent_states = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LunarAscent/SimulationOutput/benchmarks/benchmark_1_states.dat');
ascent_dependent = load('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LunarAscent/SimulationOutput/benchmarks/benchmark_1_dependent_variables.dat');

figure
plot3(ascent_states(:,2),ascent_states(:,3),ascent_states(:,4))
axis equal

figure
for i=1:3
    subplot(2,2,i)
    plot(ascent_dependent(:,1),ascent_dependent(:,i+1))
    hold on
end



