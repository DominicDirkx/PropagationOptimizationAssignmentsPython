clc
clear all
close all

for i=1:13
    
    ascent_states = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LunarAscent/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_states.dat'));
    ascent_dependent = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LunarAscent/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_dependent_variables.dat'));
    
    figure(1)
    plot3(ascent_states(:,2),ascent_states(:,3),ascent_states(:,4))
    axis equal
    hold on
    
    figure(2)
    for j=1:3
        subplot(2,2,j)
        plot(ascent_dependent(:,1),ascent_dependent(:,j+1))
        hold on
    end
    
    
    transfer_states = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LowThrust/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_states.dat'));
    transfer_dependent = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/LowThrust/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_dependent_variables.dat'));
    
    figure(3)
    plot3(transfer_states(:,2),transfer_states(:,3),transfer_states(:,4))
    axis equal
    hold on
    
    figure(4)
    for j=1:3
        subplot(2,2,j)
        plot(transfer_dependent(:,1)-transfer_dependent(1,1),transfer_dependent(:,j+1))
        hold on
    end
    subplot(2,2,4)
    
    plot(transfer_states(:,1)-transfer_states(1,1),transfer_states(:,8))
    hold on
    
    entry_states = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/ShapeOptimization/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_states.dat'));
    entry_dependent = load(strcat('/home/dominic/Software/PropagationOptimization2022/PropagationOptimizationAssignmentsPython/ShapeOptimization/SimulationOutput/benchmarks_',num2str(i-1),'/benchmark_1_dependent_variables.dat'));
    
    figure(5)
    plot3(entry_states(:,2),entry_states(:,3),entry_states(:,4))
    axis equal
    hold on
    
    figure(6)
    for j=1:2
        subplot(2,2,j)
        plot(entry_dependent(:,1)-entry_dependent(1,1),entry_dependent(:,j+1))
        hold on
    end
    
    
end



