%% ---[ Simulator: Random Electron Scattering (single or multiple particles) ]---
%   Monte Carlo Modeling
%
%   Author: Ragini Bakshi, Jan 2021

set(0,'DefaultFigureWindowStyle','docked')
set(0, 'defaultaxesfontsize', 20)
set(0, 'defaultaxesfontname', 'Times New Roman')
set(0, 'DefaultLineLineWidth',2);

clear all
close all

x = 0; % electron x position
v = 0; % electron velocity
t = 0; % time elapsed

F = 1;               % Force experienced
m = 9.10938215e-31;  % electron mass
%m=1;
np = 1;              
dt = 1;              % 1 sec time intervals
steps = 100;
tau = 10;            % est mean time between collisions

% multiple particles
v = zeros(np,1);
x = zeros(np,1);

re = 0; % rebound behaviour: make 0 for testing

for i = 2:steps
    t(i) = t(i-1) + dt;
    v(:,i) = v(:, i-1) + F/m*dt;
    x(:,i) = x(:, i-1) + v(:, i-1)*dt + F/m*dt^2/2;
    
    if np == 1
        p = 0.05;
    else
        p = 1-exp(dt/tau);
    end
    
    r = rand(np,1) < p;
    v(r,i) = re*v(r,i);
    average_v(i,:) = mean(v,2)
    
    subplot(3,1,1)
    plot(t,v,'-');
    hold on
    subplot(3,1,1)
    plot(t,average_v,'g*');
    hold off
    xlabel('time')
    ylabel('avg velocity')
    title(['Average Drift Velocity:' num2str(average_v(i))])
    
    subplot(3,1,2)
    plot(x(1,:),v(1,:),'r-');
    hold on
    subplot(3,1,2)
    if np == 1
        plot(x(1,:),average_v(:,1),'g*');
    end
    hold off
    xlabel('x (position)')
    ylabel('avg velocity')    

    subplot(3,1,3)
    plot(t,x,'-');
    xlabel('time')
    ylabel('x (position)')
    
    pause(0.01)
end

display('Average V')
mean(v)

