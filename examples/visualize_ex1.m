% Script demonstrating how to load and visualize example 1.
% Example 1 is a quasidynamic earthquake cycle simulation for a spring slider
% with linear elastic off-fault material.
%
% Run ex1 from the repository root (./source/main examples/ex1.in), then run
% this script from the examples directory.

addpath('../matlab/visualizePetsc')

% output prefix (outputDir in ex1.in), relative to this directory
sourceDir = '../data/ex1_';

% context (domain.txt, ..., data_context.h5) and the fault time series (data_1D.h5)
d = loadContext_hdf5(sourceDir);
d = loadData1D_hdf5(d, sourceDir);

%% make plots
yr = 3.1536e7; % seconds per year

figure(1),clf

% plot shear stress
subplot(2,1,1)
plot(d.time./yr, d.tau)
xlabel('time (years)'),ylabel('\tau (MPa)')
title('Spring Slider')

% plot slip velocity
subplot(2,1,2)
semilogy(d.time./yr, d.slipVel)
xlabel('time (years)'),ylabel('V (m/s)'),ylim([1e-14 10])

% phase plot: slip velocity vs shear stress
figure(2),clf
semilogx(d.slipVel, d.tau) % phase plot
hold on
semilogx(d.slipVel(1), d.tau(1), 'g*') % initial condition
semilogx(d.slipVel(end), d.tau(end), 'r*') % final condition
xlabel('slip velocity (m/s)'),ylabel('shear stress (MPa)')
legend('simulation','initial condition','final condition','Location','Northwest')
