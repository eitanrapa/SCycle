% Script demonstrating how to load and visualize example 2.
% Example 2 is a quasidynamic earthquake cycle simulation with a vertical
% strike-slip fault, and linear elastic off-fault material.
%
% Run ex2 from the repository root (./source/main examples/ex2.in), then run
% this script from the examples directory.

addpath('../matlab/visualizePetsc')

% output prefix (outputDir in ex2.in), relative to this directory
sourceDir = '../data/ex2_';

% context (domain.txt, ..., data_context.h5): d.y and d.z are [Nz, Ny],
% so d.z(:,1) is depth along the fault
d = loadContext_hdf5(sourceDir);

% fault time series (data_1D.h5): d.tau, d.slipVel, d.slip, d.psi are [Nz, Nt]
d = loadData1D_hdf5(d, sourceDir);

%% plot results
z = d.z(:,1);

% friction parameters and normal stress
figure(1),clf
subplot(1,2,1)
plot(d.fault.sNEff, z)
set(gca,'YDir','reverse')
xlabel('\sigma_N (MPa)'),ylabel('depth (km)'),title('effective normal stress')

subplot(1,2,2)
plot(d.fault.a, z, 'k-'),hold on
plot(d.fault.a - d.fault.b, z, 'r-')
plot(d.fault.b, z, 'b-')
plot([0 0],[0 z(end)],'-','Color',[1 1 1].*0.6)
set(gca,'YDir','reverse')
legend('a','a-b','b'),ylabel('depth (km)'),title('rate-and-state parameters')

% shear stress
figure(2),clf
tauStride = ceil(length(d.time)/100);
h3 = plot(d.tau(:,1:tauStride:end), z, 'c-');
hold on
h1 = plot(d.tau(:,1), z, '.-','Color',[0,128,0]./255,'Linewidth',1);
h2 = plot(d.tau(:,end), z, 'r.-','Linewidth',1);
set(gca,'YDir','reverse')
title(sprintf('%i: %.9e',length(d.time),d.time(end)))
ylabel('depth (km)'),xlabel('\tau (MPa)')
grid on, grid minor
legend([h1 h2 h3(1)],{'initial \tau','final \tau','intermediate values of \tau'},'Location','Southeast')

% slip velocity
[map,~] = createDivColormap(-14,-9,1,100); % colormap highlighting the loading velocity
figure(3),clf
pcolor(1:size(d.slipVel,2), z', log10(abs(d.slipVel)))
shading flat
colormap(map)
hcb = colorbar; ylabel(hcb,'log V (m/s)')
set(hcb,'YTick',-14:1:1)
caxis([-14 1]),ylim([0 d.dom.Lz])
set(gca,'YDir','reverse')
title('slip velocity'),ylabel('depth (km)'),xlabel('model step count')
grid on, grid minor

% phase plot: depth-integrated slip velocity vs depth-integrated shear stress
p = trapz(z, d.slipVel);
t = trapz(z, d.tau);

figure(4),clf
semilogx(p,t) % phase plot
hold on
semilogx(p(1),t(1),'g*') % initial condition
semilogx(p(end),t(end),'r*') % final condition
xlabel('integrated slip velocity (m/s km)'),ylabel('integrated shear stress (MPa km)')
legend('simulation','initial condition','final condition','Location','Northwest')
