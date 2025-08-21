% Attenuation mechanism maps! 
% DH 8/2025
%
% This code calculates J1 and J2 for the
% linearized backstress model, the
% high-temperature background (HTB) of the premelt
% model and the peak of the premelt model for
% defined ranges of frequency (f), temperature
% (T), grain size (d), and background tectonic 
% stress (sig) using the VBRc. Then their complex 
% compliances are summed to obtain their combined 
% response, which is stored in J1_tot, J2_tot, 
% invQ_tot, and G_eff_tot. 


addpath('/home/diedehein/Documents/GitHub/VBR')             %add path to VBR
% addpath('/home/diedehein/Documents/MATLAB/ZachEs_Doodads')  %add path to Zach's mindex function (optional)
close all; clear all; clc
vbr_init

%% set variable ranges and constants
Resolution = 61;

%set frequencies to calculate over 
f = logspace(-13,2,Resolution); % [Hz] 
VBR.in.SV.f = f;

%  build 3D grid of temeprature, grain size, stress
T   = linspace(1000, 1600, Resolution) + 273; %K, temperatures
d   = logspace(1,    5,    Resolution)      ; %um, grain sizes
sig = logspace(-1,   3,    Resolution)      ; %MPa, differential stress

[VBR.in.SV.dg_um, VBR.in.SV.T_K, VBR.in.SV.sig_MPa] = meshgrid(d, T, sig); %creates an array with dimensions in order of T, d, and then sig

sz = size(VBR.in.SV.T_K);

% set constants
VBR.in.SV.phi = 0 * ones(sz);
VBR.in.SV.P_GPa = 3.5 * ones(sz); % pressure [GPa]
VBR.in.SV.Tsolidus_K = -5.104.*VBR.in.SV.P_GPa.^2 + 132.899.*VBR.in.SV.P_GPa + 1120.661 +273; %K , solidus temperature (for premelt model) from Hirschmann, 2000
VBR.in.SV.rho = 3300 * ones(sz); %kgm-3, density

VBR_in = VBR.in;

%% run VBR for each model
% Run VBR with linearized backstress model and save
VBR.in.elastic.methods_list = {'anharmonic'}; % set methods list
VBR.in.viscous.methods_list = {'HZK2011'};
VBR.in.anelastic.methods_list = {'backstress_linear'};
VBR.in.elastic.anharmonic = Params_Elastic('anharmonic'); 
VBR.in.elastic.anharmonic.temperature_scaling = 'isaak';
VBR.in.elastic.anharmonic.pressure_scaling = 'abramson';
[VBR_linBackstress] = VBR_spine(VBR); % run VBR

% Run VBR with peak of premelt model only
clearvars -except VBR_in f T d sig sz VBR_linBackstress % clear previous VBR
VBR.in = VBR_in; % restore state variables
VBR.in.anelastic.methods_list = {'xfit_premelt'}; % set methods list
VBR.in.anelastic.xfit_premelt.A_B=0; % set amplitude of HTB to zero
[VBR_peak] = VBR_spine(VBR); % run VBR

% Run VBR with HTB of premelt model only
clearvars -except VBR_in f T d sig sz VBR_linBackstress VBR_peak % clear previous VBR
VBR.in = VBR_in; %restore state variables
VBR.in.anelastic.methods_list = {'xfit_premelt'}; % set methods list
VBR.in.anelastic.xfit_premelt.Ap_fac_1=0; % set amplitudes of high-frequency peak to zero
VBR.in.anelastic.xfit_premelt.Ap_fac_2=0;
VBR.in.anelastic.xfit_premelt.Ap_fac_3=0;
[VBR_HTB] = VBR_spine(VBR); % run VBR

% Run VBR with steady-state dislocation creep only 
% [TBD!]

%% Combine complex compliances and calculate attenuation and effective modulus for all

% find unrelaxed modulus
Ju = (1/VBR_linBackstress.out.elastic.anharmonic.Gu); %GPa, unrelaxed shear compliance, is equal for all three models

% % summing complex compliances
% J_tot = VBR.out.anelastic.xfit_premelt.J1 + 1i.*VBR.out.anelastic.xfit_premelt.J2 + VBR.out.anelastic.backstress_linear.J1 + 1i.*VBR.out.anelastic.backstress_linear.J2 - Ju;
% J1_tot = real(J_tot);
% J2_tot = imag(J_tot);
% invQ_tot = J2_tot./J1_tot;

% or alternatively, sum J1 and J2s to get invQ and Geff (produces same results as above)
J1_tot =  VBR_HTB.out.anelastic.xfit_premelt.J1 + VBR_peak.out.anelastic.xfit_premelt.J1 + VBR_linBackstress.out.anelastic.backstress_linear.J1 - 2*Ju; %subtracting Ju twice as it is incorporated in J1 of both the HTB, peak, and backstress models%
J2_tot = (VBR_HTB.out.anelastic.xfit_premelt.J2 + VBR_peak.out.anelastic.xfit_premelt.J2 + VBR_linBackstress.out.anelastic.backstress_linear.J2);
invQ_tot = J2_tot./J1_tot; %calculate attenuation, Q-1
G_tot = 1./(J1_tot+1i.*J2_tot); %calculate combined complex shear modulus
G_eff_tot = abs(G_tot); %calculate effective shear modulus

%% plotting spectra of J1, J2, Geff, and invQ at arbitrary conditions
% for troubleshooting purposes, plot J1, J2,
% G_eff, and invQ as a function of frequency for a
% given temperature, grain size, and differential stress

%find indices of conditions that you want to plot
T_target   = 1500+273; %K, target temperature
d_target   = 1000; %um, target grain size
sig_target = 1;    %MPa, target stress

% print selected target conditions to screen
T_slice   = T(mindex(T, T_target)) - 273
d_slice   = d(mindex(d, d_target))
sig_slice = sig(mindex(sig, sig_target))

% find indices of target conditions
iT = mindex(T,T_target);
id = mindex(d,d_target);
isig = mindex(sig,sig_target);

% setup figure
fig = figure('Position',[1500  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

% J1
subplot(4,1,1)
hold on
loglog(VBR.in.SV.f, squeeze(J1_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J1(iT,id,isig,:)),'r-')
loglog(VBR.in.SV.f, squeeze(VBR_peak.out.anelastic.xfit_premelt.J1(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J1(iT,id,isig,:)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
title(['T= ' num2str(T_slice) '^oC, d=' num2str(round(d_slice)) ' um, stress=' num2str(sig_slice) ' MPa'])
xlabel('Frequency (Hz)')
ylabel('J_1')
legend('Total', 'Pre-melt, HTB', 'Premelt, peak', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(VBR.in.SV.f, squeeze(J2_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J2(iT,id,isig,:)),'r-')
loglog(VBR.in.SV.f, squeeze(VBR_peak.out.anelastic.xfit_premelt.J2(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J2(iT,id,isig,:)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Frequency (Hz)')
ylabel('J_2')

% G_eff
subplot(4,1,3)
hold on
loglog(VBR.in.SV.f, squeeze(G_eff_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR_HTB.out.anelastic.xfit_premelt.M(iT,id,isig,:)),'r-')
loglog(VBR.in.SV.f, squeeze(VBR_peak.out.anelastic.xfit_premelt.M(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.M(iT,id,isig,:)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Frequency (Hz)')
ylabel({'Effective shear'; 'modulus (GPa)'})

% invQ
subplot(4,1,4)
hold on
loglog(VBR.in.SV.f, squeeze(invQ_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR_HTB.out.anelastic.xfit_premelt.Qinv(iT,id,isig,:)),'r-')
loglog(VBR.in.SV.f, squeeze(VBR_peak.out.anelastic.xfit_premelt.Qinv(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.Qinv(iT,id,isig,:)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Frequency (Hz)')
ylabel('Attenuation, Q^{-1}')
 
%% Map 1: invQ and Geff as a function of stress and grain size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find indices of conditions that you want to plot
T_target = 1400+273; %K, target temperature
F_target = 1/1e2; %Hz, target frequency
T_slice =  T(mindex(T, T_target)) - 273;
F_slice =  f(mindex(f, F_target));
iT = mindex(T,T_target);
iF = mindex(f,F_target);

% figure out which mechanism predicts the biggest 
% J2 (i.e., the dominant anelastic mechanism) 
% and biggest J1 (i.e., the dominant modulus
% relaxation mechanism) and create a grid of
% indices (1, 2, or 3 for the HTB, peak, and
% linearized backstress model, respectively) as a
% function of grain size and stress.
for i = 1:length(d)
    for j = 1:length(sig)
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J2(iT,i,j,iF) VBR_peak.out.anelastic.xfit_premelt.J2(iT,i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.J2(iT,i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J2 (i.e., is the dominant anelastic mechanism)
        Qmech(i,j) = ind;
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J1(iT,i,j,iF) VBR_peak.out.anelastic.xfit_premelt.J1(iT,i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.J1(iT,i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J1 (i.e., is the dominant relaxation mechanism)
        Mmech(i,j) = ind;
    end
end

% set up figure
fig = figure('Position',[200  600    1200   600]);
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 
clf
levels = [0 1 2 3]; %set Q/Mmech values that colored fields should fall between
invQcontours = 10.^[-5:1:0]; %
Geffcontours = [10 20 24 26 27 30 40 50 60 65 66 67 68 69 70 71];%GPa

% Dissipation mechanism map
ax1 = subplot(1,2,1);
contourf(ax1,sig,d./1000,Qmech,levels,'linestyle','none');hold on %Puts the colors on the map 
[C,h] = contour(ax1,sig,d./1000,squeeze(invQ_tot(iT,:,:,iF)),invQcontours,'k','linewidth',2,"LabelFormat","%0.1e"); %overlies contours of attenuation on the map
clim(ax1,[1 3])
clabel(C,h);

% aesthetics
axis square
title({'{\it Q}^{-1} and largest {\it J}_2 mech.,', ['{\it T} = ' sprintf('%.0f',T_slice) '^{\circ}C, {\it f} = ' sprintf('%.0e',F_slice) ' Hz']})
xlabel('Stress (MPa)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')

% Relaxation mechanism map
subplot(1,2,2)
ax2 = subplot(1,2,2);
contourf(ax2,sig,d./1000,Mmech,levels,'linestyle','none');hold on
[C,h] = contour(ax2,sig,d./1000,squeeze(G_eff_tot(iT,:,:,iF)./1e9),Geffcontours,'k','linewidth',2,"LabelFormat","%0.0f");
clim(ax2,[1 3])
clabel(C,h);

%aesthetics
axis square
title({'{\it G}_{eff} (GPa) and largest {\it J}_1 mech.', ['{\it T} = ' sprintf('%.0f',T_slice) '^{\circ}C, {\it f} = ' sprintf('%.0e',F_slice) ' Hz']})
xlabel('Stress (MPa)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
colormap(summer(3))
set(gcf,'color','w');


%% Slice through map at given grainsize
d_target = 10000; %um, target grain size
d_slice =  d(mindex(d, d_target));
id = mindex(d,d_target);

yline(ax1, d_slice./1000,'LineWidth',3)
yline(ax2, d_slice./1000,'LineWidth',3)

fig = figure('Position',[1500  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

% J1
subplot(4,1,1)
hold on
loglog(sig, squeeze(J1_tot(iT,id,:,iF)),'LineWidth',3,'Color','k')
loglog(sig, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J1(iT,id,:,iF)),'r-')
loglog(sig, squeeze(VBR_peak.out.anelastic.xfit_premelt.J1(iT,id,:,iF)),'r--')
loglog(sig, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J1(iT,id,:,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
title(['Slice at d = ' num2str(round(d_slice)) ' um, T= ' num2str(T_slice) '^oC, f=' num2str(F_slice) ' Hz'])
xlabel('Stress (MPa)')
ylabel('J_1')
legend('Total', 'Pre-melt, HTB', 'Premelt, peak', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(sig, squeeze(J2_tot(iT,id,:,iF)),'LineWidth',3,'Color','k')
loglog(sig, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J2(iT,id,:,iF)),'r-')
loglog(sig, squeeze(VBR_peak.out.anelastic.xfit_premelt.J2(iT,id,:,iF)),'r--')
loglog(sig, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J2(iT,id,:,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel('J_2')

% G_eff
subplot(4,1,3)
hold on
loglog(sig, squeeze(G_eff_tot(iT,id,:,iF))./1e9,'LineWidth',3,'Color','k')
loglog(sig, squeeze(VBR_HTB.out.anelastic.xfit_premelt.M(iT,id,:,iF))./1e9,'r-')
loglog(sig, squeeze(VBR_peak.out.anelastic.xfit_premelt.M(iT,id,:,iF))./1e9,'r--')
loglog(sig, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.M(iT,id,:,iF))./1e9,'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel({'Effective shear'; 'modulus (GPa)'})

% invQ
subplot(4,1,4)
hold on
loglog(sig, squeeze(invQ_tot(iT,id,:,iF)),'LineWidth',3,'Color','k')
loglog(sig, squeeze(VBR_HTB.out.anelastic.xfit_premelt.Qinv(iT,id,:,iF)),'r-')
loglog(sig, squeeze(VBR_peak.out.anelastic.xfit_premelt.Qinv(iT,id,:,iF)),'r--')
loglog(sig, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.Qinv(iT,id,:,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel('Q^{-1}')

%% Slice through map at given stress
sig_target = 1; %MPa, target stress
sig_slice =  sig(mindex(sig, sig_target));
isig = mindex(sig,sig_target);

xline(ax1, sig_slice,'LineWidth',3)
xline(ax2, sig_slice,'LineWidth',3)

fig = figure('Position',[750  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

% J1
subplot(4,1,1)
hold on
loglog(d, squeeze(J1_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J1(iT,:,isig,iF)),'r-')
loglog(d, squeeze(VBR_peak.out.anelastic.xfit_premelt.J1(iT,:,isig,iF)),'r--')
loglog(d, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J1(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
title(['Slice at \sigma = ' num2str(round(sig_slice)) ' MPa, T= ' num2str(T_slice) '^oC, f=' num2str(F_slice) ' Hz'])
xlabel('Stress (MPa)')
ylabel('J_1')
legend('Total', 'Pre-melt, HTB', 'Premelt, peak', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(d, squeeze(J2_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J2(iT,:,isig,iF)),'r-')
loglog(d, squeeze(VBR_peak.out.anelastic.xfit_premelt.J2(iT,:,isig,iF)),'r--')
loglog(d, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J2(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel('J_2')

% G_eff
subplot(4,1,3)
hold on
loglog(d, squeeze(G_eff_tot(iT,:,isig,iF))./1e9,'LineWidth',3,'Color','k')
loglog(d, squeeze(VBR_HTB.out.anelastic.xfit_premelt.M(iT,:,isig,iF))./1e9,'r-')
loglog(d, squeeze(VBR_peak.out.anelastic.xfit_premelt.M(iT,:,isig,iF))./1e9,'r--')
loglog(d, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.M(iT,:,isig,iF))./1e9,'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel({'Effective shear'; 'modulus (GPa)'})

% invQ
subplot(4,1,4)
hold on
loglog(d, squeeze(invQ_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d, squeeze(VBR_HTB.out.anelastic.xfit_premelt.Qinv(iT,:,isig,iF)),'r-')
loglog(d, squeeze(VBR_peak.out.anelastic.xfit_premelt.Qinv(iT,:,isig,iF)),'r--')
loglog(d, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.Qinv(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Stress (MPa)')
ylabel('Q^{-1}')


%% Map 2: invQ and Geff as a function of temperature and grain size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find indices of conditions that you want to plot
sig_target = 10; %MPa, target stress
F_target = 1/100; %Hz, target frequency
sig_slice =  sig(mindex(sig, sig_target));
F_slice =  f(mindex(f, F_target));
isig = mindex(sig,sig_target);
iF = mindex(f,F_target);

% as before, but now for grain size and
% temperature
for i = 1:length(T)
    for j = 1:length(d)
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J2(i,j,isig,iF) VBR_peak.out.anelastic.xfit_premelt.J2(i,j,isig,iF) VBR_linBackstress.out.anelastic.backstress_linear.J2(i,j,isig,iF)]); %find out premelt (ind=1) or backstress (ind=2) etc. mechanism predicts larger J2 (i.e., is the dominant anelastic mechanism)
        Qmech2(i,j) = ind;
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J1(i,j,isig,iF) VBR_peak.out.anelastic.xfit_premelt.J1(i,j,isig,iF) VBR_linBackstress.out.anelastic.backstress_linear.J1(i,j,isig,iF)]); %find out premelt (ind=1) or backstress (ind=2) etc. mechanism predicts larger J1 (i.e., is the dominant relaxation mechanism)
        Mmech2(i,j) = ind;
    end
end

fig = figure('Position',[200  0    1200   600]);
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 
clf
levels = [0 1 2 3]; %set Q/Mmech values that colored fields should fall between
invQcontours = 10.^[-5:1:0]; %
Geffcontours = [10 20 30 40 50 60 65 66 67 68 69 70 71];%GPa

% Dissipation mechanism map
ax1 = subplot(1,2,1);
contourf(ax1,T-273,d./1000,Qmech2',levels,'linestyle','none');hold on
[C,h] = contour(ax1,T-273,d./1000,squeeze(invQ_tot(: ,:,isig,iF))',invQcontours,'k','linewidth',2,"LabelFormat","%0.1e");
clim(ax1,[1 3])
clabel(C,h);

%aesthetics
axis square
title({'{\it Q}^{-1} and largest {\it J}_2 mech.,', ['{\it \sigma} = ' sprintf('%.0f',sig_slice) ' MPa, {\it f} = ' sprintf('%.0e',F_slice) ' Hz']})
xlabel('T (^{\circ}C)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')

% Relaxation mechanism map
subplot(1,2,2)
ax2 = subplot(1,2,2);
contourf(ax2,T-273,d./1000,Mmech2',levels,'linestyle','none');hold on
[C,h] = contour(ax2,T-273,d./1000,squeeze(G_eff_tot(:,:,isig,iF)./1e9)',Geffcontours,'k','linewidth',2,"LabelFormat","%0.0f");
clim(ax2,[1 3])
clabel(C,h);

%aesthetics
axis square
title({'{\it G}_{eff} (GPa) and largest {\it J}_1 mech.', ['{\it \sigma} = ' sprintf('%.0f',sig_slice) ' MPa, {\it f} = ' sprintf('%.0e',F_slice) ' Hz']})
xlabel('T (^{\circ}C)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
colormap(summer(3))
set(gcf,'color','w');

%% Slice through map at given grainsize
d_target = 1000; %um, target grain size
d_slice =  d(mindex(d, d_target));
id = mindex(d,d_target);

yline(ax1, d_slice./1000,'LineWidth',3)
yline(ax2, d_slice./1000,'LineWidth',3)

fig = figure('Position',[1500  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

% J1
subplot(4,1,1)
hold on
loglog(T-273, squeeze(J1_tot(:,id,isig,iF)),'LineWidth',3,'Color','k')
loglog(T-273, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J1(:,id,isig,iF)),'r-')
loglog(T-273, squeeze(VBR_peak.out.anelastic.xfit_premelt.J1(:,id,isig,iF)),'r--')
loglog(T-273, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J1(:,id,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
title(['Slice at d = ' num2str(round(d_slice./1000,3)) ' mm, \sigma = ' num2str(sig_slice) ' MPa, f =' num2str(F_slice) ' Hz'])
xlabel('Temperautre (^{\circ}C)')
ylabel('J_1')
legend('Total', 'Pre-melt, HTB', 'Premelt, peak', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(T-273, squeeze(J2_tot(:,id,isig,iF)),'LineWidth',3,'Color','k')
loglog(T-273, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J2(:,id,isig,iF)),'r-')
loglog(T-273, squeeze(VBR_peak.out.anelastic.xfit_premelt.J2(:,id,isig,iF)),'r--')
loglog(T-273, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J2(:,id,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Temperature (^{\circ}C)')
ylabel('J_2')

% G_eff
subplot(4,1,3)
hold on
loglog(T-273, squeeze(G_eff_tot(:,id,isig,iF))./1e9,'LineWidth',3,'Color','k')
loglog(T-273, squeeze(VBR_HTB.out.anelastic.xfit_premelt.M(:,id,isig,iF))./1e9,'r-')
loglog(T-273, squeeze(VBR_peak.out.anelastic.xfit_premelt.M(:,id,isig,iF))./1e9,'r--')
loglog(T-273, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.M(:,id,isig,iF))./1e9,'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Temperature (^{\circ}C)')
ylabel({'Effective shear'; 'modulus (GPa)'})

% invQ
subplot(4,1,4)
hold on
loglog(T-273, squeeze(invQ_tot(:,id,isig,iF)),'LineWidth',3,'Color','k')
loglog(T-273, squeeze(VBR_HTB.out.anelastic.xfit_premelt.Qinv(:,id,isig,iF)),'r-')
loglog(T-273, squeeze(VBR_peak.out.anelastic.xfit_premelt.Qinv(:,id,isig,iF)),'r--')
loglog(T-273, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.Qinv(:,id,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Temperature (^{\circ}C)')
ylabel('Q^{-1}')

%% Slice through map at given Temperature
T_target = 1200 + 273; %K, target temperature
T_slice =  T(mindex(T, T_target));
iT = mindex(T,T_target);

xline(ax1, T_slice-273,'LineWidth',3)
xline(ax2, T_slice-273,'LineWidth',3)

fig = figure('Position',[750  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

% J1
subplot(4,1,1)
hold on
loglog(d./1000, squeeze(J1_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d./1000, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J1(iT,:,isig,iF)),'r-')
loglog(d./1000, squeeze(VBR_peak.out.anelastic.xfit_premelt.J1(iT,:,isig,iF)),'r--')
loglog(d./1000, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J1(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
title(['Slice at \sigma = ' num2str(round(sig_slice)) ' MPa, T= ' num2str(T_slice) '^oC, f=' num2str(F_slice) ' Hz'])
xlabel('Grain size (mm)')
ylabel('J_1')
legend('Total', 'Pre-melt, HTB', 'Premelt, peak', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(d./1000, squeeze(J2_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d./1000, squeeze(VBR_HTB.out.anelastic.xfit_premelt.J2(iT,:,isig,iF)),'r-')
loglog(d./1000, squeeze(VBR_peak.out.anelastic.xfit_premelt.J2(iT,:,isig,iF)),'r--')
loglog(d./1000, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.J2(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Grain size (mm)')
ylabel('J_2')

% G_eff
subplot(4,1,3)
hold on
loglog(d./1000, squeeze(G_eff_tot(iT,:,isig,iF))./1e9,'LineWidth',3,'Color','k')
loglog(d./1000, squeeze(VBR_HTB.out.anelastic.xfit_premelt.M(iT,:,isig,iF))./1e9,'r-')
loglog(d./1000, squeeze(VBR_peak.out.anelastic.xfit_premelt.M(iT,:,isig,iF))./1e9,'r--')
loglog(d./1000, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.M(iT,:,isig,iF))./1e9,'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Grain size (mm)')
ylabel({'Effective shear'; 'modulus (GPa)'})

% invQ
subplot(4,1,4)
hold on
loglog(d./1000, squeeze(invQ_tot(iT,:,isig,iF)),'LineWidth',3,'Color','k')
loglog(d./1000, squeeze(VBR_HTB.out.anelastic.xfit_premelt.Qinv(iT,:,isig,iF)),'r-')
loglog(d./1000, squeeze(VBR_peak.out.anelastic.xfit_premelt.Qinv(iT,:,isig,iF)),'r--')
loglog(d./1000, squeeze(VBR_linBackstress.out.anelastic.backstress_linear.Qinv(iT,:,isig,iF)),'b-')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Grain size (mm)')
ylabel('Q^{-1}')


%% functions
function [ minX_ind ] = mindex( X,a )
% [ minX_ind ] = mindex( X,a )
%   simple function to return the index of the minimum point in vector X
%   
%   if a second argument is given, the function outputs the index of the
%   point in X closest to the water level, "a"
%   basically just outputs the second output of the "min" function,
%   without giving you the magnitude of the minimum value.
% 
%   N.B. can be used as a zero finder if a==0
%
%   Intended for use when calling the value in one vector corresponding to
%   the minimum value in X - i.e more efficient than the clunkier:
%       Y(find(X==min(X)) or, more often, Y(find((X-a)==min(X-a)))
%   Instead, can now use
%       Y(mindex(X)) or Y(mindex(X,a))
%   See also maxdex.m

if nargin<2
    [~,minX_ind] = min(X);
else
    [~,minX_ind] = min(abs(X-a));
end

end