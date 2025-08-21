% attenuation mechanism maps! 
% DH 7/2025

close all; clear all; clc
vbr_init

%% set methods list
VBR.in.elastic.methods_list = {'anharmonic'};
VBR.in.viscous.methods_list = {'HZK2011'};
VBR.in.anelastic.methods_list = {'xfit_premelt','backstress_linear'};
VBR.in.elastic.anharmonic = Params_Elastic('anharmonic'); 
VBR.in.elastic.anharmonic.temperature_scaling = 'isaak';
VBR.in.elastic.anharmonic.pressure_scaling = 'abramson';

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
%% run VBR
[VBR] = VBR_spine(VBR);

%% Combine complex compliances and calculate attenuation and effective modulus for all

% find unrelaxed modulus
Ju = (1/VBR.out.elastic.anharmonic.Gu); %GPa, unrelaxed shear compliance

% % summing complex compliances
% J_tot = VBR.out.anelastic.xfit_premelt.J1 + 1i.*VBR.out.anelastic.xfit_premelt.J2 + VBR.out.anelastic.backstress_linear.J1 + 1i.*VBR.out.anelastic.backstress_linear.J2 - Ju;
% J1_tot = real(J_tot);
% J2_tot = imag(J_tot);
% invQ_tot = J2_tot./J1_tot;

% or alternatively, sum J1 and J2s to get invQ and Geff (produces same results)
J1_tot = VBR.out.anelastic.xfit_premelt.J1 + VBR.out.anelastic.backstress_linear.J1 - Ju; %subtracting Ju once as it is incorporated in J1 of both the premelt and backstress models%
J2_tot = (VBR.out.anelastic.xfit_premelt.J2 + VBR.out.anelastic.backstress_linear.J2);
invQ_tot = J2_tot./J1_tot;
G_tot = 1./(J1_tot+1i.*J2_tot); 
G_eff_tot = abs(G_tot);

%% plotting a single spectrum at arbitrary conditions

%find indices of conditions that you want to plot
T_target   = 1300+273; %K, target temperature
d_target   = 1000; %um, target grain size
sig_target = 1;    %MPa, target stress

% print selected target conditions to screen
T_slice =     T(mindex(T, T_target)) - 273
d_slice =     d(mindex(d, d_target))
sig_slice = sig(mindex(sig, sig_target))

iT = mindex(T,T_target);
id = mindex(d,d_target);
isig = mindex(sig,sig_target);

fig = figure('Position',[1500  800    750   1000]);
set(gcf,'color','w');
set(fig,'DefaultAxesFontSize',15,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.03 .02],'DefaultAxesLineWidth',2)
fig.Renderer = 'painters'; %ensures that figure is rendered as vector graphics 

%J1
subplot(4,1,1)
hold on
loglog(VBR.in.SV.f, squeeze(J1_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.xfit_premelt.J1(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.backstress_linear.J1(iT,id,isig,:)),'b--')
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
legend('Total', 'Pre-melt', 'Linear backstress', 'Location','northeast')

% J2
subplot(4,1,2)
hold on
loglog(VBR.in.SV.f, squeeze(J2_tot(iT,id,isig,:)),'LineWidth',3,'Color','k')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.xfit_premelt.J2(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.backstress_linear.J2(iT,id,isig,:)),'b--')
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
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.xfit_premelt.M(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.backstress_linear.M(iT,id,isig,:)),'b--')
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
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.xfit_premelt.Qinv(iT,id,isig,:)),'r--')
loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.backstress_linear.Qinv(iT,id,isig,:)),'b--')
set(gca,"YScale",'log')
set(gca,"XScale",'log')
ax = gca;
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
box on %enables box
set(gca, 'layer','top') %plots the axes on top of the data so data don't overlap with axes

%annotate figure
xlabel('Frequency (Hz)')
ylabel('Q^{-1}')

 
%% plotting maps

colors = ['r','b']; %for now, red for pre-melt and blue for linear backstress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% map 1: stress, grain size at fixed T and f %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find indices of conditions that you want to plot
T_target = 1500+273; %K, target temperature
F_target = 1/100; %Hz, target frequency
T_slice =  T(mindex(T, T_target)) - 273;
F_slice =  f(mindex(f, F_target));
iT = mindex(T,T_target);
iF = mindex(f,F_target);

% figure out which mechanism predicts the biggest J2 (i.e., the dominant anelastic mechanism) 
% and biggest J1 (i.e., the dominant modulus relaxation mechanism) 
for i = 1:length(d)
    for j = 1:length(sig)
        [~,ind] = max([VBR.out.anelastic.xfit_premelt.J2(iT,i,j,iF) VBR.out.anelastic.backstress_linear.J2(iT,i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J2 (i.e., is the dominant anelastic mechanism)
        Qmech(i,j) = ind;
        [~,ind] = max([VBR.out.anelastic.xfit_premelt.J1(iT,i,j,iF) VBR.out.anelastic.backstress_linear.J1(iT,i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J1 (i.e., is the dominant relaxation mechanism)
        Mmech(i,j) = ind;
    end
end

fig = figure('Position',[200  800    800   400]);
clf
colormap(summer(2))

% Dissipation mechanism map
ax1 = subplot(1,2,1);
contourf(ax1,sig,d./1000,Qmech,[1:1.5:2],'linestyle','none');hold on %Puts the colors on the map 
[C,h] = contour(ax1,sig,d./1000,squeeze(invQ_tot(iT,:,:,iF)),10.^[-2:0.5:0],'k','linewidth',2,"LabelFormat","%0.2f"); %overlies contours of attenuation on the map
clim(ax1,[1 2])
clabel(C,h);

axis square
title('Q^{-1} contours and mechanism with largest J_2')
% title(['T= ' num2str(T_slice) '^oC, period=' num2str(round(1/F_slice)) ' s'])
xlabel('Stress (MPa)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
text(0.03,0.92, ['{\it T} = ' sprintf('%.0f',T_slice) ' ^{\circ}C'],'Units','normalized','FontSize',15) %annotates variable T at a position of 10% from the left vertical axis, and 90% from the bottom horizontal axis
text(0.03,0.80, ['{\it f} = ' sprintf('%.0e',F_slice) ' Hz'],'Units','normalized','FontSize',15) %annotates variable T at a position of 10% from the left vertical axis, and 90% from the bottom horizontal axis
colormap(summer(2))

% Relaxation mechanism map
subplot(1,2,2)
ax2 = subplot(1,2,2);
contourf(ax2,sig,d./1000,Mmech,[1:1.5:2],'linestyle','none');hold on
[C,h] = contour(ax2,sig,d./1000,squeeze(G_eff_tot(iT,:,:,iF)./1e9),10.^[-2:0.1:2],'k','linewidth',2,"LabelFormat","%0.0f");
clim(ax2,[1 2])
clabel(C,h);
axis square
title('G_{eff} contours and mechanism with largest J_1')
xlabel('Stress (MPa)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
colormap(summer(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% map 2: Temperature and grain size at fixed stress and f %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find indices of conditions that you want to plot
sig_target = 10; %K, target temperature
F_target = 1/100; %Hz, target frequency
sig_slice =  sig(mindex(sig, sig_target));
F_slice =  f(mindex(f, F_target));
iT = mindex(T,T_target);
iSig = mindex(sig,sig_target);

% figure out which mechanism predicts the biggest J2 (i.e., the dominant anelastic mechanism) 
% and biggest J1 (i.e., the dominant modulus relaxation mechanism) 
for i = 1:length(T)
    for j = 1:length(d)
        [~,ind] = max([VBR.out.anelastic.xfit_premelt.J2(i,j,iSig,iF) VBR.out.anelastic.backstress_linear.J2(i,j,iSig,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J2 (i.e., is the dominant anelastic mechanism)
        Qmech2(i,j) = ind;
        [~,ind] = max([VBR.out.anelastic.xfit_premelt.J1(i,j,iSig,iF) VBR.out.anelastic.backstress_linear.J1(i,j,iSig,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J1 (i.e., is the dominant relaxation mechanism)
        Mmech2(i,j) = ind;
    end
end

fig = figure('Position',[335  100    780   400]);
clf
colormap(summer(2))

% Dissipation mechanism map
ax1 = subplot(1,2,1);
contourf(ax1,T,d./1000,Qmech2,[1:1.5:2],'linestyle','none');hold on
[C,h] = contour(ax1,T  ,d./1000,squeeze(invQ_tot(: ,:,iSig,iF)),10.^[-2:0.5:0],'k','linewidth',2,"LabelFormat","%0.2f");
clim(ax1,[1 2])
clabel(C,h);

axis square
title('Q^{-1} contours and mechanism with largest J_2')
xlabel('T (K)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
text(0.03,0.92, ['{\it \sigma} = ' sprintf('%.0f',sig_slice) ' MPa'],'Units','normalized','FontSize',15) %annotates variable T at a position of 10% from the left vertical axis, and 90% from the bottom horizontal axis
text(0.03,0.80, ['{\it f} = ' sprintf('%.0e',F_slice) ' Hz'],'Units','normalized','FontSize',15) %annotates variable T at a position of 10% from the left vertical axis, and 90% from the bottom horizontal axis
colormap(summer(2))

% Relaxation mechanism map
subplot(1,2,2)
ax2 = subplot(1,2,2);
contourf(ax2,T,d./1000,Mmech2,[1:1.5:2],'linestyle','none');hold on
[C,h] = contour(ax2,T,d./1000,squeeze(G_eff_tot(:,:,iSig,iF)./1e9),10.^[-2:0.1:2],'k','linewidth',2,"LabelFormat","%0.0f");
clim(ax2,[1 2])
clabel(C,h);
axis square
title('G_{eff} contours and mechanism with largest J_1')
xlabel('T (K)')
ylabel('Grain size (mm)')
set(gca,'XScale','log')
set(gca,'YScale','log')
colormap(summer(2))

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