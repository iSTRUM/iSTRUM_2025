clear; clc; close all;

flist = dir('LON_LAT_Q_DEPTH_VS*.mat');
Depth_Range = [100 200]
for ijk = 1:length(flist)

fname = flist(ijk).name;
load(fname)
vs = zzz(:,5);
Q = zzz(:,3);
idx = find(zzz(:,4) > min(Depth_Range) & zzz(:,4) < max(Depth_Range));
scatter(vs,Q,50,'filled','MarkerFaceAlpha',0.5);
hold on
Attribution = extractBetween(fname,'LON_LAT_Q_DEPTH_VS_','.mat');
legentry{ijk} = Attribution{1};
end
legend(legentry,'Location','northwest')
ylim([0 600])
set(gca,'fontsize',18,'fontweight','bold')
grid on; box on;
xlabel('Vs(km/s)')
ylabel('Q')
title(['Ensemble of Models at Depth Between ' num2str(min(Depth_Range)) '-' num2str(max(Depth_Range)) 'km' ])
