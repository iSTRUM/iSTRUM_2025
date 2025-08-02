addpath(getenv('vbrdir'))
    vbr_init();

%%  write method list %%
  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.anelastic.methods_list={'eburgers_psp';};


  VBR.in.anelastic.eburgers_psp=Params_Anelastic('eburgers_psp');

  % use the single sample background only fit:
  VBR.in.anelastic.eburgers_psp.eBurgerMethod='bg_peak'; % 'bg_only' or 'bg_peak' or 's6585_bg_only'
  VBR.in.anelastic.eburgers_psp.bg_peak.DeltaB = 0;

  % frequencies to calculate at
  VBR.in.SV.f = 1./logspace(-2,4,100);

  %% Define the Thermodynamic State %%
  VBR.in.SV.T_K=700:50:1200;
  VBR.in.SV.T_K=VBR.in.SV.T_K+273;
  sz=size(VBR.in.SV.T_K); % temperature [K]

  % remaining state variables
  VBR.in.SV.dg_um=3.1*ones(sz);
  VBR.in.SV.P_GPa = 0.2 * ones(sz); % pressure [GPa]
  VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
  VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
  VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction

  %% call VBR_spine %%
  [VBR] = VBR_spine(VBR) ;
