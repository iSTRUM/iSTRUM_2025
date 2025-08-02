function VBR = backstress_shear_mod()

    addpath(getenv('vbrdir'))
    vbr_init();

    VBR.in.anelastic.methods_list = {'backstress_linear'};
    VBR.in.elastic.methods_list = {'anharmonic'};
    VBR.in.elastic.anharmonic = Params_Elastic('anharmonic');
    VBR.in.elastic.anharmonic.temperature_scaling = 'isaak';
    VBR.in.elastic.anharmonic.pressure_scaling = 'abramson';

    % set state variables
    VBR.in.SV.T_K = [1300, 1400, 1500] + 273;
    sz = size(VBR.in.SV.T_K);
    VBR.in.SV.sig_MPa = full_nd(3., sz);
    VBR.in.SV.dg_um = full_nd(0.001 * 1e6, sz);

    % following are needed for anharmonic calculation
    VBR.in.SV.P_GPa = full_nd(5., sz);
    VBR.in.SV.rho = full_nd(3300, sz);
    VBR.in.SV.f = logspace(-8, 0, 500);%[0.001, 0.01];

    % calculations
    VBR = VBR_spine(VBR);

    % plotting
    Qinv = VBR.out.anelastic.backstress_linear.Qinv;
    Vs = VBR.out.anelastic.backstress_linear.V / 1e3;

    fsize = 18;
    figure
    subplot(4,1,1)
    loglog(VBR.in.SV.f, VBR.out.anelastic.backstress_linear.Qinv(1,1,:), 'k','displayname', 'Qinv from G')
    hold on
    loglog(VBR.in.SV.f, VBR.out.anelastic.backstress_linear.Qinv_E(1,1,:), '--r','displayname', 'Qinv from E')
    legend()
    ylabel('Qinv', 'fontsize',fsize)
    xlabel('f [Hz]', 'fontsize',fsize)

    subplot(4,1,2)
    rat = VBR.out.anelastic.backstress_linear.Qinv(1,1,:) ./ VBR.out.anelastic.backstress_linear.Qinv_E(1,1,:);
    semilogx(VBR.in.SV.f,rat, 'k')
    ylabel('Qinv G / Qinv E', 'fontsize',fsize)
    xlabel('f [Hz]', 'fontsize',fsize)

    subplot(4,1,3)
    semilogx(VBR.in.SV.f, VBR.out.anelastic.backstress_linear.M(1,1,:)/1e9, 'k','displayname', 'Qinv from G')
    hold on
    semilogx(VBR.in.SV.f, VBR.out.anelastic.backstress_linear.E(1,1,:)/1e9, '--r','displayname', 'Qinv from E')
    ylabel('M [GPa]', 'fontsize',fsize)
    xlabel('f [Hz]', 'fontsize', fsize)


    subplot(4,1,4)
    rat = VBR.out.anelastic.backstress_linear.M ./ VBR.out.anelastic.backstress_linear.E;
    plot(VBR.in.SV.f, rat(1,1,:), 'k','displayname', 'Qinv from G')
    ylabel('M G / M E', 'fontsize',fsize)
    xlabel('f [Hz]', 'fontsize',fsize)

    for iax = 1:4
        subplot(4,1,iax)
        set(gca, 'fontsize', fsize)
    end

end