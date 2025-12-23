function out = hdo_reactor_v4p1_pH2coupled(feed, op, sel, h2, thermo, comp)
% HDO reactor model v4.2 (pH2 coupled + gas mixture props + optional MT + optional T-profile)
% New in v4.2:
%   - op.mt.enabled: intraparticle diffusion effectiveness factor eta(Thiele)
%   - op.rx_mode: "isothermal"(default) or "adiabatic" (simple energy-coupled along bed)
%   - op.Nseg: axial segments for tau-integration when mt or adiabatic is enabled
%
% Backward-compatible: if op.mt.enabled=false AND op.rx_mode="isothermal", results match v4.1.

% ------------------------------
% Defaults / setup
% ------------------------------
out.warnings = {};

if ~isfield(op,'T_feed_K'); op.T_feed_K = 298.15; end
if ~isfield(op,'Ptot_MPa') && isfield(op,'p_MPa')
    op.Ptot_MPa = op.p_MPa;
    out.warnings{end+1} = "op.Ptot_MPa not found; using op.p_MPa as total pressure Ptot.";
end

% --- NEW defaults (v4.2) ---
if ~isfield(op,'rx_mode'); op.rx_mode = "isothermal"; end          % "isothermal" | "adiabatic"
if ~isfield(op,'Nseg');    op.Nseg = 200; end                     % segments for tau integration
if ~isfield(op,'mt');      op.mt = struct(); end                  % mass transfer (intraparticle)
if ~isfield(op.mt,'enabled');  op.mt.enabled = false; end
if ~isfield(op.mt,'Rp_m');      op.mt.Rp_m = 0.7e-3; end           % pellet radius [m], guess
if ~isfield(op.mt,'Deff_m2s');  op.mt.Deff_m2s = 1.0e-9; end       % effective diffusivity [m2/s], guess
% (ADD) eta override defaults
if ~isfield(op.mt,'override');  op.mt.override = false; end
if ~isfield(op.mt,'eta_fixed'); op.mt.eta_fixed = 1.0; end

rx_mode = string(op.rx_mode);
Nseg    = max(10, round(op.Nseg));

if nargin < 4 || isempty(h2); h2 = struct(); end
if ~isfield(h2,'lambda_H2');     h2.lambda_H2 = 20; end
if ~isfield(h2,'purge_frac');    h2.purge_frac = 0.02; end
if ~isfield(h2,'preheat_basis'); h2.preheat_basis = "makeup"; end

if nargin < 5 || isempty(thermo); thermo = struct(); end
if ~isfield(thermo,'cp_oil_kJ_per_kgK'); thermo.cp_oil_kJ_per_kgK = 2.0; end

% Default gas cp (kJ/kg-K) and k (-)
if ~isfield(thermo,'cp_H2_kJ_per_kgK');   thermo.cp_H2_kJ_per_kgK = 14.3; end
if ~isfield(thermo,'cp_CO2_kJ_per_kgK');  thermo.cp_CO2_kJ_per_kgK = 0.90; end
if ~isfield(thermo,'cp_CO_kJ_per_kgK');   thermo.cp_CO_kJ_per_kgK  = 1.04; end
if ~isfield(thermo,'cp_C3H8_kJ_per_kgK'); thermo.cp_C3H8_kJ_per_kgK = 1.67; end

if ~isfield(thermo,'k_H2');   thermo.k_H2   = 1.41; end
if ~isfield(thermo,'k_CO2');  thermo.k_CO2  = 1.30; end
if ~isfield(thermo,'k_CO');   thermo.k_CO   = 1.40; end
if ~isfield(thermo,'k_C3H8'); thermo.k_C3H8 = 1.13; end

dH_fields = {'dH1_kJ_per_molTAG','dH2_kJ_per_molTAG','dH3_kJ_per_molTAG'};
has_all_dH = all(isfield(thermo, dH_fields));
if ~isfield(thermo,'enabled')
    thermo.enabled = has_all_dH;
end
if thermo.enabled && ~has_all_dH
    error("thermo.enabled=true requires thermo.dH1/dH2/dH3 (kJ/mol TAG).");
end

if nargin < 6 || isempty(comp); comp = struct(); end
if ~isfield(comp,'enabled'); comp.enabled = false; end
if ~isfield(comp,'eta');     comp.eta = 0.75; end
if ~isfield(comp,'T_in_K');  comp.T_in_K = 313.15; end

% ------------------------------
% Constants
% ------------------------------
R = 8.314; % J/mol-K

% Molecular weights [kg/mol]
MW.TAG  = 885.45e-3;
MW.H2   = 2.016e-3;
MW.CO   = 28.01e-3;
MW.CO2  = 44.01e-3;
MW.C3H8 = 44.10e-3;
MW.H2O  = 18.015e-3;
MW.C18  = 254.50e-3;
MW.C17  = 240.47e-3;

% ------------------------------
% Input checks
% ------------------------------

if ~isfield(feed,'mdot_oil_kgph'); error("feed.mdot_oil_kgph is required."); end
if ~isfield(op,'T_K'); error("op.T_K is required."); end
if ~isfield(op,'Ptot_MPa'); error("op.Ptot_MPa (total pressure) is required."); end
if ~isfield(op,'LHSV_per_h'); error("op.LHSV_per_h is required."); end

if ~isfield(sel,'s1_HDO') || ~isfield(sel,'s2_DCO2') || ~isfield(sel,'s3_DCO')
    error("sel.s1_HDO, sel.s2_DCO2, sel.s3_DCO are required.");
end

s1 = sel.s1_HDO; s2 = sel.s2_DCO2; s3 = sel.s3_DCO;
if abs((s1+s2+s3)-1) > 1e-6
    error("Selectivities must sum to 1. Current sum = %.6f", s1+s2+s3);
end

f = h2.purge_frac;
if f <= 0
    f = 1e-4;
    out.warnings{end+1} = "purge_frac was <=0; set to 1e-4 to avoid inert accumulation singularity.";
end
if f >= 1
    error("purge_frac must be < 1.");
end
if h2.lambda_H2 < 1
    error("lambda_H2 must be >= 1.");
end

if rx_mode ~= "isothermal" && ~thermo.enabled
    out.warnings{end+1} = "rx_mode is not isothermal but thermo is disabled -> forcing rx_mode='isothermal' (no T-profile).";
    rx_mode = "isothermal";
end

% ------------------------------
% Feed molar flow
% ------------------------------
out.mdot_oil_kgph = feed.mdot_oil_kgph;
out.nTAG_in_molph = out.mdot_oil_kgph / MW.TAG;

% ------------------------------
% Space time from LHSV
% ------------------------------
out.LHSV_per_h = op.LHSV_per_h;
out.tau_h = 1 / op.LHSV_per_h;

% ------------------------------
% Effective H2 stoichiometry per TAG reacted
% ------------------------------
out.nuH2_eff = 15*s1 + 6*s2 + 9*s3; % mol H2 / mol TAG reacted

% ------------------------------
% Gas production per TAG reacted (recycled gas set)
% ------------------------------
alpha_CO2  = 3*s2;
alpha_CO   = 3*s3;
alpha_C3H8 = 1;
alpha_inert = alpha_CO2 + alpha_CO + alpha_C3H8;

% ------------------------------
% Coupled recycle composition -> yH2 -> pH2
% ------------------------------
out.lambda_H2 = h2.lambda_H2;
out.purge_frac = f;

den = out.lambda_H2*out.nuH2_eff + ((1-f)/f)*alpha_inert;
out.yH2_in = (out.lambda_H2*out.nuH2_eff) / den;

out.Ptot_MPa = op.Ptot_MPa;
out.pH2_MPa = out.yH2_in * out.Ptot_MPa;

if out.yH2_in < 0.2
    out.warnings{end+1} = sprintf("Low yH2_in=%.3f -> pH2 may be too low. Increase purge or lambda, or add gas purification.", out.yH2_in);
end

% ------------------------------
% Precompute recycle gas composition (independent of extent)
% ------------------------------
yH2   = out.yH2_in;
yCO2  = (((1-f)/f)*alpha_CO2)  / den;
yCO   = (((1-f)/f)*alpha_CO)   / den;
yC3H8 = (((1-f)/f)*alpha_C3H8) / den;

yvec = [yH2, yCO2, yCO, yC3H8];
MWvec = [MW.H2, MW.CO2, MW.CO, MW.C3H8];
cpvec = [thermo.cp_H2_kJ_per_kgK, thermo.cp_CO2_kJ_per_kgK, thermo.cp_CO_kJ_per_kgK, thermo.cp_C3H8_kJ_per_kgK];
kvec  = [thermo.k_H2, thermo.k_CO2, thermo.k_CO, thermo.k_C3H8];

MWmix = sum(yvec .* MWvec);                    % kg/mol
wvec  = (yvec .* MWvec) / max(MWmix, 1e-30);   % mass fractions
cp_mix = sum(wvec .* cpvec);
k_mix  = sum(yvec .* kvec);

% ------------------------------
% Kinetics (base) + optional MT + optional T-profile
% k(T,pH2) = 9.2*exp(-17.8e3/(R*T))*(pH2)^0.77 [1/h]
% ------------------------------
k_kin0_per_h = 9.2 * exp(-(17.8e3)/(R*op.T_K)) * (out.pH2_MPa^0.77);

% intraparticle effectiveness at inlet (for reference)
eta0 = 1.0;
if op.mt.enabled
    if op.mt.override
        eta0 = op.mt.eta_fixed;
    else
        k_s0 = k_kin0_per_h / 3600; % 1/s  (FIX: k_kin0_per_h 사용)
        phi0 = op.mt.Rp_m * sqrt(max(k_s0,1e-30)/max(op.mt.Deff_m2s,1e-30));
        if phi0 > 1e-8
            eta0 = (3/(phi0^2))*(phi0*coth(phi0) - 1);
        end
    end
    eta0 = max(0.05, min(1.0, eta0));
end

out.k_per_h        = k_kin0_per_h;            % backward compat
out.eta0           = eta0;
out.k_eff0_per_h   = eta0 * k_kin0_per_h;
out.X_kin          = 1 - exp(-out.k_eff0_per_h * out.tau_h);

% If no MT and isothermal, keep exact closed-form (v4.1 behavior)
if (rx_mode == "isothermal") && (~op.mt.enabled)
    out.X_final = min([out.X_kin, 1]);
    out.T_out_K = op.T_K;
    out.eta_avg = 1.0;
else
    % tau-integration with optional T update
    dtau = out.tau_h / Nseg;

    % effective heat of reaction per mol TAG reacted (kJ/mol, negative exo)
    dH_eff_kJ_per_molTAG = NaN;
    if thermo.enabled
        dH_eff_kJ_per_molTAG = s1*thermo.dH1_kJ_per_molTAG + s2*thermo.dH2_kJ_per_molTAG + s3*thermo.dH3_kJ_per_molTAG;
    end

    % adiabatic coupling requires a heat-capacity flow Cdot (kJ/h-K)
    % We make it mildly self-consistent by iterating Cdot using a conversion guess.
    X_guess = out.X_kin;
    nTAG_in = out.nTAG_in_molph;

    nGas_in_per_molTAGrxn = out.lambda_H2*out.nuH2_eff + ((1-f)/f)*alpha_inert; % mol gas in / mol TAG reacted

    n_iter = 1;
    if rx_mode == "adiabatic"; n_iter = 2; end

    X_new = X_guess;
    T_out = op.T_K;
    eta_sum = 0;

    for it = 1:n_iter
        % estimate gas flow for Cp damping
        nTAG_rxn_guess = X_guess * nTAG_in;
        nGas_in_guess  = nGas_in_per_molTAGrxn * nTAG_rxn_guess; % mol/h
        mGas_in_guess_kgph = nGas_in_guess * MWmix;              % kg/h

        Cdot_kJphK = out.mdot_oil_kgph*thermo.cp_oil_kJ_per_kgK + mGas_in_guess_kgph*cp_mix;
        Cdot_kJphK = max(Cdot_kJphK, 1e-30);

        % integrate
        X = 0.0;
        T = op.T_K;
        eta_sum = 0.0;

        for i = 1:Nseg
            k_kin = 9.2 * exp(-(17.8e3)/(R*T)) * (out.pH2_MPa^0.77);

            eta = 1.0;
            if op.mt.enabled
                if op.mt.override
                    eta = op.mt.eta_fixed;
                else
                    k_s = k_kin / 3600; % 1/s
                    phi = op.mt.Rp_m * sqrt(max(k_s,1e-30)/max(op.mt.Deff_m2s,1e-30));
                    if phi > 1e-8
                        eta = (3/(phi^2))*(phi*coth(phi) - 1);
                    end
                end
                eta = max(0.05, min(1.0, eta));
            end


            k_eff = eta * k_kin;

            dX = k_eff * (1 - X) * dtau;
            dX = max(0, min(1-X, dX));
            X  = X + dX;

            eta_sum = eta_sum + eta;

            if rx_mode == "adiabatic"
                % heat released in this slice (kJ/h): (-dH_eff)*d(nTAG_rxn)
                % d(nTAG_rxn) = nTAG_in * dX (mol/h)
                Qslice_kJph = (-dH_eff_kJ_per_molTAG) * (nTAG_in * dX);
                dT = Qslice_kJph / Cdot_kJphK;
                T = T + dT;
            end
        end

        X_new = min(max(X,0),1);
        T_out = T;

        X_guess = X_new;
    end

    out.X_final = X_new;
    out.T_out_K = T_out;
    out.eta_avg = eta_sum / Nseg;
end

% ------------------------------
% Extent and liquid products
% ------------------------------
out.nTAG_rxn_molph = out.X_final * out.nTAG_in_molph;
out.nTAG_out_molph = out.nTAG_in_molph - out.nTAG_rxn_molph;

out.nH2_cons_molph = out.nuH2_eff * out.nTAG_rxn_molph;

% gas produced
out.nCO2_prod_molph  = alpha_CO2  * out.nTAG_rxn_molph;
out.nCO_prod_molph   = alpha_CO   * out.nTAG_rxn_molph;
out.nC3H8_prod_molph = alpha_C3H8 * out.nTAG_rxn_molph;

% inlet inerts (steady-state accumulation with purge only)
out.nCO2_in_molph  = ((1-f)/f) * out.nCO2_prod_molph;
out.nCO_in_molph   = ((1-f)/f) * out.nCO_prod_molph;
out.nC3H8_in_molph = ((1-f)/f) * out.nC3H8_prod_molph;
out.nInert_in_molph = out.nCO2_in_molph + out.nCO_in_molph + out.nC3H8_in_molph;

% H2 inlet by lambda
out.nH2_in_molph = out.lambda_H2 * out.nH2_cons_molph;

% outlet gas
out.nH2_out_molph   = out.nH2_in_molph - out.nH2_cons_molph;
out.nCO2_out_molph  = out.nCO2_in_molph  + out.nCO2_prod_molph;
out.nCO_out_molph   = out.nCO_in_molph   + out.nCO_prod_molph;
out.nC3H8_out_molph = out.nC3H8_in_molph + out.nC3H8_prod_molph;

out.nGas_in_molph  = out.nH2_in_molph  + out.nInert_in_molph;
out.nGas_out_molph = out.nH2_out_molph + out.nCO2_out_molph + out.nCO_out_molph + out.nC3H8_out_molph;

% purge/recycle component-wise
out.nH2_purge_molph   = f * out.nH2_out_molph;
out.nCO2_purge_molph  = f * out.nCO2_out_molph;
out.nCO_purge_molph   = f * out.nCO_out_molph;
out.nC3H8_purge_molph = f * out.nC3H8_out_molph;

out.nH2_recycle_molph   = (1-f) * out.nH2_out_molph;
out.nCO2_recycle_molph  = (1-f) * out.nCO2_out_molph;
out.nCO_recycle_molph   = (1-f) * out.nCO_out_molph;
out.nC3H8_recycle_molph = (1-f) * out.nC3H8_out_molph;

% H2 makeup to close loop
out.nH2_makeup_molph = out.nH2_cons_molph + out.nH2_purge_molph;

out.yH2_out = out.nH2_out_molph / max(out.nGas_out_molph, 1e-30);
out.recycle_to_makeup = out.nH2_recycle_molph / max(out.nH2_makeup_molph, 1e-30);

% hydrocarbon products (proxy)
out.nC18_out_molph = 3*s1 * out.nTAG_rxn_molph;
out.nC17_out_molph = 3*(s2+s3) * out.nTAG_rxn_molph;
out.nH2O_out_molph = (6*s1 + 3*s3) * out.nTAG_rxn_molph;

% ------------------------------
% Mass flows (kg/h)
% ------------------------------
out.mTAG_out_kgph = out.nTAG_out_molph * MW.TAG;

out.mH2_in_kgph      = out.nH2_in_molph * MW.H2;
out.mH2_makeup_kgph  = out.nH2_makeup_molph * MW.H2;
out.mH2_cons_kgph    = out.nH2_cons_molph * MW.H2;

out.mCO2_in_kgph     = out.nCO2_in_molph * MW.CO2;
out.mCO_in_kgph      = out.nCO_in_molph  * MW.CO;
out.mC3H8_in_kgph    = out.nC3H8_in_molph * MW.C3H8;

out.mGas_in_kgph = out.mH2_in_kgph + out.mCO2_in_kgph + out.mCO_in_kgph + out.mC3H8_in_kgph;

out.mC18_out_kgph = out.nC18_out_molph * MW.C18;
out.mC17_out_kgph = out.nC17_out_molph * MW.C17;
out.m_jetproxy_kgph = out.mC18_out_kgph + out.mC17_out_kgph;

out.mH2O_out_kgph = out.nH2O_out_molph * MW.H2O;

% ------------------------------
% Gas mixture properties (as before; now using precomputed composition too)
% ------------------------------
out.gas_yH2   = yH2;
out.gas_yCO2  = yCO2;
out.gas_yCO   = yCO;
out.gas_yC3H8 = yC3H8;

out.gas_wH2   = wvec(1);
out.gas_wCO2  = wvec(2);
out.gas_wCO   = wvec(3);
out.gas_wC3H8 = wvec(4);

out.cp_gas_mix_kJ_per_kgK = cp_mix;
out.k_gas_mix = k_mix;

% ------------------------------
% ENERGY (Level-1)
% ------------------------------
out.energy_enabled = thermo.enabled;
out.Q_rxn_kW = NaN;
out.Q_remove_kW = NaN;
out.Q_preheat_kW = NaN;
out.deltaT_ad_K = NaN;
out.deltaT_profile_K = NaN;
out.W_comp_kW = 0;

if thermo.enabled
    % reaction heat
    xi1 = s1 * out.nTAG_rxn_molph;
    xi2 = s2 * out.nTAG_rxn_molph;
    xi3 = s3 * out.nTAG_rxn_molph;

    Qrxn_kJph = xi1*thermo.dH1_kJ_per_molTAG + xi2*thermo.dH2_kJ_per_molTAG + xi3*thermo.dH3_kJ_per_molTAG;
    out.Q_rxn_kW = -Qrxn_kJph/3600;   % (+) released

    if rx_mode == "isothermal"
        out.Q_remove_kW = out.Q_rxn_kW;  % isothermal duty
    else
        out.Q_remove_kW = 0;            % adiabatic model
    end

    dT = op.T_K - op.T_feed_K;

    % Preheat basis
    if h2.preheat_basis == "total_gas"
        mGas_preheat = out.mGas_in_kgph;
        cpGas = out.cp_gas_mix_kJ_per_kgK;
    else
        mGas_preheat = out.mH2_makeup_kgph;
        cpGas = thermo.cp_H2_kJ_per_kgK;
    end

    Qpreheat_kJph = out.mdot_oil_kgph*thermo.cp_oil_kJ_per_kgK*dT + mGas_preheat*cpGas*dT;
    out.Q_preheat_kW = Qpreheat_kJph/3600;

    % Adiabatic deltaT estimate (same as v4.1)
    Cdot_kJphK = out.mdot_oil_kgph*thermo.cp_oil_kJ_per_kgK + out.mGas_in_kgph*out.cp_gas_mix_kJ_per_kgK;
    out.deltaT_ad_K = out.Q_rxn_kW*3600 / max(Cdot_kJphK, 1e-30);

    out.deltaT_profile_K = out.T_out_K - op.T_K;
else
    out.warnings{end+1} = "Energy disabled: provide thermo.dH1/dH2/dH3 to compute Q_rxn/Q_preheat.";
end

% ------------------------------
% Compressor power (unchanged)
% ------------------------------
if comp.enabled
    if ~isfield(comp,'P_in_bar') || ~isfield(comp,'P_out_bar')
        error("comp.enabled=true requires comp.P_in_bar and comp.P_out_bar.");
    end
    if comp.P_out_bar <= comp.P_in_bar
        out.warnings{end+1} = "Compressor PR <= 1: W_comp set to 0. Check pressures.";
        out.W_comp_kW = 0;
    else
        mdot_gas_kgps = out.mGas_in_kgph / 3600; % kg/s
        cp_JpkgK = out.cp_gas_mix_kJ_per_kgK * 1000;

        if isfield(comp,'k_override')
            k_use = comp.k_override;
        else
            k_use = out.k_gas_mix;
        end

        PR = comp.P_out_bar / comp.P_in_bar;

        out.W_comp_kW = (mdot_gas_kgps*cp_JpkgK*comp.T_in_K/comp.eta) * (PR^((k_use-1)/k_use) - 1) / 1000;
        out.k_comp_used = k_use;
    end
else
    out.W_comp_kW = 0;
end

end
