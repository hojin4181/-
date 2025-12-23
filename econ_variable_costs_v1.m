function E = econ_variable_costs_v1(feed, op, outHDO, SUM, P, varargin)
% Variable OPEX-only economics (raw materials + energy only)
% Ignores: CAPEX, separation costs, catalyst costs, labor (per assignment)

% ===== opts optional parsing =====
if nargin >= 6 && ~isempty(varargin{1})
    opts = varargin{1};
else
    opts = struct();
end

if ~isfield(opts,'heat_source');      opts.heat_source = "ng"; end           % "ng" | "electric"
if ~isfield(opts,'cooling_model');    opts.cooling_model = "chiller"; end    % "chiller" | "cw" | "none"
if ~isfield(opts,'include_revenue');  opts.include_revenue = true; end
if ~isfield(opts,'SAF_premium_mult'); opts.SAF_premium_mult = 1.0; end

heat_source      = string(opts.heat_source);
cooling_model    = string(opts.cooling_model);
include_revenue  = logical(opts.include_revenue);
SAF_premium_mult = double(opts.SAF_premium_mult);

E = struct();

% --------- Base flows (per hour) ---------
m_oil_in = feed.mdot_oil_kgph;
m_jet    = SUM.m_jet_out_HI_kgph;
mH2      = SUM.H2_makeup_total_kgph;

W_comp_kW = SUM.W_comp_total_kW;
Q_heat_kW = SUM.Q_heat_total_kW;
Q_cool_kW = SUM.Q_cool_total_kW;

% --------- IMPORTANT: isothermal reactor cooling duty ---------
rx_mode = "isothermal";
if isfield(op,'rx_mode')
    rx_mode = string(op.rx_mode);
end

% If HDO is isothermal, exothermic heat must be removed by utilities.
if rx_mode == "isothermal"
    if isfield(outHDO,'Q_rxn_kW')
        Q_cool_kW = Q_cool_kW + max(0, outHDO.Q_rxn_kW);
    end
end

% --------- Costs (KRW/h) ---------
E.C_oil_krwph = m_oil_in * P.soy_oil_krw_per_kg;
E.C_H2_krwph  = mH2      * P.H2_krw_per_kg;

% Electricity (compressor)
E.C_elec_comp_krwph = W_comp_kW * P.elec_krw_per_kWh;

% Heating utility
switch heat_source
    case "electric"
        E.C_heat_krwph = (Q_heat_kW / P.eheater_eta) * P.elec_krw_per_kWh;
        E.heat_basis = "electric";
    otherwise % "ng"
        E.C_heat_krwph = (Q_heat_kW / P.boiler_eta) * P.ng_krw_per_kWh_th;
        E.heat_basis = "ng_boiler";
end

% Cooling utility
switch cooling_model
    case "none"
        E.C_cool_krwph = 0;
    case "cw"
        % Cooling water: energy-only basis -> 0 (or add pump power if modeled)
        E.C_cool_krwph = 0;
    otherwise % "chiller"
        E.C_cool_krwph = (Q_cool_kW / P.COP_cool) * P.elec_krw_per_kWh;
end

E.C_total_krwph = E.C_oil_krwph + E.C_H2_krwph + E.C_elec_comp_krwph + E.C_heat_krwph + E.C_cool_krwph;

% --------- Unit costs ---------
E.OPEX_var_krw_per_kgjet = E.C_total_krwph / max(1e-9, m_jet);
E.OPEX_var_krw_per_Ljet  = E.OPEX_var_krw_per_kgjet * P.jet_density_kg_per_L;

% --------- Daily totals ---------
E.C_total_krw_per_day = E.C_total_krwph * 24;
E.jet_kg_per_day      = m_jet * 24;

% --------- Optional revenue compare ---------
if include_revenue
    P_sell = P.jet_krw_per_kg * SAF_premium_mult;
    E.P_sell_krw_per_kg = P_sell;

    E.Rev_krwph = m_jet * P_sell;
    E.GM_krwph  = E.Rev_krwph - E.C_total_krwph;

    E.GM_krw_per_day = E.GM_krwph * 24;
else
    E.P_sell_krw_per_kg = NaN;
    E.Rev_krwph = NaN;
    E.GM_krwph  = NaN;
    E.GM_krw_per_day = NaN;
end

% --------- Print report ---------
fprintf("\n=== ECON (Variable OPEX only) ===\n");
fprintf("Feed(oil): %.2f kg/h @ %.0f KRW/kg => %.0f KRW/h\n", m_oil_in, P.soy_oil_krw_per_kg, E.C_oil_krwph);
fprintf("H2 make-up: %.3f kg/h @ %.0f KRW/kg => %.0f KRW/h\n", mH2, P.H2_krw_per_kg, E.C_H2_krwph);
fprintf("Elec (comp): %.1f kW @ %.1f KRW/kWh => %.0f KRW/h\n", W_comp_kW, P.elec_krw_per_kWh, E.C_elec_comp_krwph);
fprintf("Heat (%s): %.1f kW => %.0f KRW/h\n", E.heat_basis, Q_heat_kW, E.C_heat_krwph);
fprintf("Cool (%s): %.1f kW => %.0f KRW/h\n", cooling_model, Q_cool_kW, E.C_cool_krwph);

fprintf("TOTAL: %.0f KRW/h | Jet out: %.2f kg/h\n", E.C_total_krwph, m_jet);
fprintf("OPEX_var: %.0f KRW/kg-jet (%.0f KRW/L @ rho=%.2f)\n", E.OPEX_var_krw_per_kgjet, E.OPEX_var_krw_per_Ljet, P.jet_density_kg_per_L);

if include_revenue
    fprintf("Bench sell price: %.0f KRW/kg (mult=%.2f) => GM: %.0f KRW/h (%.0f KRW/day)\n", ...
        E.P_sell_krw_per_kg, SAF_premium_mult, E.GM_krwph, E.GM_krw_per_day);
end
end
