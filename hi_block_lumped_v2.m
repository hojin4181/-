function out = hi_block_lumped_v2(m_jet_in_kgph, opHI, cat, h2HI, thermoHI, compHI)
% HI block (lumped) + utilities (H2 + energy + compression)
%
% Inputs
%   m_jet_in_kgph : HDO jet-range feed to HI (e.g., C17+C18 proxy) [kg/h]
%   opHI.T_C      : HI temperature [°C]
%   opHI.pH2_MPa  : hydrogen partial pressure axis [MPa] (surrogate input)
%   opHI.use_LHSV : true/false
%   opHI.LHSV_per_h : if use_LHSV=true
%   opHI.Vcat_cm3 : catalyst bed volume [cm3]
%   opHI.rho_kg_per_L : liquid density [kg/L]
%
%   cat.activity  : activity multiplier (default 1.0)
%
%   h2HI (optional):
%     h2HI.mode = "once_through" | "recycle"
%     h2HI.H2_per_kg_feed = kg H2 / kg liquid feed (co-feed basis)
%     h2HI.purge_frac = purge fraction (recycle mode)
%     h2HI.nuH2_kg_per_kg_light = effective H2 consumption per kg light (default 0)
%
%   thermoHI (optional):
%     thermoHI.enabled = true/false
%     thermoHI.T_in_K, thermoHI.T_H2_in_K
%     thermoHI.cp_liq_kJ_per_kgK, thermoHI.cp_H2_kJ_per_kgK
%     thermoHI.dH_iso_kJ_per_kgIso, thermoHI.dH_crk_kJ_per_kgLight (default 0)
%
%   compHI (optional):
%     compHI.enabled = true/false
%     compHI.P_in_bar, compHI.P_out_bar, compHI.eta, compHI.T_in_K
%     compHI.cp_H2_kJ_per_kgK, compHI.k_H2
%
% Outputs:
%   기존 v1 출력 + out.mH2_in_kgph, out.mH2_cons_kgph, out.mH2_makeup_kgph
%                + out.Q_preheat_HI_kW (±), out.Q_rxn_HI_kW, out.W_comp_HI_kW

% -------- defaults --------
if ~isfield(opHI,'rho_kg_per_L'); opHI.rho_kg_per_L = 0.80; end
if ~isfield(opHI,'Vcat_cm3');     opHI.Vcat_cm3 = 50; end
if ~isfield(opHI,'use_LHSV');     opHI.use_LHSV = true; end
if ~isfield(cat,'activity');      cat.activity = 1.0; end

if nargin < 4 || isempty(h2HI); h2HI = struct(); end
if ~isfield(h2HI,'mode');                   h2HI.mode = "once_through"; end
if ~isfield(h2HI,'H2_per_kg_feed');         h2HI.H2_per_kg_feed = 0.0; end
if ~isfield(h2HI,'purge_frac');             h2HI.purge_frac = 0.02; end
if ~isfield(h2HI,'nuH2_kg_per_kg_light');   h2HI.nuH2_kg_per_kg_light = 0.0; end

if nargin < 5 || isempty(thermoHI); thermoHI = struct(); end
if ~isfield(thermoHI,'enabled');            thermoHI.enabled = false; end
if ~isfield(thermoHI,'T_in_K');             thermoHI.T_in_K = 298.15; end
if ~isfield(thermoHI,'T_H2_in_K');          thermoHI.T_H2_in_K = 298.15; end
if ~isfield(thermoHI,'cp_liq_kJ_per_kgK');  thermoHI.cp_liq_kJ_per_kgK = 2.0; end
if ~isfield(thermoHI,'cp_H2_kJ_per_kgK');   thermoHI.cp_H2_kJ_per_kgK  = 14.3; end
if ~isfield(thermoHI,'dH_iso_kJ_per_kgIso');   thermoHI.dH_iso_kJ_per_kgIso = 0.0; end
if ~isfield(thermoHI,'dH_crk_kJ_per_kgLight'); thermoHI.dH_crk_kJ_per_kgLight = 0.0; end

if nargin < 6 || isempty(compHI); compHI = struct(); end
if ~isfield(compHI,'enabled');        compHI.enabled = false; end
if ~isfield(compHI,'P_in_bar');       compHI.P_in_bar = 20; end
if ~isfield(compHI,'P_out_bar');      compHI.P_out_bar = 10; end
if ~isfield(compHI,'eta');            compHI.eta = 0.75; end
if ~isfield(compHI,'T_in_K');         compHI.T_in_K = 313.15; end
if ~isfield(compHI,'cp_H2_kJ_per_kgK'); compHI.cp_H2_kJ_per_kgK = 14.3; end
if ~isfield(compHI,'k_H2');           compHI.k_H2 = 1.405; end

T = opHI.T_C;
P = opHI.pH2_MPa;

% -------- contact time tau --------
Vcat_L = opHI.Vcat_cm3 / 1000; % cm3 -> L
if opHI.use_LHSV
    tau_h = 1 / opHI.LHSV_per_h;
else
    if ~isfield(opHI,'mdot_liq_kgph')
        error("opHI.use_LHSV=false이면 opHI.mdot_liq_kgph 필요합니다.");
    end
    Vdot_Lph = opHI.mdot_liq_kgph / opHI.rho_kg_per_L;
    tau_h = Vcat_L / Vdot_Lph;
end

% Reference contact time from paper condition: LHSV=0.5 -> tau_ref=2 h
tau_ref = 2.0;

% ======================================================================
% 1) k_crk(T,P) from anchors (Fig.10 surrogate)
% ======================================================================
Tvec = [300 320 350];
fL_T = [0.06 0.12 0.43];          % at P=1 MPa
Pvec = [1 3 5];
fL_P = [0.12 0.10 0.10];          % at T=320°C

kcrk_T = -log(1 - fL_T) / tau_ref;
kcrk_P = -log(1 - fL_P) / tau_ref;

kcrk_ref = kcrk_T(2); % T=320 at P=1
kcrk_Ti = interp1(Tvec, kcrk_T, T, 'linear', 'extrap');
kcrk_Pi = interp1(Pvec, kcrk_P, P, 'linear', 'extrap');
k_crk = (kcrk_Ti * kcrk_Pi / kcrk_ref) * cat.activity;

% ======================================================================
% 2) k_iso(T,P) from n-alkane content anchors (Fig.10 surrogate)
% ======================================================================
wn_T = [0.43 0.28 0.19];          % n-alkane content vs T at P=1
wn_P = [0.28 0.50 0.57];          % n-alkane content vs P at T=320

kiso_T = -log(wn_T) / tau_ref;
kiso_P = -log(wn_P) / tau_ref;

kiso_ref = kiso_T(2); % T=320 at P=1
kiso_Ti = interp1(Tvec, kiso_T, T, 'linear', 'extrap');
kiso_Pi = interp1(Pvec, kiso_P, P, 'linear', 'extrap');
k_iso = (kiso_Ti * kiso_Pi / kiso_ref) * cat.activity;

% ======================================================================
% 3) Predict yields at tau
% ======================================================================
f_light = 1 - exp(-k_crk * tau_h);
f_light = max(0, min(0.95, f_light));  % clamp

m_light = m_jet_in_kgph * f_light;
m_jet_out = m_jet_in_kgph - m_light;

w_n_jet = exp(-k_iso * tau_h);
w_n_jet = max(0.02, min(0.98, w_n_jet));
w_iso_jet = 1 - w_n_jet;

m_n = m_jet_out * w_n_jet;
m_iso = m_jet_out * w_iso_jet;

% ======================================================================
% 4) Split iso into mono vs multi (Fig.12 surrogate)
% ======================================================================
mono_T = [0.53 0.50 0.63];
multi_T = [0.08 0.19 0.21];
ms_T = multi_T ./ (mono_T + multi_T);

mono_P = [0.50 0.47 0.44];
multi_P = [0.19 0.12 0.05];
ms_P = multi_P ./ (mono_P + multi_P);

ms_Ti = interp1(Tvec, ms_T, T, 'linear', 'extrap');
ms_Pi = interp1(Pvec, ms_P, P, 'linear', 'extrap');

ms_ref = ms_T(2); % T=320, P=1
multi_share = (ms_Ti * ms_Pi / ms_ref);
multi_share = max(0.02, min(0.60, multi_share));

m_multi = m_iso * multi_share;
m_mono  = m_iso - m_multi;

% -------- H2 accounting (co-feed + effective consumption) --------
mH2_in   = h2HI.H2_per_kg_feed * m_jet_in_kgph;               % kg/h
mH2_cons = h2HI.nuH2_kg_per_kg_light * m_light;               % kg/h (optional)

mode = string(h2HI.mode);
if mode == "once_through"
    % separated and vented: makeup equals feed
    mH2_makeup = mH2_in;
    mH2_purge  = mH2_in;     % 의미상 배출(once-through)
elseif mode == "recycle"
    % simple purge model: purge loss proportional to reactor inlet H2
    f = max(1e-4, h2HI.purge_frac);
    mH2_purge  = f * mH2_in;
    mH2_makeup = mH2_cons + mH2_purge;
else
    error("h2HI.mode must be 'once_through' or 'recycle'");
end

% -------- Energy accounting --------
T_HI_K = opHI.T_C + 273.15;
Q_preheat_HI_kW = 0;
Q_rxn_HI_kW     = 0;

if thermoHI.enabled
    Q_liq_kW = (m_jet_in_kgph/3600) * thermoHI.cp_liq_kJ_per_kgK * (T_HI_K - thermoHI.T_in_K);
    Q_H2_kW  = (mH2_in/3600)        * thermoHI.cp_H2_kJ_per_kgK  * (T_HI_K - thermoHI.T_H2_in_K);
    Q_preheat_HI_kW = Q_liq_kW + Q_H2_kW;

    % effective reaction heat (tunable); keep 0 if unknown
    Q_rxn_HI_kW = (m_iso*thermoHI.dH_iso_kJ_per_kgIso + m_light*thermoHI.dH_crk_kJ_per_kgLight) / 3600;
end

% -------- Compressor power (fresh H2 makeup only) --------
W_comp_HI_kW = 0;
if compHI.enabled
    PR = compHI.P_out_bar / compHI.P_in_bar;
    if PR > 1.0001
        mdot = mH2_makeup/3600;           % kg/s
        cp   = compHI.cp_H2_kJ_per_kgK;   % kJ/kg-K
        k    = compHI.k_H2;
        W_comp_HI_kW = mdot * cp * compHI.T_in_K * (PR^((k-1)/k) - 1) / compHI.eta;
    else
        W_comp_HI_kW = 0; % no compression needed
    end
end

% -------- pack outputs --------
out.tau_h = tau_h;

out.k_crk_per_h = k_crk;
out.k_iso_per_h = k_iso;

out.f_light = f_light;
out.m_light_kgph = m_light;
out.m_jet_out_kgph = m_jet_out;

out.w_n_jet = w_n_jet;
out.w_iso_jet = w_iso_jet;

out.m_n_kgph = m_n;
out.m_iso_kgph = m_iso;

out.multi_share_in_iso = multi_share;
out.m_mono_kgph = m_mono;
out.m_multi_kgph = m_multi;

% utilities
out.mH2_in_kgph     = mH2_in;
out.mH2_cons_kgph   = mH2_cons;
out.mH2_purge_kgph  = mH2_purge;
out.mH2_makeup_kgph = mH2_makeup;

out.Q_preheat_HI_kW = Q_preheat_HI_kW; % (+) heating, (-) cooling
out.Q_rxn_HI_kW     = Q_rxn_HI_kW;     % sign per user convention
out.W_comp_HI_kW    = W_comp_HI_kW;
end
