function [outHDO, outHI, SUM] = run_hdo_hi_and_summary(feed, op, sel, h2, thermo, comp, opHI, cat, h2HI, thermoHI, compHI)
% Runs HDO -> HI and returns integrated KPI summary.

% --- 1) HDO ---
outHDO = hdo_reactor_v4p1_pH2coupled(feed, op, sel, h2, thermo, comp);

% --- 2) HI inlet from HDO ---
m_jet_in = outHDO.m_jetproxy_kgph;

% HI inlet temperature linking
% If HDO is isothermal at op.T_K, assume HI liquid feed enters at that temperature.
if isfield(outHDO,'T_out_K')
    thermoHI.T_in_K = outHDO.T_out_K;
else
    thermoHI.T_in_K = op.T_K;
end

% If H2 comes from compressor outlet at compHI.T_in_K then it is heated to HI temp
if ~isfield(thermoHI,'T_H2_in_K')
    thermoHI.T_H2_in_K = compHI.T_in_K;
end

% --- 3) HI with utilities ---
outHI = hi_block_lumped_v2(m_jet_in, opHI, cat, h2HI, thermoHI, compHI);

% --- 4) Summary KPI ---
SUM = struct();

% mass
SUM.m_oil_in_kgph     = feed.mdot_oil_kgph;
SUM.m_jetproxy_HDO_kgph = outHDO.m_jetproxy_kgph;
SUM.m_jet_out_HI_kgph = outHI.m_jet_out_kgph;
SUM.m_light_HI_kgph   = outHI.m_light_kgph;

SUM.yield_jetproxy_vs_oil = outHDO.m_jetproxy_kgph / max(1e-9, feed.mdot_oil_kgph);
SUM.yield_jet_final_vs_oil = outHI.m_jet_out_kgph / max(1e-9, feed.mdot_oil_kgph);

% composition
SUM.w_n_jet = outHI.w_n_jet;
SUM.w_iso_jet = outHI.w_iso_jet;
SUM.multi_share_in_iso = outHI.multi_share_in_iso;

% hydrogen
SUM.H2_makeup_HDO_kgph = outHDO.mH2_makeup_kgph;
SUM.H2_makeup_HI_kgph  = outHI.mH2_makeup_kgph;
SUM.H2_makeup_total_kgph = SUM.H2_makeup_HDO_kgph + SUM.H2_makeup_HI_kgph;

% power (kW)
SUM.W_comp_HDO_kW = outHDO.W_comp_kW;
SUM.W_comp_HI_kW  = outHI.W_comp_HI_kW;
SUM.W_comp_total_kW = SUM.W_comp_HDO_kW + SUM.W_comp_HI_kW;

% heat duties: split into heating vs cooling for clarity
[SUM.Q_heat_HDO_kW, SUM.Q_cool_HDO_kW] = split_heat(outHDO.Q_preheat_kW);
[SUM.Q_heat_HI_kW,  SUM.Q_cool_HI_kW]  = split_heat(outHI.Q_preheat_HI_kW);

SUM.Q_heat_total_kW = SUM.Q_heat_HDO_kW + SUM.Q_heat_HI_kW;
SUM.Q_cool_total_kW = SUM.Q_cool_HDO_kW + SUM.Q_cool_HI_kW;

% reaction heat
SUM.Q_rxn_HDO_kW = outHDO.Q_rxn_kW;
SUM.Q_rxn_HI_kW  = outHI.Q_rxn_HI_kW;
SUM.Q_rxn_total_kW = SUM.Q_rxn_HDO_kW + SUM.Q_rxn_HI_kW;

% --- 5) Print compact report ---
fprintf("\n=== HDO ===\n");
fprintf("yH2_in=%.3f, pH2=%.2f MPa (Ptot=%.1f MPa), X=%.4f\n", outHDO.yH2_in, outHDO.pH2_MPa, outHDO.Ptot_MPa, outHDO.X_final);
fprintf("jet-proxy=%.2f kg/h | H2_makeup=%.3f kg/h | W_comp=%.1f kW\n", outHDO.m_jetproxy_kgph, outHDO.mH2_makeup_kgph, outHDO.W_comp_kW);
fprintf("Q_rxn=%.1f kW | Q_preheat=%.1f kW\n", outHDO.Q_rxn_kW, outHDO.Q_preheat_kW);
fprintf("T_out=%.2f K (%.2f C) | eta_avg=%.3f\n", outHDO.T_out_K, outHDO.T_out_K-273.15, outHDO.eta_avg);

fprintf("\n=== HI ===\n");
fprintf("jet_in=%.2f kg/h -> jet_out=%.2f kg/h, light=%.2f kg/h (light=%.1f%%)\n", ...
    m_jet_in, outHI.m_jet_out_kgph, outHI.m_light_kgph, 100*outHI.m_light_kgph/max(1e-9,m_jet_in));
fprintf("w_n=%.2f | w_iso=%.2f | multi_share_in_iso=%.2f\n", outHI.w_n_jet, outHI.w_iso_jet, outHI.multi_share_in_iso);
fprintf("H2_in=%.3f kg/h | H2_makeup=%.3f kg/h | W_comp=%.1f kW\n", outHI.mH2_in_kgph, outHI.mH2_makeup_kgph, outHI.W_comp_HI_kW);
fprintf("Q_preheat_HI=%.1f kW (+:heating, -:cooling) | Q_rxn_HI=%.1f kW\n", outHI.Q_preheat_HI_kW, outHI.Q_rxn_HI_kW);

fprintf("\n=== TOTAL (HDO+HI) ===\n");
fprintf("oil_in=%.2f kg/h | jet_final=%.2f kg/h (yield=%.3f)\n", SUM.m_oil_in_kgph, SUM.m_jet_out_HI_kgph, SUM.yield_jet_final_vs_oil);
fprintf("H2_makeup_total=%.3f kg/h | W_comp_total=%.1f kW\n", SUM.H2_makeup_total_kgph, SUM.W_comp_total_kW);
fprintf("Q_heat_total=%.1f kW | Q_cool_total=%.1f kW | Q_rxn_total=%.1f kW\n", ...
    SUM.Q_heat_total_kW, SUM.Q_cool_total_kW, SUM.Q_rxn_total_kW);
end

function [Q_heat, Q_cool] = split_heat(Q)
% returns positive heating duty and positive cooling duty
Q_heat = max(Q, 0);
Q_cool = max(-Q, 0);
end
