function P = econ_prices_baseline_2025KR()
% Baseline variable-price set (KRW basis), sourced from public references.
% NOTE: Use as a report baseline; adjust to your scenario/contract.

P = struct();

% ---- FX (KRW per USD) ----
P.usdkrw = 1474.96;  % Korea Customs exchange rate (example value)

% ---- Feedstock: Soybean oil, crude degummed (USD/mt) ----
P.soy_oil_usd_per_t = 1126;   % World Bank Pink Sheet (Nov 2025, USD/mt)
P.soy_oil_krw_per_kg = (P.soy_oil_usd_per_t/1000) * P.usdkrw;

% ---- Hydrogen (KRW/kg) ----
P.H2_krw_per_kg = 10285;  % H2nbiz average (often retail-like). Adjust if industrial.

% ---- Electricity (KRW/kWh) ----
% Example: industrial high-voltage A option I, spring/fall energy charge
P.elec_krw_per_kWh = 101.1;

% ---- Natural gas (KRW/MJ) -> (KRW/kWh_th) ----
P.ng_krw_per_MJ = 16.2467;          % KOGAS wholesale (example)
P.ng_krw_per_kWh_th = P.ng_krw_per_MJ * 3.6;

% ---- Utility efficiencies / COP ----
P.boiler_eta = 0.90;    % thermal efficiency for NG boiler (assumption)
P.eheater_eta = 0.98;   % electric heater efficiency (assumption)
P.COP_cool = 5.0;       % chiller COP (assumption; CW cooling could be ~0 cost)

% ---- Jet fuel benchmark (optional revenue compare) ----
P.jet_usd_per_bbl = 86.88;   % IATA global average (example)
P.jet_density_kg_per_L = 0.80;  % typical Jet-A1 density (assumption)
P.bbl_to_L = 158.987;

kg_per_bbl = P.bbl_to_L * P.jet_density_kg_per_L;
P.jet_krw_per_kg = (P.jet_usd_per_bbl / kg_per_bbl) * P.usdkrw;
end
