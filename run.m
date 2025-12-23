% --- (대표님 기존) HDO inputs ---
feed.mdot_oil_kgph = 10000/24;

op.T_K = 340 + 273.15;
op.Ptot_MPa = 10;
op.LHSV_per_h = 0.5;
op.T_feed_K = 298.15;

do_eta_sweep = true;   % 민감도 실행 여부
do_plot      = true;   % 그래프 출력 여부
quiet_sweep  = true;   % 스윕 구간 콘솔출력 억제(추천)

% ===== 여기(=HDO op 블록)에 추가 =====
op.rx_mode = "isothermal";   % 기본(기존과 동일)
%op.rx_mode = "adiabatic";  % 단열 프로파일(에너지 연동)

%op.mt.enabled = false;

op.mt.enabled = true;
op.mt.Rp_m = 0.7e-3;        % 촉매 pellet 반경 [m] (예: 0.7 mm)
op.mt.Deff_m2s = 1.0e-9;    % 유효확산계수 [m^2/s] (예: 1e-9)
op.Nseg = 200;              % (권장) 분할 수: mt 켜면 적분 방식
op.mt.override = false;   % 기본 실행은 계산모드(=Rp, Deff로 eta 계산)
op.mt.eta_fixed = 1.0;    % override=true일 때만 의미


sel.s1_HDO = 0.4;
sel.s2_DCO2 = 0.3;
sel.s3_DCO  = 0.3;

h2.lambda_H2 = 20;
h2.purge_frac = 0.02;
h2.preheat_basis = "makeup";

thermo.enabled = true;
thermo.dH1_kJ_per_molTAG = -1200;
thermo.dH2_kJ_per_molTAG = -800;
thermo.dH3_kJ_per_molTAG = -900;
thermo.cp_oil_kJ_per_kgK = 2.0;

comp.enabled = true;
comp.P_in_bar = 20;
comp.P_out_bar = 100;
comp.eta = 0.75;
comp.T_in_K = 313.15;

% --- HI operating ---
opHI.T_C = 320;
opHI.pH2_MPa = 1;           % surrogate axis
opHI.use_LHSV = true;
opHI.LHSV_per_h = 0.5;
opHI.Vcat_cm3 = 50;

cat.activity = 1.0;

% --- HI hydrogen model (대표님이 민감도 대상) ---
h2HI.mode = "once_through";         % 우선 단순화(공급=배출)
h2HI.H2_per_kg_feed = 0.02;         % 예: 액상 대비 2 wt% co-feed (가정값)
h2HI.purge_frac = 0.02;             % recycle 모드에서만 의미
h2HI.nuH2_kg_per_kg_light = 0.0;    % 근거 없으면 0, 민감도만 제시 권장

% --- HI thermo (유입온도는 HDO 출구로 자동 연결됨) ---
thermoHI.enabled = true;
thermoHI.cp_liq_kJ_per_kgK = 2.0;
thermoHI.cp_H2_kJ_per_kgK  = 14.3;
thermoHI.T_H2_in_K = 298.15;        % H2 storage temp 가정
thermoHI.dH_iso_kJ_per_kgIso   = 0; % 근거 없으면 0
thermoHI.dH_crk_kJ_per_kgLight = 0; % 근거 없으면 0

% --- HI compressor (HI가 1 MPa=10 bar면, 별도 압축이 필요 없을 수도 있음) ---
compHI.enabled = false;             % 필요 시 true로
compHI.P_in_bar = 20;
compHI.P_out_bar = 10;
compHI.eta = 0.75;
compHI.T_in_K = 313.15;

% --- Run ---
% ===== BASE CASE (no override) =====
op_base = op;
op_base.mt.override = false;    % 자동 계산 모드
op_base.mt.eta_fixed = 1.0;     % 의미 없음(안전용)
[outHDO, outHI, SUM] = run_hdo_hi_and_summary(feed, op_base, sel, h2, thermo, comp, opHI, cat, h2HI, thermoHI, compHI);


P = econ_prices_baseline_2025KR();
opts = struct();
opts.heat_source = "ng";          % "ng" or "electric"
opts.cooling_model = "chiller";   % "chiller" | "cw" | "none"
opts.include_revenue = true;
opts.SAF_premium_mult = 1.0;      % SAF 프리미엄 반영하려면 1.2~4.2 등으로 민감도

E = econ_variable_costs_v1(feed, op_base, outHDO, SUM, P, opts);

% ===== eta sweep (sensitivity) =====
if do_eta_sweep
    etas = linspace(0.6, 1.0, 9);

    X     = zeros(size(etas));
    yield = zeros(size(etas));
    opex  = zeros(size(etas));

    for i = 1:numel(etas)
        op_i = op_base;                 % 베이스로부터 복사(중요)
        op_i.mt.enabled   = true;
        op_i.mt.override  = true;       % 여기서만 "강제지정"
        op_i.mt.eta_fixed = etas(i);

        if quiet_sweep
            evalc('[outHDO_i, outHI_i, SUM_i] = run_hdo_hi_and_summary(feed, op_i, sel, h2, thermo, comp, opHI, cat, h2HI, thermoHI, compHI);');
            evalc('E_i = econ_variable_costs_v1(feed, op_i, outHDO_i, SUM_i, P, opts);');
        else
            [outHDO_i, outHI_i, SUM_i] = run_hdo_hi_and_summary(feed, op_i, sel, h2, thermo, comp, opHI, cat, h2HI, thermoHI, compHI);
            E_i = econ_variable_costs_v1(feed, op_i, outHDO_i, SUM_i, P, opts);
        end

        X(i)     = outHDO_i.X_final;
        yield(i) = SUM_i.yield_jet_final_vs_oil;
        opex(i)  = E_i.OPEX_var_krw_per_kgjet;
    end

    % 표(보고서용)
    T_eta = table(etas(:), X(:), yield(:), opex(:), ...
        'VariableNames', {'eta','X_HDO','yield_jet_vs_oil','OPEX_KRW_per_kgjet'});
    disp(T_eta);

    % 그래프(보고서용)
    if do_plot
        figure; plot(etas, X, '-o'); xlabel('\eta'); ylabel('X_{HDO}');
        figure; plot(etas, yield, '-o'); xlabel('\eta'); ylabel('Jet yield (kg/kg-oil)');
        figure; plot(etas, opex, '-o'); xlabel('\eta'); ylabel('OPEX (KRW/kg-jet)');
    end
end

% ===== SAF premium multiplier sweep (BASE economics only) =====
do_premium_sweep = true;

if do_premium_sweep
    % base outputs
    m_jet_kgph = SUM.m_jet_out_HI_kgph;

    % 아래 줄은 변수명에 맞게 택1
    C_total_krwph = E.C_total_krwph;                 % (대표님이 E로 쓰는 경우)
    OPEX_krw_per_kg = E.OPEX_var_krw_per_kgjet;      % (대표님이 E로 쓰는 경우)
    % C_total_krwph = E_base.C_total_krwph;          % (E_base로 쓰는 경우)
    % OPEX_krw_per_kg = E_base.OPEX_var_krw_per_kgjet;

    P0 = P.jet_krw_per_kg; % benchmark jet price (KRW/kg)

    % 손익분기 multiplier (판매가가 원가와 같아지는 지점)
    mult_BE = OPEX_krw_per_kg / max(1e-9, P0);
    fprintf("\n=== Premium sweep (BASE) ===\n");
    fprintf("Break-even: mult_BE=%.2f (P_sell=%.0f KRW/kg)\n", mult_BE, OPEX_krw_per_kg);

    % 스윕 구간(보고서용으로 1.0~4.2 중심 + 여유구간)
    mults = [1.0 1.2 1.5 2.0 2.5 3.0 3.2 3.4 3.6 4.0 4.2];

    P_sell = P0 .* mults;                           % KRW/kg
    Rev_krwph = m_jet_kgph .* P_sell;               % KRW/h
    GM_krwph  = Rev_krwph - C_total_krwph;          % KRW/h
    GM_day    = GM_krwph * 24;                      % KRW/day
    margin_per_kg = P_sell - OPEX_krw_per_kg;       % KRW/kg-jet

    T_prem = table(mults(:), P_sell(:), margin_per_kg(:), GM_krwph(:), GM_day(:), ...
        'VariableNames', {'mult','P_sell_KRW_per_kg','Margin_KRW_per_kg','GM_KRW_per_h','GM_KRW_per_day'});
    disp(T_prem);

    if do_plot
        figure;
        plot(mults, GM_krwph, '-o'); grid on;
        yline(0,'--');
        xlabel('SAF premium multiplier'); ylabel('Gross margin (KRW/h)');

        figure;
        plot(mults, margin_per_kg, '-o'); grid on;
        yline(0,'--');
        xlabel('SAF premium multiplier'); ylabel('Unit margin (KRW/kg-jet)');
    end
end
