%% Dual Cycle Analysis
clc; clear; close all;

 
% Constants
P1 = 100;       % kPa (initial pressure)
T1 = 300;       % K (initial temperature)
k = 1.4;        % Specific heat ratio (gamma)
rp = 1.7;       % Pressure ratio (P3/P2)
cutoff_ratio_percent = 5;  % 5% of stroke for cut-off
Cv = 0.718;     % kJ/kg·K (specific heat atconstant volume)
Cp = 1.005;     % kJ/kg·K (specific heat at constant pressure)
 
% Initialize results
results = struct();
resultIndex = 1;
 
% Calculate for rc = 14 to 18
for rc = 14:18
    % 1. Isentropic compression (1 -> 2)
    T2 = T1 * (rc ^ (k - 1));
    P2 = P1 * (rc ^ k);
    
    % 2. Constant volume heataddition (2 -> 3)
    P3 = rp * P2;
    T3 = T2 * rp;
    
    % 3. Constant pressure heat addition (3 -> 4)
    rho = 1 + (cutoff_ratio_percent / 100) * (rc - 1);
    T4 = T3 * rho;
    
    % Skip if T4 exceeds 2500 K
    if T4 > 2500
        continue;
    end
    
    % Heat calculations
    Qv = Cv * (T3 - T2);
    Qp = Cp * (T4 - T3);
    Qin = Qv + Qp;
    
    % 4. Isentropic expansion (4 -> 5)
    V4_V5 = rho / rc;
    T5 = T4 * (V4_V5 ^ (k - 1));
    
    % Heat rejection
    Qout = Cv * (T5 - T1);
    
    % Efficiency
    eta = 1 - (Qout / Qin);
    
    % Store results (using proper roundingsyntax)
    results.rc(resultIndex) = rc;
    results.T2(resultIndex) = round(T2*100)/100;
    results.P2(resultIndex) = round(P2*100)/100;
    results.T3(resultIndex) = round(T3*100)/100;
    results.P3(resultIndex) = round(P3*100)/100;
    results.T4(resultIndex) = round(T4*100)/100;
    results.rho(resultIndex) = round(rho*1000)/1000;
    results.T5(resultIndex) = round(T5*100)/100;
    results.Qv(resultIndex) = round(Qv*100)/100;
    results.Qp(resultIndex) = round(Qp*100)/100;
    results.Qin(resultIndex) = round(Qin*100)/100;
    results.Qout(resultIndex) = round(Qout*100)/100;
    results.Efficiency(resultIndex) = round(eta*10000)/100;
    
    resultIndex =resultIndex + 1;
end
 
% Display results
disp('Dual Cycle Analysis (rc= 14 to 18, T4 ≤ 2500 K)');
disp('================================================================================');
disp(' rc  |  T2(K) | P2 (kPa) | T3 (K) | P3 (kPa) | T4 (K) |  ρ   | T5 (K) | Qin   |  Qout  | η (%) ');
disp('================================================================================');
 
for i = 1:length(results.rc)
    fprintf('%3d  | %7.2f |%8.2f | %6.2f | %8.2f | %6.2f | %4.3f | %6.2f | %6.2f | %6.2f | %5.2f\n',results.rc(i), results.T2(i), results.P2(i), results.T3(i),results.P3(i), results.T4(i), results.rho(i), results.T5(i),results.Qin(i), results.Qout(i), results.Efficiency(i));
end

 
% Plot efficiency vs compression ratio
figure;
plot(results.rc, results.Efficiency, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
title('Dual Cycle Efficiencyvs Compression Ratio', 'FontSize', 12);
xlabel('Compression Ratio (rc)', 'FontSize', 10);
ylabel('Efficiency (%)', 'FontSize', 10);
grid on;
 
% Add data labels
for i = 1:length(results.rc)
    text(results.rc(i), results.Efficiency(i), ...
        sprintf('%.2f%%', results.Efficiency(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', 8);
end