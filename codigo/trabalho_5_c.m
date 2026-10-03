clc;
clear;
close all;

%% =============================
% 1. COMPONENTES (APENAS VISUALIZAÇÃO)
%% =============================

% DNA
data = load('espectros_c/DNA.ttt');
lambda = data(:,1);
mu_DNA = data(:,2);
mu_DNA_c = mu_DNA / max(mu_DNA);

figure;
plot(lambda, mu_DNA_c, 'LineWidth', 2);
title('DNA'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% Proteínas
data = load('espectros_c/Proteinas.ttt');
mu_prot = data(:,2);
mu_prot_c = mu_prot / max(mu_prot);

figure;
plot(lambda, mu_prot_c, 'LineWidth', 2);
title('Proteínas'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% HbO2
data = load('espectros_c/HbO.ttt');
mu_HbO2 = data(:,2);
mu_HbO2_c = mu_HbO2 / max(mu_HbO2);

figure;
plot(lambda, mu_HbO2_c, 'LineWidth', 2);
title('HbO2'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% Hb
data = load('espectros_c/Hb.ttt');
mu_Hb = data(:,2);
mu_Hb_c = mu_Hb / max(mu_Hb);

figure;
plot(lambda, mu_Hb_c, 'LineWidth', 2);
title('Hb'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% Lípidos
data = load('espectros_c/Lipids.ttt');
mu_lip = data(:,2);
mu_lip_c = mu_lip / max(mu_lip);

figure;
plot(lambda, mu_lip_c, 'LineWidth', 2);
title('Lípidos'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% Água
data = load('espectros_c/H2O.ttt');
mu_H2O = data(:,2);
mu_H2O_c = mu_H2O / max(mu_H2O);

figure;
plot(lambda, mu_H2O_c, 'LineWidth', 2);
title('Água'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;

% Melanina + Lipofuscina
mel  = load('espectros_c/Melanina.ttt');
lipo = load('espectros_c/Lipofuscina.ttt');

mu_mel  = mel(:,2);
mu_lipo = lipo(:,2);

mu_mel_c  = mu_mel  / max(mu_mel);
mu_lipo_c = mu_lipo / max(mu_lipo);

figure;
plot(lambda, mu_mel_c,  'b', 'LineWidth', 2); hold on;
plot(lambda, mu_lipo_c, 'r', 'LineWidth', 2);
legend('Melanina', 'Lipofuscina');
title('Pigmentos'); xlabel('\lambda (nm)'); ylabel('\mu_a (A.U.)'); grid on;


%% =============================
% 2. CÁLCULO DOS 10 ESPECTROS
%% =============================

N = 10;
d = 0.05; % cm

for i = 1:N
    Rt_data = load(sprintf('espectros_c/Rt%d.ttt', i));
    Tt_data = load(sprintf('espectros_c/Tt%d.ttt', i));

    lambda  = Rt_data(:,1);
    Rt      = 0.01 .* Rt_data(:,2);
    Tt      = 0.01 .* Tt_data(:,2);

    mu_a(:,i) = (1 - (Tt + Rt)) / d;
end


%% =============================
% 3. MÉDIA + DESVIO PADRÃO
%% =============================

mu_mean = mean(mu_a, 2);
mu_sd   = std(mu_a, 0, 2);

figure;
stdshade(mu_a', 0.2, 'b', lambda', 1);
xlabel('\lambda (nm)');
ylabel('\mu_a (cm^{-1})');
title('Espectro médio com desvio padrão');
grid on;


%% =============================
% 4. RECONSTRUÇÃO
%% =============================

% Matriz X com espetros normalizados pelo máximo de cada componente
% para equilibrar magnitudes e não penalizar nenhum componente no ajuste.
% Ordem: DNA, Proteinas, HbO2, Hb, H2O, Lipidos, Melanina, Lipofuscina
X = [mu_DNA_c, mu_prot_c, mu_HbO2_c, mu_Hb_c, ...
     mu_H2O_c, mu_lip_c,  mu_mel_c,  mu_lipo_c];

% Vetor alvo: espetro médio experimental (unidades reais, cm^-1)
y = mu_mean;

% -------------------------------------------------------------------------
% Passo 1: correr lsqnonneg sem lípidos para estimar sum(w_outros)
% -------------------------------------------------------------------------
X_sem_lip = [mu_DNA_c, mu_prot_c, mu_HbO2_c, mu_Hb_c, ...
             mu_H2O_c, mu_mel_c,  mu_lipo_c];

w_resto_temp  = lsqnonneg(X_sem_lip, y);
soma_w_outros = sum(w_resto_temp);

% -------------------------------------------------------------------------
% Passo 2: estimar w_lip_min proporcional à concentração mínima do artigo
% No artigo, lípidos = 0.03% no tecido com cancro (mínimo reportado)
% Os restantes componentes (sem água) somam 23%
% → w_lip_min = (0.03 / 23) * soma_w_outros
% -------------------------------------------------------------------------
conc_lip_min = 0.03; % % (concentração mínima do artigo, tecido com cancro)
w_lip_min    = (conc_lip_min / 23) * soma_w_outros;

fprintf('w_lip_min estimado: %.6f\n', w_lip_min);

% -------------------------------------------------------------------------
% Passo 3: subtrair contribuição mínima dos lípidos ao espetro
% experimental e correr lsqnonneg nos restantes 7 componentes
% -------------------------------------------------------------------------
y_corr   = y - w_lip_min * mu_lip_c;
w_resto  = lsqnonneg(X_sem_lip, y_corr);

% Vetor completo de pesos (lípidos na posição 6)
% Ordem: DNA, Prot, HbO2, Hb, H2O, Lip, Mel, Lipo
w = [w_resto(1:5); w_lip_min; w_resto(6:7)];

% Espetro reconstruído
mu_recon = X * w;


%% =============================
% 5. GRÁFICOS DA RECONSTRUÇÃO
%% =============================

% Comparação em unidades reais
figure;
plot(lambda, y,        'b',  'LineWidth', 2); hold on;
plot(lambda, mu_recon, 'r--','LineWidth', 2);
legend('Experimental médio', 'Reconstruído');
xlabel('\lambda (nm)');
ylabel('\mu_a (cm^{-1})');
title('Reconstrução do espectro');
grid on;

% Comparação normalizada (apenas visual)
figure;
plot(lambda, y/max(y),               'b',  'LineWidth', 2); hold on;
plot(lambda, mu_recon/max(mu_recon), 'r--','LineWidth', 2);
legend('Experimental norm.', 'Reconstruído norm.');
xlabel('\lambda (nm)');
ylabel('\mu_a (A.U.)');
title('Comparação normalizada');
grid on;


%% =============================
% 6. CÁLCULO DAS CONCENTRAÇÕES
%% =============================
% Os pesos w foram estimados livremente pelo lsqlin com lb >= 0.
% A concentração da água é fixada a 77% (valor de referência bibliográfico).
% Os restantes 23% são distribuídos pelos outros 7 componentes
% proporcionalmente aos seus pesos.

nomes = {'DNA','Proteínas','HbO2','Hb','Água','Lípidos','Melanina','Lipofuscina'};

idx_agua   = 5;                % posição da água na matriz X
idx_outros = [1, 2, 3, 4, 6, 7, 8]; % restantes componentes

w_outros    = w(idx_outros);
conc_outros = (w_outros / sum(w_outros)) * 23;

conc_all = zeros(8,1);
conc_all(idx_agua)   = 77;
conc_all(idx_outros) = conc_outros;

% Tabela de resultados
fprintf('\n====== CONCENTRAÇÕES DOS COMPONENTES ======\n');
fprintf('%-15s  %10s  %12s\n', 'Componente', 'Peso (w)', 'Conc. (%)');
fprintf('%s\n', repmat('-', 1, 42));
for k = 1:8
    fprintf('%-15s  %10.4f  %12.4f\n', nomes{k}, w(k), conc_all(k));
end
fprintf('%s\n', repmat('-', 1, 42));
fprintf('%-15s  %10s  %12.4f\n', 'TOTAL', '', sum(conc_all));

% Gráfico de barras
figure;
bar(conc_all, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', nomes, 'XTickLabelRotation', 30);
ylabel('Concentração (%)');
title('Concentração dos componentes no tecido');
grid on;


%% =============================
% 7. DIFERENÇA ESPETRAL (EXPERIMENTAL vs RECONSTRUÍDO)
%% =============================

% Diferença absoluta ponto a ponto
diff_abs = abs(y - mu_recon);

% Diferença percentual relativa ao experimental
% (epsilon para evitar divisão por zero)
epsilon   = 1e-10;
diff_perc = (diff_abs ./ (abs(y) + epsilon)) * 100;

diff_perc_media = mean(diff_perc);
diff_perc_max   = max(diff_perc);

fprintf('\n====== DIFERENÇA ESPETRAL ======\n');
fprintf('Diferença percentual média (200-1000 nm): %.4f %%\n', diff_perc_media);
fprintf('Diferença percentual máxima:              %.4f %%\n', diff_perc_max);

% Gráfico da diferença absoluta
figure;
plot(lambda, diff_abs, 'k', 'LineWidth', 2);
xlabel('\lambda (nm)');
ylabel('|\mu_{a,exp} - \mu_{a,recon}| (cm^{-1})');
title('Diferença absoluta: Experimental - Reconstruído');
grid on;

% Gráfico da diferença percentual
figure;
plot(lambda, diff_perc, 'm', 'LineWidth', 2);
xlabel('\lambda (nm)');
ylabel('Diferença relativa (%)');
title('Diferença percentual: Experimental vs Reconstruído');
grid on;
yline(diff_perc_media, 'r--', 'LineWidth', 1.5, ...
      'Label', sprintf('Média: %.2f%%', diff_perc_media));