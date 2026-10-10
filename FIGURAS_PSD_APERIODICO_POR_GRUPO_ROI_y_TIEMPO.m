%% SCRIPT: PSD COMPARATIVO (GRUPO 1 VS GRUPO 2) POR ROI
clearvars; close all; clc;

% --- 1. CONFIGURACIÓN DE RUTAS ---
baseDir = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\EEG_MATLAB';
outputBase = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\APERIODICO';
saveDir = fullfile(outputBase, 'GRAFICAS_PSD_POR_ROI');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Carga de resultados de Excel
masterFileName = fullfile(outputBase, 'DEF_FOOOF_limpio_colap.xlsx');
masterData = readtable(masterFileName);

% Colores y etiquetas [cite: 2026-03-30]
col_g1 = [1.00, 0.55, 0.41]; % Placebo (Salmón)
col_g2 = [0.25, 0.88, 0.82]; % Real (Turquesa)
p_names = {'Delta','Theta','Alfa','Beta','Gamma'};
centros = [2, 6, 10, 22, 35];
rois_labels = {'Fp1-F7', 'Fp2-F4', 'C3-P3', 'C4-P4', 'Midline', 'T3-T5', 'T4-T6', 'Occipital'};
tiempos_nombres = {'Baseline', 'tSMS 1', 'tSMS 2', 'tSMS 3', 'Post Tarea'};

%% --- 2. CARGA DE DATOS (PSD POR CANAL) ---
fprintf('🚀 Iniciando carga de espectros por ROI...\n');
all_data = []; 
freqs_real = [];

for g = 1:2
    grupoPath = fullfile(baseDir, num2str(g));
    sujetos = dir(fullfile(grupoPath, 'S*')); 
    for s = 1:length(sujetos)
        sujPath = fullfile(grupoPath, sujetos(s).name);
        archivos = {dir(fullfile(sujPath, '*pre_OA.mat')), ...
                    dir(fullfile(sujPath, '*tSMS.mat')), ...
                    dir(fullfile(sujPath, '*post.mat'))};
        for a = 1:3
            if ~isempty(archivos{a})
                dataStruct = load(fullfile(archivos{a}.folder, archivos{a}.name));
                if isfield(dataStruct, 'EEG')
                    EEG = dataStruct.EEG;
                    [pxx, f] = pwelch(EEG.data', hanning(EEG.srate*2), [], EEG.srate*2, EEG.srate);
                    psd_per_roi = pxx'; 
                    freqs_real = f;
                    
                    t_vector = [];
                    if a == 1, t_vector = 1;
                    elseif a == 2, t_vector = [2, 3, 4];
                    elseif a == 3, t_vector = 5; end
                    
                    for t = t_vector
                        for r = 1:8
                            all_data = [all_data; g, s, t, r, psd_per_roi(r, :)];
                        end
                    end
                end
            end
        end
    end
end

%% --- 3. GENERACIÓN DE FIGURAS (UNA POR ROI) ---
f_mask = freqs_real >= 3.5 & freqs_real <= 40;
f_plot_model = linspace(3.5, 40, 200);

for r = 1:8
    % Crear figura para el ROI actual
    h = figure('Name', ['PSD_ROI_' num2str(r)], 'Color', 'w', 'Position', [100 50 900 950]);
    
    for t = 1:5
        for g = 1:2
            % Calcular posición en el subplot (5 filas, 2 columnas)
            % Columna 1: Placebo | Columna 2: Real
            pos = (t-1)*2 + g;
            subplot(5, 2, pos); hold on; grid on;
            
            % Filtrar datos
            idx_raw = (all_data(:,1) == g & all_data(:,3) == t & all_data(:,4) == r);
            idx_excel = (masterData.Grupo == g & masterData.Tiempo == t & masterData.ROI == r);
            
            if any(idx_raw) && any(idx_excel)
                % Promedio señal bruta
                mean_raw_log = log10(mean(all_data(idx_raw, 5:end), 1, 'omitnan'));
                
                % Parámetros FOOOF
                m_off = mean(masterData.Offset(idx_excel), 'omitnan');
                m_exp = mean(masterData.Exponent(idx_excel), 'omitnan');
                
                % Alineación y modelo
                model_ap_vals = m_off - log10(freqs_real(f_mask).^m_exp);
                offset_correction = mean(model_ap_vals) - mean(mean_raw_log(f_mask));
                mean_aligned = mean_raw_log + offset_correction;
                
                % Plot original (gris)
                plot(freqs_real(f_mask), mean_aligned(f_mask), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
                
                % Plot modelo aperiódico (azul punteado)
                ap_fit = m_off - log10(f_plot_model.^m_exp);
                plot(f_plot_model, ap_fit, '--b', 'LineWidth', 1);
                
                % Plot modelo completo (Color según grupo)
                col = if_func(g==1, col_g1, col_g2);
                full_fit = ap_fit;
                for j = 1:5
                    m_amp = mean(masterData.(['Peak_' p_names{j}])(idx_excel), 'omitnan');
                    if m_amp > 0
                        full_fit = full_fit + m_amp * exp(-(f_plot_model - centros(j)).^2 / (2 * 1.5^2));
                    end
                end
                plot(f_plot_model, full_fit, 'Color', col, 'LineWidth', 2);
            end
            
            % Estética
            ylim([-1.2 0.7]); xlim([3.5 40]);
            if t == 1, title(if_func(g==1, 'PLACEBO (Grupo 1)', 'REAL (Grupo 2)'), 'FontSize', 12); end
            if g == 1, ylabel(tiempos_nombres{t}, 'FontWeight', 'bold', 'FontSize', 10); end
            if t < 5, xticklabels([]); end
            if g == 2, yticklabels([]); end
        end
    end
    
    % Título general de la figura (ROI)
    annotation('textarrow',[0.5 0.5],[0.95 0.95],'string', ['Análisis PSD + FOOOF - ROI: ' rois_labels{r}], ...
        'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','FontWeight','bold','FontSize',14);
    
    % Guardar
    saveas(h, fullfile(saveDir, ['PSD_ROI_' num2str(r) '_' rois_labels{r} '.png']));
end

fprintf('✅ Proceso completado. 8 figuras generadas en: %s\n', saveDir);

% Función auxiliar para evitar conflictos
function out = if_func(cond, valT, valF)
    if cond, out = valT; else, out = valF; end
end