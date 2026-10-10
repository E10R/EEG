%% SCRIPT DE CONECTIVIDAD TOTAL: INTRA, INTER E INFORME POR HOJAS
basePath = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\EEG_MATLAB';
outPath  = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\DEF_PLV_APERIODICO';
picosPath = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\APERIODICO_DEF\Presencia_picos.xlsx';

if ~exist(outPath, 'dir'), mkdir(outPath); end

% 1. Configuración de Grupos
grupos = zeros(42,1);
grupos([2,4,6,8,11,12,14,19,20,21,23,26,28,29,33,34,35,37,39,42]) = 1; % Placebo
grupos([1,3,5,7,9,10,13,15,17,18,22,24,25,27,30,31,32,36,38,40,41]) = 2; % Real

% 2. ROIs y Bandas
clusters.ROI1 = {'Fp1', 'F7'};
clusters.ROI2 = {'Fp2', 'F4', 'F8'};
clusters.ROI3 = {'C3', 'P3'};
clusters.ROI4 = {'C4', 'P4'};
clusters.ROI5 = {'Fz', 'Cz', 'Pz'};
clusters.ROI6 = {'T3', 'T5'};
clusters.ROI7 = {'T4', 'T6'};
clusters.ROI8 = {'O1', 'O2', 'Oz'};

bandas = struct('Delta',[0.5 3.5], 'Theta',[4 7.5], 'Alfa',[8 12.5], 'Beta',[13 29.5], 'Gamma',[30 40]);
bNames = fieldnames(bandas);
cNames = fieldnames(clusters);

% Mapeo de tiempos para Excel y Picos
mapTiempos = struct('Pre', 1, 'tSMS_1', 2, 'tSMS_2', 3, 'tSMS_3', 4, 'Post', 5);
tiempos_etiquetas = struct('Pre','Tiempo 1', 'tSMS_1','Tiempo 2', 'tSMS_2','Tiempo 3', 'tSMS_3','Tiempo 4', 'Post','Tiempo 5');

% 3. Carga de Picos
T_picos = readtable(picosPath);
T_picos.Properties.VariableNames = strtrim(T_picos.Properties.VariableNames); 

% Acumuladores para el Informe N
n_report_G1 = struct(); 
n_report_G2 = struct(); 

% 4. Bucle Principal de Sujetos
for s = 1:42
    if grupos(s) == 0, continue; end
    subjStr = sprintf('S%02d', s);
    gNum = grupos(s);
    subjFullPath = fullfile(basePath, num2str(gNum), subjStr);
    targetDir = fullfile(outPath, num2str(gNum));
    
    if ~exist(subjFullPath, 'dir'), continue; end
    
    Sujeto_PLV_Matrices = struct(); 
    sujeto_table = table(); 
    
    cond_files = {'1_pre_OA', '3_tSMS', '4_post'};
    cond_labels = {'Pre', 'tSMS', 'Post'};

    for k = 1:length(cond_files)
        fileFull = fullfile(subjFullPath, sprintf('EEG_preprocess_%s_%s.mat', subjStr, cond_files{k}));
        if ~exist(fileFull, 'file'), continue; end
        
        d = load(fileFull); v = fieldnames(d); cleanEEG = d.(v{1});
        is_tSMS = strcmp(cond_labels{k}, 'tSMS');
        n_segs = is_tSMS * 3 + ~is_tSMS * 1;
        pnts_seg = floor(size(cleanEEG.data, 2) / n_segs);

        for i_seg = 1:n_segs
            cond_actual = cond_labels{k};
            if is_tSMS, cond_actual = sprintf('tSMS_%d', i_seg); end
            t_pico = mapTiempos.(cond_actual);
            
            idx_range = ((i_seg-1)*pnts_seg + 1) : (i_seg*pnts_seg);
            current_data = cleanEEG.data(:, idx_range);
            
            % Fila base para el Excel Intra
            row_data = table({subjStr}, {cond_actual}, 'VariableNames', {'Sujeto', 'Condicion'});

            for b = 1:numel(bNames)
                bN = bNames{b};
                lim = bandas.(bN);
                
                % Filtrado
                tempEEG_data = eegfiltfft(double(current_data), cleanEEG.srate, lim(1), lim(2));
                fasi_tutti = angle(hilbert(tempEEG_data.'))';
                segnali_medi = zeros(8, size(tempEEG_data, 2));
                picos_ok = zeros(8, 1);

                % --- INTRA (WITHIN) ---
                for c = 1:8
                    cN = cNames{c};
                    col_within = ['Within_' bN '_' cN];
                    idx_p = (T_picos.Sujeto == s) & (T_picos.ROI == c) & (T_picos.Tiempo == t_pico);
                    
                    if any(idx_p) && T_picos.(['Peak_' bN])(find(idx_p,1)) == 1
                        picos_ok(c) = 1;
                        ch_idx = find(ismember({cleanEEG.chanlocs.labels}, clusters.(cN)));
                        if ~isempty(ch_idx)
                            segnali_medi(c, :) = mean(tempEEG_data(ch_idx, :), 1);
                            if length(ch_idx) > 1
                                m_w = ones(length(ch_idx));
                                for x = 1:length(ch_idx), for y = (x+1):length(ch_idx)
                                    val = abs(mean(exp(1i*(fasi_tutti(ch_idx(x),:)-fasi_tutti(ch_idx(y),:)))));
                                    m_w(x,y)=val; m_w(y,x)=val;
                                end; end
                                row_data.(col_within) = mean(m_w(~eye(size(m_w))));
                            else, row_data.(col_within) = 1; end
                        end
                    else
                        row_data.(col_within) = NaN; % Asegura que no salga vacío
                    end
                end

                % --- INTER (BETWEEN) ---
                m_b = nan(8);
                fN = matlab.lang.makeValidName([cond_actual '_' bN]);
                for c1 = 1:8, for c2 = (c1+1):8
                    if picos_ok(c1) == 1 && picos_ok(c2) == 1
                        df_b = angle(hilbert(segnali_medi(c1,:))) - angle(hilbert(segnali_medi(c2,:)));
                        m_b(c1, c2) = abs(mean(exp(1i * df_b)));
                        m_b(c2, c1) = m_b(c1, c2);
                        
                        % Sumar al conteo por grupo
                        if gNum == 1
                            if ~isfield(n_report_G1, fN), n_report_G1.(fN) = zeros(8); end
                            n_report_G1.(fN)(c1, c2) = n_report_G1.(fN)(c1, c2) + 1;
                        else
                            if ~isfield(n_report_G2, fN), n_report_G2.(fN) = zeros(8); end
                            n_report_G2.(fN)(c1, c2) = n_report_G2.(fN)(c1, c2) + 1;
                        end
                    end
                end; end
                Sujeto_PLV_Matrices.(cond_actual).(bN).Between.Matrice = m_b;
            end
            if isempty(sujeto_table), sujeto_table = row_data; else, sujeto_table = [sujeto_table; row_data]; end
        end
    end
    
    if ~isempty(sujeto_table)
        if ~exist(targetDir, 'dir'), mkdir(targetDir); end
        save(fullfile(targetDir, [subjStr '_Conectividad_Matrices.mat']), 'Sujeto_PLV_Matrices');
        writetable(sujeto_table, fullfile(targetDir, [subjStr '_Conectividad_PLV.xlsx']));
    end
end

% 5. GENERAR INFORME N FINAL (ORGANIZADO POR HOJAS Y TIEMPOS)
informePath = fullfile(outPath, 'Informe_N_Sujetos_POR_GRUPO.xlsx');
if exist(informePath, 'file'), delete(informePath); end

cond_lista = fieldnames(mapTiempos);
for b = 1:numel(bNames)
    bN = bNames{b};
    hoja_banda = table();
    
    for t = 1:length(cond_lista)
        cLabel = cond_lista{t};
        fN = matlab.lang.makeValidName([cLabel '_' bN]);
        
        m1 = zeros(8); if isfield(n_report_G1, fN), m1 = n_report_G1.(fN); end
        m2 = zeros(8); if isfield(n_report_G2, fN), m2 = n_report_G2.(fN); end
        
        for r1 = 1:8, for r2 = (r1+1):8
            if (m1(r1,r2) + m2(r1,r2)) > 0
                nueva_fila = table({tiempos_etiquetas.(cLabel)}, {['ROI' num2str(r1)]}, {['ROI' num2str(r2)]}, ...
                    m1(r1,r2), m2(r1,r2), 'VariableNames', ...
                    {'Tiempo', 'ROI_A', 'ROI_B', 'N_Grupo_1_Placebo', 'N_Grupo_2_Real'});
                hoja_banda = [hoja_banda; nueva_fila];
            end
        end; end
    end
    if ~isempty(hoja_banda)
        writetable(hoja_banda, informePath, 'Sheet', bN);
    end
end
fprintf('\n>>> PROCESO COMPLETADO.\n');