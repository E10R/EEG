% 1. Rutas de entrada y salida
basePath = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\EEG_MATLAB';
outPath  = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\DEF_CONECTIVIDAD';

if ~exist(outPath, 'dir'), mkdir(outPath); end

% 2. Configuración de Grupos
grupos = zeros(42,1);
grupos([2,4,6,8,11,12,14,19,20,21,23,26,28,29,33,34,35,37,39,42]) = 1;
grupos([1,3,5,7,9,10,13,15,17,18,22,24,25,27,30,31,32,36,38,40,41]) = 2;

% 3. Definición de Clusters y Bandas
clusters.Cluster1 = {'Fp1', 'F7'};
clusters.Cluster2 = {'Fp2', 'F4', 'F8'};
clusters.Cluster3 = {'C3', 'P3'};
clusters.Cluster4 = {'C4', 'P4'};
clusters.Cluster5 = {'Fz', 'Cz', 'Pz'};
clusters.Cluster6 = {'T3', 'T5'};
clusters.Cluster7 = {'T4', 'T6'};
clusters.Cluster8 = {'O1', 'O2', 'Oz'};

bandas = struct('Delta',[0.5 3.5], 'Theta',[4 7.5], 'Alpha',[8 12.5], 'Beta',[13 29.5], 'Gamma',[30 40 ]);
bNames = fieldnames(bandas);
cNames = fieldnames(clusters);

condiciones_archivos = {'1_pre_OA', '3_tSMS', '4_post'};
condiciones_nombres  = {'Pre', 'tSMS', 'Post'};

% 4. Bucle Principal
for s = 1:42
    if grupos(s) == 0, continue; end
    
    subjStr = sprintf('S%02d', s);
    grupoStr = num2str(grupos(s));
    subjFullPath = fullfile(basePath, grupoStr, subjStr);
    targetGroupDir = fullfile(outPath, grupoStr);
    
    if ~exist(targetGroupDir, 'dir'), mkdir(targetGroupDir); end
    if ~exist(subjFullPath, 'dir'), continue; end

    Sujeto_PLV_Matrices = struct(); 
    sujeto_table = table(); 
    archivos_procesados = 0;

    for cond = 1:length(condiciones_archivos)
        fileName = sprintf('EEG_preprocess_%s_%s.mat', subjStr, condiciones_archivos{cond});
        fileFull = fullfile(subjFullPath, fileName);
        
        if exist(fileFull, 'file')
            try
                temp_data = load(fileFull);
                varnames = fieldnames(temp_data);
                cleanEEG = temp_data.(varnames{1}); 

                if ~isstruct(cleanEEG) || ~isfield(cleanEEG, 'chanlocs')
                    fprintf('!! Error de formato en %s. Saltando...\n', fileName);
                    continue;
                end

                % --- LÓGICA DE SEGMENTACIÓN ---
                % Si es tSMS, dividimos en 3. Si no, 1 solo bloque.
                is_tSMS = strcmp(condiciones_nombres{cond}, 'tSMS');
                if is_tSMS
                    n_segs = 3;
                    longitud_total = size(cleanEEG.data, 2);
                    pnts_seg = floor(longitud_total / 3);
                else
                    n_segs = 1;
                end

                for i_seg = 1:n_segs
                    % Definir nombre de la sub-condición
                    if is_tSMS
                        cond_actual = sprintf('tSMS_%d', i_seg);
                        idx_range = ( (i_seg-1)*pnts_seg + 1 ) : (i_seg*pnts_seg);
                        current_data = cleanEEG.data(:, idx_range);
                    else
                        cond_actual = condiciones_nombres{cond};
                        current_data = cleanEEG.data;
                    end

                    fprintf('Procesando %s (Grupo %s) -> %s\n', subjStr, grupoStr, cond_actual);
                    row_data = table({subjStr}, {cond_actual}, 'VariableNames', {'Sujeto', 'Condicion'});

                    % Procesamiento por Bandas
                    for b = 1:numel(bNames)
                        lim = bandas.(bNames{b});
                        
                        % Filtrado (Usando current_data del segmento)
                        % Creamos un dummy temporal para usar pop_eegfiltnew si es necesario
                        dummyEEG = eeg_emptyset();
                        dummyEEG.data = current_data;
                        dummyEEG.srate = cleanEEG.srate;
                        dummyEEG.chanlocs = cleanEEG.chanlocs;
                        dummyEEG.nbchan = size(current_data, 1);
                        dummyEEG.pnts = size(current_data, 2);
                        dummyEEG.trials = 1;
                        dummyEEG = eeg_checkset(dummyEEG);
                        
                        tempEEG = pop_eegfiltnew(dummyEEG, lim(1), lim(2));
                        fasi_tutti = angle(hilbert(tempEEG.data.'))';
                        segnali_medi = zeros(numel(cNames), size(tempEEG.data, 2));

                        % PLV Within
                        for c = 1:numel(cNames)
                            c_name = cNames{c};
                            ch_labels = {cleanEEG.chanlocs.labels};
                            ch_idx = find(ismember(ch_labels, clusters.(c_name)));
                            
                            if length(ch_idx) > 1
                                m_w = ones(length(ch_idx), length(ch_idx)); 
                                for x = 1:length(ch_idx)
                                    for y = (x+1):length(ch_idx)
                                        df = fasi_tutti(ch_idx(x), :) - fasi_tutti(ch_idx(y), :);
                                        pv = abs(mean(exp(1i * df)));
                                        m_w(x,y) = pv; m_w(y,x) = pv;
                                    end
                                end
                                Sujeto_PLV_Matrices.(cond_actual).(bNames{b}).Within.(c_name) = m_w;
                                mask = ~eye(size(m_w));
                                row_data.(['Within_' bNames{b} '_' c_name]) = mean(m_w(mask));
                            end
                            if ~isempty(ch_idx)
                                segnali_medi(c, :) = mean(tempEEG.data(ch_idx, :), 1);
                            end
                        end

                        % PLV Between
                        fasi_cl = angle(hilbert(segnali_medi.'))';
                        m_b = ones(numel(cNames), numel(cNames));
                        for c1 = 1:numel(cNames)
                            for c2 = (c1+1):numel(cNames)
                                df_b = fasi_cl(c1, :) - fasi_cl(c2, :);
                                val_b = abs(mean(exp(1i * df_b)));
                                m_b(c1, c2) = val_b; m_b(c2, c1) = val_b;
                            end
                        end
                        Sujeto_PLV_Matrices.(cond_actual).(bNames{b}).Between.Matrice = m_b;
                    end
                    
                    if isempty(sujeto_table), sujeto_table = row_data; else, sujeto_table = [sujeto_table; row_data]; end
                end
                archivos_procesados = archivos_procesados + 1;
            catch ME
                fprintf('!! Error en %s: %s\n', fileName, ME.message);
            end
        end
    end

    if archivos_procesados > 0
        save(fullfile(targetGroupDir, [subjStr '_Conectividad_Matrices.mat']), 'Sujeto_PLV_Matrices');
        writetable(sujeto_table, fullfile(targetGroupDir, [subjStr '_Conectividad_PLV.xlsx']));
    end
end
fprintf('\n¡Todo listo! tSMS dividido en 3 partes.\n');