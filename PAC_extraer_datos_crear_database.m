% =========================================================================
% SCRIPT: PAC MAESTRO - VERSIÓN ROBUSTA (IGNORA ARCHIVOS FALTANTES)
% =========================================================================
clear; close all; clc;

% --- CONFIGURACIÓN DE RUTAS ---
root_data = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\EEG_MATLAB';
output_root = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\DEF_PAC';
if ~exist(output_root, 'dir'), mkdir(output_root); end

% --- TUS 8 CLÚSTERES ESPECÍFICOS ---
clusters.Cluster_1 = {'Fp1', 'F7'};
clusters.Cluster_2 = {'Fp2', 'F4', 'F8'};
clusters.Cluster_3 = {'C3', 'P3'};
clusters.Cluster_4 = {'C4', 'P4'};
clusters.Cluster_5 = {'Fz', 'Cz', 'Pz'};
clusters.Cluster_6 = {'T3', 'T5'};
clusters.Cluster_7 = {'T4', 'T6'};
clusters.Cluster_8 = {'O1', 'Oz', 'O2'};

cNames = fieldnames(clusters);

% --- PARÁMETROS PAC (LÍMITES FIJOS) ---
band_limits.Theta = [4 7.5]; 
band_limits.Alpha = [8 12.5];
band_limits.Beta  = [13 29.5]; 
band_limits.Gamma = [30 40];

phase_bands = {'Theta', 'Alpha'}; 
amp_bands   = {'Beta', 'Gamma'};
num_bins    = 18; 
phase_bins  = linspace(-pi, pi, num_bins+1);

% Acumulador Principal
master_data = table();

% --- BUCLE DE GRUPOS (1 y 2) ---
for g = 1:2
    group_path = fullfile(root_data, num2str(g));
    if ~exist(group_path, 'dir'), continue; end
    subj_dirs = dir(fullfile(group_path, 'S*'));
    
    for s = 1:numel(subj_dirs)
        subj_folder = subj_dirs(s).name;
        subj_id = str2double(regexp(subj_folder, '\d+', 'match'));
        subj_path = fullfile(group_path, subj_folder);
        
        % Definición de archivos
        cond_files = {
            sprintf('EEG_preprocess_%s_1_pre_OA.mat', subj_folder), 1;
            sprintf('EEG_preprocess_%s_3_tSMS.mat', subj_folder), 2; 
            sprintf('EEG_preprocess_%s_4_post.mat', subj_folder), 5
        };
        
        for f = 1:size(cond_files, 1)
            file_name = cond_files{f,1};
            base_time_code = cond_files{f,2};
            full_file_path = fullfile(subj_path, file_name);
            
            % --- SI EL ARCHIVO NO EXISTE, SALTAMOS ---
            if ~exist(full_file_path, 'file')
                fprintf('AVISO: Falta archivo %s. Saltando...\n', file_name);
                continue; 
            end
            
            try
                fprintf('Procesando Sujeto %d | Archivo: %s\n', subj_id, file_name);
                % Carga limpia
                data_in = load(full_file_path);
                if isfield(data_in, 'cleanEEG'), EEG = data_in.cleanEEG; 
                elseif isfield(data_in, 'EEG'), EEG = data_in.EEG; 
                else, error('No se encontró estructura EEG'); end
                
                % Guardar etiquetas de canales antes de procesar
                chan_labels = {EEG.chanlocs.labels};
                
                n_segments = 1;
                if base_time_code == 2, n_segments = 3; end
                L = size(EEG.data, 2);
                seg_len = floor(L / n_segments);
                
                for seg = 1:n_segments
                    current_time_code = base_time_code + (n_segments > 1)*(seg - 1);
                    
                    data_chunk = double(EEG.data(:, (1+(seg-1)*seg_len):(seg*seg_len)));
                    
                    row_res = struct();
                    row_res.Sujeto = subj_id;
                    row_res.Grupo = g;
                    row_res.Tiempo = current_time_code;
                    
                    % Filtrado
                    filt_an = struct();
                    all_b = unique([phase_bands, amp_bands]);
                    for b = 1:numel(all_b)
                        bn = all_b{b};
                        tmp = EEG; tmp.data = data_chunk; tmp.pnts = size(data_chunk,2);
                        tmp_f = pop_eegfiltnew(tmp, band_limits.(bn)(1), band_limits.(bn)(2), [], 0, [], 0);
                        filt_an.(bn) = hilbert(tmp_f.data')';
                    end
                    
                    % Cálculo PAC
                    for c = 1:numel(cNames)
                        % Usamos chan_labels que guardamos antes
                        ch_idx = find(ismember(chan_labels, clusters.(cNames{c})));
                        if isempty(ch_idx), continue; end
                        
                        for pb = 1:numel(phase_bands)
                            for ab = 1:numel(amp_bands)
                                p_band = phase_bands{pb}; a_band = amp_bands{ab};
                                p_vals = angle(filt_an.(p_band)(ch_idx, :));
                                a_vals = abs(filt_an.(a_band)(ch_idx, :));
                                
                                mi_list = zeros(length(ch_idx), 1);
                                for e = 1:length(ch_idx)
                                    m_amp = zeros(1, num_bins);
                                    for bi = 1:num_bins
                                        m = p_vals(e,:) >= phase_bins(bi) & p_vals(e,:) < phase_bins(bi+1);
                                        if any(m), m_amp(bi) = mean(a_vals(e, m)); end
                                    end
                                    if sum(m_amp) > 0
                                        P = m_amp/sum(m_amp); P(P==0) = eps;
                                        mi_list(e) = (log(num_bins) - (-sum(P.*log(P)))) / log(num_bins);
                                    end
                                end
                                col = sprintf('%s_%s_mod_%s', cNames{c}, p_band, a_band);
                                row_res.(col) = mean(mi_list);
                            end
                        end
                    end
                    master_data = [master_data; struct2table(row_res)];
                end
                clear EEG data_in filt_an data_chunk;
            catch ME
                fprintf('ERROR en Sujeto %d, archivo %s: %s\n', subj_id, file_name, ME.message);
                continue; % Salta al siguiente archivo si algo falla internamente
            end
        end
    end
end

% --- GUARDADO FINAL ---
writetable(master_data, fullfile(output_root, 'Master_PAC_Fijo_Final.xlsx'));
fprintf('\n=== PROCESO COMPLETADO: Tabla generada ===\n');