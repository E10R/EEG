%% =============================================================
%  SOLO FOOOF POR ELECTRODO - VERSIÓN Q1 (CAR + PARÁMETROS ROBUSTOS)
%  Analiza 19 can. × 20 tiempos, exporta mediciones brutas
%% =============================================================
clear; close all; clc;

%% --- 1. EEGLAB ---
ruta_eeglab = 'C:\Users\evadi\OneDrive\Escritorio\eeglab_current';
if exist(ruta_eeglab,'dir')
    addpath(ruta_eeglab); eeglab nogui;
else
    error('No se encuentra EEGLAB.');
end

%% --- 2. PYTHON ---
pyenv('Version','C:\Users\evadi\AppData\Local\Programs\Python\Python310\python.exe');
currentDir = fileparts(mfilename('fullpath'));
if count(py.sys.path,currentDir)==0
    insert(py.sys.path,int32(0),currentDir);
end

%% --- 3. RUTAS ---
baseDir   = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\EEG_MATLAB';
outputDir = 'C:\Users\evadi\OneDrive\Escritorio\SUPERDEF_ESTE\Base_datos_def\ESTOS_DEFINITIVOS\SUPERDEFINITIVOS\DEF_APERIODICO_RE-REF';
if ~exist(outputDir,'dir'), mkdir(outputDir); end
masterFile = fullfile(outputDir,'Master_FOOOF_Bruto.xlsx');

%% --- 4. MAPEO DE ELECTRODOS ---
electrodeMap = struct( ...
    'Fp1',1,'Fp2',2,'F4',3,'F7',4,'F8',5,'Fz',6,...
    'C3',7,'C4',8,'Cz',9,'T3',10,'T4',11,'T5',12,'T6',13,...
    'P3',14,'P4',15,'Pz',16,'O1',17,'O2',18,'Oz',19);
electrodeNames = fieldnames(electrodeMap);
nElectrodes = numel(electrodeNames);

%% --- 5. PARAMETROS OPTIMIZADOS Q1 ---
bandas = {'Delta',[0.5 3.5];'Theta',[4 7.5];'Alfa',[8 12.5];
          'Beta',[13 29.5];'Gamma',[30 40]};
f_range = [0.5 40];

% [min_peak_width, max_peak_width, max_n_peaks, min_peak_height, peak_threshold]
settings_list = [1.0, 8.0, 5.0, 0.1, 2.0]; 

segmentDuration = 120;

%% --- 6. TABLA BASE ---
vars = {'Sujeto','Grupo','Electrodo','Electrodo_Num','Tiempo',...
        'Offset','Exponent','Peak_Delta','Peak_Theta','Peak_Alfa',...
        'Peak_Beta','Peak_Gamma','R2'};
types = {'double','double','string','double','double',...
          'double','double','double','double','double','double','double','double'};
masterData = table('Size',[0,numel(vars)],'VariableTypes',types,'VariableNames',vars);
errorLog = {};

%% --- 7. LOOP PRINCIPAL ---
fprintf('\n========== INICIO FOOOF BRUTO (CAR) ==========\n');
grupoNombre = {'Placebo','Real'};
for grupo=1:2
    grupoDir = fullfile(baseDir,num2str(grupo));
    if ~exist(grupoDir,'dir'), warning('Grupo %d no encontrado.',grupo); continue; end
    sujetos = dir(fullfile(grupoDir,'S*'));
    fprintf('\n📁 Grupo %d (%s): %d sujetos\n',grupo,grupoNombre{grupo},numel(sujetos));
    for s=1:numel(sujetos)
        subj = sujetos(s).name;
        subjNum = str2double(regexp(subj,'\d+','match','once'));
        subjPath = fullfile(grupoDir,subj);
        fprintf(' └─ %s:',subj);
        
        archivos = {
            sprintf('EEG_preprocess_%s_1_pre_OA.mat',subj),'pre',1;
            sprintf('EEG_preprocess_%s_3_tSMS.mat',subj),'stim',2;
            sprintf('EEG_preprocess_%s_4_post.mat',subj),'post',17};
        filas=0;
        
        for a=1:3
            fPath = fullfile(subjPath,archivos{a,1});
            if ~exist(fPath,'file'), continue; end
            load(fPath,'EEG'); EEG=eeg_checkset(EEG);
            
            % APLICAR CAR (Rereferencia) ANTES DE PROCESAR
            EEG = pop_reref(EEG, []); 
            
            fs=EEG.srate; tipo=archivos{a,2};
            switch tipo
                case 'pre'
                    segmentos={EEG}; tiempos=1;
                case 'stim'
                    segmentos=segmentarEEG(EEG,15,segmentDuration); tiempos=2:16;
                case 'post'
                    EEG.data=EEG.data(:,1:min(size(EEG.data,2),fs*8*60));
                    segmentos=segmentarEEG(EEG,4,segmentDuration); tiempos=17:20;
            end
            
            for seg=1:numel(segmentos)
                if isempty(segmentos{seg}), continue; end
                EEGs=segmentos{seg}; tiempo=tiempos(seg);
                for e=1:nElectrodes
                    elec=electrodeNames{e};
                    elecNum=electrodeMap.(elec);
                    ch=find(strcmpi({EEGs.chanlocs.labels},elec),1);
                    if isempty(ch), continue; end
                    
                    % PSD con Welch
                    [p_spec,f_axis]=pwelch(EEGs.data(ch,:),hanning(fs*2),fs,fs*2,fs);
                    
                    try
                        % Python Bridge
                        py.fooof_bridge.fit_fooof(py.numpy.array(f_axis'),...
                            py.numpy.array(p_spec'),f_range,py.list(settings_list));
                        load('bridge_result.mat','ap_params','peaks','r_squared');
                        
                        peakVals=extractPeaksByBand(peaks,bandas);
                        
                        row={subjNum,grupo,string(elec),elecNum,tiempo,...
                             ap_params(1),ap_params(2),...
                             peakVals.Delta,peakVals.Theta,peakVals.Alfa,...
                             peakVals.Beta,peakVals.Gamma,r_squared};
                        masterData=[masterData;row];
                        filas=filas+1;
                    catch MEfoo
                        errorLog{end+1}=sprintf('S%d_%s_t%d: %s',subjNum,elec,tiempo,MEfoo.message);
                    end
                end
            end
        end
        fprintf(' %d filas\n',filas);
    end
end

%% --- 8. GUARDAR RESULTADOS ---
writetable(masterData,masterFile);
if ~isempty(errorLog)
    fid=fopen(fullfile(outputDir,'ErrorLog.txt'),'w');
    fprintf(fid,'%s\n',errorLog{:});
    fclose(fid);
end
fprintf('\n✅ PROCESO COMPLETADO. Guardado en: %s\n',masterFile);

%% --- FUNCIONES AUXILIARES ---
function segmentos=segmentarEEG(EEG,nSeg,duracion)
    fs=EEG.srate; N=size(EEG.data,2);
    pts=floor(N/nSeg); segmentos=cell(1,nSeg);
    for i=1:nSeg
        ini=(i-1)*pts+1; fin=min(i*pts,N);
        if ini>=N, segmentos{i}=[]; continue; end
        try
            segmentos{i}=pop_select(EEG,'point',[ini fin]);
        catch, segmentos{i}=[]; end
    end
end

function p=extractPeaksByBand(peaks,bands)
    for b=1:size(bands,1)
        n=bands{b,1}; r=bands{b,2};
        if isempty(peaks)
            val=0;
        else
            idx=peaks(:,1)>=r(1)&peaks(:,1)<=r(2);
            if any(idx), val=max(peaks(idx,2)); else, val=0; end
        end
        p.(n)=val;
    end
end