%START EEGLAB
clear all
close all
clc


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
eeglab;

Fs=512;

folder_path = 'C:\Users\jimen\OneDrive\Escritorio\HNPT\Datos_Eva\';
output_path = 'C:\Users\jimen\OneDrive\Escritorio\HNPT\Datos_Eva\'; % Update with the path to your output folder

columns=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"];

%LOAD EEG

letra = ["A" "B"];
% letra = ["A"];
test = ["OA" "OC" "tSMS" "POST"];

%%
for expIdx= 1:length(letra)
    modo = char(letra(expIdx));
    sujetos = {dir(strcat(folder_path,modo)).name};
    sujetos = sujetos(3:end);
    % for subIdx = 1:length(sujetos)
    for subIdx = 2:3    
        paciente = char(sujetos(subIdx));
        output_subject_path = fullfile(output_path,modo, paciente);
        files = dir(fullfile(strcat(folder_path,modo, "\",paciente), '*.edf'));
        column=1;
        for setIdx = 1:length(files)
            currently = msgbox(strcat(modo, " Sujeto ", paciente, " Archivo ", files(setIdx).name), "Subject window");

            eeg_filepath=strcat(files(setIdx).folder,'\',files(setIdx).name);
            EEG.etc.eeglabvers = '2021.1'; % this tracks which version of EEGLAB is being used, you may ignore it
            EEG = pop_biosig(eeg_filepath);
            EEG.setname=strcat('EDF_',files(setIdx).name(1:end-4));
            eeg_12=EEG;
            disp(eeg_filepath)
            channelToRemove = {'ECG1+'}; % You can specify the channel by name or number, e.g., {'FCz'} or {'1'}
            EEG = pop_select(EEG, 'nochannel', channelToRemove);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG= preprocess_EEG(EEG,strcat(files(setIdx).folder,'\'),'C:\Users\jimen\OneDrive\Escritorio\HNPT\EEG_dolor\EEG\eeglab2021.1\plugins\dipfit\standard_BEM\elec\standard_1020.elc');
            save(strcat(files(setIdx).folder,'\','EEG_preprocess_',files(setIdx).name(1:end-4),'.mat'),'EEG')
            EEG = pop_saveset( EEG, 'filename',strcat(files(setIdx).name(1:end-4),'.set'),'filepath',strcat(files(setIdx).folder,'\'));


            asr_latency = [EEG.event.latency].';
            asr_duration = [EEG.event.duration].';
            length_eeg_before_asr= length(EEG.etc.clean_sample_mask); %number of samples in EEG before applying ASR
            length_eeg_after_asr= EEG.pnts; %number of samples in EEG after applying ASR
            samples_asr=1:length_eeg_before_asr;

            samples_asr=samples_asr.*EEG.etc.clean_sample_mask;

            j=0;
            for i=1:length(samples_asr)
                if samples_asr(i)~=0
                    j=j+1;
                    samples_asr(i)=j;
                end
            end


            num_channels = EEG.nbchan; % Number of EEG channels
            segment_length = 2 * 60; % Length of each segment in seconds

            positions_segments_before_asr=1:(segment_length*Fs):length_eeg_before_asr;
            positions_segments_after_asr = zeros(1,length(positions_segments_before_asr));
            for n=1:num_channels
                % Loop through each position
                for i = 1:length(positions_segments_before_asr)
                    pos = positions_segments_before_asr(i);
                    % Check if the value at the position is 0
                    if samples_asr(pos) == 0
                        % Search for the closest value with a smaller index
                        for j = pos-1:-1:1
                            if samples_asr(j) ~= 0
                                positions_segments_after_asr(i) = samples_asr(j);
                                break;
                            end
                        end
                    else
                        positions_segments_after_asr(i) = samples_asr(pos);
                    end
                end
            end
            % Preallocate a cell array to store the PSD and frequency vectors
            num_positions=length(positions_segments_after_asr);
            psd_values = cell(num_positions, 1);
            psd_matrix = [];
            frequencies = cell(num_positions, 1);
            bands = cell(num_positions, 1);


            if contains(files(setIdx).name, 'OA') || contains(files(setIdx).name, 'OC')
                length_interval=ones(1,EEG.nbchan).*length(EEG.data);

                [psd, freq, band_power_table] = calculate_eeg_psd(EEG.data,{EEG.chanlocs.labels},length_interval, Fs);
                plot(psd)
                for chIdx = 1: EEG.nbchan
                    % Create table for Excel
                    psd_table = array2table(cell2mat(psd(chIdx)), 'VariableNames', {files(setIdx).name(end-5:end-4)});
                    if contains(files(setIdx).name, 'OA')
                        % Create the new column with the label "Frequencies" and the associated data
                        freq_table = array2table(freq, 'VariableNames', {'Frequencies'});

                        % Concatenate the new row at the beginning of the table
                        psd_table = [freq_table, psd_table];
                    end
                    % Sheet name for the channel
                    sheet_name = EEG.chanlocs(chIdx).labels;

                    % Write to Excel
                    output_file = fullfile(output_subject_path, sprintf('%s.xlsx', paciente));
                    writetable(psd_table, output_file, 'Sheet', sheet_name,'Range',columns(column)+string(2) );
                end
                if contains(files(setIdx).name, 'OA')
                    column=column+2;
                else
                    column=column+1;
                end
            else

                for pos = 1:num_positions-1
                    if pos == 1
                        length_interval=ones(1,EEG.nbchan).*length(EEG.data(:,1:positions_segments_after_asr(pos+1)));
                        [psd, freq, band_power_table] = calculate_eeg_psd(EEG.data(:,1:positions_segments_after_asr(pos+1)),{EEG.chanlocs.labels}, length_interval,Fs);
                    elseif pos == num_positions-1
                        length_interval=ones(1,EEG.nbchan).*length(EEG.data(:,positions_segments_after_asr(pos):end));
                        [psd, freq, band_power_table] = calculate_eeg_psd(EEG.data(:,positions_segments_after_asr(pos):end),{EEG.chanlocs.labels}, length_interval, Fs);

                    else
                        length_interval=ones(1,EEG.nbchan).*length(EEG.data(:,positions_segments_after_asr(pos):positions_segments_after_asr(pos+1)));
                        [psd, freq, band_power_table] = calculate_eeg_psd(EEG.data(:,positions_segments_after_asr(pos):positions_segments_after_asr(pos+1)),{EEG.chanlocs.labels}, length_interval,Fs);

                    end
                    psd_matrix = [psd_matrix, psd];
                    psd_values{pos} = psd;
                    frequencies{pos} = freq;
                    bands{pos} = band_power_table;
                    segments=strcat(files(setIdx).name(end-7:end-4),'_', string(1:num_positions));

                    for chIdx = 1: EEG.nbchan
                        % Initialize structure to store PSD data for all channels
                        all_psd_data = struct();

                        % Create table for Excel
                        psd_table = array2table(cell2mat(psd(chIdx)), 'VariableNames', segments(pos));


                        % Sheet name for the channel
                        sheet_name = EEG.chanlocs(chIdx).labels;

                        % Write to Excel
                        output_file = fullfile(output_subject_path, sprintf('%s.xlsx', paciente));
                        writetable(psd_table, output_file, 'Sheet', sheet_name,'Range', columns(column)+string(2) );
                    end
                    column=column+1;
                end

            end
            close all
            close(currently);
        end


    end
end


%% FUNCTIONS

function [EEG] = preprocess_EEG(EEG, pathtosave,standard_eeg)
%START EEGLAB
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% eeglab;
% if endsWith(filename,'.vhdr')
%     EEG = pop_loadbv(filename); %carga el eeg vhdr
% elseif endsWith(filename,'.edf')
%     EEG = pop_readedf(filename); %carga el eeg vhdr
% end
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); %saves the dataset in EEGLAB memory
EEG = eeg_checkset( EEG ); %check the dataset consistency

%CHANNELS

EEG = pop_chanedit(EEG, 'lookup', standard_eeg);
EEG = eeg_checkset( EEG );

%FILTER 40-05

EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'plotfreqz',1);
EEG = pop_eegfiltnew(EEG, 'hicutoff',40,'plotfreqz',1);
EEG.comments = pop_comments(EEG.comments,'','Filters applied, locutoff 0.5, hicutoff 40',1);
EEG = eeg_checkset(EEG);

%GUARDO ESTE DATASET PARA USARLO LUEGO

direc1 = pathtosave;
archivo1 = char('EEG_PRE');
if isfile(strcat(direc1,archivo1))
    delete(strcat(direc1,archivo1,'.set'))
    EEG = pop_saveset( EEG,archivo1,direc1);
else
    EEG = pop_saveset( EEG,archivo1,direc1);
    EEG = eeg_checkset( EEG );
end

%ASR

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','channels_ignore',{'ECG1+'},'WindowCriterionTolerances',[-Inf 7] );

EEG.comments = pop_comments(EEG.comments,'','ASR applied 1',1);
EEG = eeg_checkset(EEG);

%CREO ARRAY CHANNELS PARA INTERPOLAR
PostASR = [];
for u = 1:length(EEG.chanlocs)
    PostASR(u) = EEG.chanlocs(u).urchan;
end

ChansPostASR = [];
j = 1;
for h = 1:length({EEG.chanlocs.labels})
    if j > length(PostASR)
        ChansPostASR(h) = 0;
        break;
    end
    if h == PostASR(j)
        ChansPostASR(h) = PostASR(j);
        %ChansPostASR(h) = h;
        j = j + 1;
    else
        ChansPostASR(h) = 0;
    end
end

Inter = find(ChansPostASR==0);
InterNO = find(ChansPostASR~=0);

%CARGO EL DATASET SIN ELIMINAR CANALES
EEG = pop_loadset('EEG_PRE.set',direc1);

%ASR SIN ELIMINAR CANALES
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
% EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
EEG.comments = pop_comments(EEG.comments,'','ASR applied 2',1);
EEG = eeg_checkset(EEG);

%ICA
InterNO = find(ChansPostASR~=0);
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on','chanind',InterNO);
EEG.comments = pop_comments(EEG.comments,'','Runica applied',1);
EEG = eeg_checkset( EEG );

%ELIMINAR COMPONENTES ICA
EEG = pop_iclabel(EEG, 'default');
IcaRem = zeros([1 length(EEG.etc.ic_classification.ICLabel.classifications)]);
for k = 1:length(IcaRem)
    %                 if (EEG.etc.ic_classification.ICLabel.classifications(k,1)<0.04) %brain
    %                     IcaRem(k) = 1;
    if (EEG.etc.ic_classification.ICLabel.classifications(k,2)>0.1) %muscle
        IcaRem(k) = 1;
    elseif (EEG.etc.ic_classification.ICLabel.classifications(k,3)>0.1) %eye
        IcaRem(k) = 1;
    elseif (EEG.etc.ic_classification.ICLabel.classifications(k,4)>0.1) %heart
        IcaRem(k) = 1;
    elseif (EEG.etc.ic_classification.ICLabel.classifications(k,5)>0.1) %line noise
        IcaRem(k) = 1;
    elseif (EEG.etc.ic_classification.ICLabel.classifications(k,6)>0.1) %channel noise
        IcaRem(k) = 1;
    end
end
IcaRemBad = find(IcaRem==1);

EEG = pop_subcomp( EEG, IcaRemBad, 0);
EEG.comments = pop_comments(EEG.comments,'','ICA components removed',1);
EEG = eeg_checkset( EEG );

%INTERPOLACIÓN

EEG = pop_interp(EEG, Inter,'spherical');
EEG.comments = pop_comments(EEG.comments,'','Interpolation applied',1);
EEG = eeg_checkset( EEG );


end

function [psd, freq, band_power_table] = calculate_eeg_psd(EEG, channel_names, length_interval, fs)
% signal: EEG signal with n channels (each column represents a channel)
% fs: sampling frequency

% Define frequency bands
bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
band_limits = [0.5 4; 4 8; 8 12; 12 30; 30 100];

% Compute PSD for each channel
n_channels = length(EEG(:,1));
psd = cell(1, n_channels);
freq = cell(1, n_channels);
band_power_table = cell(n_channels * length(bands), 3); % Initialize table

row = 1; % Row index for band_power_table
for i = 1:n_channels
    % Example parameters
    NFFT = fs / 0.5; % FFT length to get 0.5 Hz resolution
    window = hamming(NFFT); % Window function
    noverlap = 0; % Overlap between segments
    [p, freq] = pwelch(EEG.data(i,:), window, noverlap, NFFT, fs);
    psd{i}= [p;length_interval(i)];

    for j = 1:length(bands)
        % Get band limits
        band_start = band_limits(j, 1);
        band_end = band_limits(j, 2);

        % Compute band power
        band_power = trapz(freq(freq >= band_start & freq <= band_end), psd{i}(freq >= band_start & freq <= band_end));

        % Store in table
        band_power_table{row, 1} = channel_names{i};
        band_power_table{row, 2} = bands{j};
        band_power_table{row, 3} = band_power;
        row = row + 1;
    end
end

freq=[freq;0];
% Convert cell array to table for better readability
band_power_table = cell2table(band_power_table, 'VariableNames', {'Channel', 'Band', 'Power'});
end

