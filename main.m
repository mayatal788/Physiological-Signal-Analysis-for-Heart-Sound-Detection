clear all;
clc;

% load data 
all_data_3 = readmatrix('ecg_3.txt');
data_scope = all_data_3(: , 1);
data_ecg = all_data_3(: , 2);
fs = 500;

% Define segments for Rest and Exercise
data_rest_scope = data_scope(1:20 * fs);
data_exercise_scope = data_scope(21 * fs:31 * fs);

% --- Execute Detection ---
[s1_idx_rest, s2_idx_rest] = detect_heart_sounds(data_rest_scope, fs);
[s1_idx_exercise, s2_idx_exercise] = detect_heart_sounds(data_exercise_scope, fs);

%% Heart sound identification algorithm function
function [s1_idx, s2_idx] = detect_heart_sounds(signal, fs)
    [C, L] = wavedec(signal, 5, 'db6'); %wavelet denoising type db6, level 5
    
    % Reconstructing signal using d2 and d3
    d3 = wrcoef('d', C, L, 'db6', 3); 
    d2 = wrcoef('d', C, L, 'db6', 2);
    %d4 = wrcoef('d', C, L, 'db6', 4);
    %d5 = wrcoef('d', C, L, 'db6', 5);
    signal_preprocessed = d2 + d3; % chosen as it gives the best results 
   
    N = 32;    % Window size
    hop = 16;  % Step size
    sig_norm = signal_preprocessed / max(abs(signal_preprocessed)); %normalize to [-1,1] range 
    
    E_hs = [];
    for i = 1:hop:(length(sig_norm)-N) %extraction of heart sound envelope (Shannon Energy)  
        window = sig_norm(i:i+N-1);
        energy = -mean(abs(window).^3 .* log(abs(window).^3 + eps)); %3rd-order Shannon Energy formula 
        E_hs = [E_hs, energy];
    end
   
    P_ha = (E_hs - mean(E_hs)) / std(E_hs);  % Standardization
    
    Th = prctile(P_ha, 55); %thresholding of sound - adaptive and amplitude based ( 55 percentile )
    time_gates = P_ha > Th; % creating time gates of s1 and s2

   %% noise removal 
    % Removing gates that are too short to be physiological heart sounds
    for repeat = 1:2 
        diff_g = diff([0 time_gates 0]);
        starts = find(diff_g == 1);
        ends   = find(diff_g == -1) - 1;
        
        % Minimum gate duration constraint 
        min_gate_ms  = 90 - (repeat-1)*10; % starting with min 90 ms then allowing 80 ms 
        min_gate_len = ceil((min_gate_ms/1000)*fs/hop);  % convert minimum gate duration from milliseconds to number of envelope points
        
        for k = 1:length(starts) % deleting time gates shorter than the minimum lenght 
            if (ends(k) - starts(k) + 1) < min_gate_len
                time_gates(starts(k):ends(k)) = 0;
            end
        end
    end
    
    
    % Merge gates separated by small gaps 
    diff_g = diff([0 time_gates 0]);
    starts = find(diff_g == 1);
    ends   = find(diff_g == -1) - 1;
    
    merge_gap_ms  = 15; % merging gates with gaps smaller than 15 ms                              
    merge_gap_len = ceil((merge_gap_ms/1000)*fs/hop); 
    
    for k = 1:length(starts)-1
        gap_len = starts(k+1) - ends(k) - 1;
        if gap_len > 0 && gap_len <= merge_gap_len
            time_gates(ends(k)+1 : starts(k+1)-1) = 1;
        end
    end
    
    %Peak Marking within Gates 
    diff_gates = diff([0, time_gates, 0]);
    starts = find(diff_gates == 1);   
    ends   = find(diff_gates == -1);  
    
    cand_env_idx = zeros(1, numel(starts));
    for k = 1:numel(starts)
        seg = P_ha(starts(k):ends(k)-1);
        [~, localMaxIdx] = max(seg);
        cand_env_idx(k) = starts(k) + localMaxIdx - 1; 
    end
    
    % Map envelope indices back to signal sample indices
    cand_sample_idx = (cand_env_idx - 1) * hop + round(N/2); 
    cand_sample_idx = cand_sample_idx(cand_sample_idx >= 1 & cand_sample_idx <= length(sig_norm));

    % Refinement (Remove too-close peaks)
    min_sep_ms = 80 ; % minimum separetion between peaks 
    min_sep_samples = round((min_sep_ms/1000) * fs);
    cand_sample_idx = sort(unique(cand_sample_idx));
    keep = true(size(cand_sample_idx));
    for k = 2:numel(cand_sample_idx)
        if cand_sample_idx(k) - cand_sample_idx(k-1) < min_sep_samples
            env_k   = min(max(round(cand_sample_idx(k)/hop)+1, 1), length(P_ha));
            env_km1 = min(max(round(cand_sample_idx(k-1)/hop)+1, 1), length(P_ha));
            if P_ha(env_k) > P_ha(env_km1), keep(k-1) = false; else, keep(k) = false; end % keeping the peak with a higher envelope value 
        end
    end
    cand_sample_idx = cand_sample_idx(keep); % filtering the peaks list and keeping only the ones we wanted 

    % Physiological Classification - verifies that a sufficient number of candidate heart sound events were detected
    %if numel(cand_sample_idx) < 3
        %s1_idx = []; s2_idx = []; return;
    %end
    
    % Calculate intervals between consecutive events
    diffs = diff(cand_sample_idx) / fs; 
    labels = strings(1, numel(cand_sample_idx));
   
    [~, k_long] = max(diffs); % finding the longest time space idx between peaks ( the diastole )
    labels(k_long) = "S2";% peak before diastole is s2
    labels(k_long+1) = "S1";% peak after diastole is s1 
    for i = k_long+2:numel(labels) % moving forward - if the former peak was s1 the current is s2 and vice versa 
        if labels(i-1) == "S1", labels(i) = "S2"; else, labels(i) = "S1"; end
    end
    for i = k_long-1:-1:1 % moving backwords - if next peak is s1 current is s2 and vice versa 
        if labels(i+1) == "S1", labels(i) = "S2"; else, labels(i) = "S1"; end
    end


    s1_idx = cand_sample_idx(labels == "S1");
    s2_idx = cand_sample_idx(labels == "S2");

    % Plotting
    figure;
    subplot(2,1,1);
    plot(sig_norm); hold on;
    plot(s1_idx, sig_norm(round(s1_idx)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    plot(s2_idx, sig_norm(round(s2_idx)), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    title('Preprocessed Signal with Classified S1 (Green) and S2 (Red)');
    
    subplot(2,1,2);
    plot(P_ha); hold on;
    plot(time_gates * max(P_ha), 'r--');
    title('Shannon Energy Envelope and Time Gates');
    legend('Envelope', 'Time Gates');
end


%% --- Plot ECG + Heart Sounds together (Rest vs Post Exercise) ---

% Create matching ECG segments (same indices as scope segments)
data_rest_ecg     = data_ecg(1:20*fs);
data_exercise_ecg = data_ecg(21*fs:31*fs);

% Time vectors (seconds) for each segment
t_rest     = (0:length(data_rest_scope)-1)/fs;
t_exercise = (0:length(data_exercise_scope)-1)/fs;

% Normalize/scale
rest_hs = data_rest_scope;
ex_hs   = data_exercise_scope;

% Convert detected sample indices to time (sec)
tS1_rest = (s1_idx_rest-1)/fs;
tS2_rest = (s2_idx_rest-1)/fs;

tS1_ex   = (s1_idx_exercise-1)/fs;
tS2_ex   = (s2_idx_exercise-1)/fs;

% -------- REST FIGURE --------
figure('Name','REST: ECG + Heart Sounds','Color','w');

subplot(2,1,1);
plot(t_rest, data_rest_ecg, 'k'); grid on; hold on;
title('REST - ECG ');
xlabel('Time [s]'); ylabel('ECG [a.u.]');

% Mark S1/S2 on ECG
for k = 1:numel(tS1_rest), xline(tS1_rest(k),'g--','S1','LabelVerticalAlignment','bottom'); end
for k = 1:numel(tS2_rest), xline(tS2_rest(k),'r--','S2','LabelVerticalAlignment','bottom'); end
legend('ECG','Location','best');

subplot(2,1,2);
plot(t_rest, rest_hs, 'b'); grid on; hold on;
title('REST - Heart Sounds');
xlabel('Time [s]'); ylabel('HS [a.u.]');

for k = 1:numel(tS1_rest), xline(tS1_rest(k),'g--','S1','LabelVerticalAlignment','bottom'); end
for k = 1:numel(tS2_rest), xline(tS2_rest(k),'r--','S2','LabelVerticalAlignment','bottom'); end
legend('Heart Sounds','Location','best');

% -------- EXERCISE FIGURE --------
figure('Name','EXERCISE: ECG + Heart Sounds','Color','w');

subplot(2,1,1);
plot(t_exercise, data_exercise_ecg, 'k'); grid on; hold on;
title('POST EXERCISE - ECG ');
xlabel('Time [s]'); ylabel('ECG [a.u.]');

for k = 1:numel(tS1_ex), xline(tS1_ex(k),'g--','S1','LabelVerticalAlignment','bottom'); end
for k = 1:numel(tS2_ex), xline(tS2_ex(k),'r--','S2','LabelVerticalAlignment','bottom'); end
legend('ECG','Location','best');

subplot(2,1,2);
plot(t_exercise, ex_hs, 'b'); grid on; hold on;
title('POST EXERCISE - Heart Sounds ');
xlabel('Time [s]'); ylabel('HS [a.u.]');

for k = 1:numel(tS1_ex), xline(tS1_ex(k),'g--','S1','LabelVerticalAlignment','bottom'); end
for k = 1:numel(tS2_ex), xline(tS2_ex(k),'r--','S2','LabelVerticalAlignment','bottom'); end
legend('Heart Sounds','Location','best');
