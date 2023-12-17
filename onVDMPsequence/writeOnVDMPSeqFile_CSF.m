%% onVDMP_3T
% Create a sequence file for an onVDMP protocol according to:
% 

% Jianpan Huang 2023
% jianpanhuang@outlook.com

% author name for sequence file
author = 'Jianpan Huang';

%% get correct path
script_fp = fullfile(getPulseqCESTRootDir, 'seq-generation');

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 2*60           ; % 60 bp pairs (120 pulses)
defs.tp            = 7.5e-3         ; % pulse duration [s]
defs.td            = eps            ; % interpulse delay [s], t_delay = 0 ms
defs.Trec          = 2.4            ; % recovery time [s]
defs.Trec_M0       = 12             ; % recovery time before M0 [s]
defs.M0_offset     = -300           ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [0]; % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
                         defs.td ;  % saturation time [s]
defs.B0            = 3                ; % B0 [T]
defs.seq_id_string = 'onVDMP_CSF_3T'; % unique seq id

defs.B1pa          = 3.1;  % mean sat pulse b1 [uT]
defs.spoiling      = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = fullfile(script_fp, strcat(defs.seq_id_string,'.seq')); % filename

%% scanner limits 
% see pulseq doc for more ino
% init sequence
scanner_lims = getScannerLimits();
seq = SequenceSBB(scanner_lims);
gamma_hz  = seq.sys.gamma*1e-6;                  % for H [Hz/uT]
%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = defs.B1pa*gamma_rad*defs.tp; % flip angle of sat pulse

% create pulseq saturation pulse object
% satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', defs.tp, 'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse = makeSaturationPulseFromCWPE('block', defs.B1pa, defs.tp, defs.td, scanner_lims);
% resample pulse for reduced file size and io time
satPulse      = resamplePulseForRLE(satPulse, 1000); 

[B1cwpe,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,0,gamma_hz);
defs.B1cwpe = B1cwpe;


%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*gamma_hz*defs.B0;
phase_cycle = [0, pi, 0, pi]; % Set the alternate phase (+pi), added by Jianpan Huang
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == defs.M0_offset*gamma_hz*defs.B0
        if defs.Trec_M0 > 0
            seq.addBlock(mr.makeDelay(defs.Trec_M0));
        end
    else
        if defs.Trec > 0
            seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
        end
    end
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:defs.n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        phase_ind = mod(np,4); % Get the pulse index accordting to the pulse order, added by Jianpan Huang
        if phase_ind == 0
            phase_ind = 4;
        end
        satPulse.phaseOffset = phase_cycle(phase_ind); % Set phase offset according to the pulse index, added by Jianpan Huang
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(defs.td)); % add delay
        end
    end
    if defs.spoiling % spoiling before readout
        seq.addSpoilerGradients();
    end
    seq.addPseudoADCBlock(); % readout trigger event
end

%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

% %% plot
% seq.plotSaturationPhase();

