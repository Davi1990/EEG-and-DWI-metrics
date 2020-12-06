%% TMS-EEG preprocessing and artifact correction
% Script prepared by Mouhsin Shafi, MD/PhD
% Email: mouhsin.shafi@gmail.com, mshafi@bidmc.harvard.edu
%% 
% %========================================================================
% %Step 1: Load Data, Upload Channels, Delete non-EEG Channels, Save
% %========================================================================

clear all; close all; clc;
%homefolder = 'C:\Users\mshafi\Dropbox\MATLAB';
%datafolder = 'C:\Users\mshafi\Desktop\TMS-EEG Processing';
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
cd(datafolder);

Files = dir('*.vhdr');
gate_channel = [];
channelsremoved = [];

cd(homefolder);

 %-------------------------------------------------------------------------
for i = 1 : length(Files)
    %Load the file
    [ALLEEG, EEG, ~, ~] = eeglab;
    FileName = Files(i).name;
    cd(datafolder);
    thisfilename = FileName;
    Ind1 = find(FileName == '.');
    basefilename = thisfilename(1:Ind1-1);
    EEG = pop_loadbv([datafolder '/'], thisfilename, [],[]);
%    EEG = pop_loadbv([datafolder '\'], thisfilename, [],[]);
    EEG.setname= [basefilename '.set'];
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',EEG.setname,'gui','off');
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw
    counter = 1;
    
    %Now define channel locations
    EEG=pop_chanedit(EEG, 'lookup', ...
        '/Users/davidemomi/Documents/MATLAB/toolbox/eeglab14_1_0b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');

    %Save resulting dataset
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.setname = [basefilename '_S01.set'];
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
    clear EEG ALLEEG CURRENTSET
    
end



%%
clear all; close all; clc;
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab 
EEG1 = pop_loadset( '52SC_V4_SP_M1_1_S01.set', datafolder); % read in the dataset
EEG2 = pop_loadset( '52SC_V4_SP_M1_2_S01.set', datafolder); % read in the dataset
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG1); % copy it to ALLEEG


pop_mergeset(EEG1, EEG2)
basefilename = EEG1.filename;
fs_indices = strfind(basefilename,'.set'); 
basefilename2 = basefilename(1:fs_indices-1)
EEG.setname = [basefilename2 '_Merged.set'];

pop_saveset(ans, EEG.setname, datafolder)


%% 
% %========================================================================
% %Step 2: Epoch data & baseline correction
% %========================================================================

clear all; close all; clc;
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';

%     Epoch durations; note that using longer periods here to ensure that
%     time-frequency analysis is accurate in the -1000 to +1000 time period
%     for frequencies above 2-3 Hz
epoch_begin = -1; %Duration of time before pulse to begin epoch, in seconds
epoch_end = 2; %Duration of time after pulse to end epoch, in seconds
baseremovebegin = -900; %Beginning of time period of baseline correction
baseremoveend = -100; %End of time period of baseline correction
% eventmarker = 'R 15';
eventmarker = 'S128';

cd(datafolder);
Files = dir('*_Merged.set');
cd(homefolder);

% -------------------------------------------------------------------------
for i = 1 : length(Files)
    
    counter = 1;
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    % Load the file
    FileName = Files(i).name;
    cd(datafolder);
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    EEG = eeg_checkset( EEG );
    
    Ind1 = strfind(FileName, 'Merged.set');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    % Epoch data
    newfilename = [basefilename '_S02_Ep.set'];
    EEG = pop_epoch( EEG, { eventmarker }, [epoch_begin  epoch_end], 'newname', newfilename, 'epochinfo', 'yes');
    EEG.datahistory{counter} = ['Epochs extracted from ' num2str(epoch_begin) ...
        's to ' num2str(epoch_end) 's'];
    counter = counter+1;
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    % Baseline correct
    EEG = pop_rmbase( EEG, [baseremovebegin baseremoveend]);
    EEG.datahistory{counter} = ['Baseline removed from ' num2str(baseremovebegin) ...
                                'ms to ' num2str(baseremoveend) 'ms'];
    counter = counter+1;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    EEG.setname = newfilename;
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
    ALLEEG = pop_delset( ALLEEG, [1] );
    clear EEG ALLEEG CURRENTSET
end


%%
% %=========================================================================
% %                  Step 3: Reject bad channels
% %========================================================================

clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
bandlow = 1; %Lower band for visualization bandpass filter
bandhigh = 50; %Upper band for visualization bandpass filter
stoplow = 57; %Lower band for visualization bandstop filter
stophigh = 63; %Upper band for visualization bandstop filter
zeropad = 100; %Time, in ms, to zero-pad for visualization
filtord = 4; %Filter order
cd(datafolder);
Files = dir('*_S02_Ep.set');

for i = 1 : length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    cd(datafolder);
    FileName = Files(i).name;    %Bandstop (notch) filter
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S02_Ep.set');
    basefilename = FileName(1:Ind1-1);
    
    %For visualization only, zeropad to target ms and then 
    % bandpass filter to <50 Hz
    %Original unfiltered data stored in EEG.origdata
    EEG.origdata = EEG.data;
    ind = find(EEG.times>=0, 1, 'first');
    ind2 = find(EEG.times<=zeropad, 1, 'last');
    EEG.data(:,ind:ind2,:) = 0;   

    % Filter data for visualization only ---------------------------------
    EEG = tesa_filtbutter( EEG, stoplow, stophigh, filtord, 'bandstop' );
    %Bandpass filter
    EEG = tesa_filtbutter( EEG, bandlow, bandhigh, filtord, 'bandpass' );
    
    %     Now select channels to delete
    disp('Identify channels that need to be deleted by right-clicking on the channel name.');
    disp('Then close the EEG viewer and type "dbcont"');
    pop_VisEd(EEG);
    keyboard
    
    numbadchans = 0;
    for count=1:size(EEG.chanlocs, 2)
        if EEG.chanlocs(1,count).badchan == 1
            numbadchans = numbadchans+1;
            EEG.BadEl{numbadchans} = EEG.chanlocs(1,count).labels;
        end
    end
    
    %Put unfiltered data back in the EEG.data
    EEG.data = EEG.origdata;
    EEG.origdata = [];
        
    %And now delete bad channels
    if (numbadchans>0)
        EEG=pop_select(EEG,'nochannel',EEG.BadEl); % delete bad channels
        EEG = eeg_checkset( EEG );
        EEG.data=double(EEG.data);
    else
        EEG.BadEl = {};
    end
    
    thisfile = [basefilename '_S03_elec.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    clear ALLEEG EEG numbadchans;
end

%%
% %========================================================================
% %Step 4: Zero-pad early TMS artifact
% %========================================================================
% Changes from prior version: Using TESA for cutting data
% Determines voltages at 10, 15 and 20 ms after TMS pulse, as well as time
% at which voltage is less than the desired limit. Then asks user to define
% end of zero-pad period. Note that uses a default value defined in
% cuttimes if user does not provide a value

clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
cuttimes = [-2 14]; %VAR - times to cut out (zero-pad)
voltlimit = 150; %VAR - voltage limit for zero padding

cd(datafolder);
Files = dir('*_S03_elec.set');

timevolt = zeros(length(Files),4);
for i = 1 : length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    
    %Find max mean (across epochs) voltage across all channels at 10, 15 
    % and 20ms, and determine time at which max voltage < limit
    EEG.meandata = mean(EEG.data,3);
    tempdata = max(abs(EEG.meandata));
    ind1 = find(EEG.times>=0,1,'first') + 3; %deliberately choose point 3ms after pulse
    ind2 = find(tempdata(1,ind1:size(tempdata,2))<voltlimit,1,'first') + ind1 - 1;
    EEG.timevolt.threshvalue = voltlimit;
    EEG.timevolt.threshold = EEG.times(ind2);
    timevolt(i,1) = EEG.timevolt.threshold;
    EEG.timevolt.t10 = tempdata(1,find(EEG.times>=10,1,'first'));
    timevolt(i,2) = EEG.timevolt.t10;
    EEG.timevolt.t15 = tempdata(1,find(EEG.times>=15,1,'first'));
    timevolt(i,3) = EEG.timevolt.t15;
    EEG.timevolt.t20 = tempdata(1,find(EEG.times>=20,1,'first'));
    timevolt(i,4) = EEG.timevolt.t20;
    
%     %Display that information, and then prompt user for input for final
%     %zero-pad time
%     fprintf('\nFor subject %i the voltage at time 10 is %f, time 15 is %f, time 20 is %f. \n', ...
%         i, EEG.timevolt.t10, EEG.timevolt.t15, EEG.timevolt.t20);
%     fprintf('The first timepoint at which the voltage is less than the threshold of %i is %i.\n\n', ...
%         EEG.timevolt.threshvalue, EEG.timevolt.threshold);
    EEG = pop_saveset( EEG, 'filename', FileName, 'filepath',datafolder);
    clear ALLEEG EEG 
end
    
    disp('tThresh V10 V15 20');
    disp(timevolt);
    prompt = {'Enter last time to zero-pad:'};
    dlg_title = 'Zero pad end time'; num_lines = 1; defaultans = {num2str(cuttimes(1,2))};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    cuttimes(1,2) = str2num(answer{1,1});
    
for i = 1 : length(Files)
    
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    
    Ind1 = strfind(FileName, '_S03_');
    basefilename = FileName(1:Ind1-1);
    %Now remove the data in the desired interval
    EEG = tesa_removedata(EEG, cuttimes);
    thisfile = [basefilename '_S04_zero.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    clear ALLEEG EEG 
end

cd(homefolder);
%%
% % ========================================================================
% Step 5: Remove bad epochs - No changes made
% ========================================================================
clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
sinchan_prob_thresh = 3.5; %VAR: Single-channel threshold for rejection based
                            % on probability; default use 3.5
allchan_prob_thresh = 3; %VAR: All-channel threshold for rejection based
                            % on probability; default use 3
sinchan_kurt_thresh = 5; %VAR: Single-channel threshold for rejection based
                            % on kurtosis; default use 5
allchan_kurt_thresh = 3; %VAR: All-channel threshold for rejection based
                            % on kurtosis; default use 3
voltage_reject = 1; %VAR: Whether to reject based on a voltage threshold;
                        % Use 1 for yes, 0 otherwise
voltage_thresh =  100; %VAR: Voltage threshold for epoch rejection
ChansExclude = {'EOG1', 'EOG2' 'Fp1', 'Fpz', 'Fp2', 'AF7', 'AF8'}; %VAR: Channels to exclude from voltage thresholding
BeginTimeInclude = -0.5; %VAR: Beginning of time period to do voltage thresholding on
EndTimeInclude = 1; %VAR: End of time period to do voltage thresholding on
BeginTimeExclude = 0; %VAR: Beginning of time period to exclude from voltage thresholding
EndTimeExclude = 0.05; %VAR: End of time period to exclude from voltage thresholding (in seconds)

cd(datafolder);
Files = dir('*_S04_zero.set');

for i = 1 : length(Files)
  

    FileName = Files(i).name;
    clear EEG ALLEEG CURRENTSET ALLCOM count2 chanind;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S04_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
% ------------------------------------------------------------------------
%           5a) Copy original data into new variable, then bandpass and
%           average reference data. Note that bandpassing and rereferecing 
%           is strictly for visualization to help determine which epochs to 
%           delete, and the original unfiltered data is replaced at the end
%           of this step
% ------------------------------------------------------------------------ 
    EEG.origdata = EEG.data;
    filtord = 4;
    stoplow = 57;
    stophigh = 63;
    bandlow = 1;    %Lower edge of bandpass filter
    bandhigh = 50;   %Upper edge of bandpass filter
    % -------------- Filtering
    % Uses fourth-order Butterworth filter
    % Backwards and forwards filtering
    %Bandstop (notch) filter
    EEG = tesa_filtbutter( EEG, stoplow, stophigh, filtord, 'bandstop' );
    %Bandpass filter
    EEG = tesa_filtbutter( EEG, bandlow, bandhigh, filtord, 'bandpass' );
    % -------------- Average Referencing
    EEG = pop_reref(EEG, []);
    

% ------------------------------------------------------------------------
%           5b) Tag trials based on amplitude, probability, and kurtosis  
% ------------------------------------------------------------------------ 
    %Use EEGLAB function to calcuate Kurtosis and tag bad channels
    search_array = [];
    
    % Now identify abnormal trials based on probability and kurtosis
    EEG = pop_jointprob(EEG,1,[1:size(EEG.chanlocs,2)],sinchan_prob_thresh,allchan_prob_thresh,1,0);
    EEG = pop_rejkurt(EEG,1,[1:size(EEG.chanlocs,2)],sinchan_kurt_thresh,allchan_kurt_thresh,1,0);
    
    % Below will also identify trials with activity greater than the
    % voltage threshold, but will exclude eye channels
    if (voltage_reject==1)
        for count = 1 : EEG.nbchan
            electrodes{count} = EEG.chanlocs(count).labels;
            search_array = [search_array count];
        end
        
        for count2 = 1 : length(ChansExclude);
            tempchan = find(strcmp(electrodes, ChansExclude{count2}));
            if ~isempty(tempchan) 
                chanind(count2) = tempchan;
                search_array = search_array(search_array~=chanind(count2));
            end
        end
        

        EEG = pop_eegthresh(EEG,1,search_array,-voltage_thresh,voltage_thresh,...
            [BeginTimeInclude EndTimeExclude],[BeginTimeExclude EndTimeInclude],1,0);
    end
            
% ------------------------------------------------------------------------
%  5c)  Scroll, visualize tagged trials, manually tag trials with lots of 
%       muscle and other artifacts, untag what is not noisy, finalize
%       trials to delete, and save resulting dataset
% ------------------------------------------------------------------------ 
    eeglab redraw;
    disp('Visually select trials with lots of muscle and other artifacts, using pop_rejmenu');
    disp('Update marked trials, but do NOT delete. Then type dbcont when done');
    pop_rejmenu(EEG,1);
    keyboard;

%             Run pop_rejmenu from the GUI here, update marked trials, but
%             do NOT delete yet (needs to be saved)

    EEG.badtr = union(union(union(find(EEG.reject.rejmanual>0), ...
        find(EEG.reject.rejjp>0)), union(find(EEG.reject.rejkurt>0), ...
        find(EEG.reject.rejthresh>0))),find(EEG.reject.rejconst>0));
    EEG.setname = [basefilename '_marked.set'];
    EEG = pop_saveset( EEG, 'filename',[basefilename '_S04b_marked.set'],...
            'filepath',datafolder);

% ------------------------------------------------------------------------
%   5d) Delete Bad Trials
% ------------------------------------------------------------------------
            
            EEG.data = EEG.origdata;  %Copy original unfiltered data back
            EEG.origdata = [];
            EEG = pop_rejepoch( EEG, EEG.badtr ,0); % EEGLAB function to delete all tagged trials
            EEG = eeg_checkset( EEG );
            EEG.data=double(EEG.data);

    thisfile = [basefilename '_S05_ClEp.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    clear ALLEEG EEG
end
cd(homefolder);

%%
% %========================================================================
% %Step 6: First fICA run and removing TMS pulse/muscle artifact
% %========================================================================
clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
DoPCA = 1; %VAR: Whether to do dimensionality reduction w PCA prior to ICA: 0 for no, 1 for yes
CalculateDimensions = 0; %VAR: Whether to calculate number of dimensions to reduce data to prior to ICA
PercentVar = 99; %VAR: Percent of variance to explain in PCA reduction
MinComp = 30; %VAR: Minimum number of components to include if calculating variance
PCAdimensions = 60; %VAR: Number of dimensions to reduce to if doing PCA and not calculating
plottime = [-100 300]; %VAR: Display window for TESA
MuscleThresh = 8; %VAR: Threshold to start using for detection of artifacts
MuscleWin = [20 35]; %VAR: Window to determine muscle artifact; will neeed
                        %to be adjusted if more than 15s taken out

cd(datafolder);

Files = dir('*_S05_ClEp.set');

for i = 1 : length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    
    Ind1 = strfind(FileName, '_S05_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    % Run ICA
    if DoPCA == 1
        if CalculateDimensions == 1
            PCAdimensions = [];
            PCAdimensions = fcn_EstimateNrICAComp(EEG, PercentVar, MinComp);
            str = ['For file ' basefilename ', ' num2str(PCAdimensions) 'components are necessary to keep 99%% of the variance in the data'];
            disp(str);
            EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
        else
            EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
        end
        EEG.ica1_dimensions = PCAdimensions;
    else
        EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');
    end
    
    %First sort components by percent of variance explained
    [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
    
    thisfile = [basefilename '_S05b_fICA' '.set'];
    EEG = pop_saveset( EEG, 'filename',thisfile,'filepath',datafolder);
    EEG = eeg_checkset( EEG );
    
    %Open Component activations in eegplot window
	  pop_eegplot(EEG, 0, 1, 0, [], 'winlength', 3, 'dispchans', 5);

    %Section below uses TESA for artifact selection, focusing only on TMS
    %pulse artifact. 
    EEG = tesa_compselect( EEG,'compCheck','on','comps',[],'figSize','medium', ...
        'plotTimeX', plottime,'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh', ...
        MuscleThresh,'tmsMuscleWin', MuscleWin,'tmsMuscleFeedback','on', ...
        'blink','off','move','off','muscle','off','elecNoise','off');
    EEG.filepath = []; %Removes filepath so it is not stored for later
    
    %Save filename for after component marking
    thisfile = [basefilename '_S05c_ICAmarked.set'];
    EEG.setname= thisfile; EEG.filepath = datafolder;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    
    thisfile = [basefilename '_S06_ICA1.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    clear ALLEEG EEG
    
end
cd(homefolder);

%%
% %========================================================================
% %Step 7: Interpolate missing data, then filter and average reference
% %========================================================================
% % ========================================================================
clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';

% datafolder = strcat('C:\Users\pbouche1\Desktop\66_SB\VS', num2str(j),'\Merged');
interp_method = 'linear'; %VAR: Interpolation method; use either this or cubic
cubwindow = [-25 35]; %VAR: IF using cubic interpolation, window to use
% interp_method = 'cubic'; %VAR: Interpolation method; alternative is linear
cuttimes = [-2 14];
stoplow=57; % VAR - BAND STOP LOWER EDGE
stophigh=63; % VAR - BAND STOP UPPER EDGE
bandlow=1; % VAR - BANDPASS LOWER EDGE
bandhigh=100; % VAR - BANDPASS UPPER EDGE
filtord = 4; %VAR - FILTER ORDER

Rereference = 1; % VAR - Set to 0 for no rereferencing, 1 for average reference

Subepoch = [-0.5 1]; %VAR: Subepoch, in seconds, to extract prior to ICA

cd(datafolder);
Files = dir('*_S06_ICA1.set');

for i = 1 : length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S06_ICA1.set');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    
%     EEG = tesa_removedata(EEG, cuttimes);
    %Interpolate missing data
    XX={'R 7', 'R 15'};
    EEG.tmscut.cutEvent=XX;
    EEG = tesa_interpdata(EEG, interp_method, cubwindow);
    
    %Bandstop (notch) filter
    EEG = tesa_filtbutter( EEG, stoplow, stophigh, filtord, 'bandstop' );
    
    %Bandpass filter
    EEG = tesa_filtbutter( EEG, bandlow, bandhigh, filtord, 'bandpass' );
    
    %And average reference (if desired)%     
    if Rereference == 1                 %LEAVE UNCOMMENTED
        EEG = pop_reref(EEG, []);
    end
    
    %Now extract subepoch, getting rid of edges, to do second round ICA
    EEG = pop_select( EEG,'time', Subepoch);
    
    thisfile = [basefilename '_S07_interpfilteravref.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
    clear ALLEEG EEG
end


%%
% %=======================================================================
% %Step 8: Load Clean Data, Run ICA 
% %======================================================================== 

clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
cd(datafolder);
DoPCA = 1; %VAR: Whether to do dimensionality reduction w PCA prior to ICA: 0 for no, 1 for yes
CalculateDimensions = 0; %VAR: Whether to calculate number of dimensions to reduce data to prior to ICA
PercentVar = 99; %VAR: If calculating dimensions, percent of variance to explain in PCA reduction
MinComp = 57; %Minimum number of components to include if calculating variance
PCAdimensions = 57; %VAR: Number of dimensions to reduce to if doing PCA but not calculating dimensions
ICAmethod = 'fast ICA'; %VAR: ICA type
% ICAmethod = 'infomax';

Files = dir('*_S07_interpfilteravref.set'); 
for i = 1 : length(Files)
    FileName = Files(i).name;
    
    clear EEG ALLEEG CURRENTSET ALLCOM;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S07_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw

    % -------------------------------------------------------------
    %               10a)Perform ICA - fastica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'fast ICA')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp(EEG, PercentVar, MinComp);
                str = ['For file ' basefilename ', ' num2str(PCAdimensions) 'components are necessary to keep 99%% of the variance in the data'];
                disp(str);
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            else
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            end
            EEG.ica2_dimensions = PCAdimensions;
        else
            EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_S08_fICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_S08_fICA.set'];
        end
        
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    % -------------------------------------------------------------
    %               10b)Perform ICA - runica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'infomax')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp_RunPCA_Var(EEG, PercentVar);
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            else
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            end
        else
            EEG = pop_runica(pop_runica(EEG, 'extended',1));
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_S08_rICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_S08_rICA.set'];
        end
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    string = ['Finished running ICA on ' FileName];
    disp(string);
    
end
      
% ************************************************************************
%                               END
% ************************************************************************

%% 
% % =======================================================================
% % Step 9: Second round of component selection using TESA
% % ======================================================================== 
clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
cd(datafolder);

Files = dir('*_S08_fICA57.set'); 
for i = 1 : length(Files)
    FileName = Files(i).name;
    
    clear EEG ALLEEG CURRENTSET ALLCOM;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S08_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    %First sort components by percent of variance explained
    [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);

    %Open Component activations in eegplot window
	  pop_eegplot(EEG, 0, 1, 0, [], 'winlength', 5, 'dispchans', 5);
    
    %Create filename for after component marking
    thisfile = [basefilename '_S08b_marked.set'];
    EEG.setname= thisfile; EEG.filepath = datafolder;

    %Now call TESA for automatic component selection
    %Does TMS pulse, blink, lateral eye movements, muscle and electrode
    firsttime = EEG.tmscut(1,1).cutTimesTMS(2) + 1;
    EEG = tesa_compselect( EEG,'compCheck','on','comps',[],'figSize','medium','plotTimeX',[-100 300],'plotFreqX',[1 100], ...
        'tmsMuscle','on','tmsMuscleThresh', 6,'tmsMuscleWin',[firsttime firsttime+20],'tmsMuscleFeedback','off',...
        'blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off',...
        'move','on','moveThresh',2, 'moveElecs',{'F7','F8'},'moveFeedback','off', ...
        'muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off',...
        'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    EEG.filepath = []; %Removes filepath so it is not stored for later   
    thisfile = [basefilename '_S09_ICA2.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
end

%%
% % =======================================================================
% %  Step 10: Low-Pass Filter at 50Hz, then interpolate missing channels
% % =======================================================================

clear all; close all; clc; eeglab;
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
UpperEdge = 50; % VAR - LOWPASS EDGE

cd(datafolder);
EEG = pop_loadset('filename', 'BV_ElectrodeTemplate.set', 'filepath', datafolder);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
Files = dir('*_S09_ICA2.set');
for ii = 1 : length(Files)
    FileName = Files(ii).name;
    Ind1 = strfind(FileName, 'S09_');
    basefilename = FileName(1:Ind1-2);
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

   % Low-pass filter the data
    Fs=EEG.srate;ord=4;
    [xall,yall]=butter(ord,UpperEdge/(Fs/2),'low');
    
    for trial = 1:size(EEG.data,3)
        for ch=1:size(EEG.data,1)
            EEG.data(ch,:,trial) = double(filtfilt(xall,yall, double(EEG.data(ch,:,trial))));
        end
    end
    
    %And interpolate missing channels
    EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
%     
    %And delete EOG leads
%     EEG = pop_select( EEG,'nochannel',{'EOG1' 'EOG2'});
    EEG = pop_select( EEG,'nochannel',{'EOG1' 'EOG2' 'Iz' 'Fp1' 'Fpz'});
    
    EEG.setname = [basefilename '_S10_Final'];
    EEG = pop_saveset( EEG, 'filename', EEG.setname, 'filepath', datafolder);
    ALLEEG = pop_delset( ALLEEG, [2] );
    [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
end

%% 
% =======================================================================
%  Step 11: GMFA Analysis: Calculate & Plot GMFA and peaks
% =======================================================================

clear all; close all; clc; 
homefolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
datafolder = '/Volumes/WD_Elements/BROAD/EEG/M1_IPL/Subjects/52/M1/visit_2';
cd(datafolder);
Linetime = 0; %Time of pulse
Graph1_start = -200; %VAR Time at which to begin topoplot
Graph1_end = 400; %VAR Time at which to end topoplot
TOI_start = 20; %VAR Time of interest start for peakfinding
TOI_end = 400; %VAR Time of interest end for peakfinding
threshval = 2; %VAR Number of standard deviations above the mean to be a 
                % minimum peak after stimulation; note that this does not
                % need to be defined at all, but probably should be, at
                % least as the max baseline voltage
selval = 3; %VAR Criteria for labeling something a peak; does not need to 
            %be set manually
first_basetime = -450; %VAR Beginning of official baseline period
last_basetime = -50; %VAR End of official baseline period
FStitle = 16; %VAR Figure title font size
FS = 14; %VAR Figure axes font size

Files = dir('*_S10_Final.set'); 
for i = 1 : length(Files)
    FileName = Files(i).name;
    clear EEG ALLEEG CURRENTSET ALLCOM;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S10_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    %Define variables used in script
    graphstart_bin = find((EEG.times>=Graph1_start), 1, 'first');
    graphend_bin = find((EEG.times>=Graph1_end),1, 'first');
    firstbasebin = find((EEG.times>=first_basetime),1,'first');
    lastbasebin = find((EEG.times<=last_basetime),1,'last');
    pulse1_latency = EEG.event(1,1).latency;
    TOIstart_bin = find((EEG.times>=TOI_start), 1, 'first');
    TOIend_bin = find((EEG.times>=TOI_end),1, 'first');
    First_peak_time = ceil(EEG.tmscut.cutTimesTMS(2));
    
    % GMFA calculations
    %----------------------------------------------------------------------
    %Do Standard GMFA analysis here using TESA
    EEG = pop_tesa_tepextract( EEG, 'GMFA'); %Standard GMFA analysis
    
    % Calculate the normalized GMFA, which is GMFA divided by the mean GMFA
    % in the baseline period
    AvgbaseGMFA=mean(EEG.GMFA.R1.tseries(firstbasebin:lastbasebin));
    EEG.GMFA.norm.tseries = EEG.GMFA.R1.tseries/AvgbaseGMFA;
    EEG.GMFA.norm.time = EEG.GMFA.R1.time;
    
    %Now detect peaks automatically (note that this is NOT a TESA function,
    %but rather uses the PeakFinder script, and is my own analysis)
%     %----------------------------------------------------------------------
    %First define the threshold - threshval SD above mean baseline value
    tempdata = [];
    tempdata = EEG.GMFA.R1.tseries;
    peakLoc = []; peakMag = []; peaks=[]; 
    mean_baseval = mean(tempdata(firstbasebin:lastbasebin));
    std_baseval = std(tempdata(firstbasebin:lastbasebin));
    maxbaseval = max(tempdata(firstbasebin:lastbasebin));
    thresh = max(((std_baseval * threshval) + mean_baseval),maxbaseval); 
   
    %And selection criteria for something to be defined a peak
    sel = (max(tempdata(firstbasebin:lastbasebin))-min(tempdata(firstbasebin:lastbasebin)))/selval;
    
    %Now find and save the peaks
    [peakLoc, peakMag] = peakfinder(tempdata, sel, thresh);
    peakLoc=EEG.GMFA.R1.time(peakLoc);
    peaks(:,1)=peakLoc';
    peaks(:,2)=peakMag';
    
    %Delete peaks before time region of interest
    numpeaks = size(peaks,1);
    for count=numpeaks:-1:1
        if ((peaks(count,1)<TOI_start))
            peaks(count,:) = [];
        end
    end
    
    %And Delete peaks after time region of interest
    numpeaks = size(peaks,1);
    for count=numpeaks:-1:1
        if ((peaks(count,1))>TOI_end)
            peaks(count,:) = [];
        end
    end
    
    EEG.GMFA.R1.peaks(:,1)= peaks(:,1);
    EEG.GMFA.R1.peaks(:,2)= peaks(:,2);
    
    %Plot nGMFP
    figure;plot(EEG.times(graphstart_bin:graphend_bin), EEG.GMFA.norm.tseries(1,graphstart_bin:graphend_bin)); hold on;
    titlename = [basefilename ' Normalized GMFP'];
    title(titlename);
    ymax = max(EEG.GMFA.norm.tseries);
    axis([Graph1_start Graph1_end 0 (ymax + (ymax/10))]);
    line([Linetime, Linetime], [0 ymax], 'Color', 'r'); hold off;  
    saveas(gcf,titlename);
    saveas(gcf,titlename, 'png');
    
%     %Save the resulting figure
    % exportfig(gcf, titlename,'format', 'png', 'Color', 'cmyk', ...
    %   'Resolution', 600, 'FontMode', 'scaled', 'Bounds', 'tight', ...
    %  'LineMode','scaled');
    % clf;
    
    %Plot ERPtopo
    topotitle = [basefilename ' ERPTopo'];
    figure; pop_timtopo(EEG, [Graph1_start Graph1_end], peaks(:,1));
    gtext(topotitle, 'fontsize', 12);
    saveas(gcf,topotitle);
    saveas(gcf,topotitle, 'png');
    
%     %Save the resulting figure
    %exportfig(gcf, topotitle,'format', 'png', 'Color', 'cmyk', ...
     %   'Resolution', 600, 'FontMode', 'scaled', 'FontSize', 10, ...
     %  'FontSizeMin', 8, 'FontSizeMax', 12, 'Bounds', 'tight', ...
     % 'LineMode','scaled', 'Width', 7);
     %clf;

    
    %And save the resulting file
    EEG.setname = [basefilename '_S11_GMFA.set'];
    EEG = pop_saveset( EEG, 'filename', EEG.setname, 'filepath', datafolder);
    close all;
    eeglab redraw;   
end
