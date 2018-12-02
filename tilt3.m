%----------------------------------------------------------------------
%                Action-Tilt Experiment (Beta) - 2018-11-06
%----------------------------------------------------------------------

% Developed by:
% Raymond MacNeil 
% University of British Columbia
% Department of Psychology | Vision Lab
% 2136 West Mall
% Vancouver, BC Canada V6T 1Z4
% T: 604-827-5763
% E: raymond.macneil@psych.ubc.ca

% This experiment was created using the Matlab Psychtoolbox 3.0.14 
% Flavor: beta - Corresponds to SVN Revision 8715 
% For more info visit: https://github.com/Psychtoolbox-3/Psychtoolbox-3

% This is experiment requires participants to perform
% an orientation discimination task. It has a 2x4x6 factorial design. 
% Factors are: 1) Illusion Type (RFI or STI); 2) Inducer Orientation (15 deg 
% CW or CCW); 3) Target Orientation (8 deg CW or CCW); and 4) SOAs of
% 1.6, 0.8, 0.4, 0.2, 0.1, 0 seconds.
% The inducing context of the RFI or STI is presented for 200ms
% Modeled after experiment described in Corbett et al. (2009) 
% doi:10.1016/j.visres.2008.09.020

%----------------------------------------------------------------------
%                          Exp. Launch Dialog
%----------------------------------------------------------------------
prompt = {'Please Enter the Subject ID'};
blank = {''};

options.Resize='off';
options.WindowStyle='normal';
options.Interpreter='none';

answer = inputdlg(prompt, 'Action Tilt', [1 50], blank, options);
if isempty(answer)
   return
end

%----------------------------------------------------------------------
%                   Setup Pyschtoolbox & Key Defaults
%----------------------------------------------------------------------

% clear the screen
sca;
close all;
clearvars;

% Hide the cursor
HideCursor;

% load some default settings for Psychtoolbox
PsychDefaultSetup(2)

% reaarange the random number generator
rng(sum(100 * clock));

% debug mode, skip SyncTests
Screen('Preference', 'SkipSyncTests', 1);              

% get screen numbers
screens = Screen('Screens');

% if there is an external displayed, then we want to draw to this
screenNumber = max(screens);

% define white, black, and grey
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;

% open a window
[windowExp, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

% measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', windowExp);

% get window size in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', windowExp);

% get center screen coordinates
[xCenter, yCenter] = RectCenter(windowRect);

% get the refresh rate of our screen 
% the relationship between the two is: ifi = 1 / hertz
hertz = FrameRate(windowExp);

% get color depth of the pixel in bits
pixelSize = Screen('PixelSize', windowExp);

% get the maximum coded luminance level (this should be 1)
maxLum = Screen('ColorRange', windowExp);
    


% flip to clear
Screen('Flip', windowExp);

%----------------------------------------------------------------------
%                        Experimental Setup
%----------------------------------------------------------------------
% Generate a condition matrix representing all possible trials
% Illusion type, SOA, inducercw[1]orccw[-1], target(cw[1] or ccw[-1])
% only the STI is included for the beta version, development in progress 
% Col2; 1=1.6s; 2=0.8s; 3=0.4s; 4=0.2s, 5=0.1s, 6=0s, 

sti_table=[1   1   1   1; 
           1   1   1  -1;
           1   1  -1   1;
           1   1  -1  -1;
           1   2   1   1; 
           1   2   1  -1;
           1   2  -1   1;
           1   2  -1  -1;
           1   3   1   1;
           1   3   1  -1;
           1   3  -1   1;
           1   3  -1  -1;
           1   4   1   1;
           1   4   1  -1;
           1   4  -1   1;
           1   4  -1  -1; 
           1   5   1   1;
           1   5   1  -1;
           1   5  -1   1;
           1   5  -1  -1;
           1   6   1   1;
           1   6   1  -1;
           1   6  -1   1;
           1   6  -1  -1];  
       
rfi_table=[2   1   1   1; 
           2   1   1  -1;
           2   1  -1   1;
           2   1  -1  -1;
           2   2   1   1; 
           2   2   1  -1;
           2   2  -1   1;
           2   2  -1  -1;
           2   3   1   1;
           2   3   1  -1;
           2   3  -1   1;
           2   3  -1  -1;
           2   4   1   1;
           2   4   1  -1;
           2   4  -1   1;
           2   4  -1  -1; 
           2   5   1   1;
           2   5   1  -1;
           2   5  -1   1;
           2   5  -1  -1;
           2   6   1   1;
           2   6   1  -1;
           2   6  -1   1;
           2   6  -1  -1];

exp_matrix = [sti_table; rfi_table]       
       
% define number of trials
nTrials = size(exp_matrix, 1);

% re-order the 64 trials randomly for x, and generate four blocks
[m,n] = size(exp_matrix);
idx = randperm(m);
ExpTrialSeq = exp_matrix;
ExpTrialSeq(idx,:) = exp_matrix(:,:);

% number of training trials
nTrain = 16;

% compile training trials from a random subset of CondTable 
TrainTrials = [sti_table(randperm(nTrials/2, nTrain/2),:); 
               rfi_table(randperm(nTrials/2, nTrain/2),:)];
TrainTrials = Shuffle(TrainTrials);           
           
% append to top of random trial sequence matrix already generated
FullExpSeq = [TrainTrials; ExpTrialSeq]; 
l = transpose([1:length(FullExpSeq(:,1))]); %#ok<*NBRAK>

% preallocate matrix for the illusion condition schedule
illus_sched = false(length(l), 2);


% determine row indices from FullExpSeq corresponding with the randomly
% generated sequencing of STI and RFI trials
illus_vect = FullExpSeq(:,1);
[sti_row_idx] = find(illus_vect==1);
[rfi_row_idx] = find(illus_vect==2);

% define row indices of the STI illus_sched column (col 1) accordingly
illus_sched(sti_row_idx, 1) = true;
% define row indices of the RFI illus_sched column (col 2) accordingly
illus_sched(rfi_row_idx, 2) = true;

% preallocate results matrices 
keyResp = NaN * ones(length(l), 1);
response_time = NaN * ones(length(l), 1);
corr = strings(length(l),1);

% define primary grating parameters
width = 300;  % width of circular grating in pixels
height = 300; % height of circular grating in pixels
radius = [150, 30]; % by assigning a value to radius a circular aperture is generated
angle = [15, 8]; % orientation angles for inducing and target stimuli
phase = 0; % phase of the sin wave grating
freq = 15/180; 
contrastPreMultiplicator = 0.50;
contrast = 2;
backgroundColorOffset = [0, 0, 0, 0];
auxParameters = [phase, freq, contrast, 0;
                 phase+90, freq, contrast, 0];

% define primary RFI parameters   
frameRead = imread('rfisq.jpg');
yCord = 12;
smooth = 2;
theta = 6;
% Generating an anti-alised (smooth) line to form the rod requires that
% we call the Screen('DrawLines') function. Therefore, we have to
% define the coordinates of each end of the line that forms an angle
% theta with the vertical of our screen (i.e. the y-axis). Recall
% that the trigonometric function tan(theta) can be expressed 
% as the ratio: opposite/adjacent (of a given right angle triangle); If we 
% imagine that the middle of our screen is the point of origin on a  
% cartesian plane (0,0), then we define y-coordinate corresponding to
% half the length of our line; so we simply use y/tan(theta) to get the x 
% coordinate; note that the y-coordinate in this case represents the 
% length of the adjacent side of a right angle triangle that forms one of
% its vertices at the point of origin. Our line--well half of it--is the 
% hypotenuse! (In MatLab, the default is for trigonometric functions to 
% operate on radians, so you have to convert from degrees.) We could have 
% similarly approached this problem by directly defining the length of our 
% line instead of a more 'arbitrary' value. Then after dividing that by 
% two, we can use the the trig identities sin and cos to get absolute 
% values for the x and y coordinates, respectively. The matrix below
% for linCords gives you the coordinate positions for each end of the line
xCord = yCord/tan(pi/2-(theta*(pi/180))); 
lineCords = [-xCord xCord; yCord -yCord]; 
ExpLineCords = [NaN NaN; yCord -yCord];


% create the grating patterns
[induceID, STIinduceRect] = CreateProceduralSquareWaveGrating(windowExp,...             
        width, height, backgroundColorOffset, radius(1), contrastPreMultiplicator);  %#ok<*ASGLU>
[targetID, STItargetRect] = CreateProceduralSquareWaveGrating(windowExp,... 
        (width/5), (height/5), backgroundColorOffset, radius(2), contrastPreMultiplicator);  
RFIinduceRect = [-radius(1) -radius(1) radius(1) radius(1)];    
frame_image = Screen('MakeTexture', windowExp, frameRead);    
% define placement location for the grating patterns
STIinduceRect = CenterRectOnPoint(STIinduceRect, xCenter, yCenter);  
STItargetRect = CenterRectOnPoint(STItargetRect, xCenter, yCenter);
RFIinduceRect = CenterRectOnPoint(RFIinduceRect, xCenter, yCenter); 

% Set the pixel width for stimuli outlines
penWidthPixels = 4;

% Preparing and displaying the welcome screen
% A text size of 26 pixels is readable on most screens:
Screen('TextSize', windowExp, 26);
    
% This is our intro text. The '\n' sequence creates a line-feed:
InstructText = ['In this experiment you are asked to judge and report on the orientation of two stimulus types.\n' ...
                '                                             \n' ...
                'Each trial will consist of TWO brief stimulus presentations: 1) A context 2) A target.\n' ...                                                        
                '                                                \n' ...
                'Your job is to judge and report on the orientation of the SMALLER circular grating pattern\n' ...                      
                '                                                \n' ...
                '            Press the ''Z'' KEY if the target grating is oriented counter clock-wise (to the left)\n' ...
                'Press the ''/'' KEY if the target grating is oriented clock-wise (to the right)\n' ...
                '                                                \n' ...
                'You will begin with ' num2str(nTrain) ' training trials\n' ...
                '(Press any key to initiate the training.)\n'];
         
% Draw 'InstructionsText', centered in the display window:
DrawFormattedText(windowExp, InstructText, 'center', 'center');

% Show the drawn text at next display refresh cycle:
Screen('Flip', windowExp)
    
% Wait for key stroke. This will first make sure all keys are
% released, then wait for a keypress and release:
KbWait([], 3);
    
%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Timing of inter-trial interval in seconds and frames
InterTrlIntSecs = 0.2;
InterTrlIntFrames = round(InterTrlIntSecs / ifi);

% Presentation time for the inducing stimulus in seconds and frames
presTimeSecs = 0.2;
presTimeFrames = round(presTimeSecs / ifi);

% SOA ref vectors in seconds and frames 
% SOA vector represents the stimulus onset aychrony of the inducing
% stimulus and the target stimulus
soaRefSecs = [1.6 0.8 0.4 0.2 0.1 0.0]; 
soaRefFrames = round(soaRefSecs / ifi);


% Numer of frames to wait before re-drawing
waitframes = 1;

%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% Define the keyboard keys that are listened for. We will use the left
% and right arrow keys as response keys for the task and the escape key as
% an exit/reset key
escapeKey = KbName('ESCAPE');
leftKey = KbName('z');
rightKey = KbName('/?');

%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------
    
    % initiate loop through the radomized training and exp trial sequence
    k = 0
    while (k < length(l))
        k = k + 1;
        % get the trial number
        trialNo = k; 
        % trial ID (illusion type)
        trialID = FullExpSeq(k, 1);  
        % determine the SOA (in frames) for the given trial
        soaTrial = soaRefFrames(FullExpSeq(k, 2)); 
        % get the inducer tilt angle ref
        induceTilt = FullExpSeq(k, 3) * angle(1);   
        % get the target tilt angle ref
        targetTilt = FullExpSeq(k, 4) * angle(2);
        % get rod target drawing coordinates
        ExpLineCords(1,:) = FullExpSeq(k, 4) * lineCords(1,:);
        
        % Create a stimulus schedule matrix to determine when inducing
        % stimulus and target stimulus will be shown
        if soaTrial == 0
            targetTimeFrames = 0;
            target_on_frames = (soaTrial+1):(soaTrial...
                +(presTimeFrames)+targetTimeFrames);
        elseif soaTrial ~= 0
            targetTimeFrames = 12;
            target_on_frames = (presTimeFrames+soaTrial+1):(soaTrial...
                +presTimeFrames+targetTimeFrames);
        end
   
        % Row 1 - inducing stimulus
        % Row 2 - target stimulus
        stim_sched_len = presTimeFrames+soaTrial+targetTimeFrames;        
        stim_sched = zeros(2,stim_sched_len);
        
        % Set inducing "on" frames
        stim_sched(1,1:presTimeFrames) = 1; % set first (presTimeFrames) on 
  
        % Set target "on" frames
        stim_sched(2,target_on_frames) = 1;
        
        % Change the blend function to draw an antialiased fixation point
        % in the centre of the screen   
        Screen('BlendFunction', windowExp, 'GL_SRC_ALPHA',...
            'GL_ONE_MINUS_SRC_ALPHA');
        
        if k==nTrain+1  % before the first test trial
            DrawFormattedText(windowExp, 'Are you ready for the experiment?\n(Press any key to start experiment)',...
                'center', 'center', black);
           vbl = Screen('Flip', windowExp);
            KbWait([], 3); 
        end
     
        % Flip again to sync us to the vertical retrace at the same time as
        % drawing our fixation point
        Screen('DrawDots', windowExp, [xCenter; yCenter], 10, black, [], 2);
        vbl = Screen('Flip', windowExp);
        
        % now we present the inter-tial interval
        for frame = 1:InterTrlIntFrames - 1
            % Draw the fixation point
            Screen('DrawDots', windowExp, [xCenter; yCenter], 10, black, [], 2);
            % Flip to the screen
            vbl = Screen('Flip', windowExp, vbl + (waitframes - 0.5) * ifi);
        end
        
         vbl = Screen('Flip', windowExp, vbl + (waitframes - 0.5) * ifi);
         
         for f = 1:stim_sched_len
             if stim_sched(1,f) && illus_sched(k, 1) % draw sti inducing stim
                 % Set the right blend function for drawing the inducing stim
                 Screen('BlendFunction', windowExp, 'GL_ONE', 'GL_ZERO');
                 Screen('DrawTexture', windowExp, induceID, [], STIinduceRect,...
                     induceTilt, [], [], [], [], [], auxParameters(1,:));
                 Screen('FrameOval', windowExp, black, STIinduceRect, penWidthPixels);
             end
             if stim_sched(2,f) && illus_sched(k, 1) % draw sti target stim
                 Screen('BlendFunction', windowExp, 'GL_ONE', 'GL_ZERO');
                 Screen('DrawTexture', windowExp, targetID, [], STItargetRect,...
                     targetTilt, [], [], [], [], [], auxParameters(2,:));
                 Screen('FrameOval', windowExp, black, STItargetRect, penWidthPixels);
             end
             if stim_sched(1,f) && illus_sched(k, 2) % draw rfi inducing stim
                 % Set the right blend function for drawing the inducing stim
                 Screen('BlendFunction', windowExp, 'GL_ONE', 'GL_ZERO');
                 Screen('DrawTexture', windowExp, frame_image, [],... 
                     RFIinduceRect, induceTilt, [], [], [], [], [], []);
             end
             if stim_sched(2,f) && illus_sched(k, 2)
                 Screen('BlendFunction', windowExp, 'GL_SRC_ALPHA',...
                     'GL_ONE_MINUS_SRC_ALPHA');
                 Screen('DrawLines', windowExp, ExpLineCords, penWidthPixels,...
                     black, [xCenter yCenter], 2, []);
             end
             vbl = Screen('Flip', windowExp, vbl + (waitframes - 0.5) * ifi);
         end
                  
    Screen('BlendFunction', windowExp, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    % Flip again to sync us to the vertical retrace 
    Screen('DrawDots', windowExp, [xCenter; yCenter], 10, black, [], 2);
    vbl = Screen('Flip', windowExp);
 
    % response loop
    timeStart = vbl;
    respToBeMade = true;
    while respToBeMade
        [keyIsDown, secs, keyCode] = KbCheck;
        checkAbort = GetSecs - timeStart;
        
        if keyCode(escapeKey)
            ShowCursor;
            sca;
            return
       
        elseif keyCode(leftKey)    
            response = -1;
            response_time(k) = 1000 .* (GetSecs - timeStart);
            respToBeMade = false;
            
        elseif keyCode(rightKey)
            response = 1;
            response_time(k) = 1000 .* (GetSecs - timeStart);
            respToBeMade = false;
        
        elseif checkAbort >= 4
            response = 2;
            response_time(k) = -1;
            respToBeMade = false;
        
        end    
             
    end % end of response loop
    
    vbl = Screen('Flip', windowExp, vbl + (waitframes - 0.5) * ifi);
    response = response * 8;
    check = 16;
    
    if response == targetTilt && trialNo <= nTrain
            DrawFormattedText(windowExp, 'Correct!', 'center',...
                'center', [0, 1, 0]);
            vbl = Screen('Flip', windowExp);
            WaitSecs(1);
            
        elseif response == -targetTilt && trialNo <= nTrain
            DrawFormattedText(windowExp, 'Incorrect!', 'center',...
                'center', [1, 0, 0]);
            vbl = Screen('Flip', windowExp);
            WaitSecs(1);
            
        elseif abs(response) == check 
            DrawFormattedText(windowExp, 'Timeout! Trial aborted...',... 
                'center', 'center', [0, 0, 0]);
            vbl = Screen('Flip', windowExp);
            k = k - 1
            WaitSecs(1);
        
        else
            DrawFormattedText(windowExp, 'Ready...', 'center',...
                'center', [0, 0, 0]);
            vbl = Screen('Flip', windowExp);
            WaitSecs(1);
            
    end % end of feedback display loop 
    
    vbl = Screen('Flip', windowExp, vbl + (waitframes - 0.5) * ifi);
        if k ~= 0
           keyResp(k) = response; 
        end
    checkAbort = 0;
    end % end of experimental loop
    
%----------------------------------------------------------------------
%                   Data Pre-Processing and Output
%----------------------------------------------------------------------
    
    FullExpSeq(:,3) = FullExpSeq(:,3) * 15;
    FullExpSeq(:,4) = FullExpSeq(:,4) * 8;
    FullExpSeq = [l, FullExpSeq, keyResp, response_time];
    for i = 1:length(l)
        if FullExpSeq(i,5) == FullExpSeq(i,6)
           corr(i) = "correct";
        elseif (FullExpSeq(i,5) * -1) == FullExpSeq(i,6)
           corr(i) = "incorrect";
        elseif isnan(FullExpSeq(i,5))
           corr(i) = "timeout";
        end
    end    
    
    format bank;          
    idx_subset_1 = find(FullExpSeq(:,4) <0 & FullExpSeq(:,5) <0);
    idx_subset_2 = find(FullExpSeq(:,4) >0 & FullExpSeq(:,5) >0);
    congr_trials_midx = sort([idx_subset_1; idx_subset_2]);
    incon_trials_midx = setdiff(l, congr_trials_midx);
    catg_congruency = strings(40,1);
    catg_congruency(congr_trials_midx) = "congruent";
    catg_congruency(incon_trials_midx) = "incongruent";
    FullExpSeq = FullExpSeq(:, 2:5);
    trial_num = l;
    data_atilt = table(trial_num, FullExpSeq, catg_congruency,...
        response_time, corr)   
    ShowCursor;
    sca;
 
        