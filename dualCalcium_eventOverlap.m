function dualCalcium_eventOverlap ( Hz )

% script that analyses ROIs for two channels at the same time
% from excel files showing changes to mean gray values 
% will have one ROI for the green channel and one ROI for the red
% channel, so each graph will be composed of two traces

% Hz = sampling rate (frames/second) of calcium imaging

% choose folder where all data is stored
myBigFolder = '2.temp excel data';
myBigFolderInfoPre = dir(myBigFolder);
myBigFolderInfo = myBigFolderInfoPre(4:end);

% run analysis for data in each subfolder
for kk = 1:length(myBigFolderInfo)

    % clear variables
    myFolder = [];
    fPattern_green = [];
    fPattern_red = [];
    files_green = [];
    files_red = [];
    trend_traces = [];
    ROI_normTracesAll = [];
    ROI_normTracesAll_events = [];

    % choose subfolder for each experiment
    myFolder = myBigFolderInfo(kk).name;

    % load matrices from excel data corresponding to green and red fluorescence
    fPattern_green = fullfile( myBigFolder,myFolder, '*GREEN.xlsx' );
    fPattern_red = fullfile( myBigFolder, myFolder, '*RED.xlsx' );
    files_green = dir( fPattern_green );
    files_red = dir( fPattern_red );
    % for dir to work, need to be on the MATLAB folder as current folder

    for noROI = 1:length( files_green )
        ROIdata_greenAll{noROI} = xlsread( files_green(noROI).name );
        ROIdata_redAll{noROI} = xlsread( files_red(noROI).name );
    end

    % open for loop so that each ROI is done individually
    for h = 1:length( ROIdata_greenAll )

        % firstly clear repeated variables
        ROI_fluoData = [];
        trend_traces = [];
        ROIdata_green = [];
        ROIdata_red = [];
        dF_overF = [];
        time = [];

        %load data for new ROI
        ROIdata_green = ROIdata_greenAll{h};
        ROIdata_red = ROIdata_redAll{h};

        % convert timepoints (as frames - in first colum) to time in sec
        time = ( ROIdata_green( :,1 ) ./ Hz ) / 60;
        ROI_green = ROIdata_green( :,2:end );
        ROI_red = ROIdata_red( :,2:end );

        % make matrix that includes fluorescence changes from OPC (GCaMP)
        % and axonal (RGECO) ROIs
        % first column is OPC GCaMP (green channel) and second is axonal
        % RGECO (red channel)
        ROI_fluoData(:,1) = ROI_green(:,1);
        ROI_fluoData(:,2) = ROI_red(:,2);
        
        % :::::::::::::::::::::::::: MULTIPLE ACTIVITY :::::::::::::::::::::::::: %

        % Calculating dF/F traces

        % Concentrating analysis on green vs red channels
        for GreRed = 1:size( ROI_fluoData,2 )

            ROIno_fluoData = [];
            ROIno_fluoData = ROI_fluoData(:,GreRed);

            % Getting dF/F of for each timepoint

            % setting size of moving window to calculate baseline
            % fluorescence
            window = 6;

            for t_one = 1:size( ROI_fluoData,1 )
                % Calculating F0 as the minimum value from the previous n
                % timepoints in the trace

                if t_one <= window

                    fZero = min(ROIno_fluoData(1:t_one));
                else
                    fZero = min(ROIno_fluoData(t_one-window:t_one-1));
                    deltaF = ROIno_fluoData(t_one) - fZero;

                    dF_overF(t_one,GreRed) = deltaF./fZero;
                end

            end
        end


        % Separately (but in same figure) plot the OPC and neuron traces:

        figure(1)
        subplot(2,1,1)
        plot( time(1:end-1), dF_overF(1:end-1,1), "color", [0.6510 0.6510 0.6510])
        title(files_green(h).name)
        ylim([-0.1 1])
        subplot(2,1,2)
        plot( time(1:end-1), dF_overF(1:end-1,2), "color", [0.8000 0.0549 0.1922])
        ylim([-0.1 0.5])

 
        % ::::::::::::::::::::::::::: EVENT DETECTION ::::::::::::::::::::::::::: %

        % Set threshold for calcium signal
        threshold_GCaMP = max(1.*std(dF_overF(:,1)),0.3);
        threshold_RGECO_gad1b = max(4.*std(dF_overF(:,2)),0.08);
        threshold_RGECO_vglut2a = max(3.*std(dF_overF(:,2)),0.08); 

        %:: Event isolation using "findpeaks" :: %

        % GCaMP ROI
        [~, xLocG, widthG, ~] = findpeaks(dF_overF(:,1), time,...
            'MinPeakDistance',0.05,'MinPeakHeight', threshold_GCaMP);

        % RGECO ROI
        [~, xLocR, widthR, ~] = findpeaks(dF_overF(:,2), time,...
            'MinPeakDistance',0.05,'MinPeakHeight', threshold_RGECO_gad1b);

        % Determine axis units for plotting rectangles indicating detected
        % signal

        % starting corner of rectangle

        % starting point will be the half-width location from peak
        XstartP_G = ( xLocG - ( widthG ) ) - 0.05 ;
        XstartP_R = ( xLocR - ( widthR./2 ) );
        
        % starting point for the rectangle indicating time after the OPC event
        XstartP_G_after = XstartP_G + widthG + 0.05;

        % Plot ROI traces showing event coincidence
        
        figure(2)
        
        % plot GCaMP rectangles
        for j=1:length(xLocG)
            rectangle('Position', [XstartP_G(j),-0.1,(widthG(j)+0.05),1.2],...
                'EdgeColor',[0.9020 0.9020 0.9020],'FaceColor',[0.9020 0.9020 0.9020]);
            rectangle('Position', [XstartP_G_after(j),-0.1,(widthG(j)+0.05),1.2],...
                'EdgeColor',[0.9294 0.8667 0.6039],'FaceColor',[0.9294 0.8667 0.6039]);
        end
        % plot RGECO rectangles
        for jj=1:length(xLocR)
            rectangle('Position', [XstartP_R(jj),-0.1,(widthR(jj)),1.2],...
                'EdgeColor',[0.9216 0.5961 0.6588],'FaceColor',[0.9216 0.5961 0.6588]);
        end
        hold on
        % overlay activity traces
        
        % plot GCaMP trace
        plot( time, dF_overF(:,1), "color", "black")
        title(files_green(h).name)
        ylim([-0.1 1])
        xlim([0 1])

        % plot RGECO trace
        plot( time, dF_overF(:,2), "color", [0.8000 0.0549 0.1922])
        ylim([-0.1 0.5])
        xlim([0 1])
        hold off

    end
    

end



