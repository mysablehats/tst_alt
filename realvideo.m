function [labellabel, allskel3] = realvideo(gases, arc_conn, simvar, runvars )
global VERBOSE
VERBOSE = false;
realtimevideo = false;
justwrite = false;

chunksize = 60;
chunk.chunk = zeros(20,3,chunksize);
chunk.timers = zeros(1,chunksize, 'uint64');
chunk.times = zeros(1,chunksize);
chunk.counter = 0;

allskel3 = [];

% Set-up webcam video input
[vid, ~] = startkinect();

for i=1:length(runvars)
    switch runvars{i}
        case 'realtime'
            realtimevideo = 1;
        case 'justwrite'
            justwrite = 1;
    end
end

if realtimevideo
    % Open figure
    hFigure=figure();
    % Define frame rate
    NumberFrameDisplayPerSecond=10;
    % set up timer object
    TimerData=timer('TimerFcn', {@FrameRateDisplay,vid, chunk, gases, arc_conn, simvar, justwrite},'Period',1/NumberFrameDisplayPerSecond,'ExecutionMode','fixedRate','BusyMode','drop');
end
% Start video object
start(vid);

if realtimevideo
    start(TimerData)
    uiwait(hFigure);
end

% Clean up everything
if exist('TimerData', 'var')
    
    %     %pause(1)
    %     while(TimerData.Running)
    %     end
    stop(TimerData);
    %waitfor(TimerData);
    delete(TimerData);
    stop(vid);
    pause(1)
    delete(vid);
    % clear persistent variables
    clear vid
    clear TimerData
    clear functions;
    %clear
elseif ~realtimevideo
    labellabel = [];
    while(isempty(labellabel))
        [chunk, labellabel] = FrameRateDisplay([], [],vid, chunk, gases, arc_conn, simvar, justwrite);        
    end
    stop(vid);
    pause(1)
    delete(vid);
    % clear persistent variables
    clear functions;
    %clear
end

end

% This function is called by the timer to display one frame of the figure
function [chunk, labellabel] = FrameRateDisplay(obj, event,vid,chunk, gases, arc_conn, simvar, justwrite)
persistent IM; % im not sure this is necessary
labellabel = [];
try
    trigger(vid);
    [IM,~,metaData]=getdata(vid,1,'uint8');
    [skelskel, chunk ] = readskeleton(metaData, chunk);
    if chunk.counter > size(chunk.chunk,3)
     
        allskel3 = generate_skel_online(chunk.chunk);
        save(savefilesave2('onlineclass', simvar.env),'allskel3','arc_conn', 'simvar','chunk')
        if justwrite        
            %error('I finished.')
            labellabel = 'hi';
        else
            labellabel = online_classifier(gases,allskel3, arc_conn, simvar);
        end
    end
    makeimage(IM, skelskel)
catch ME
    ME.getReport
    return
end
end
function [skelskel, chunk]= readskeleton(metaData, chunk)
skelskel = [];
if any(metaData.IsSkeletonTracked)==1
    dbgmsg(strcat('Tracked: ',num2str(sum(metaData.IsSkeletonTracked)),' skeletons.'),0)
    for i = 1:length(metaData.IsSkeletonTracked)
        if metaData.IsSkeletonTracked(i)==1
            dbgmsg('Reached inside of loop',0)
            skelskel =  coordshift(skeldraw_(metaData.JointWorldCoordinates(:,:,i),false));
            chunk.chunk(:,:,2:end) = chunk.chunk(:,:,1:end-1);
            chunk.chunk(:,:,1) = metaData.JointWorldCoordinates(:,:,i);
            chunk.counter = chunk.counter +1;
            if chunk.counter>1
                chunk.times(chunk.counter-1) = toc(chunk.timers(chunk.counter-1));
            end
            chunk.timers(chunk.counter) = tic;
        end
    end
end
end
function makeimage(IM, skelskel)
persistent handlesRaw;
persistent handlesPlot;
persistent handlesmyskel;
persistent myaxes;
if isempty(handlesRaw)
    % if first execution, we create the figure objects
    %subplot(2,1,1);
    handlesRaw=imagesc(IM);
    title('CurrentImage');
    hold on
    % Plot first value
    %Values=mean(IM(:));
    %subplot(2,2,2);
    %handlesPlot=plot(Values);
    %title('Average of Frame');
    %xlabel('Frame number');
    %ylabel('Average value (au)');
    
    %my skeleton
    sampleskel = [0.0697    0.1773    1.6761;
        0.0756    0.2420    1.6839;
        0.0678    0.5732    1.6773;
        0.0010    0.7354    1.5891;
        -0.0791    0.4813    1.7441;
        -0.1515    0.3129    1.4867;
        -0.1649    0.2954    1.2563;
        -0.1067    0.2954    1.2395;
        0.2255    0.4464    1.5866;
        0.2237    0.2958    1.4024;
        -0.0567    0.2926    1.2936;
        -0.1113    0.3175    1.2855;
        0.0002    0.1035    1.7097;
        0.0094   -0.3763    1.6717;
        -0.1009   -0.6928    1.6735;
        -0.1548   -0.7310    1.5964;
        0.1391    0.0965    1.6379;
        0.1254   -0.3289    1.6740;
        0.1750   -0.4106    1.2409;
        0.2620   -0.3947    1.1648];
    
    handlesmyskel = plot3(sampleskel(1,:),-sampleskel(2,:), sampleskel(3,:));
    set(handlesmyskel, 'LineWidth',4, 'Color', 'y');
    %hold off
    %subplot(2,1,2);
    %handlesmyskel=plot([]);
    %    try
    %        [~ , handlesmyskel] = skeldraw(sampleskel,true);
    %    catch
    %        disp('cant initialize axes handle')
    %    end
    myaxes = gca; %get(handlesmyskel,'Parent');
    %set(myaxes, 'color', 'none')
else
    % We only update what is needed
    set(handlesRaw,'CData',IM);
    %Value=mean(IM(:));
    %OldValues=get(handlesPlot,'YData');
    %set(handlesPlot,'YData',[OldValues Value]);
    %%%
    if exist('skelskel','var')&&~isempty(skelskel)
        %plot3(skelskel(1,:),skelskel(2,:), skelskel(3,:))
        %disp('reachedplot')
        set(handlesmyskel, 'XData',skelskel(1,:),'YData',skelskel(2,:),'ZData', skelskel(3,:))
        %set(handlesmyskel,'CData',skelskel)
        set(myaxes,'XLim', [0 640]);
        set(myaxes,'YLim', [0 480]);
        %         set(myaxes,'ZLim', [-0 5]);
        view(0,90);
    end
end
end
function b =  coordshift(a)
b = zeros(size(a));
b(1,:) = a(1,:)*240 +320;
b(2,:) = -a(2,:)*240 +240;
b(3,:) = a(3,:)*240;
end