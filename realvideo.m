function [sksks, allskel3] = realvideo(outstruct, arc_conn, simvar, realtimevideo)
global VERBOSE
VERBOSE = false;
% persistent vid
% if exist('vid','var')
%     closepreview(vid)
% end
% Define frame rate
NumberFrameDisplayPerSecond=10;

% Open figure
hFigure=figure();

% Set-up webcam video input
try
   [vid, src] = startkinect();
    
catch  ME
    ME.getReport
    disp('Couldnt open kinect')
    try
        vid = imaqfind; %in case i am already aquiring
    catch  ME2
                 ME2.getReport
        try
            % For macs.
            % this is dumb
            vid = videoinput('macvideo', 1);
        catch  ME3
                 ME3.getReport
            errordlg('No webcam available');
        end
    end
end

if realtimevideo
    % set up timer object
    TimerData=timer('TimerFcn', {@FrameRateDisplay,vid,outstruct, arc_conn, simvar},'Period',1/NumberFrameDisplayPerSecond,'ExecutionMode','fixedRate','BusyMode','drop');
end
% Start video and timer object
try %maybe I have already started vid... or failed to stop it? 
start(vid);
catch
    stop(vid);
    start(vid);
end
%preview(vid)
%try
if realtimevideo
    start(TimerData)
    uiwait(hFigure);
end
%catch
%     stop(TimerData);
%     delete(TimerData);
%     stop(vid);
%     delete(vid);
%     % clear persistent variables
%     clear functions;
% end
% We go on until the figure is closed



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
    sksks = [];
    chunksize = 500;
    chunk.chunk = zeros(20,3,chunksize);
    chunk.timers = zeros(1,chunksize, 'uint64');
    chunk.times = zeros(1,chunksize);
    chunk.counter = 0;
    while(isempty(sksks))
        [sksks, allskel3, chunk] = FrameRateDisplay([], [],vid,outstruct, arc_conn, simvar, chunk);
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


function [skeldata, allskel3, chunk] = FrameRateDisplay(obj, event,vid,outstruct, arc_conn, simvar, chunk)
persistent IM; % im not sure this is necessary
skeldata = [];
allskel3 = [];
try
trigger(vid);
[IM,~,metaData]=getdata(vid,1,'uint8');

[skelskel, skeldata, allskel3, chunk ] = readskeleton(metaData, outstruct, arc_conn, simvar, chunk);
makeimage(IM, skelskel)
catch ME
    ME.getReport
    return
end

end
function [skelskel, labellabel, allskel3, chunk ]= readskeleton(metaData, outstruct, arc_conn, simvar, chunk)
skelskel = [];
labellabel = [];
allskel3 = [];

%try
if any(metaData.IsSkeletonTracked)==1
    disp(strcat('Tracked: ',num2str(sum(metaData.IsSkeletonTracked)),' skeletons.'))
%    dbgmsg(metaData.IsSkeletonTracked,0)
    for i = 1:length(metaData.IsSkeletonTracked)
        if metaData.IsSkeletonTracked(i)==1
            %disp(metaData.JointWorldCoordinates(:,:,i))
            %try
                dbgmsg('Reached inside of loop',0)
                skelskel =  coordshift(skeldraw_(metaData.JointWorldCoordinates(:,:,i),false));
                chunk.chunk(:,:,2:end) = chunk.chunk(:,:,1:end-1);
                chunk.chunk(:,:,1) = metaData.JointWorldCoordinates(:,:,i);                
                chunk.counter = chunk.counter +1;
                if chunk.counter>1
                    chunk.times(chunk.counter-1) = toc(chunk.timers(chunk.counter-1)); 
                end
                chunk.timers(chunk.counter) = tic;
                if chunk.counter > size(chunk.chunk,3)
                    allskel3 = generate_skel_online(chunk.chunk);
                    save(savefilesave2('onlineclass', simvar.env), 'outstruct','allskel3','arc_conn', 'simvar','chunk')
                    labellabel = online_classifier(outstruct,allskel3, arc_conn, simvar);
                end
           % catch ME
%                 ME.getReport
%                 size(chunk(:,:,1))
%                 size(metaData.JointWorldCoordinates(:,:,i))
%                 find(chunk==0)
%                 disp(metaData.JointWorldCoordinates(:,:,i))
%                 error('Something fishy happened') %why error catch if you are goint to break the program??
%             %end
        end
    end
end
%catch
%    disp('Can''t draw! :/')
%end
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