%%%start kinect
function [vid, src] = startkinect()
global VERBOSE
try
    vid = videoinput('kinect',2,'Depth_640x480');
    % Set parameters for video
    % Acquire only one frame each time
    set(vid,'FramesPerTrigger',1);
    % Go on forever until stopped
    set(vid,'TriggerRepeat',Inf);
    % Get a grayscale image
    set(vid,'ReturnedColorSpace','grayscale');
    triggerconfig(vid, 'Manual');
    %set(vid,'TrackingMode','Skeleton')
catch  ME
    if VERBOSE
    ME.getReport
    warning('Couldnt open kinect. Maybe it is already running?')
    end
    try
        vid = imaqfind; %in case i am already aquiring
    catch  ME2
        warning('No Kinect. trying for a normal camera????')
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
src = getselectedsource(vid);
src.TrackingMode = 'Skeleton';
src.CameraElevationAngle = 0;
