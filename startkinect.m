%%%start kinect
function [vid, src] = startkinect()
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
src = getselectedsource(vid);
src.TrackingMode = 'Skeleton';
src.CameraElevationAngle = 0