function [conformstruc, skelldef] = conformskel(varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbgmsg('Applies preconditioning functions to both training and validation datasets',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = false;
skelldef = struct();
if isempty(varargin)||strcmp(varargin{1},'test')
    return
else        
    if isstruct(varargin{1})
        outputisstruc = true;
        conformstruc = varargin{1};
        data_train = conformstruc.train.data;
        data_val = conformstruc.val.data;
        data_ytrain = conformstruc.train.y ;
        data_yval = conformstruc.val.y;
        
        if isfield(conformstruc,'awk')
            awk = conformstruc.awk;
        else
            awk = generate_awk(data_val);
            dbgmsg('awk not defined. considering all joints as having equal importance.',1)
        end
        lindx = 2;
        singleskelset = false;
    else
        outputisstruc = false;
        data_train = varargin{1};
        if nargin>2
            data_val = varargin{2};            
            if isnumeric(varargin{3})
                awk = varargin{3};
                lindx = 4;
            elseif ischar(varargin{3})
                awk = generate_awk(data_val);
                lindx = 3;
                dbgmsg('awk not defined. considering all joints as having equal importance.',1)
            else
                error('Strange input. I don''t know what to do with it. ')
            end
            singleskelset = false;
        else
            awk = generate_awk(data_train);
            lindx = 2;
            singleskelset = true;
        end 
        %awk = generate_awk;
        data_ytrain = []; %% I need to change this if I plan on increasing the size of the data, such as with the mirror func
        data_yval = []; %%same        
    end
    %%% initiallize variables to make skelldef
    killdim = [];
    skelldef.realkilldim = [];
    skelldef.length = size(data_train,1);
    skelldef.elementorder = 1:skelldef.length;
    
    switch skelldef.length
        case {75,60,45}
            skelldef.awk.pos = repmat(awk(setdiff(1:skelldef.length/3,killdim)),3,1);
            skelldef.awk.vel = [];
            skelldef.novel = true;
        case {150,120,90}
            skelldef.awk.pos = repmat(awk(setdiff(1:skelldef.length/6,killdim)),3,1);
            skelldef.awk.vel = repmat(awk(setdiff(1:skelldef.length/6,killdim)),3,1);
            skelldef.novel = false;
    end
    
    %%%
    skelldef.bodyparts = genbodyparts(skelldef.length);
    
    %%%errorkarling
    if ~singleskelset&&(size(data_val,1)~=skelldef.length)
        error('data_train and data_val must have the same length!!!')
    end

    % creates the function handle cell array
    conformations = {};
    killdim = [];
    
    for i =lindx:length(varargin)
        switch varargin{i}
            case 'test'
                test = true;
            case 'highhips'
                conformations = [conformations, {@highhips}]; %#ok<*AGROW>
            case 'nohips'
                conformations = [conformations, {@centerhips}];
                killdim = [killdim, skelldef.bodyparts.hip_center];
            case 'normal'
                conformations = [conformations, {@normalize}];
                %dbgmsg('Unimplemented normalization: ', varargin{i} ,true);
            case 'mirrorx'
                conformations = [conformations, {@mirrorx}];
            case 'mirrory'
                conformations = [conformations, {@mirrory}];
            case 'mirrorz'
                conformations = [conformations, {@mirrorz}];
            case 'mahal'
                %conformations = [conformations, {@mahal}];
                dbgmsg('Unimplemented normalization: ', varargin{i} ,true);
            case 'norotate'
                conformations = [conformations, {@norotatehips}];
                dbgmsg('WARNING, the normalization: ' , varargin{i},' is performing poorly, it should not be used.', true);
            case 'norotatehips'
                conformations = [conformations, {@norotatehips2}];
                %dbgmsg('WARNING, the normalization: ' , varargin{i},' is performing poorly, it should not be used.', true);
            case 'norotateshoulders'
                conformations = [conformations, {@norotateshoulders}];
                dbgmsg('WARNING, the normalization: ' , varargin{i},' is performing poorly, it should not be used.', true);
            case 'notorax'
                conformations = [conformations, {@centertorax}];
                dbgmsg('WARNING, the normalization: ' , varargin{i},' is performing poorly, it should not be used.', true);
                killdim = [killdim, skelldef.bodyparts.TORSO];
            case 'nofeet'
                conformations = [conformations, {@nofeet}]; %not sure i need this...
                killdim = [killdim, skelldef.bodyparts.RIGHT_FOOT, skelldef.bodyparts.LEFT_FOOT];
            case 'nohands'
                dbgmsg('WARNING, the normalization: ' , varargin{i},' is performing poorly, it should not be used.', true);
                killdim = [killdim, skelldef.bodyparts.RIGHT_HAND, skelldef.bodyparts.LEFT_HAND];
            case 'axial'
                %conformations = [conformations, {@axial}];
                dbgmsg('Unimplemented normalization: ', varargin{i} ,true);
            case 'addnoise'
                conformations = [conformations, {@abnormalize}];
            case 'spherical'
                conformations = [conformations, {@to_spherical}];
            case 'intostick'
                conformations = [conformations, {@intostick}];
                killdim = [4:(skelldef.length/6) (skelldef.length/6+4):(skelldef.length/3) ];
            case 'intostick2'
                conformations = [conformations, {@intostick2}];
                killdim = [4:(skelldef.length/6) (skelldef.length/6+4):(skelldef.length/3) ];
            otherwise
                dbgmsg('ATTENTION: Unimplemented normalization/ typo.',varargin{i},true);
        end
    end
    
    % execute them for training and validation sets
    if ~isempty(skelldef.bodyparts)
        for i = 1:length(conformations)
            func = conformations{i};
            dbgmsg('Applying normalization: ', varargin{i+lindx-1},true);
            if isequal(func, @mirrorx)||isequal(func,@mirrory)||isequal(func, @mirrorz)
                data_trainmirror = data_train;
                data_ytrainmirror = data_ytrain;
                if ~singleskelset
                    data_valmirror = data_val;
                    data_yvalmirror = data_yval;
                end
            else
                data_trainmirror = [];                
                data_ytrainmirror = [];
                if ~singleskelset
                    data_valmirror = [];
                    data_yvalmirror = [];
                end
            end
            
            if isequal(func,@normalize)
                %%% must go through whole dataset!
                %%% if there is ever another function that requires this,
                %%% then I should probably use a switch - if that works...
                allskels = makefatskel(data_train);
                %%% calculating the magic constants for our data
                if skelldef.novel
                    vectdata_pos = reshape(allskels(1:skelldef.length/3,:,:),1,[]);
                    skelldef.pos_std = std(vectdata_pos);
                    skelldef.pos_mean=mean(vectdata_pos);
                    skelldef.vel_std = [];
                    skelldef.vel_mean=[];
                else
                    vectdata_pos = reshape(allskels(1:skelldef.length/6,:,:),1,[]);
                    skelldef.pos_std = std(vectdata_pos);
                    skelldef.pos_mean=mean(vectdata_pos);
                    vectdata_vel = reshape(allskels((skelldef.length/6+1):end,:,:),1,[]);
                    skelldef.vel_std = std(vectdata_vel);
                    skelldef.vel_mean=mean(vectdata_vel);
                end
                %                 case {@mirrorx,@mirrory,@mirrorz}
                %                     data_trainmirror = data_train;
                %                     data_valmirror = data_val;
            end
            for j = 1:size(data_train,2)
                [tdskel,skelldef.hh] = makefatskel(data_train(:,j));
                data_train(:,j) = makethinskel(func(tdskel, skelldef));
            end
            data_train = [data_train data_trainmirror];
            data_ytrain = [data_ytrain data_ytrainmirror];
            if ~singleskelset
                for j = 1:size(data_val,2)
                    [tdskel,skelldef.hh] = makefatskel(data_val(:,j));
                    data_val(:,j) = makethinskel(func(tdskel, skelldef));
                end
                data_val = [data_val data_valmirror];
                data_yval = [data_yval data_yvalmirror];
            end
        end
    end
    % squeeze them accordingly?
    if ~test
        whattokill = reshape(1:skelldef.length,skelldef.length/3,3);
        realkilldim = whattokill(killdim,:);
        conform_train = data_train(setdiff(1:skelldef.length,realkilldim),:); %sorry for the in-liners..
        if ~singleskelset
            conform_val = data_val(setdiff(1:skelldef.length,realkilldim),:);
        end
        skelldef.elementorder = skelldef.elementorder(setdiff(1:skelldef.length,realkilldim));
    else
        conform_train = data_train;
        if ~singleskelset
            conform_val = data_val;
        end
    end
    if outputisstruc
        conformstruc.train.data = conform_train;
        conformstruc.val.data = conform_val;
        conformstruc.train.y = data_ytrain;
        conformstruc.val.y = data_yval;
    else
        conformstruc = conform_train;
    end
end
skelldef.realkilldim = realkilldim;
[skelldef.pos, skelldef.vel] = generateidx(skelldef.length, skelldef);

end
function newskel = centerhips(tdskel, skelldef)
bod = skelldef.bodyparts;
if isempty(bod.hip_center)&&skelldef.hh==30
    %%% then this possibly it is the 15 joint skeleton
    hip = (tdskel(bod.LEFT_HIP,:) + tdskel(bod.RIGHT_HIP,:))/2;
else
    hip = tdskel(bod.hip_center,:);
end
if skelldef.novel
    hips = repmat(hip,skelldef.hh,1); % this is so that we dont subtract the velocities
else
    hips = [repmat(hip,skelldef.hh/2,1);zeros(skelldef.hh/2,3)]; % this is so that we dont subtract the velocities
end
newskel = tdskel - hips;
end
function newskel = highhips(tdskel, skelldef)
bod = skelldef.bodyparts;
if isempty(bod.hip_center)&&skelldef.hh==30
    %%% then this possibly it is the 15 joint skeleton
    hip = (tdskel(bod.LEFT_HIP,:) + tdskel(bod.RIGHT_HIP,:))/2;
else
    hip = tdskel(bod.hip_center,:);
end
hip(1,2) = 0; %%
if skelldef.novel
    hips = repmat(hip,skelldef.hh,1);
else
    hips = [repmat(hip,skelldef.hh/2,1);zeros(skelldef.hh/2,3)]; % this is so that we dont subtract the velocities
end
newskel = tdskel - hips;
end
function newskel = centertorax(tdskel, skelldef)
bod = skelldef.bodyparts;
if skelldef.novel
    torax = repmat(tdskel(bod.TORSO,:),skelldef.hh,1);
else
    torax = [repmat(tdskel(bod.TORSO,:),skelldef.hh/2,1);zeros(skelldef.hh/2,3)];
end
newskel = tdskel - torax;
end
function tdskel = norotatehips(tdskel, skelldef)

bod = skelldef.bodyparts;
rvec = tdskel(bod.RIGHT_HIP,:)-tdskel(bod.LEFT_HIP,:);
rotmat = vecRotMat(rvec/norm(rvec),[1 0 0 ]); % arbitrary direction. hope I am not doing anything too stupid
for i = 1:skelldef.hh
    tdskel(i,:) = (rotmat*tdskel(i,:)')';
end
end
function tdskel = norotatehips2(tdskel, skelldef)
a = detectrotation(tdskel, skelldef, 'hips');
tdskel = rotskel(tdskel,0,a,0);
end
function [angle, vectvect] = detectrotation(skel, skelldef, joint)
%angle = 0;
switch joint
    case 'hips'
        vectvect = skel(skelldef.bodyparts.LEFT_HIP,:) - skel(skelldef.bodyparts.RIGHT_HIP,:);
    case 'torax'
        vectvect = skel(skelldef.bodyparts.LEFT_SHOULDER,:) - skel(skelldef.bodyparts.RIGHT_SHOULDER,:);
end
[angle, ~, ~] = cart2pol(vectvect(1), vectvect(3), vectvect(2));
end
function tdskel = norotateshoulders(tdskel, skelldef)
bod = skelldef.bodyparts;
rvec = tdskel(bod.RIGHT_SHOULDER,:)-tdskel(bod.LEFT_SHOULDER,:);
rotmat = vecRotMat(rvec/norm(rvec),[1 0 0 ]); % arbitrary direction. hope I am not doing anything too stupid
for i = 1:skelldef.hh
    tdskel(i,:) = (rotmat*tdskel(i,:)')';
end
end
function newskel = abnormalize(tdskel, ~)
newskel = tdskel + rand(size(tdskel));
end
function tdskel = nofeet(tdskel, skelldef)
bod = skelldef.bodyparts;
sizeofnans = size(tdskel([bod.RIGHT_FOOT, bod.LEFT_FOOT],:));
tdskel([bod.RIGHT_FOOT, bod.LEFT_FOOT],:) = NaN(sizeofnans);
end
function tdskel = mirrorx(tdskel, ~)
tdskel(:,1) = -tdskel(:,1);
end
function tdskel = mirrory(tdskel, ~)
tdskel(:,2) = -tdskel(:,2);
end
function tdskel = mirrorz(tdskel, ~)
tdskel(:,3) = -tdskel(:,3);
end
function tdskel = normalize(tdskel, skelldef)
if skelldef.novel
    for i = 1:skelldef.hh
        tdskel(i,:) = (tdskel(i,:) - skelldef.pos_mean)/skelldef.pos_std; %- skelldef.pos_mean
    end
else
    for i = 1:skelldef.hh/2
        tdskel(i,:) = (tdskel(i,:) - skelldef.pos_mean)/skelldef.pos_std; %- skelldef.pos_mean
    end
    for i = (skelldef.hh/2+1):skelldef.hh
        tdskel(i,:) = (tdskel(i,:) - skelldef.vel_mean)/skelldef.vel_std;%- skelldef.vel_mean
    end
end
end
function newskel = to_spherical(tdskel, ~)
newskel = zeros(size(tdskel));
[newskel(:,1),newskel(:,2),newskel(:,3)] = cart2sph(tdskel(:,1),tdskel(:,2),tdskel(:,3));
end
function zeroskel = intostick(tdskel, skelldef)
bod = skelldef.bodyparts;
UCI = [bod.HEAD bod.NECK bod.LEFT_SHOULDER  bod.RIGHT_SHOULDER bod.LEFT_ELBOW bod.RIGHT_ELBOW bod.LEFT_HAND bod.RIGHT_HAND  ];
uppercentroid = mean(tdskel(UCI ,:));
uppercentroidvel  = mean(tdskel(UCI+skelldef.hh/2,:));
middlecentroid = tdskel(bod.TORSO,:);
middlecentroidvel = tdskel(bod.TORSO+skelldef.hh/2,:);
LCI = [bod.LEFT_FOOT bod.RIGHT_FOOT bod.LEFT_KNEE bod.RIGHT_KNEE bod.LEFT_HIP bod.RIGHT_HIP];
lowercentroid =mean(tdskel(LCI,:));
lowercentroidvel =mean(tdskel(LCI+skelldef.hh/2,:));
zeroskel = zeros(size(tdskel));
zeroskel(1:3,:) = [uppercentroid;middlecentroid;lowercentroid];
zeroskel((skelldef.hh/2+1):(skelldef.hh/2+3),:) = [uppercentroidvel;middlecentroidvel;lowercentroidvel];
end
function zeroskel = intostick2(tdskel, skelldef)
bod = skelldef.bodyparts;
UCI = [bod.HEAD bod.NECK bod.LEFT_SHOULDER  bod.RIGHT_SHOULDER bod.LEFT_ELBOW bod.RIGHT_ELBOW  ];
uppercentroid = mean(tdskel(UCI ,:));
uppercentroidvel  = mean(tdskel(UCI+skelldef.hh/2,:));
middlecentroid = tdskel(bod.TORSO,:);
middlecentroidvel = tdskel(bod.TORSO+skelldef.hh/2,:);
LCI = [ bod.LEFT_KNEE bod.RIGHT_KNEE bod.LEFT_HIP bod.RIGHT_HIP];
lowercentroid =mean(tdskel(LCI,:));
lowercentroidvel =mean(tdskel(LCI+skelldef.hh/2,:));
zeroskel = zeros(size(tdskel));
zeroskel(1:3,:) = [uppercentroid;middlecentroid;lowercentroid];
zeroskel((skelldef.hh/2+1):(skelldef.hh/2+3),:) = [uppercentroidvel;middlecentroidvel;lowercentroidvel];
end
function bodyparts = genbodyparts(lenlen)
bodyparts = struct();
switch lenlen
    case {150,75}
        bodyparts.SpineBase = 1;
        bodyparts.SpineMid = 2;
        bodyparts.Neck = 3;
        bodyparts.Head = 4;
        bodyparts.ShoulderLeft = 5;
        bodyparts.ElbowLeft = 6;
        bodyparts.WristLeft = 7;
        bodyparts.HandLeft = 8;
        bodyparts.ShoulderRight = 9;
        bodyparts.ElbowRight = 10;
        bodyparts.WristRight = 11;
        bodyparts.HandRight = 12;
        bodyparts.HipLeft = 13;
        bodyparts.KneeLeft = 14;
        bodyparts.AnkleLeft = 15;
        bodyparts.FootLeft = 16;
        bodyparts.HipRight = 17;
        bodyparts.KneeRight = 18;
        bodyparts.AnkleRight = 19;
        bodyparts.FootRight = 20;
        bodyparts.SpineShoulder = 21;
        bodyparts.HandTipLeft = 22;
        bodyparts.ThumbLeft = 23;
        bodyparts.HandTipRight = 24;
        bodyparts.ThumbRight = 25;
        %%% synonyms
        bodyparts.hip_center = bodyparts.SpineBase;
        bodyparts.NECK = bodyparts.Neck;
        bodyparts.RIGHT_HIP = bodyparts.HipRight;
        bodyparts.LEFT_HIP = bodyparts.HipLeft;        
        bodyparts.LEFT_SHOULDER = bodyparts.ShoulderLeft;
        bodyparts.RIGHT_SHOULDER = bodyparts.ShoulderRight;        
        bodyparts.LEFT_ELBOW = bodyparts.WristRight;
        bodyparts.RIGHT_ELBOW = bodyparts.WristLeft;
        bodyparts.LEFT_KNEE =  bodyparts.KneeLeft;
        bodyparts.RIGHT_KNEE = bodyparts.KneeRight;        
        bodyparts.RIGHT_FOOT =  [bodyparts.AnkleRight,	 bodyparts.FootRight];
        bodyparts.LEFT_FOOT =  [bodyparts.AnkleLeft,	 bodyparts.FootLeft];
        bodyparts.HEAD	=	 bodyparts.Head;
        bodyparts.TORSO = bodyparts.SpineShoulder;
        bodyparts.RIGHT_HAND = [bodyparts.HandTipRight bodyparts.ThumbRight];
        bodyparts.LEFT_HAND = [bodyparts.HandTipLeft bodyparts.ThumbLeft];        
    case {120,60}
        bodyparts.hip_center = 1;
        bodyparts.spine = 2;
        bodyparts.shoulder_center = 3;
        bodyparts.head = 4;
        bodyparts.shoulder_left = 5;
        bodyparts.elbow_left = 6;
        bodyparts.wrist_left = 7;
        bodyparts.hand_left = 8;
        bodyparts.shoulder_right = 9;
        bodyparts.elbow_right = 10;
        bodyparts.wrist_right = 11;
        bodyparts.hand_right = 12;
        bodyparts.hip_left = 13;
        bodyparts.knee_left = 14;
        bodyparts.ankle_left = 15;
        bodyparts.foot_left = 16;
        bodyparts.hip_right = 17;
        bodyparts.knee_right = 18;
        bodyparts.ankle_right = 19;
        bodyparts.foot_right = 20;
        %%% synonyms        
        bodyparts.RIGHT_HIP = bodyparts.hip_right;
        bodyparts.LEFT_HIP = bodyparts.hip_left;        
        bodyparts.LEFT_SHOULDER = bodyparts.shoulder_left;
        bodyparts.RIGHT_SHOULDER = bodyparts.shoulder_right;        
        bodyparts.RIGHT_FOOT =  [bodyparts.ankle_right,	 bodyparts.foot_right];
        bodyparts.LEFT_FOOT =  [bodyparts.ankle_left,	 bodyparts.foot_left];
        bodyparts.HEAD	=	 bodyparts.head;
        bodyparts.TORSO = bodyparts.shoulder_center;
        bodyparts.LEFT_HAND = [bodyparts.wrist_left, bodyparts.hand_left];
        bodyparts.RIGHT_HAND = [bodyparts.wrist_right, bodyparts.hand_right]; 
    case {90,45}
        bodyparts.HEAD = 1;
        bodyparts.NECK = 2;
        bodyparts.TORSO = 3;
        bodyparts.LEFT_SHOULDER = 4;
        bodyparts.LEFT_ELBOW = 5;
        bodyparts.RIGHT_SHOULDER = 6;
        bodyparts.RIGHT_ELBOW = 7;
        bodyparts.LEFT_HIP = 8;
        bodyparts.LEFT_KNEE = 9;
        bodyparts.RIGHT_HIP = 10;
        bodyparts.RIGHT_KNEE = 11;
        bodyparts.LEFT_HAND = 12;
        bodyparts.RIGHT_HAND = 13;
        bodyparts.LEFT_FOOT = 14;
        bodyparts.RIGHT_FOOT = 15;
        %%%
        bodyparts.hip_center = [];    
    otherwise
        dbgmsg('No idea from this size from what type of skeleton this is. I will assume it is a randomstick.')
        return
end
end
function awk = generate_awk(data_val)
%% Awk definition:
important = 1;%0.1;
relevant = 1;%0.03;
minor = 1;%0.005;

awk = [...
    important;...   %1    hips
    important;...   %2    abdomen
    important;...   %3    neck or something
    relevant;...    %4    tip of the head
    important;...   %5    right shoulder
    relevant;...    %6    right also shoulder or elbow
    relevant;...    %7    right elbow maybe
    relevant;...    %8    right hand
    important;...   %9    left part of shoulder
    relevant;...    %10   left something maybe elbow
    relevant;...    %11   left maybe elbow
    relevant;...    %12   left hand
    important;...   %13   left hip
    relevant;...    %14   left knee
    minor;...       %15   left part of foot
    minor;...       %16   left tip of foot
    important;...   %17   right hip %important because hips dont lie
    relevant;...    %18   right knee
    minor;...       %19   right part of foot
    minor;...       %20   right tip of foot
    important;...   %21   middle of upper torax
    minor;...       %22   right some part of the hand
    minor;...       %23   right some other part of the hand
    minor;...       %24   left some part of the hand
    minor];         %25   left some other part of the hand

if size(awk,1)*6~=size(data_val,1)
    if floor(size(data_val,1)/2)==size(data_val,1)/2
        awk = ones(size(data_val,1)/6,1);
    else
        awk = ones(size(data_val,1)/3,1);
    end
    dbgmsg('Must update awk for a skeleton this size.',1)
end


end
function R = vecRotMat(f,t)
%% Purpose:
%Commonly, it is desired to have a rotation matrix which will rotate one
%unit vector, f,  into another unit vector, t. It is desired to
%find R(f,t) such that R(f,t)*f = t.
%
%This program, vecRotMat is the most
%efficient way to accomplish this task. It uses no square roots or
%trigonometric functions as they are very computationally expensive.
%It is derived from the work performed by Moller and Hughes, which have
%suggested that this method is the faster than any previous transformation
%matrix methods tested.
%
%
%% Inputs:
%f                      [N x 3]                         N number of vectors
%                                                       in which to
%                                                       transform into
%                                                       vector t.
%
%t                      [N x 3]                         N number of vectors
%                                                       in which it is
%                                                       desired to rotate
%                                                       f.
%
%
%% Outputs:
%R                      [3 x 3 x N]                     N number of
%                                                       rotation matrices
%
%% Source:
% Moller,T. Hughes, F. "Efficiently Building a Matrix to Rotate One
% Vector to Another", 1999. http://www.acm.org/jgt/papers/MollerHughes99
%
%% Created By:
% Darin C. Koblick (C) 07/17/2012
% Darin C. Koblick     04/22/2014       Updated when lines are close to
%                                       parallel by checking
%% ---------------------- Begin Code Sequence -----------------------------
%It is assumed that both inputs are in vector format N x 3
dim3 = 2;
%Declare function handles for multi-dim operations
normMD = @(x,y) sqrt(sum(x.^2,y));
anyMD  = @(x) any(x(:));
% Inputs Need to be in Unit Vector Format
if anyMD(single(normMD(f,dim3)) ~= single(1)) || anyMD(single(normMD(t,dim3)) ~= single(1))
    error('Input Vectors Must Be Unit Vectors');
end
%Pre-Allocate the 3-D transformation matrix
R = NaN(3,3,size(f,1));

v = permute(cross(f,t,dim3),[3 2 1]);
c = permute(dot(f,t,dim3),[3 2 1]);
h = (1-c)./dot(v,v,dim3);

idx  = abs(c) > 1-1e-13;
%If f and t are not parallel, use the following computation
if any(~idx)
    %For any vector u, the rotation matrix is found from:
    R(:,:,~idx) = ...
        [c(:,:,~idx) + h(:,:,~idx).*v(:,1,~idx).^2,h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)-v(:,3,~idx),h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)+v(:,2,~idx); ...
        h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)+v(:,3,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,2,~idx).^2,h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)-v(:,1,~idx); ...
        h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)-v(:,2,~idx),h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)+v(:,1,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,3,~idx).^2];
end
%If f and t are close to parallel, use the following computation
if any(idx)
    f = permute(f,[3 2 1]);
    t = permute(t,[3 2 1]);
    p = zeros(size(f));
    iidx = abs(f(:,1,:)) <= abs(f(:,2,:)) & abs(f(:,1,:)) < abs(f(:,3,:));
    if any(iidx & idx)
        p(:,1,iidx & idx) = 1;
    end
    iidx = abs(f(:,2,:)) < abs(f(:,1,:)) & abs(f(:,2,:)) <= abs(f(:,3,:));
    if any(iidx & idx)
        p(:,2,iidx & idx) = 1;
    end
    iidx = abs(f(:,3,:)) <= abs(f(:,1,:)) & abs(f(:,3,:)) < abs(f(:,2,:));
    if any(iidx & idx)
        p(:,3,iidx & idx) = 1;
    end
    u = p(:,:,idx)-f(:,:,idx);
    v = p(:,:,idx)-t(:,:,idx);
    rt1 = -2./dot(u,u,dim3);
    rt2 = -2./dot(v,v,dim3);
    rt3 = 4.*dot(u,v,dim3)./(dot(u,u,dim3).*dot(v,v,dim3));
    R11 = 1 + rt1.*u(:,1,:).*u(:,1,:)+rt2.*v(:,1,:).*v(:,1,:)+rt3.*v(:,1,:).*u(:,1,:);
    R12 = rt1.*u(:,1,:).*u(:,2,:)+rt2.*v(:,1,:).*v(:,2,:)+rt3.*v(:,1,:).*u(:,2,:);
    R13 = rt1.*u(:,1,:).*u(:,3,:)+rt2.*v(:,1,:).*v(:,3,:)+rt3.*v(:,1,:).*u(:,3,:);
    R21 = rt1.*u(:,2,:).*u(:,1,:)+rt2.*v(:,2,:).*v(:,1,:)+rt3.*v(:,2,:).*u(:,1,:);
    R22 = 1 + rt1.*u(:,2,:).*u(:,2,:)+rt2.*v(:,2,:).*v(:,2,:)+rt3.*v(:,2,:).*u(:,2,:);
    R23 = rt1.*u(:,2,:).*u(:,3,:)+rt2.*v(:,2,:).*v(:,3,:)+rt3.*v(:,2,:).*u(:,3,:);
    R31 = rt1.*u(:,3,:).*u(:,1,:)+rt2.*v(:,3,:).*v(:,1,:)+rt3.*v(:,3,:).*u(:,1,:);
    R32 = rt1.*u(:,3,:).*u(:,2,:)+rt2.*v(:,3,:).*v(:,2,:)+rt3.*v(:,3,:).*u(:,2,:);
    R33 = 1 + rt1.*u(:,3,:).*u(:,3,:)+rt2.*v(:,3,:).*v(:,3,:)+rt3.*v(:,3,:).*u(:,3,:);
    R(:,:,idx) = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
end
end