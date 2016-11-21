function [ss, metrics, gases] = starter_sc(ss, allconn)
gases = [];
%% starter_sc
% This is the main function to run the chained classifier, label and
% generate confusion matrices and recall, precision and F1 values for the
% skeleton classifier of activities using an architecture of chained neural
% gases on skeleton activities data (the STS V2 Dataset). This work is an
% attempt to implement Parisi, 2015's paper.

%%
global VERBOSE LOGIT
VERBOSE = true;
LOGIT = true;
dbgmsg('ENTERING MAIN LOOP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% making metrics structure

metrics = struct('confusions',[],'conffig',[],'outparams',[]);

arq_connect = baq(allconn);

%ss = makess(length(arq_connect));


%%% end of gas structures region


%% Gas-chain Classifier
% This part executes the chain of interlinked gases. Each iteration is one
% gas, and currently it works as follows:
% 1. Function setinput() chooses based on the input defined in allconn
% 2. Run either a Growing When Required (gwr) or Growing Neural Gas (GNG)
% on this data
% 3. Generate matrix of best matching units to be used by the next gas
% architecture

dbgmsg('Starting chain structure for GWR and GNG for nodes:',num2str(labindex),0)
dbgmsg('###Using multilayer GWR and GNG ###',0)

for j = 1:length(arq_connect)
    [ss.gas(j), ss.train] = gas_method(ss, ss.train,'train', arq_connect(j),j, size(ss.train.data,1)); % I had to separate it to debug it.
    metrics(j).outparams = ss.gas(j).outparams;
    if size(ss.val.data,1) > 0
<<<<<<< HEAD
        [ss, ss.val ]= gas_method(ss, ss.val,'val', arq_connect(j),j, size(ss.train.data,1));
=======
        [~, ss.val ]= gas_method(ss, ss.val,'val', arq_connect(j),j, size(ss.train.data,1));
>>>>>>> refs/remotes/origin/alice_variant
    end
end


%% Gas Outcomes
% This should be made into a function... It is a nice graph to perhaps
% debug the gas...
if arq_connect(j).params.PLOTIT
    figure
    for j = 1:length(arq_connect)        
        subplot (1,length(arq_connect),j)
        hist(ss.gas(j).outparams.graph.errorvect)
        title((ss.gas(j).name))
    end
end
%% Labelling
% The current labelling procedure for both the validation and training datasets. As of this moment I label all the gases
% to see how adding each part increases the overall performance of the
% structure, but since this is slow, the variable whatIlabel can be changed
% to contain only the last gas.
%
% The labelling procedure is simple. It basically picks the label of the
% closest point and assigns to that. In a sense the gas can be seen as a
% dimensional (as opposed to temporal) filter, encoding prototypical
% action-lets.
dbgmsg('Labelling',num2str(labindex),0)

whatIlabel = 1:length(ss.gas); %change this series for only the last value to label only the last gas

%%
% Specific part on what I want to label
for j = whatIlabel
    dbgmsg('Applying labels for gas: ''',ss.gas(j).name,''' (', num2str(j),') for process:',num2str(labindex),0)
    %%% hack to fix non existing validation set in online mode
    if ~isfield(ss.val,'gas')||j>length(ss.val.gas)
        [ss.train.gas(j).class, ~, ss.gas(j).nodesl ] = labeller(ss.gas(j).nodes, ss.train.gas(j).bestmatchbyindex,  ss.train.gas(j).bestmatchbyindex(1,1), ss.train.gas(j).inputs.input, ss.train.gas(j).y);
    else
        %standard line of code before hack.
        [ss.train.gas(j).class, ss.val.gas(j).class, ss.gas(j).nodesl ] = labeller(ss.gas(j).nodes, ss.train.gas(j).bestmatchbyindex,  ss.val.gas(j).bestmatchbyindex, ss.train.gas(j).inputs.input, ss.train.gas(j).y);
    end
    %%%% I dont understand what I did, so I will code this again.
    %%% why did I write .bestmatch when it should be nodes??? what was I thinnking [ss.train.gas(j).class, ss.val.gas(j).class] = untitled6(ss.gas(j).bestmatch, ss.train.gas(j).inputs.input,ss.val.gas(j).inputs.input, y_train);
    
end

%% Displaying multiple confusion matrices for GWR and GNG for nodes
% This part creates the matrices that can later be shown with the
% plotconfusion() function.

ss.figset = {}; %% you should clear the set first if you want to rebuild them
dbgmsg('Displaying multiple confusion matrices for GWR and GNG for nodes:',num2str(labindex),0)

for j = whatIlabel
    [~,metrics(j).confusions.train,~,~] = confusion(ss.train.gas(j).y,ss.train.gas(j).class);
    ss.gas(j).fig.train = {ss.train.gas(j).y,   ss.train.gas(j).class,strcat(ss.gas(j).name,ss.gas(j).method,'T')};
    
    if isfield(ss.val, 'data')&&~isempty(ss.val.data)%&&~isempty(ss.val.data)
        [~,metrics(j).confusions.val,~,~] = confusion(ss.val.gas(j).y,ss.val.gas(j).class);      
        
        dbgmsg(ss.gas(j).name,' Confusion matrix on this validation set:',writedownmatrix(metrics(j).confusions.val),0)
        ss.gas(j).fig.val =   {ss.val.gas(j).y,     ss.val.gas(j).class,  strcat(ss.gas(j).name,ss.gas(j).method,'V')};
        %ss.figset = [ss.figset, ss.gas(j).fig.val, ss.gas(j).fig.train];
        %%%
    else
        metrics(j).confusions.val = [];
        ss.gas(j).fig.val = [];
    end
    metrics(j).conffig = ss.gas(j).fig;
end

%% Actual display of the confusion matrices:
metitems = [];
for j = whatIlabel
    if arq_connect(j).params.PLOTIT
        metitems = [metitems j*arq_connect(j).params.PLOTIT];
    end
end
if ~isempty(metitems)
    figure
    plotconf(metrics(metitems))
end
%plotconf(ss.figset{:})
if isfield(ss.gas(end).fig, 'val')&&~isempty(ss.gas(end).fig.val)
    figure
    plotconfusion(ss.gas(end).fig.val{:})
end
end
function stringmatrix = writedownmatrix(matrix)
stringmatrix = '';
for i =1:size(matrix,1)
    for j = 1:size(matrix,2)
        stringmatrix = strcat(stringmatrix,'(', num2str(i),',',num2str(j),'):', num2str(matrix(i,j)),'\t');
    end
end
end