function data = online_classifier(data, allskel, arq_connect, simvar)
global VERBOSE LOGIT
VERBOSE = 1; % true;
LOGIT = false;
%a = struct();

%runs conformactions

[allskel,~] = conformactions(allskel, [], simvar.prefilter); %%% should be split in 2

%extractdata
[data.test, ~] = extractdata(allskel, simvar.activity_type, simvar.labels_names,simvar.extract{:}); %%% I dont know if I need this, or if it makes sense to use this function
%conformskel
[data.test, ~ ] = conformskel(data.test, simvar.preconditions{:}, 'test', simvar.paramsZ.skelldef);
%loop from l.65 of starter_sc -> make it into a function???
for j = 1:length(arq_connect)
    [~, data.test] = gas_method(data, data.test, 'vot', arq_connect(j),j,size(data.train.data,1));
end


for i = 1:length(arq_connect)
    %labeller
    %for i gases... like above
    [data.test.gas(i).class, ~] = find(data.gas(i).nodesl(:, data.test.gas(i).bestmatchbyindex));
end
if isfield(simvar,'labels_names')
    %size(simvar.labels_names)
else
    disp('no label names defined!')
end
try
    %size(arq_connect)
    for i = 1:length(arq_connect)
        a = []; for j = 1:length(data.test.gas(i).class), a = ([ a simvar.labels_names{data.test.gas(i).class(j)} '-']);end, disp(a), clear a
    end
    %pause(1)
catch
    disp('something went wrong')
end
%disp('')
if simvar.paramsZ.PLOTIT
    figure
    plot_act(data, data.test.data, simvar.paramsZ.skelldef) %% data.test!!!!
end

end