function data = online_classifier(gases, allskel, arq_connect, simvar)
global VERBOSE LOGIT
VERBOSE = 0; % true;
LOGIT = false;
[data.test,~,~] = all3(allskel, simvar);
whatIlabel = 1:length(gases);
[data.gas, data.test] = sameclassfunc(gases, data.test, 'test', whatIlabel, arq_connect);
numskells = 4;
showoutcomes(simvar,data.test, numskells)
if simvar.paramsZ.PLOTIT
    figure
    plot_act(data.gas, data.test, 5, simvar.paramsZ.skelldef, numskells) 
end
end