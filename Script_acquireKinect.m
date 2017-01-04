clear all;
close all;
env = aa_environment;
classfdata = loadfileload('realclassifier',env); 
realclass = classfdata.realclass;
if isfield(realclass, 'gases')
    realvideo(realclass.gases, realclass.allconn, realclass.simvar,1)
else
    warning('Old file??')
    realvideo(realclass.outstruct.gas, realclass.allconn, realclass.simvar,0) %%% should work, right?
end

%realvideo(classfdata.outstruct, baq(classfdata.pallconn), classfdata.simvar,0); 