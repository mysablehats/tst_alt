function acquirekinect(ax)
env = aa_environment;
classfdata = loadfileload('realclassifier',env); 
realclass = classfdata.realclass;
if isfield(realclass, 'gases')
    realvideo(ax, realclass.gases, realclass.allconn, realclass.simvar,1)
else
    warning('Old file??')
    realvideo(ax, realclass.outstruct.gas, realclass.allconn, realclass.simvar,0) %%% should work, right?
end

%realvideo(classfdata.outstruct, baq(classfdata.pallconn), classfdata.simvar,0); 