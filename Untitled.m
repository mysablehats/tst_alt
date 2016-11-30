clear all;
%close all;
env = aa_environment;
onlineclassdata = loadfileload('onlineclass',env);

classfdata = loadfileload('realclassifier',env);
%realclass = classfdata.realclass;
%realvideo(classfdata.realclass.outstruct, classfdata.realclass.allconn, classfdata.realclass.simvar,0)
%tic
onlineclassdata.allskel3 = generate_skel_online(onlineclassdata.chunk.chunk);
%classfdata.realclass.simvar.paramsZ.PLOTIT = 0;
if classfdata.realclass.simvar.paramsZ.PLOTIT
    figure
    skeldraw(onlineclassdata.allskel3.skel(:,:,1:3),'ts');
    figure
    skeldraw(makefatskel(onlineclassdata.outstruct.train.data(1:45,1:3)),'ts');
end
labellabel = online_classifier(classfdata.realclass.outstruct,onlineclassdata.allskel3, classfdata.realclass.allconn, classfdata.realclass.simvar);
%toc