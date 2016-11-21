clear all;
close all;
load('..\\onlineclass.mat'); 
allskel3 = generate_skel_online(chunk.chunk);
skeldraw(allskel3.skel(:,:,1:3),'ts');
figure
skeldraw(makefatskel(outstruct.train.data(1:45,1:3)),'ts');
labellabel = online_classifier(outstruct,allskel3, arc_conn, simvar);