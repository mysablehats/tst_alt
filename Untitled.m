clear all;
%close all;
env = aa_environment;
onlineclassdata = loadfileload('onlineclass',env); 
onlineclassdata.allskel3 = generate_skel_online(onlineclassdata.chunk.chunk);
figure
skeldraw(onlineclassdata.allskel3.skel(:,:,1:3),'ts');
figure
skeldraw(makefatskel(onlineclassdata.outstruct.train.data(1:45,1:3)),'ts');
labellabel = online_classifier(onlineclassdata.outstruct,onlineclassdata.allskel3, onlineclassdata.arc_conn, onlineclassdata.simvar);