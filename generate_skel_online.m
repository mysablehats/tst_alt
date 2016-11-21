function skel = generate_skel_online(chunk)
chunk15 = translate_from(chunk);
skel.skel = chunk15(:,:,1:9);
skel.vel = diff(chunk15,1,3);
skel.act_type = '';
end