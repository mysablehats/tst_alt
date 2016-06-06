rotateskels = zeros(25,3,100);
domdom = linspace(0,2*pi);
detectedrot = zeros(1,100);
for i = 1:100
    rotateskels(:,:,i) = rotskel(protoskel,0,domdom(i),0);
    detectedrot(i) = detectrotation(rotateskels(:,:,i), skelldef, 'hips');
end
skeldraw(rotateskels);
figure 
plot(domdom,detectedrot)
figure
hold on
for i = 1:100
    skeldraw(rotskel(rotateskels(:,:,i),0,detectedrot(i),0),'f');
end
hold off
%load('protoskel.mat')
%[~, skelldef] = conformskel(makethinskel(protoskel),'norotateshoulders');
%a= detectrotation(protoskel2, skelldef, 'hips');skeldraw(rotskel(protoskel2,0,a,0),'f');