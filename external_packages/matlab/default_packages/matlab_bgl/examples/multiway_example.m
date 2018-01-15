load ../graphs/minnesota.mat

n = size(A,1);
k1 = 75;
k2 = 50;
%vs = ceil(size(A,1)*rand(1,k));
start1 = 1;
start2 = 800;
vs1 = start1:start1+k1;
vs2 = start2:start2+k2;
vs = [vs1 vs2];
C = approx_multiway_cut(A,vs);

gplot(triu(A),xy,':');
hold on;
gplot(triu(C),xy,'r-');
plot(xy(vs,1),xy(vs,2),'.');
hold off;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'Box','off');