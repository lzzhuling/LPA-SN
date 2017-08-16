%% Linearized Proximal Algorithm with semismooth Newton iteration (LPA-SN)
%%
%% Main Function
dim=2;
num=100;
R=0.3;
P0=-0.5+rand(2,10);
PP=-0.5+rand(2,100);
DD=randistance(P0,PP,R);
[X0,rmsd] = LPAsolverSN(P0,PP,DD,dim,num,R);
scatter(1:30,rmsd);
set(gca,'yscale','log');