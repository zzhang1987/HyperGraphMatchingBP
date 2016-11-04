function [tend, id, res] = runFGM(Fname)
load(Fname);
global footpath;
%delete(gcp)
%parpool
footpath = '/Users/zhen/HungarianBP/';
footpath = strcat(footpath,'/fgm');
rng(45678);



rmpath(footpath);
rmpath(genpath([footpath '/src']));
rmpath(genpath([footpath '/lib']));


addpath(footpath);
addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));
prSet(1);

NofNodes = prod(size(GT));

NT1 = size(Triplets,1);
NT2 = size(NTriplets,1);

I = repmat(Triplets(:,1), 1, NT2);
J = repmat(Triplets(:,2), 1, NT2);
K = repmat(Triplets(:,3), 1, NT2);

I1 = repmat(NTriplets(:,1), 1, NT1)';
J1 = repmat(NTriplets(:,2), 1, NT1)';
K1 = repmat(NTriplets(:,3), 1, NT1)';

idx1 = I * NofNodes + I1;
idx2 = J * NofNodes + J1;
idx3 = K * NofNodes + K1;

indH3=[idx1(:), idx2(:), idx3(:)];
valH3=Similarity(:);




% added by Quynh Nguyen
%remove duplicated tuples
indH3 = sort(indH3, 2);
[indH3 id1 id2] = unique(indH3, 'rows');
valH3 = valH3(id1);
%remove duplicated indices: (i,i,i), (i,i,j), (i,j,i), (i,j,j), etc
t1 = indH3(:, 1) - indH3(:, 2);
t2 = indH3(:, 1) - indH3(:, 3);
t3 = indH3(:, 2) - indH3(:, 3);
t4 = (t1 == 0) + (t2 == 0) + (t3 == 0);
indH3 = indH3(t4 == 0, :);
valH3 = valH3(t4 == 0);
% upperbound the number of nonzeros
maxL = min(5*10^6, length(valH3));
[v id] = sort(valH3, 'descend');
% id = randperm(length(valH3));
id = id(1:maxL);
valH3 = valH3(id);
indH3 = indH3(id, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Super-symmetrize the 3rd order tensor: if (i,j,k) with i#j#k is a nonzero
% entry of the tensor, then all of its six permutations should also be the entries
% NOTE that this step is important for the current implementation of our algorithms !!!
% Thus, if a tuple (i,j,k) for i#j#k is stored in 'indH3' below, then all of its
% permutations should also be stored in 'indH3' and 'valH3'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps = perms([1 2 3]);
Nt3 = length(valH3);
valH3 = [valH3; valH3; valH3; valH3; valH3; valH3];
old_indH3 = indH3;
indH3 = [];
for i = 1:size(ps, 1)
    indH3((i-1)*Nt3+1:i*Nt3, :) = old_indH3(:, ps(i, :));
end

% sort tuples in ascending order
dim = NofNodes * NofNodes;
uid = indH3(:, 1)*dim*dim + indH3(:, 2)*dim + indH3(:, 3);
[v id] = sort(uid);
valH3 = valH3(id);
indH3 = indH3(id, :);
if exist('KP', 'var')
    idx1 = repmat((1:NofNodes)', [1 NofNodes]) - 1;
    idx2 = repmat((1:NofNodes), [NofNodes 1]) - 1;
    indH1 = idx1 * NofNodes + idx2;
    indH1 = int32(indH1(:));
    valH1 = KP(:);
end

if exist('KQ', 'var')
    E1 = Edges;
    E2 = [NEdges; NEdges(:,2), NEdges(:,1)];
    NEdgesSize = size(NEdges,1)
    KQ = [KQ; KQ(:, (NEdgesSize+1):2*NEdgesSize), KQ(:, 1:NEdgesSize)] / 2;
    
    E1 = [Edges; Edges(:,2), Edges(:,1)];
    
    [G1,H1] = gphEg2IncA(E1' + 1, NofNodes);
    [G2,H2] = gphEg2IncA(E2' + 1, NofNodes);
    gph1.G=G1;
    gph1.H=H1;
    gph2.G=G2;
    gph2.H=H2;
end
gphs{1}=gph1;
gphs{2}=gph2;
[pars, algs] = gmPar(2);
tstart = tic;
asgFgmD = fgmD(KP, KQ, ones(NofNodes,NofNodes), gphs, [], pars{9}{:});
tend = toc(tstart);

[~, id] = max(asgFgmD.X');
id = id - 1;
res = asgFgmD.obj;

end



function [G, H, U, V] = gphEg2IncA(Eg, n)
% Obtain incidence matrix for asymmetric edges.
%
% Input
%   Eg      -  graph edge, 2 x (2m)
%   n       -  #nodes
%
% Output
%   G       -  node-edge incidence (starting), n x (2m)
%   H       -  node-edge incidence (ending), n x (2m)
%   U       -  incidence component, n x m
%   V       -  incidence component, n x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-11-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-01-2012

% dimension
m = round(size(Eg, 2) / 2);

% incidence matrix
[U, V] = zeross(n, m);
for c = 1 : m
    U(Eg(1, c), c) = 1;
    V(Eg(2, c), c) = 1;
end

% incidence
G = [U, V];
H = [V, U];

end