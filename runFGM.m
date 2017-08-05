function [tend, id, res] = runFGM(Fname)
load(Fname);
global footpath;
%delete(gcp)
%parpool
homepath = getenv('HOME');
footpath = strcat(homepath, '/Inference/code/');
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
