function [tend, id, res] = runRRWHM(MFname)
load(MFname);

NofNodes = prod(size(GT));
% if consider Triplets
if exist('Triplets','var')
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
end

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
    I1 = repmat(E1(:,1), [1, size(E2,1)]);
    J1 = repmat(E1(:,2), [1, size(E2,1)]);
    I2 = repmat(E2(:,1)', [size(E1,1), 1]);
    J2 = repmat(E2(:,2)', [size(E1,1), 1]);
    
    idx1 = I1 * NofNodes + I2;
    idx2 = J1 * NofNodes + J2;
    
    indH2 = [idx1(:), idx2(:)];
    valH2 = KQ(:)/2;
    
    indH2 = [indH2; indH2(:,2), indH2(:,1)];
    valH2 = [valH2; valH2];
    
end

if ~exist('indH1','var')
    indH1 = [];
    valH1 = [];
end

if ~exist('indH2', 'var')
    indH2 = [];
    valH2 = [];
end

if isempty(indH1) && isempty(valH1)
  indH1=int32(zeros(0,1));
  valH1=zeros(0,1);
end
if isempty(indH2) && isempty(valH2)
  indH2=int32(zeros(0,2));
  valH2=zeros(0,1);
end
if isempty(indH3) && isempty(valH3)
  indH3=int32(zeros(0,3));
  valH3=zeros(0,1);
end


iterMax = 300;
c = 0.2;
X = ones(NofNodes,NofNodes);
X = X ./ sum(sum(X));
tstart = tic;
[X, nIter] = mexRRWHM(X,int32(indH1),valH1,int32(indH2),valH2,int32(indH3),valH3,iterMax,c);
tend = toc(tstart);
X = asgHun(X);

res = getMatchScore(int32(indH3), valH3, X);
[~, id] = max(X);
id = int32(id - 1);
