function  [tend, id, res] = runSM(Fname)
load(Fname);

NofNodes = prod(size(GT));

M = [];
% if consider Edges
if exist('KQ','var')
        % Added by Lee 2016-10-31
        % Failed because KQ didn't containe ALL the edges
%     Mat = sparse(NofNodes*NofNodes,NofNodes*NofNodes);
%     dialog_KP = KP.*eye(NofNodes);
%     Mat = Mat + kron(dialog_KP,eye(NofNodes));
%     E1 = Edges+1;
%     E2 = NEdges;
%     for i=1:size(E1,1)
%         for j=1:size(E2,1)
%             Mat(E2(j,1)*NofNodes+E1(i,1),E2(j,2)*NofNodes+E1(i,2)) = KQ(i,j);
%         end
%     end
%     Mat = Mat + Mat';
% 
% 
% E12 = ones(NofNodes,NofNodes);
% [L12(:,1),L12(:,2)] = find(E12);
% [group1,group2] = make_group12(L12);

    %-------------------------
    if size(P1, 1) == 2 
        P1 = P1';
        P2 = P2';
    else 
        P1 = P1;
        P2 = P2;
    end
    nP1 = size(P1,1);
    nP2 = size(P2,1);

    E12 = ones(nP1,nP2);
    n12 = nnz(E12);
    [L12(:,1) L12(:,2)] = find(E12);
    %[group1 group2] = make_group12(L12);

    E1 = ones(nP1); E2 = ones(nP2);
    [L1(:,1) L1(:,2)] = find(E1);
    [L2(:,1) L2(:,2)] = find(E2);
    G1 = P1(L1(:,1),:)-P1(L1(:,2),:);
    G2 = P2(L2(:,1),:)-P2(L2(:,2),:);

    if 0
        G1_x = reshape(G1(:,1), [nP1 nP1]);
        G1_y = reshape(G1(:,2), [nP1 nP1]);
        G2_x = reshape(G2(:,1), [nP2 nP2]);
        G2_y = reshape(G2(:,2), [nP2 nP2]);
        M = (repmat(G1_x, nP2, nP2)-kron(G2_x,ones(nP1))).^2 + (repmat(G1_y, nP2, nP2)-kron(G2_y,ones(nP1))).^2;
    else
        G1 = sqrt(G1(:,1).^2+G1(:,2).^2);
        G2 = sqrt(G2(:,1).^2+G2(:,2).^2);
        G1 = reshape(G1, [nP1 nP1]);
        G2 = reshape(G2, [nP2 nP2]);
        M = (repmat(G1, nP2, nP2)-kron(G2,ones(nP1))).^2;
    end
    M = exp(-M./0.5);
    %------------------------------
end

tstart = tic;
[X] = SM(M);
tend = toc(tstart);
X = reshape(X,[nP2 nP1]);
X = asgHun(X');

res = [];
[~, id] = max(X);
id = int32(id - 1);
