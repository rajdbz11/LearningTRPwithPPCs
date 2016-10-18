% Generate valid Precision matrices of desired structure
% Valid precision matrix => positive definite. Check for positive
% eigenvalues.

%  --------------------------- Chains -----------------------------------
L = 3+ ceil(7*rand(1));
SMat = eye(L);
SMat = SMat + [zeros(L,1),SMat(:,1:L-1)] + [zeros(1,L);SMat(1:L-1,:)];
TempMat = sprandsym(SMat);
% % -------------------------------------------------------------------------

% %  --------------------------- Loops -----------------------------------
% L = 3+ ceil(22*rand(1));
% SMat = eye(L);
% SMat = SMat + [zeros(L,1),SMat(:,1:L-1)] + [zeros(1,L);SMat(1:L-1,:)];
% SMat(1,L) = 1; SMat(L,1) = 1;
% TempMat = sprandsym(SMat);
% -------------------------------------------------------------------------



% %  --------------------------- Grids ------------------------------------
% Specify the size of the grid. Total no. of nodes in the graph is MxN
% M = 2 + ceil(8*rand(1)); % No. of rows
% N = M; % No. of columns
% 
% % Generate indices for the nodes in this graph
% IndexMat = reshape(cumsum(ones(M*N,1)),M,N);
% 
% % Now to generate the structure of the Precision Matrix
% % Start by populating the diagonal elements
% SMat = eye(M*N); 
% % Now populate the precision matrix based on the grid connectivity: 
% % up-down and left-right neigbors only.
% % Also include extra checks for the nodes at the edges of the grid.
% for ii = 1:M
%     for jj = 1:N
%         % first look for the left and right entries
%         switch jj
%             case 1
%                 SMat(IndexMat(ii,jj),IndexMat(ii,jj+1)) = 1;
%             case N
%                 SMat(IndexMat(ii,jj),IndexMat(ii,jj-1)) = 1;
%             otherwise
%                 SMat(IndexMat(ii,jj),IndexMat(ii,jj+1)) = 1;
%                 SMat(IndexMat(ii,jj),IndexMat(ii,jj-1)) = 1;
%         end
%         % Next look for the up and down entries
%         switch ii
%             case 1
%                 SMat(IndexMat(ii+1,jj),IndexMat(ii,jj)) = 1;
%             case M
%                 SMat(IndexMat(ii-1,jj),IndexMat(ii,jj)) = 1;
%             otherwise
%                 SMat(IndexMat(ii+1,jj),IndexMat(ii,jj)) = 1;
%                 SMat(IndexMat(ii-1,jj),IndexMat(ii,jj)) = 1;
%         end
%     end
% end
% L = M*N; % total no. of variables
% TempMat = sprandsym(SMat);
% % -------------------------------------------------------------------------


% %  --------------------------- Dense ------------------------------------
% L = 3+ ceil(22*rand(1));
% TempMat = sprandsym(L,0.9);
% % -------------------------------------------------------------------------



AMat    = zeros(L,L);
for ii = 1:L
    for jj = 1:L
        AMat(ii,jj) = TempMat(ii,jj);
    end
end
clear TempMat;
EVals   = eig(AMat);
c       = 2*rand(1,1); % arbitrary positive constant
AMat  = AMat + eye(L)*(c - min(EVals));

CovMat      = inv(AMat);
StdVec      = sqrt(diag(CovMat));
CheckMat    = CovMat./(StdVec*StdVec');



