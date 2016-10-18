clear
% Generate the graphical model
GeneratePMats;
N_V   = size(AMat,1); % No. of variables/nodes in the graph
BVec  = randn(N_V,1);

% Run the dynamics
lambda = 0.1; %low pass filtering for the dynamics
[~, ~, ~, ~, ~, ~, ~,Dynamics,NIter] = GaussTRP_Fn(AMat,BVec,lambda);

EdgeStruct  = Dynamics.EdgeStruct; % Get the edge data

N_E   = length(EdgeStruct); % No. of edges
dsVec = sum(abs(AMat)>0,2)-1; % No. of neighbors of each node

NodesMat    = [];

for ii = 1:N_V-1
    for jj = ii+1:N_V
        if AMat(ii,jj) ~= 0
            NodesMat    = [NodesMat; ii, jj];
        end
    end
end

% Edge participation. For each node, find out the edges it belongs to
NEdges = size(NodesMat,1);
EdgeSet = {};
for kk = 1:N_V
    temp1 = union(find(NodesMat(:,1) == kk),find(NodesMat(:,2) == kk));
    EdgeSet{kk} =  temp1(:);     
end





% First list the model paramters

% Invariances for this particular graph
JMat = AMat;

% Global parameters
alpVec  = [1 2]';
bet1Vec = [-1 1]';
bet2Vec = [1 2]';
gamVec  = [-1 -2]';
Nalp    = length(alpVec);
Nbet1   = length(bet1Vec);
Nbet2   = length(bet2Vec);
Ngam    = length(gamVec);

N_VParams = 2; % No. of node params
N_EParams = 4; % No. of edge params

% Type 1 parameters - setting true values for now 
del1_Vec = 1 - dsVec;
W1       = [1 0 0 0]';
c1       = 0;
d1       = 0;
G1       = zeros(Nalp,Nbet1,N_EParams);
G1(2,1,2) = -1;

% Type 2 parameters - setting true values for now 
del2_Vec = 1 - dsVec;
W2       = [0 0 1 0]';
c2       = 0;
d2       = 0;
G2       = zeros(Nalp,Nbet2,Ngam,N_EParams,N_EParams);
G2(1,1,1,4,2) = -1;

% Type 3 parameters - setting true values for now 
W3       = [0 0 0 0 1 0 0 0]';
c3       = 0;
d3       = 0;
G3       = zeros(Nalp,Nbet1,N_EParams);
G3(2,1,2) = 1;

% Type 4 parameters - setting true values for now 
W4       = [0 0 0 0 0 1 0 0]';
c4       = 0;
d4       = 0;
G4       = zeros(Nalp,Nbet2,Ngam,N_EParams,N_EParams);
G4(1,1,1,4,2) = 1;




% Generate dynamics using true model params 


A_model     = zeros(N_V,NIter);
B_model     = zeros(N_V,NIter);
Pss_model   = zeros(N_E,NIter);
Ptt_model   = zeros(N_E,NIter);
Vs_model    = zeros(N_E,NIter);
Vt_model    = zeros(N_E,NIter);

% Initialize the dynamics
A_model(:,1)  = Dynamics.AMat(:,1);
B_model(:,1)  = Dynamics.BMat(:,1);


for iter = 1:NIter-1 
    
    for kk = 1:N_V
        ReqEdges     = EdgeSet{kk};  % Node kk belongs to edges in ReqEdges
        NNeighs      = dsVec(kk); % No. of neighbors of each nodes
        
        temp1 = 0;
        temp2 = 0;
        
        for nn = 1:NNeighs
            EdgeId    = ReqEdges(nn);
            ReqNodes  = EdgeStruct{EdgeId}.Nodes; 
            idx_s     = find(ReqNodes == kk);
            
            if idx_s == 1
            theta_vec = [Dynamics.Pss(EdgeId,iter); Dynamics.Ptt(EdgeId,iter); ...
                Dynamics.Vs(EdgeId,iter); Dynamics.Vt(EdgeId,iter)];
            else 
            theta_vec = [Dynamics.Ptt(EdgeId,iter); Dynamics.Pss(EdgeId,iter); ...
                Dynamics.Vt(EdgeId,iter); Dynamics.Vs(EdgeId,iter)];
            end
            
            Jsu  = JMat(ReqNodes(1),ReqNodes(2));
            
            % Type 1
            Jtens1 = repmat(Jsu.^alpVec,1,Nbet1,N_EParams);
            
            thetamat1 = zeros(Nbet1,N_EParams);
            for bb = 1:Nbet1
                thetamat1(bb,:) = (theta_vec.^bet1Vec(bb))';
            end
            
            the1 = G1*0;
            for aa = 1:Nalp
                the1(aa,:,:) = thetamat1;
            end
            
            tens_prod1 = G1.*Jtens1.*the1;
            temp1 = temp1 + W1'*theta_vec + c1*Jsu + d1*Jsu^2 + sum(tens_prod1(:)); 
            
            % Type 2
            Jtens2 = repmat(Jsu.^alpVec,1,Nbet2,Ngam,N_EParams,N_EParams);
            theta_btens = G2*0;
            theta_gtens = G2*0;

            for bb = 1:Nbet2
                for ll = 1:N_EParams
                    theta_btens(:,bb,:,ll,:) = theta_vec(ll)^bet2Vec(bb);
                end
            end
            
            for gg = 1:Ngam
                for mm = 1:N_EParams
                    theta_gtens(:,:,gg,:,mm) = theta_vec(mm)^gamVec(gg);
                end
            end
            
            tens_prod2 = G2.*Jtens2.*theta_btens.*theta_gtens;
            temp2 = temp2 + W2'*theta_vec + c2*Jsu + d2*Jsu^2 + sum(tens_prod2(:)); 
            
            
            
        end
        
        A_model(kk,iter+1) = (1-lambda)*A_model(kk,iter) + ...
            lambda*(temp1 + del1_Vec(kk)*A_model(kk,iter));
        
        B_model(kk,iter+1) = (1-lambda)*B_model(kk,iter) + ...
            lambda*(temp2 + del2_Vec(kk)*B_model(kk,iter));
        
    end
    
    
    
end






























