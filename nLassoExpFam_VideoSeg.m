clear all
close all
profile off 
profile on
%%%% Generate two separate ER graphs of size N/2 with edge probability p
%%%% connect nodes between two cluter ERs using edge probability q
%restoredefaultpath
%rehash toolboxcache

[pathtothismfile,name,ext] = fileparts(mfilename('fullpath')) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph Initialisation SBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-edge probability of a cluster:
p = .1;
% Inter-edge probability between clusters
%q = 0.99;



% Number of nodes per cluster
N1 = 40;
N2 = 40;
% total number of nodes 
N=N1+N2;
RUNS = 10; 
RUNS = 1 ; 

% numer of features at each data point 



%%% read in image and construct empirical graph
patch_x = 0; 
patch_y = 0; 
patchsize = 3*(2*patch_x+1)*(2*patch_y+1); 

featurelen=2; 
d=featurelen;
CS_matrix = randn(d,patchsize); %%% random projection 


%%% instead of RP use average color channels of patch
featurelen=3; 
d=featurelen;
CS_matrix = [[ones(1,(2*patch_x+1)*(2*patch_y+1)),zeros(1,2*(2*patch_x+1)*(2*patch_y+1))]; 
             [zeros(1,(2*patch_x+1)*(2*patch_y+1)),ones(1,(2*patch_x+1)*(2*patch_y+1)),zeros(1,(2*patch_x+1)*(2*patch_y+1))];
             [zeros(1,2*(2*patch_x+1)*(2*patch_y+1)),ones(1,(2*patch_x+1)*(2*patch_y+1))]]; 
             

         
         
pic = (imread('AlexTest.png')) ; %RG_in; 
pic_dim = size(pic); 
pic  = pic(1:5:pic_dim(1), 1:5:pic_dim(2),:) ; 
pic_dim = size(pic); 
orange_level = zeros(pic_dim(1),pic_dim(2)); 
X =[]; 
tmp = zeros(3,1); 
for iter_row=1:pic_dim(1) 
  for iter_col=1:pic_dim(2) 
      tmp(:) = pic(iter_row,iter_col,:); 
     %  X = [X,tmp]; 
      orange_level(iter_row,iter_col) = sum(tmp.*[1;0;0.0])/norm(tmp); 
  end
end
out_img = zeros(pic_dim(1),pic_dim(2)); 

maxval = max(max(orange_level)) ; 

for iter_row=1:pic_dim(1) 
  for iter_col=1:pic_dim(2) 
     % tmp(:) = pic(iter_row,iter_col,:); 
      val = orange_level(iter_row,iter_col); 
      
      if (val > 0.9*maxval)  
       out_img(iter_row,iter_col) = val; 
      end
  end
end




raw_dim= size(pic); 

figure(1)
imshow(pic)

dim_len = size(pic); 
y_size = dim_len(1); 
x_size = dim_len(2); 

number_pixels = dim_len(1)*dim_len(2); 
N = number_pixels ; 

%feature_matrix = reshape (RGB,dim_len(1)*dim_len(2),dim_len(3)) ; 

feature_matrix = zeros(number_pixels,featurelen); 


for iter_pixel=1:number_pixels
    [row,col] = ind2sub(dim_len(1:2),iter_pixel) ; 
    tmp(:) = pic(row,col,:); 
    feature_matrix(iter_pixel,:) = tmp; 
end



% construct weighted adj. matrix 

W = sparse(dim_len(1),dim_len(2)) ; 

[rows, cols]= ind2sub(dim_len(1:2),1:number_pixels); 

row_col_mtx = [rows',cols'];   % the ith row contains row and col idx of pixel i 
max_degree = 10 ; 
Neighbours = sparse(N,max_degree); % ith row stores indcies of neighbors of pixel i 
degree_list = zeros(N,1); %ith row stores degree of pixel i 
edge_cntr=0; 
edge_idx = sparse(number_pixels,number_pixels); 
nodesinedge = sparse(number_pixels*3,2); 

nr_rows = pic_dim(1) ; 
nr_cols = pic_dim(2); 
iter_cntrs=1 ; 
for iter_row=1:nr_rows 
    for iter_col=1:nr_cols 
        node_idx=((iter_col-1)*nr_rows)+(iter_row-1)+1; 
        
       neighb1=((iter_col-1)*nr_rows)+rem(iter_row-1+1,nr_rows)+1; 
       nodesinedge(iter_cntrs,1) = node_idx ; 
       nodesinedge(iter_cntrs,2) = neighb1 ; 
       iter_cntrs=iter_cntrs+1; 
 
        neighb1=(rem(iter_col-1+1,nr_cols)*nr_rows)+rem(iter_row-1+1,nr_rows)+1; 
        nodesinedge(iter_cntrs,1) = node_idx ; 
       nodesinedge(iter_cntrs,2) = neighb1 ; 
       iter_cntrs=iter_cntrs+1; 
        
        neighb1=(rem(iter_col-1+1,nr_cols)*nr_rows)+rem(iter_row-1,nr_rows)+1; 
        nodesinedge(iter_cntrs,1) = node_idx ; 
       nodesinedge(iter_cntrs,2) = neighb1 ; 
       iter_cntrs=iter_cntrs+1; 
        
        
    end
end

        
% 
% 
% for i=1:number_pixels
%     for j=1:number_pixels
%         if j ~= i 
%            row_i = row_col_mtx(i,1) ; 
%            col_i = row_col_mtx(i,2); 
%            row_j = row_col_mtx(j,1); 
%            col_j = row_col_mtx(j,2); 
%            
%            
%            if abs(row_i - row_j)< 2 
%                if abs(col_i-col_j)<2 
%                    W(i,j) =1 ; 
%                    Neighbours(i,degree_list(i)+1)= j ;   % undirected neighours of pixel i 
%                    degree_list(i) = degree_list(i)+1 ; 
%                     if (i < j) 
%                         edge_cntr=edge_cntr+1; 
%                         edge_idx(i,j) = edge_cntr ; 
%                         nodesinedge = [nodesinedge;[i,j]] ; 
%                    end
%                end
%                
%            end
%         end
%         
%      end
%    
% end


%%%%% 
% generate feature vectors for each data point 
%%%%%
feature_matrix = (normalize(feature_matrix)) ; 

X_mtx = feature_matrix'; 

cols_samples =[];   % collects the col idx of each sampled node
rows_samples =[];   % collects the row idx of each sampled node



%%%%%%%%%%%%%%%%%%%
% create sampling set for foreground 
%%%%%%%%%%%%%%%%%
max_red = max(feature_matrix(:,1)) ; 
samplingset_fg = find(feature_matrix(:,1)>max_red*0.99);  
rand_idx = randperm(length(samplingset_fg)); 
samplingset_fg = samplingset_fg(rand_idx(1:floor(length(samplingset_fg)/2))); 

%%%%%%%%%%%%%%%%%%%%
% create sampling set for background 
%%%%%%%%%%%%%%%%%%%%%%

%idx = dbscan(feature_matrix,10,40);
samplingset_bg = find(feature_matrix(:,1)<max_red*0.6);  
%samplingset_bg = find(idx==1); 
rand_idx = randperm(length(samplingset_bg)); 
samplingset_bg = samplingset_bg(rand_idx(1:floor(length(samplingset_bg)/100))); 

%%%%%%%%%%%%%

samplingset = [samplingset_bg;samplingset_fg];

%D = triu(G); 
%Adjac = triu(W,1) ; 
%A_undirected = Adjac+Adjac' ; 
%degrees = sum(A_undirected,1); 
%inv_degrees = 1./degrees';

%%%% create weighted incidence matrix 
%G = digraph(triu(W,1)) ;
G = graph(nodesinedge(:,1),nodesinedge(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D = sparse(incidence(G)) ;
%D_block= kron(D,eye(d,d)) ; 

[nrnodes, nredges] = size(D); 
edge_weights = zeros(nredges,1); 
Lambda_vec = sparse(1./(sum(abs(D),1))) ; 
DtLambda= transpose(sparse(diag(Lambda_vec))*transpose(D)) ; 
%Lambda = diag(Lambda_vec) ; 
%Lambda_block = kron(Lambda,eye(d,d)) ;

Gamma_vec=sparse(1./(sum(abs(D),2)))' ;
DtGamma = transpose(D)*sparse(diag(Gamma_vec)) ; 
%Gamma = diag(Gamma_vec);  
%Gamma_block = kron(Gamma,eye(d,d)) ; 




hatx = zeros(nrnodes,1); 


lambda =1/10 ; 
lambda = 1/3; 
lambda = 100;   %%% this works find with patch_x=0,patch_y=0

labels = zeros(nrnodes,1); 
for iter_i=1:length(samplingset_bg) 
    labels(samplingset_bg(iter_i)) = 0 ; 
end
for iter_i=1:length(samplingset_fg) 
    labels(samplingset_fg(iter_i)) = 1 ; 
end

hatx = zeros(d*nrnodes,1); 
hatx = zeros(d,nrnodes)  ; 

running_average =zeros(d*nrnodes,1);
running_average =zeros(d,nrnodes);
%haty = ((1:(N-1))/(N-1))'; 
haty = zeros(nredges*d,1); 
haty = zeros(d,nredges) ; 

running_averagey = zeros(nredges*d,1); 
running_averagey = zeros(d,nredges); 

%hatxLP = zeros(N*d,1); 
hatw = zeros(d*nrnodes,1) ; 
hatw = zeros(d,nrnodes); 

newx = zeros(d,nrnodes); 
newy = zeros(d,nredges);
tic;
for iterk=1:10
    
    newx = hatx - haty*DtGamma ; 

    newx = update_x_logreg(newx,feature_matrix',samplingset,labels,Gamma_vec) ; 
    
    tildex = 2*newx-hatx; 
    
    newy = haty+0.9*tildex*DtLambda ; 
    
    haty = block_clipping(newy,d,nredges,lambda)  ;
    hatx = newx; 
    
    %% threshold weight vector nom (at 1) 
  %  weight_norms = sum(hatx.^2,1) ; 
  %  idx = find(weight_norms > 1) ; 
    
   % for iter_node=1:length(idx)
   %     idx_tmp = idx(iter_node); 
   %    hatx(:,idx_tmp) = hatx(:,idx_tmp)/weight_norms(idx_tmp) ; 
   % end
    
       
    running_average = (running_average*(iterk-1) +hatx)/iterk; 
    
    vec=(sum(running_average.*feature_matrix'));
    figure(4); 
    imagesc(reshape(vec,dim_len(1),dim_len(2)));
    hold on;
    scatter(row_col_mtx(samplingset_fg,2),row_col_mtx(samplingset_fg,1)) ; 

end
toc
runtim=toc; 

tic; 
L = kmeans(feature_matrix,2,'Replicates',5);
toc; 

RGB = pic; 
L = superpixels(pic,nrnodes);
figure(10); 
roi = true(size(L));
formask = false(size(L));
formask(row_col_mtx(samplingset_bg,1),row_col_mtx(samplingset_bg,2)) = true;
backmask = false(size(L));
BW = grabcut(pic,L,roi,samplingset_fg,samplingset_bg);
imshow(BW)
%roi(10:end-10,10:end-10,:) = true;

profile viewer
%mse_over_boundary(iter_boundary)=mse_over_boundary(iter_boundary)+norm(barw_vec-running_average)^2;


%figure(1); 
%mtx_out = reshape(running_average,[],N)' ; 
%stem(mtx_out(:,1));
%title('primal SLP')
%figure(2); 
%stem(dualSLP); 
%title('dual SLP')
%figure(2); 
%norm_w_hat=reshape(running_average,d,N) ; 
%norm_w_hat=sum(norm_w_hat.^2,1); 

%stem(norm_w_hat./sum(barw.^2,1));
%title('output nLasso');
%figure(4); 
%x_vals = sum(ratio_bound_flow,2)/RUNS ; 
%mselog = mse_over_boundary/(RUNS) ; 
%mselog=mselog/norm(barw_vec)^2 ; %norm(graphsig)^2 ; 
%stem(x_vals,mselog); 

%mtx=[x_vals mselog]; 
%mtx = flipup(mtx); 
%T = array2table(mtx,'VariableNames',{'a','b'});
%csvwrite('hingelosswoheader.csv',mtx);
%filename = sprintf('MSEoverBoundary_%s.csv',date) ;
%writetable(T,fullfile(pathtothismfile,filename));


%figure(4); 
%stem(dual);
%bound =log_bound+(1./(2*(1:K))*(hatx'*inv(Lambda)*hatx) ; 
%plot(1:N,log([log_conv log_bound])); 


function weights_out = block_clipping (weights_in,feature_dim,nr_edges,lambda) 
%%% input: weights_in vector of length featuredim*nr_datapoints, scalar
%%% lambda 
 %mtx = reshape(weights_in,feature_dim,nr_edges);
 mtx = weights_in; 
 x_norm = sqrt(sum(mtx.^2,1)) ; 

 weights_out = mtx; 
 
 idx_exceed = find(x_norm>lambda) ; 
 factor = ones(1,nr_edges); 
 
 for iter_idx=1:length(idx_exceed) 
     idx = idx_exceed(iter_idx) ; 
     factor(idx) = lambda./x_norm(idx); 
     weights_out(:,iter_idx) = mtx(:,iter_idx)*factor(idx); 
 end
 
 
 %factor = max([x_norm;lambda*ones(1,nr_edges)],[],1) ;  
 %weights_out= mtx*diag(factor);
 
% weights_out = reshape(weights_out,feature_dim*nr_edges,1) ; 
 


end 




function weights_out = update_x_logreg (weights_in,features,sampling_set,labels,Gamma_vec) 
%%% input: weights_in mtx of shape (d,N) 
%%% lambda ; labels.. shape(N,1) 
 %DiagTAU = length(sampling_set).*inv(2*diag(tau_vec));
 

 weights_out = weights_in ; % reshape(weights_in,feature_dim,nr_nodes); 
 [d N] = size(weights_in); 
 
 for iter_i =1:length(sampling_set) 
   node_idx  = sampling_set(iter_i); 
   tilde_tau = length(sampling_set)/(2*Gamma_vec(node_idx)) ; 
   vec = weights_in(:,node_idx) ; 
   feature_vec = features(:,node_idx); 
   node_y = labels(node_idx); 
   
    for iter_Newton=1:1
        predictor = vec'*feature_vec; 
        % compute gradient of loss function
        gradient = feature_vec * (node_y - sigmoid(predictor)); 
     
     % compute inverse of Hessian using "rank-1 plus identity" inversion
     % formula
     % Hessian = feature_vec * feature_vec'*(1-sigmoid(x))*sigmoid(x) + tilde_tau*eye(d,d) ; 
       fac = (1-sigmoid(predictor))*sigmoid(predictor)/tilde_tau; 
       invHessian = (1/tilde_tau)*(eye(d,d) - fac*(feature_vec * feature_vec')/(1+(1-norm(vec,2)^2*fac))) ;  
     %  update_test = (Hessian(feature_vec,vec,tilde_tau))\gradient ; %
       update = invHessian*gradient ; 
       vec = vec + update ; 
    end
   % if (norm(update)/norm(vec) > 0.5 ) 
   %     disp("error") ; 
   % end
   weights_out(:,node_idx) = vec ; 
 end 
end 



function out=sigmoid(predictor) 
    out = 1/(1+exp(-predictor)); 
end

% compute Hessian of regularized logistic loss + rho (w-w_0)^2 
% input feature vec (shape(d,1); 
%       weight_vec shape (d,1)

function H=Hessian(feature_vec,weight_vec,rho) 
    x = feature_vec'*weight_vec; 
    d = length(feature_vec); 
    H = feature_vec * feature_vec'*(1-sigmoid(x))*sigmoid(x) + rho*eye(d,d); 
end



