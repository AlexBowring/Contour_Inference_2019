function Sim_17(nSubj, SvNm, nRlz, mag, smo)
%
% Creates a 2D images of linearly increasing signal from L to R, and then applies the standardized effects Contour Inference method
% for each of the proposed options
%


%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  SvNm  = 'Normsim';  % Save name
end
if (nargin<3)
  nRlz = 5000;
end  
if exist([SvNm '.mat'], 'file')
  error('Will not overwrite sim result')
end

%------------Define parameters
% SvNm = 'LinearSig';
% nSubj  = 120_subj_124_dim;
% nRlz = 30_subj_124_dim0;
%dim     = [100 100]; 
%smo     = 3;
%mag     = 1;

thr     = mag;
rimFWHM = 2/sqrt(2*log(2));
dim     = [124 124];
dim_60  = [60 60];
dim_10  = [10 10];
stdblk  = prod(dim([1 2])/2);
stdblk_60 = prod(dim_60([1 2])/2);
stdblk_10 = prod(dim_10([1 2])/2);

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data_120_subj_124_dim = zeros([dim, nSubj]);

monte_carlo_max_resid_store_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_30_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_30_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_30_subj_124_dim = zeros([nRlz, 1]);

monte_carlo_max_resid_store_snr_std_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_120_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_60_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_30_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_30_subj_124_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_30_subj_124_dim = zeros([nRlz, 1]);

monte_carlo_max_resid_store_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_30_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_30_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_30_subj_60_dim = zeros([nRlz, 1]);

monte_carlo_max_resid_store_snr_std_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_120_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_60_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_30_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_30_subj_60_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_30_subj_60_dim = zeros([nRlz, 1]);

monte_carlo_max_resid_store_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_30_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_30_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_30_subj_10_dim = zeros([nRlz, 1]);

monte_carlo_max_resid_store_snr_std_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_120_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_60_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_max_resid_store_snr_std_30_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_min_resid_store_snr_std_30_subj_10_dim = zeros([nRlz, 1]);
monte_carlo_abs_max_resid_store_snr_std_30_subj_10_dim = zeros([nRlz, 1]);


% Creating linearly increasing signal across columns
Sig = ones(dim).*mag;
  
for t=1:nRlz
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        Noise = create_noise(wdim, 'homo', 1, smo, trnind);
        tImgs = Sig + Noise; % Creates the true image of smoothed signal + smoothed noise
        observed_data_120_subj_124_dim(:,:,i) = tImgs;
        
      end %========== Loop i (subjects)
      
      %% Getting everything for 124 x 124 dimension
      
      observed_data_60_subj_124_dim = observed_data_120_subj_124_dim(:,:,1:(nSubj/2));
      observed_data_30_subj_124_dim = observed_data_120_subj_124_dim(:,:,1:(nSubj/4));
      
      observed_mean_120_subj_124_dim = mean(observed_data_120_subj_124_dim,3);
      observed_mean_60_subj_124_dim = mean(observed_data_60_subj_124_dim,3);
      observed_mean_30_subj_124_dim = mean(observed_data_30_subj_124_dim,3);
      
      observed_std_120_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_124_dim,[prod(dim) nSubj]),stdblk),...
           dim);
      observed_std_60_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_124_dim,[prod(dim) nSubj/2]),stdblk),...
           dim);
      observed_std_30_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_124_dim,[prod(dim) nSubj/4]),stdblk),...
           dim);
       
      observed_cohen_d_120_subj_124_dim = observed_mean_120_subj_124_dim./observed_std_120_subj_124_dim;
      observed_cohen_d_60_subj_124_dim = observed_mean_60_subj_124_dim./observed_std_60_subj_124_dim;
      observed_cohen_d_30_subj_124_dim = observed_mean_30_subj_124_dim./observed_std_30_subj_124_dim;
      
      observed_cohen_d_std_120_subj_124_dim = sqrt(1+observed_cohen_d_120_subj_124_dim.^2/2);
      observed_cohen_d_std_60_subj_124_dim = sqrt(1+observed_cohen_d_60_subj_124_dim.^2/2);
      observed_cohen_d_std_30_subj_124_dim = sqrt(1+observed_cohen_d_30_subj_124_dim.^2/2);
      
      snr_resid_120_subj_124_dim     = create_resid(observed_data_120_subj_124_dim, observed_mean_120_subj_124_dim, observed_std_120_subj_124_dim, 2);
      snr_resid_60_subj_124_dim     = create_resid(observed_data_60_subj_124_dim, observed_mean_60_subj_124_dim, observed_std_60_subj_124_dim, 2);
      snr_resid_30_subj_124_dim     = create_resid(observed_data_30_subj_124_dim, observed_mean_30_subj_124_dim, observed_std_30_subj_124_dim, 2);
      
      snr_resid_std_120_subj_124_dim = reshape(...
         biasmystd(reshape(snr_resid_120_subj_124_dim,[prod(dim) nSubj]),stdblk),...
           dim);
      snr_resid_std_60_subj_124_dim = reshape(...
         biasmystd(reshape(snr_resid_60_subj_124_dim,[prod(dim) nSubj/2]),stdblk),...
           dim);
      snr_resid_std_30_subj_124_dim = reshape(...
         biasmystd(reshape(snr_resid_30_subj_124_dim,[prod(dim) nSubj/4]),stdblk),...
           dim);
      
      resid_field_120_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_124_dim - thr)./observed_cohen_d_std_120_subj_124_dim;
      abs_resid_field_120_subj_124_dim = abs(resid_field_120_subj_124_dim);
      resid_field_60_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_124_dim - thr)./observed_cohen_d_std_60_subj_124_dim;
      abs_resid_field_60_subj_124_dim = abs(resid_field_60_subj_124_dim);   
      resid_field_30_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_124_dim - thr)./observed_cohen_d_std_30_subj_124_dim;
      abs_resid_field_30_subj_124_dim = abs(resid_field_30_subj_124_dim);  
   
      
      monte_carlo_max_resid_store_120_subj_124_dim(t,1)     = max(resid_field_120_subj_124_dim(:));
      monte_carlo_min_resid_store_120_subj_124_dim(t,1)     = min(resid_field_120_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_120_subj_124_dim(t,1) = max(abs_resid_field_120_subj_124_dim(:));
      monte_carlo_max_resid_store_60_subj_124_dim(t,1)     = max(resid_field_60_subj_124_dim(:));
      monte_carlo_min_resid_store_60_subj_124_dim(t,1)     = min(resid_field_60_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_60_subj_124_dim(t,1) = max(abs_resid_field_60_subj_124_dim(:));
      monte_carlo_max_resid_store_30_subj_124_dim(t,1)     = max(resid_field_30_subj_124_dim(:));
      monte_carlo_min_resid_store_30_subj_124_dim(t,1)     = min(resid_field_30_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_30_subj_124_dim(t,1) = max(abs_resid_field_30_subj_124_dim(:));
      
      resid_field_snr_std_120_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_124_dim - thr)./snr_resid_std_120_subj_124_dim;
      abs_resid_field_snr_std_120_subj_124_dim = abs(resid_field_snr_std_120_subj_124_dim);
      resid_field_snr_std_60_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_124_dim - thr)./snr_resid_std_60_subj_124_dim;
      abs_resid_field_snr_std_60_subj_124_dim = abs(resid_field_snr_std_60_subj_124_dim);
      resid_field_snr_std_30_subj_124_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_124_dim - thr)./snr_resid_std_30_subj_124_dim;
      abs_resid_field_snr_std_30_subj_124_dim = abs(resid_field_snr_std_30_subj_124_dim);
      
      monte_carlo_max_resid_store_snr_std_120_subj_124_dim(t,1)     = max(resid_field_snr_std_120_subj_124_dim(:));
      monte_carlo_min_resid_store_snr_std_120_subj_124_dim(t,1)     = min(resid_field_snr_std_120_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_120_subj_124_dim(t,1) = max(abs_resid_field_snr_std_120_subj_124_dim(:));
      monte_carlo_max_resid_store_snr_std_60_subj_124_dim(t,1)     = max(resid_field_snr_std_60_subj_124_dim(:));
      monte_carlo_min_resid_store_snr_std_60_subj_124_dim(t,1)     = min(resid_field_snr_std_60_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_60_subj_124_dim(t,1) = max(abs_resid_field_snr_std_60_subj_124_dim(:));
      monte_carlo_max_resid_store_snr_std_30_subj_124_dim(t,1)     = max(resid_field_snr_std_30_subj_124_dim(:));
      monte_carlo_min_resid_store_snr_std_30_subj_124_dim(t,1)     = min(resid_field_snr_std_30_subj_124_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_30_subj_124_dim(t,1) = max(abs_resid_field_snr_std_30_subj_124_dim(:));
      
     
 %% Getting everything for 60 x 60 dimension
      observed_data_120_subj_60_dim = observed_data_120_subj_124_dim(1:60,1:60,1:(nSubj));
      observed_data_60_subj_60_dim = observed_data_120_subj_124_dim(1:60,1:60,1:(nSubj/2));
      observed_data_30_subj_60_dim = observed_data_120_subj_124_dim(1:60,1:60,1:(nSubj/4));
      
      observed_mean_120_subj_60_dim = mean(observed_data_120_subj_60_dim,3);
      observed_mean_60_subj_60_dim = mean(observed_data_60_subj_60_dim,3);
      observed_mean_30_subj_60_dim = mean(observed_data_30_subj_60_dim,3);

      observed_std_120_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_60_dim,[prod(dim_60) nSubj]),stdblk_60),...
           dim_60);
      observed_std_60_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_60_dim,[prod(dim_60) nSubj/2]),stdblk_60),...
           dim_60);
      observed_std_30_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_60_dim,[prod(dim_60) nSubj/4]),stdblk_60),...
           dim_60);
       
      observed_cohen_d_120_subj_60_dim = observed_mean_120_subj_60_dim./observed_std_120_subj_60_dim;
      observed_cohen_d_60_subj_60_dim = observed_mean_60_subj_60_dim./observed_std_60_subj_60_dim;
      observed_cohen_d_30_subj_60_dim = observed_mean_30_subj_60_dim./observed_std_30_subj_60_dim;
      
      observed_cohen_d_std_120_subj_60_dim = sqrt(1+observed_cohen_d_120_subj_60_dim.^2/2);
      observed_cohen_d_std_60_subj_60_dim = sqrt(1+observed_cohen_d_60_subj_60_dim.^2/2);
      observed_cohen_d_std_30_subj_60_dim = sqrt(1+observed_cohen_d_30_subj_60_dim.^2/2);
      
      snr_resid_120_subj_60_dim     = create_resid(observed_data_120_subj_60_dim, observed_mean_120_subj_60_dim, observed_std_120_subj_60_dim, 2);
      snr_resid_60_subj_60_dim     = create_resid(observed_data_60_subj_60_dim, observed_mean_60_subj_60_dim, observed_std_60_subj_60_dim, 2);
      snr_resid_30_subj_60_dim     = create_resid(observed_data_30_subj_60_dim, observed_mean_30_subj_60_dim, observed_std_30_subj_60_dim, 2);
      
      snr_resid_std_120_subj_60_dim = reshape(...
         biasmystd(reshape(snr_resid_120_subj_60_dim,[prod(dim_60) nSubj]),stdblk_60),...
           dim_60);
      snr_resid_std_60_subj_60_dim = reshape(...
         biasmystd(reshape(snr_resid_60_subj_60_dim,[prod(dim_60) nSubj/2]),stdblk_60),...
           dim_60);
      snr_resid_std_30_subj_60_dim = reshape(...
         biasmystd(reshape(snr_resid_30_subj_60_dim,[prod(dim_60) nSubj/4]),stdblk_60),...
           dim_60);
      
      resid_field_120_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_60_dim - thr)./observed_cohen_d_std_120_subj_60_dim;
      abs_resid_field_120_subj_60_dim = abs(resid_field_120_subj_60_dim);
      resid_field_60_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_60_dim - thr)./observed_cohen_d_std_60_subj_60_dim;
      abs_resid_field_60_subj_60_dim = abs(resid_field_60_subj_60_dim);   
      resid_field_30_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_60_dim - thr)./observed_cohen_d_std_30_subj_60_dim;
      abs_resid_field_30_subj_60_dim = abs(resid_field_30_subj_60_dim);  
   
      
      monte_carlo_max_resid_store_120_subj_60_dim(t,1)     = max(resid_field_120_subj_60_dim(:));
      monte_carlo_min_resid_store_120_subj_60_dim(t,1)     = min(resid_field_120_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_120_subj_60_dim(t,1) = max(abs_resid_field_120_subj_60_dim(:));
      monte_carlo_max_resid_store_60_subj_60_dim(t,1)     = max(resid_field_60_subj_60_dim(:));
      monte_carlo_min_resid_store_60_subj_60_dim(t,1)     = min(resid_field_60_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_60_subj_60_dim(t,1) = max(abs_resid_field_60_subj_60_dim(:));
      monte_carlo_max_resid_store_30_subj_60_dim(t,1)     = max(resid_field_30_subj_60_dim(:));
      monte_carlo_min_resid_store_30_subj_60_dim(t,1)     = min(resid_field_30_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_30_subj_60_dim(t,1) = max(abs_resid_field_30_subj_60_dim(:));
      
      resid_field_snr_std_120_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_60_dim - thr)./snr_resid_std_120_subj_60_dim;
      abs_resid_field_snr_std_120_subj_60_dim = abs(resid_field_snr_std_120_subj_60_dim);
      resid_field_snr_std_60_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_60_dim - thr)./snr_resid_std_60_subj_60_dim;
      abs_resid_field_snr_std_60_subj_60_dim = abs(resid_field_snr_std_60_subj_60_dim);
      resid_field_snr_std_30_subj_60_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_60_dim - thr)./snr_resid_std_30_subj_60_dim;
      abs_resid_field_snr_std_30_subj_60_dim = abs(resid_field_snr_std_30_subj_60_dim);
      
      monte_carlo_max_resid_store_snr_std_120_subj_60_dim(t,1)     = max(resid_field_snr_std_120_subj_60_dim(:));
      monte_carlo_min_resid_store_snr_std_120_subj_60_dim(t,1)     = min(resid_field_snr_std_120_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_120_subj_60_dim(t,1) = max(abs_resid_field_snr_std_120_subj_60_dim(:));
      monte_carlo_max_resid_store_snr_std_60_subj_60_dim(t,1)     = max(resid_field_snr_std_60_subj_60_dim(:));
      monte_carlo_min_resid_store_snr_std_60_subj_60_dim(t,1)     = min(resid_field_snr_std_60_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_60_subj_60_dim(t,1) = max(abs_resid_field_snr_std_60_subj_60_dim(:));
      monte_carlo_max_resid_store_snr_std_30_subj_60_dim(t,1)     = max(resid_field_snr_std_30_subj_60_dim(:));
      monte_carlo_min_resid_store_snr_std_30_subj_60_dim(t,1)     = min(resid_field_snr_std_30_subj_60_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_30_subj_60_dim(t,1) = max(abs_resid_field_snr_std_30_subj_60_dim(:));
      
      %% Getting everything for 10 x 10 dimension
      
      observed_data_120_subj_10_dim = observed_data_120_subj_124_dim(1:10,1:10,1:(nSubj));
      observed_data_60_subj_10_dim = observed_data_120_subj_124_dim(1:10,1:10,1:(nSubj/2));
      observed_data_30_subj_10_dim = observed_data_120_subj_124_dim(1:10,1:10,1:(nSubj/4));
      
      observed_mean_120_subj_10_dim = mean(observed_data_120_subj_10_dim,3);
      observed_mean_60_subj_10_dim = mean(observed_data_60_subj_10_dim,3);
      observed_mean_30_subj_10_dim = mean(observed_data_30_subj_10_dim,3);

      observed_std_120_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_10_dim,[prod(dim_10) nSubj]),stdblk_10),...
           dim_10);
      observed_std_60_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_10_dim,[prod(dim_10) nSubj/2]),stdblk_10),...
           dim_10);
      observed_std_30_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_10_dim,[prod(dim_10) nSubj/4]),stdblk_10),...
           dim_10);
       
      observed_cohen_d_120_subj_10_dim = observed_mean_120_subj_10_dim./observed_std_120_subj_10_dim;
      observed_cohen_d_60_subj_10_dim = observed_mean_60_subj_10_dim./observed_std_60_subj_10_dim;
      observed_cohen_d_30_subj_10_dim = observed_mean_30_subj_10_dim./observed_std_30_subj_10_dim;
      
      observed_cohen_d_std_120_subj_10_dim = sqrt(1+observed_cohen_d_120_subj_10_dim.^2/2);
      observed_cohen_d_std_60_subj_10_dim = sqrt(1+observed_cohen_d_60_subj_10_dim.^2/2);
      observed_cohen_d_std_30_subj_10_dim = sqrt(1+observed_cohen_d_30_subj_10_dim.^2/2);
      
      snr_resid_120_subj_10_dim     = create_resid(observed_data_120_subj_10_dim, observed_mean_120_subj_10_dim, observed_std_120_subj_10_dim, 2);
      snr_resid_60_subj_10_dim     = create_resid(observed_data_60_subj_10_dim, observed_mean_60_subj_10_dim, observed_std_60_subj_10_dim, 2);
      snr_resid_30_subj_10_dim     = create_resid(observed_data_30_subj_10_dim, observed_mean_30_subj_10_dim, observed_std_30_subj_10_dim, 2);
      
      snr_resid_std_120_subj_10_dim = reshape(...
         biasmystd(reshape(snr_resid_120_subj_10_dim,[prod(dim_10) nSubj]),stdblk_10),...
           dim_10);
      snr_resid_std_60_subj_10_dim = reshape(...
         biasmystd(reshape(snr_resid_60_subj_10_dim,[prod(dim_10) nSubj/2]),stdblk_10),...
           dim_10);
      snr_resid_std_30_subj_10_dim = reshape(...
         biasmystd(reshape(snr_resid_30_subj_10_dim,[prod(dim_10) nSubj/4]),stdblk_10),...
           dim_10);
      
      resid_field_120_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_10_dim - thr)./observed_cohen_d_std_120_subj_10_dim;
      abs_resid_field_120_subj_10_dim = abs(resid_field_120_subj_10_dim);
      resid_field_60_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_10_dim - thr)./observed_cohen_d_std_60_subj_10_dim;
      abs_resid_field_60_subj_10_dim = abs(resid_field_60_subj_10_dim);   
      resid_field_30_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_10_dim - thr)./observed_cohen_d_std_30_subj_10_dim;
      abs_resid_field_30_subj_10_dim = abs(resid_field_30_subj_10_dim);  
   
      
      monte_carlo_max_resid_store_120_subj_10_dim(t,1)     = max(resid_field_120_subj_10_dim(:));
      monte_carlo_min_resid_store_120_subj_10_dim(t,1)     = min(resid_field_120_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_120_subj_10_dim(t,1) = max(abs_resid_field_120_subj_10_dim(:));
      monte_carlo_max_resid_store_60_subj_10_dim(t,1)     = max(resid_field_60_subj_10_dim(:));
      monte_carlo_min_resid_store_60_subj_10_dim(t,1)     = min(resid_field_60_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_60_subj_10_dim(t,1) = max(abs_resid_field_60_subj_10_dim(:));
      monte_carlo_max_resid_store_30_subj_10_dim(t,1)     = max(resid_field_30_subj_10_dim(:));
      monte_carlo_min_resid_store_30_subj_10_dim(t,1)     = min(resid_field_30_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_30_subj_10_dim(t,1) = max(abs_resid_field_30_subj_10_dim(:));
      
      resid_field_snr_std_120_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_120_subj_10_dim - thr)./snr_resid_std_120_subj_10_dim;
      abs_resid_field_snr_std_120_subj_10_dim = abs(resid_field_snr_std_120_subj_10_dim);
      resid_field_snr_std_60_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_60_subj_10_dim - thr)./snr_resid_std_60_subj_10_dim;
      abs_resid_field_snr_std_60_subj_10_dim = abs(resid_field_snr_std_60_subj_10_dim);
      resid_field_snr_std_30_subj_10_dim = sqrt(nSubj)*(observed_cohen_d_30_subj_10_dim - thr)./snr_resid_std_30_subj_10_dim;
      abs_resid_field_snr_std_30_subj_10_dim = abs(resid_field_snr_std_30_subj_10_dim);
      
      monte_carlo_max_resid_store_snr_std_120_subj_10_dim(t,1)     = max(resid_field_snr_std_120_subj_10_dim(:));
      monte_carlo_min_resid_store_snr_std_120_subj_10_dim(t,1)     = min(resid_field_snr_std_120_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_120_subj_10_dim(t,1) = max(abs_resid_field_snr_std_120_subj_10_dim(:));
      monte_carlo_max_resid_store_snr_std_60_subj_10_dim(t,1)     = max(resid_field_snr_std_60_subj_10_dim(:));
      monte_carlo_min_resid_store_snr_std_60_subj_10_dim(t,1)     = min(resid_field_snr_std_60_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_60_subj_10_dim(t,1) = max(abs_resid_field_snr_std_60_subj_10_dim(:));
      monte_carlo_max_resid_store_snr_std_30_subj_10_dim(t,1)     = max(resid_field_snr_std_30_subj_10_dim(:));
      monte_carlo_min_resid_store_snr_std_30_subj_10_dim(t,1)     = min(resid_field_snr_std_30_subj_10_dim(:));
      monte_carlo_abs_max_resid_store_snr_std_30_subj_10_dim(t,1) = max(abs_resid_field_snr_std_30_subj_10_dim(:));
end

%% Saving variables

eval(['save ' SvNm ' nRlz dim smo mag '... 
      'monte_carlo_max_resid_store_120_subj_124_dim monte_carlo_min_resid_store_120_subj_124_dim monte_carlo_abs_max_resid_store_120_subj_124_dim monte_carlo_max_resid_store_snr_std_120_subj_124_dim monte_carlo_min_resid_store_snr_std_120_subj_124_dim monte_carlo_abs_max_resid_store_snr_std_120_subj_124_dim ' ...
      'monte_carlo_max_resid_store_60_subj_124_dim monte_carlo_min_resid_store_60_subj_124_dim monte_carlo_abs_max_resid_store_60_subj_124_dim monte_carlo_max_resid_store_snr_std_60_subj_124_dim monte_carlo_min_resid_store_snr_std_60_subj_124_dim monte_carlo_abs_max_resid_store_snr_std_60_subj_124_dim ' ...
      'monte_carlo_max_resid_store_30_subj_124_dim monte_carlo_min_resid_store_30_subj_124_dim monte_carlo_abs_max_resid_store_30_subj_124_dim monte_carlo_max_resid_store_snr_std_30_subj_124_dim monte_carlo_min_resid_store_snr_std_30_subj_124_dim monte_carlo_abs_max_resid_store_snr_std_30_subj_124_dim ' ...
      'monte_carlo_max_resid_store_120_subj_60_dim monte_carlo_min_resid_store_120_subj_60_dim monte_carlo_abs_max_resid_store_120_subj_60_dim monte_carlo_max_resid_store_snr_std_120_subj_60_dim monte_carlo_min_resid_store_snr_std_120_subj_60_dim monte_carlo_abs_max_resid_store_snr_std_120_subj_60_dim ' ...
      'monte_carlo_max_resid_store_60_subj_60_dim monte_carlo_min_resid_store_60_subj_60_dim monte_carlo_abs_max_resid_store_60_subj_60_dim monte_carlo_max_resid_store_snr_std_60_subj_60_dim monte_carlo_min_resid_store_snr_std_60_subj_60_dim monte_carlo_abs_max_resid_store_snr_std_60_subj_60_dim ' ...
      'monte_carlo_max_resid_store_30_subj_60_dim monte_carlo_min_resid_store_30_subj_60_dim monte_carlo_abs_max_resid_store_30_subj_60_dim monte_carlo_max_resid_store_snr_std_30_subj_60_dim monte_carlo_min_resid_store_snr_std_30_subj_60_dim monte_carlo_abs_max_resid_store_snr_std_30_subj_60_dim '...
      'monte_carlo_max_resid_store_120_subj_10_dim monte_carlo_min_resid_store_120_subj_10_dim monte_carlo_abs_max_resid_store_120_subj_10_dim monte_carlo_max_resid_store_snr_std_120_subj_10_dim monte_carlo_min_resid_store_snr_std_120_subj_10_dim monte_carlo_abs_max_resid_store_snr_std_120_subj_10_dim ' ...
      'monte_carlo_max_resid_store_60_subj_10_dim monte_carlo_min_resid_store_60_subj_10_dim monte_carlo_abs_max_resid_store_60_subj_10_dim monte_carlo_max_resid_store_snr_std_60_subj_10_dim monte_carlo_min_resid_store_snr_std_60_subj_10_dim monte_carlo_abs_max_resid_store_snr_std_60_subj_10_dim ' ...
      'monte_carlo_max_resid_store_30_subj_10_dim monte_carlo_min_resid_store_30_subj_10_dim monte_carlo_abs_max_resid_store_30_subj_10_dim monte_carlo_max_resid_store_snr_std_30_subj_10_dim monte_carlo_min_resid_store_snr_std_30_subj_10_dim monte_carlo_abs_max_resid_store_snr_std_30_subj_10_dim'])
