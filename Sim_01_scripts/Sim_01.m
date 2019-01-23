function Sim_01(nSubj,SvNm,nRlz)
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
% nSubj  = 120;
% nRlz = 300;

tau     = 1/sqrt(nSubj);
nBoot   = 5000;
dim     = [100 100]; 
smo     = 3;
mag     = 3;
rimFWHM = 2/sqrt(2*log(2)); 				 
stdblk  = prod(dim([1 2])/2);
thr     = 2;
rad     = 30;

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data  = zeros([dim nSubj]);

% This stores the vector SupG for each run
% This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
subset_success_vector_raw_80           = zeros(nRlz, 1); 
subset_success_vector_raw_90           = zeros(nRlz, 1);
subset_success_vector_raw_95           = zeros(nRlz, 1);
subset_success_vector_observed_80           = zeros(nRlz, 1); 
subset_success_vector_observed_90           = zeros(nRlz, 1);
subset_success_vector_observed_95           = zeros(nRlz, 1);
subset_success_vector_raw_80_alternate = zeros(nRlz, 1); 
subset_success_vector_raw_90_alternate = zeros(nRlz, 1);
subset_success_vector_raw_95_alternate = zeros(nRlz, 1);
subset_success_vector_observed_80_alternate = zeros(nRlz, 1); 
subset_success_vector_observed_90_alternate = zeros(nRlz, 1);
subset_success_vector_observed_95_alternate = zeros(nRlz, 1);

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

threshold_observed_80_store                  = zeros(nRlz, 1);
threshold_observed_90_store                  = zeros(nRlz, 1);
threshold_observed_95_store                  = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store                   = zeros(nBoot, nRlz);
supG_observed_store              = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_80_store                       = zeros([nRlz dim]);
upper_contour_raw_80_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_80_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_80_store                    = zeros([nRlz dim]);

lower_contour_raw_90_store                       = zeros([nRlz dim]);
upper_contour_raw_90_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_90_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_90_store                    = zeros([nRlz dim]);

lower_contour_raw_95_store                       = zeros([nRlz dim]);
upper_contour_raw_95_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_95_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_95_store                    = zeros([nRlz dim]);

lower_contour_observed_80_store                       = zeros([nRlz dim]);
upper_contour_observed_80_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_80_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_80_store                    = zeros([nRlz dim]);

lower_contour_observed_90_store                       = zeros([nRlz dim]);
upper_contour_observed_90_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_90_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_90_store                    = zeros([nRlz dim]);

lower_contour_observed_95_store                       = zeros([nRlz dim]);
upper_contour_observed_95_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_95_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_95_store                    = zeros([nRlz dim]);

supG_raw                         = zeros(nBoot,1);
supG_observed                    = zeros(nBoot,1);

% Creating linearly increasing signal across columns
Sig = create_signal(wdim, 'circle', [mag, rad], smo, trnind);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;

Sig_boundary_edges   = getBdryparams(Sig, thr);
  
for t=1:nRlz
    fprintf('.');
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        Noise = create_noise(wdim, 'homo', 1, smo, trnind);
        tImgs = Sig + Noise; % Creates the true image of smoothed signal + smoothed noise
        observed_data(:,:,i) = tImgs;
        
      end %========== Loop i (subjects)

      observed_mean = mean(observed_data,3);

      observed_std = reshape(...
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim);
       
      observed_cohen_d = observed_mean./observed_std;
      
      observed_cohen_d_std = sqrt(1+observed_cohen_d.^2/2); 
      
      observed_boundary_edges = getBdryparams(observed_cohen_d, thr);
      
      % Residuals
      unstandardized_resid = create_resid(observed_data, observed_mean, observed_std, 1);
      cohen_resid          = create_resid(observed_data, observed_mean, observed_std, 2);
            
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          % Obtaining the bootstrap std of the normal, untransformed
          % residuals
          signflips                              = randi(2,[nSubj,1])*2-3;
          unstandardized_resid_bootstrap         = unstandardized_resid*spdiags(signflips, 0, nSubj, nSubj);
          unstandardized_resid_bootstrap         = reshape(unstandardized_resid_bootstrap, [dim nSubj]);
          data_bootstrap                         = bsxfun(@plus, unstandardized_resid_bootstrap, observed_mean);
          boot_mean                              = mean(data_bootstrap,3);
          boot_std = reshape(...
            biasmystd(reshape(data_bootstrap,[prod(dim) nSubj]),stdblk),...
                dim);
          boot_cohen_d     = boot_mean./boot_std;
          boot_cohen_d_std = sqrt(1+boot_cohen_d.^2/2);
          
          % Creating the cohens d residual field and standardizing
          cohen_resid_bootstrap                  = cohen_resid*spdiags(signflips, 0, nSubj, nSubj);
          cohen_resid_bootstrap                  = reshape(cohen_resid_bootstrap, [dim nSubj]);
          cohen_resid_field                      = sum(cohen_resid_bootstrap, 3)/sqrt(nSubj);
          cohen_resid_field                      = cohen_resid_field./boot_cohen_d_std;
          
          % Calculating the maximum over the weighted interpolated true boundary edges
          true_boundary_values = getBdryvalues(cohen_resid_field, Sig_boundary_edges);
          supG_raw(k)          = max(abs(true_boundary_values)); 
          
          % Calculating the maximum over the weighted interpolated observed boundary edges
          weighted_boundary_values = getBdryvalues(cohen_resid_field, observed_boundary_edges);
          supG_observed(k)         = max(abs(weighted_boundary_values));   
      end
      
    middle_contour                = AC;
    middle_contour_volume         = sum(middle_contour(:));
    
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
       
    lower_contour_raw_80             = observed_cohen_d >= thr - supGa_raw_80*tau*observed_cohen_d_std;
    upper_contour_raw_80             = observed_cohen_d >= thr + supGa_raw_80*tau*observed_cohen_d_std;
    lower_contour_raw_80_volume_prct = sum(lower_contour_raw_80(:))/middle_contour_volume;
    upper_contour_raw_80_volume_prct = sum(upper_contour_raw_80(:))/middle_contour_volume;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_90             = observed_cohen_d >= thr - supGa_raw_90*tau*observed_cohen_d_std;
    upper_contour_raw_90             = observed_cohen_d >= thr + supGa_raw_90*tau*observed_cohen_d_std;
    lower_contour_raw_90_volume_prct = sum(lower_contour_raw_90(:))/middle_contour_volume;
    upper_contour_raw_90_volume_prct = sum(upper_contour_raw_90(:))/middle_contour_volume;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;    
    
    lower_contour_raw_95             = observed_cohen_d >= thr - supGa_raw_95*tau*observed_cohen_d_std;
    upper_contour_raw_95             = observed_cohen_d >= thr + supGa_raw_95*tau*observed_cohen_d_std;
    lower_contour_raw_95_volume_prct = sum(lower_contour_raw_95(:))/middle_contour_volume;
    upper_contour_raw_95_volume_prct = sum(upper_contour_raw_95(:))/middle_contour_volume;
    mid_on_upper_raw_95              = upper_contour_raw_95.*middle_contour;
    lower_on_mid_raw_95              = middle_contour.*lower_contour_raw_95;
    upper_subset_mid_raw_95          = upper_contour_raw_95 - mid_on_upper_raw_95;
    mid_subset_lower_raw_95          = middle_contour - lower_on_mid_raw_95;
    
    % Observed boundary
    supGa_observed_80                     = prctile(supG_observed, 80);
    supGa_observed_90                     = prctile(supG_observed, 90);
    supGa_observed_95                     = prctile(supG_observed, 95);
       
    lower_contour_observed_80             = observed_cohen_d >= thr - supGa_observed_80*tau*observed_cohen_d_std;
    upper_contour_observed_80             = observed_cohen_d >= thr + supGa_observed_80*tau*observed_cohen_d_std;
    lower_contour_observed_80_volume_prct = sum(lower_contour_observed_80(:))/middle_contour_volume;
    upper_contour_observed_80_volume_prct = sum(upper_contour_observed_80(:))/middle_contour_volume;
    mid_on_upper_observed_80              = upper_contour_observed_80.*middle_contour;
    lower_on_mid_observed_80              = middle_contour.*lower_contour_observed_80;
    upper_subset_mid_observed_80          = upper_contour_observed_80 - mid_on_upper_observed_80;
    mid_subset_lower_observed_80          = middle_contour - lower_on_mid_observed_80;
    
    lower_contour_observed_90             = observed_cohen_d >= thr - supGa_observed_90*tau*observed_cohen_d_std;
    upper_contour_observed_90             = observed_cohen_d >= thr + supGa_observed_90*tau*observed_cohen_d_std;
    lower_contour_observed_90_volume_prct = sum(lower_contour_observed_90(:))/middle_contour_volume;
    upper_contour_observed_90_volume_prct = sum(upper_contour_observed_90(:))/middle_contour_volume;
    mid_on_upper_observed_90              = upper_contour_observed_90.*middle_contour;
    lower_on_mid_observed_90              = middle_contour.*lower_contour_observed_90;
    upper_subset_mid_observed_90          = upper_contour_observed_90 - mid_on_upper_observed_90;
    mid_subset_lower_observed_90          = middle_contour - lower_on_mid_observed_90;    
    
    lower_contour_observed_95             = observed_cohen_d >= thr - supGa_observed_95*tau*observed_cohen_d_std;
    upper_contour_observed_95             = observed_cohen_d >= thr + supGa_observed_95*tau*observed_cohen_d_std;
    lower_contour_observed_95_volume_prct = sum(lower_contour_observed_95(:))/middle_contour_volume;
    upper_contour_observed_95_volume_prct = sum(upper_contour_observed_95(:))/middle_contour_volume;
    mid_on_upper_observed_95              = upper_contour_observed_95.*middle_contour;
    lower_on_mid_observed_95              = middle_contour.*lower_contour_observed_95;
    upper_subset_mid_observed_95          = upper_contour_observed_95 - mid_on_upper_observed_95;
    mid_subset_lower_observed_95          = middle_contour - lower_on_mid_observed_95;

    %
    % Storing all variables of interest
    %
    % True boundary variables
    supG_raw_store(:,t)                                    = supG_raw;
    threshold_raw_80_store(t)                              = supGa_raw_80;
    lower_contour_raw_80_store(t,:,:)                      = lower_contour_raw_80;
    upper_contour_raw_80_store(t,:,:)                      = upper_contour_raw_80;
    upper_subset_mid_raw_80_store(t,:,:)                   = upper_subset_mid_raw_80;
    mid_subset_lower_raw_80_store(t,:,:)                   = mid_subset_lower_raw_80;
    lower_contour_raw_80_volume_prct_store(t)              = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)              = upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                              = supGa_raw_90;
    lower_contour_raw_90_store(t,:,:)                      = lower_contour_raw_90;
    upper_contour_raw_90_store(t,:,:)                      = upper_contour_raw_90;
    upper_subset_mid_raw_90_store(t,:,:)                   = upper_subset_mid_raw_90;
    mid_subset_lower_raw_90_store(t,:,:)                   = mid_subset_lower_raw_90;
    lower_contour_raw_90_volume_prct_store(t)              = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)              = upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                              = supGa_raw_95;
    lower_contour_raw_95_store(t,:,:)                      = lower_contour_raw_95;
    upper_contour_raw_95_store(t,:,:)                      = upper_contour_raw_95;
    upper_subset_mid_raw_95_store(t,:,:)                   = upper_subset_mid_raw_95;
    mid_subset_lower_raw_95_store(t,:,:)                   = mid_subset_lower_raw_95;
    lower_contour_raw_95_volume_prct_store(t)              = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)              = upper_contour_raw_95_volume_prct;
    
    % Observed boundary variables
    supG_observed_store(:,t)                                    = supG_observed;
    threshold_observed_80_store(t)                              = supGa_observed_80;
    lower_contour_observed_80_store(t,:,:)                      = lower_contour_observed_80;
    upper_contour_observed_80_store(t,:,:)                      = upper_contour_observed_80;
    upper_subset_mid_observed_80_store(t,:,:)                   = upper_subset_mid_observed_80;
    mid_subset_lower_observed_80_store(t,:,:)                   = mid_subset_lower_observed_80;
    lower_contour_observed_80_volume_prct_store(t)              = lower_contour_observed_80_volume_prct;
    upper_contour_observed_80_volume_prct_store(t)              = upper_contour_observed_80_volume_prct;
 
    threshold_observed_90_store(t)                              = supGa_observed_90;
    lower_contour_observed_90_store(t,:,:)                      = lower_contour_observed_90;
    upper_contour_observed_90_store(t,:,:)                      = upper_contour_observed_90;
    upper_subset_mid_observed_90_store(t,:,:)                   = upper_subset_mid_observed_90;
    mid_subset_lower_observed_90_store(t,:,:)                   = mid_subset_lower_observed_90;
    lower_contour_observed_90_volume_prct_store(t)              = lower_contour_observed_90_volume_prct;
    upper_contour_observed_90_volume_prct_store(t)              = upper_contour_observed_90_volume_prct;

    threshold_observed_95_store(t)                              = supGa_observed_95;
    lower_contour_observed_95_store(t,:,:)                      = lower_contour_observed_95;
    upper_contour_observed_95_store(t,:,:)                      = upper_contour_observed_95;
    upper_subset_mid_observed_95_store(t,:,:)                   = upper_subset_mid_observed_95;
    mid_subset_lower_observed_95_store(t,:,:)                   = mid_subset_lower_observed_95;
    lower_contour_observed_95_volume_prct_store(t)              = lower_contour_observed_95_volume_prct;
    upper_contour_observed_95_volume_prct_store(t)              = upper_contour_observed_95_volume_prct;    
    
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the true boundary
    lower_condition_80 = thr - supGa_raw_80*tau*observed_cohen_d_std;
    upper_condition_80 = thr + supGa_raw_80*tau*observed_cohen_d_std;
    lower_condition_90 = thr - supGa_raw_90*tau*observed_cohen_d_std;
    upper_condition_90 = thr + supGa_raw_90*tau*observed_cohen_d_std;
    lower_condition_95 = thr - supGa_raw_95*tau*observed_cohen_d_std;
    upper_condition_95 = thr + supGa_raw_95*tau*observed_cohen_d_std;
    
    lower_condition_80_boundary_values = getBdryvalues(lower_condition_80, Sig_boundary_edges);
    upper_condition_80_boundary_values = getBdryvalues(upper_condition_80, Sig_boundary_edges);

    lower_condition_90_boundary_values = getBdryvalues(lower_condition_90, Sig_boundary_edges);
    upper_condition_90_boundary_values = getBdryvalues(upper_condition_90, Sig_boundary_edges);
    
    lower_condition_95_boundary_values = getBdryvalues(lower_condition_95, Sig_boundary_edges);
    upper_condition_95_boundary_values = getBdryvalues(upper_condition_95, Sig_boundary_edges);
    
    observed_cohen_d_true_boundary_values = getBdryvalues(observed_cohen_d, Sig_boundary_edges);
    
    lower_condition_80_success = observed_cohen_d_true_boundary_values < lower_condition_80_boundary_values;
    upper_condition_80_success = observed_cohen_d_true_boundary_values >= upper_condition_80_boundary_values;
    
    lower_condition_90_success = observed_cohen_d_true_boundary_values < lower_condition_90_boundary_values;
    upper_condition_90_success = observed_cohen_d_true_boundary_values >= upper_condition_90_boundary_values;
    
    lower_condition_95_success = observed_cohen_d_true_boundary_values < lower_condition_95_boundary_values;
    upper_condition_95_success = observed_cohen_d_true_boundary_values >= upper_condition_95_boundary_values;
                              
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the observed boundary
    lower_condition_80_observed = thr - supGa_observed_80*tau*observed_cohen_d_std;
    upper_condition_80_observed = thr + supGa_observed_80*tau*observed_cohen_d_std;
    lower_condition_90_observed = thr - supGa_observed_90*tau*observed_cohen_d_std;
    upper_condition_90_observed = thr + supGa_observed_90*tau*observed_cohen_d_std;
    lower_condition_95_observed = thr - supGa_observed_95*tau*observed_cohen_d_std;
    upper_condition_95_observed = thr + supGa_observed_95*tau*observed_cohen_d_std;
    
    lower_condition_80_observed_boundary_values = getBdryvalues(lower_condition_80_observed, Sig_boundary_edges);
    upper_condition_80_observed_boundary_values = getBdryvalues(upper_condition_80_observed, Sig_boundary_edges);
    
    lower_condition_90_observed_boundary_values = getBdryvalues(lower_condition_90_observed, Sig_boundary_edges);
    upper_condition_90_observed_boundary_values = getBdryvalues(upper_condition_90_observed, Sig_boundary_edges);
    
    lower_condition_95_observed_boundary_values = getBdryvalues(lower_condition_95_observed, Sig_boundary_edges);
    upper_condition_95_observed_boundary_values = getBdryvalues(upper_condition_95_observed, Sig_boundary_edges);
    
    lower_condition_80_observed_success = observed_cohen_d_true_boundary_values < lower_condition_80_observed_boundary_values;
    upper_condition_80_observed_success = observed_cohen_d_true_boundary_values >= upper_condition_80_observed_boundary_values;
    
    lower_condition_90_observed_success = observed_cohen_d_true_boundary_values < lower_condition_90_observed_boundary_values;
    upper_condition_90_observed_success = observed_cohen_d_true_boundary_values >= upper_condition_90_observed_boundary_values;
    
    lower_condition_95_observed_success = observed_cohen_d_true_boundary_values < lower_condition_95_observed_boundary_values;
    upper_condition_95_observed_success = observed_cohen_d_true_boundary_values >= upper_condition_95_observed_boundary_values;
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the true boundary in mult. bootstrap
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))==0
      subset_success_vector_raw_80(t) = 1;
      fprintf('raw nominal 80 success! \n');
    else 
      subset_success_vector_raw_80(t) = 0; 
      fprintf('raw nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))==0
      subset_success_vector_raw_90(t) = 1;
      fprintf('raw nominal 90 success! \n');
    else 
      subset_success_vector_raw_90(t) = 0; 
      fprintf('raw nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))==0
      subset_success_vector_raw_95(t) = 1;
      fprintf('raw nominal 95 success! \n');
    else 
      subset_success_vector_raw_95(t) = 0; 
      fprintf('raw nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the observed boundary in mult. bootstrap
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))==0
      subset_success_vector_observed_80(t) = 1;
      fprintf('observed nominal 80 success! \n');
    else 
      subset_success_vector_observed_80(t) = 0; 
      fprintf('observed nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))==0
      subset_success_vector_observed_90(t) = 1;
      fprintf('observed nominal 90 success! \n');
    else 
      subset_success_vector_observed_90(t) = 0; 
      fprintf('observed nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))==0
      subset_success_vector_observed_95(t) = 1;
      fprintf('observed nominal 95 success! \n');
    else 
      subset_success_vector_observed_95(t) = 0; 
      fprintf('observed nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the true boundary
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:)+sum(upper_condition_80_success)+sum(lower_condition_80_success))==0
      subset_success_vector_raw_80_alternate(t) = 1;
      fprintf('raw nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_80_alternate(t) = 0; 
      fprintf('raw nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:)+sum(upper_condition_90_success)+sum(lower_condition_90_success))==0
      subset_success_vector_raw_90_alternate(t) = 1; 
      fprintf('raw nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_90_alternate(t) = 0; 
      fprintf('raw nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:)+sum(upper_condition_95_success)+sum(lower_condition_95_success))==0
      subset_success_vector_raw_95_alternate(t) = 1; 
      fprintf('raw nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_95_alternate(t) = 0; 
      fprintf('raw nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:)+sum(upper_condition_80_observed_success)+sum(lower_condition_80_observed_success))==0
      subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('observed nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:)+sum(upper_condition_90_observed_success)+sum(lower_condition_90_observed_success))==0
      subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('observed nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:)+sum(upper_condition_95_observed_success)+sum(lower_condition_95_observed_success))==0
      subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('observed nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('observed nominal 95 alternate true boundary failure! \n');
    end 
                              
end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);

percentage_success_vector_observed_80                    = mean(subset_success_vector_observed_80, 1);
percentage_success_vector_observed_90                    = mean(subset_success_vector_observed_90, 1);
percentage_success_vector_observed_95                    = mean(subset_success_vector_observed_95, 1);

percentage_success_vector_raw_80_alternate               = mean(subset_success_vector_raw_80_alternate, 1);
percentage_success_vector_raw_90_alternate               = mean(subset_success_vector_raw_90_alternate, 1);
percentage_success_vector_raw_95_alternate               = mean(subset_success_vector_raw_95_alternate, 1);

percentage_success_vector_observed_80_alternate          = mean(subset_success_vector_observed_80_alternate, 1);
percentage_success_vector_observed_90_alternate          = mean(subset_success_vector_observed_90_alternate, 1);
percentage_success_vector_observed_95_alternate          = mean(subset_success_vector_observed_95_alternate, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_observed_80_store threshold_observed_90_store threshold_observed_95_store '...
      'lower_contour_raw_80_store lower_contour_raw_90_store lower_contour_raw_95_store lower_contour_observed_80_store lower_contour_observed_90_store lower_contour_observed_95_store '...
      'upper_contour_raw_80_store upper_contour_raw_90_store upper_contour_raw_95_store upper_contour_observed_80_store upper_contour_observed_90_store upper_contour_observed_95_store '...
      'upper_subset_mid_raw_80_store upper_subset_mid_raw_90_store upper_subset_mid_raw_95_store upper_subset_mid_observed_80_store upper_subset_mid_observed_90_store upper_subset_mid_observed_95_store '...
      'mid_subset_lower_raw_80_store mid_subset_lower_raw_90_store mid_subset_lower_raw_95_store mid_subset_lower_observed_80_store mid_subset_lower_observed_90_store mid_subset_lower_observed_95_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_observed_80 subset_success_vector_observed_90 subset_success_vector_observed_95 subset_success_vector_raw_80_alternate subset_success_vector_raw_90_alternate subset_success_vector_raw_95_alternate subset_success_vector_observed_80_alternate subset_success_vector_observed_90_alternate subset_success_vector_observed_95_alternate '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_observed_80 percentage_success_vector_observed_90 percentage_success_vector_observed_95 percentage_success_vector_raw_80_alternate percentage_success_vector_raw_90_alternate percentage_success_vector_raw_95_alternate percentage_success_vector_observed_80_alternate percentage_success_vector_observed_90_alternate percentage_success_vector_observed_95_alternate '...
      'supG_raw_store supG_observed_store '...
      'middle_contour_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_observed_80_volume_prct_store lower_contour_observed_90_volume_prct_store lower_contour_observed_95_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_observed_80_volume_prct_store upper_contour_observed_90_volume_prct_store upper_contour_observed_95_volume_prct_store'])
