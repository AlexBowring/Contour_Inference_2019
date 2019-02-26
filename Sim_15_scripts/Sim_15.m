function Sim_15(nSubj,SvNm,nRlz)
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
dim     = [100 100 100]; 
mag     = 3;
smo     = 3;
rimFWHM = 2/sqrt(2*log(2)); 				 
thr     = 2;
rad     = 30;
%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo)*ones(1,3);  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trunc_z     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(3))};
trnind      = cat(2, trunc_x, trunc_y, trunc_z);

observed_data  = zeros([prod(dim) nSubj]);

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

% This vector stores the number of violations either way
lower_condition_80_violations = zeros(nRlz, 1);
upper_condition_80_violations = zeros(nRlz, 1);
lower_condition_90_violations = zeros(nRlz, 1);
upper_condition_90_violations = zeros(nRlz, 1);
lower_condition_95_violations = zeros(nRlz, 1);
upper_condition_95_violations = zeros(nRlz, 1);

lower_condition_80_violations_observed = zeros(nRlz, 1);
upper_condition_80_violations_observed = zeros(nRlz, 1);
lower_condition_90_violations_observed = zeros(nRlz, 1);
upper_condition_90_violations_observed = zeros(nRlz, 1);
lower_condition_95_violations_observed = zeros(nRlz, 1);
upper_condition_95_violations_observed = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store                   = zeros(nBoot, nRlz);
supG_observed_store              = zeros(nBoot, nRlz);

supG_raw                         = zeros(nBoot,1);
supG_observed                    = zeros(nBoot,1);

% Creating a sphere of signal
Sig = create_signal(wdim, 'sphere', [mag, rad], smo, trnind);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC                      = Sig >= thr;
middle_contour          = AC;
middle_contour_volume   = sum(middle_contour(:));

Sig_boundary_edges   = getBdryparams(Sig, thr);
n_Sig_boundary_edges = size(getBdryvalues(Sig, Sig_boundary_edges),1);

monte_carlo_cohen_d_boundary_values = zeros(n_Sig_boundary_edges, nRlz);

for t=1:nRlz
    fprintf('.');
    observed_mean = zeros(prod(dim),1);
    observed_std  = zeros(prod(dim),1);
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        Noise = create_noise(wdim, 'homo', 1, smo, trnind);
        tImgs = Sig + Noise; % Creates the true image of smoothed signal + smoothed noise        
        tImgs = reshape(tImgs, [prod(dim), 1]);
        
        observed_data(:,i) = tImgs; 
        observed_mean = observed_mean + tImgs;
        observed_std  = observed_std + tImgs.^2;
        
      end %========== Loop i (subjects)
      
      observed_mean = observed_mean/nSubj;

      observed_std = sqrt(observed_std/nSubj - observed_mean.^2);
      
      observed_cohen_d = observed_mean./observed_std;
      
      monte_carlo_cohen_d_boundary_values(:,t) = getBdryvalues(observed_cohen_d, Sig_boundary_edges);
      
      observed_cohen_d_std        = sqrt(1+observed_cohen_d.^2/2); 
      
      observed_boundary_edges   = getBdryparams(reshape(observed_cohen_d, dim), thr);
      n_observed_boundary_edges = size(getBdryvalues(Sig, observed_boundary_edges),1);
      
      resid_boundary_values = zeros([n_Sig_boundary_edges nSubj]);
      observed_resid_boundary_values = zeros([n_observed_boundary_edges nSubj]);
      
      for i=1:nSubj
          cohen_resid = create_resid(observed_data(:,i), observed_mean, observed_std, 2); 
          standardized_cohen_resid = spdiags(1./observed_cohen_d_std, 0,prod(dim),prod(dim))*cohen_resid;
          subject_resid_field                 = reshape(standardized_cohen_resid, [dim 1]);
          resid_boundary_values(:,i)          = getBdryvalues(subject_resid_field, Sig_boundary_edges);
          observed_resid_boundary_values(:,i) = getBdryvalues(subject_resid_field, observed_boundary_edges);
      end
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;

          % True boundary
          boundary_bootstrap                = resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          boundary_resid_field              = sum(boundary_bootstrap, 2)/sqrt(nSubj);
          % Re-standardizing by bootstrap standard deviation
          boot_std                          = std(boundary_bootstrap, 0, 2);
          boundary_resid_field              = boundary_resid_field./boot_std;
          supG_raw(k)                       = max(abs(boundary_resid_field));
          supG_max_raw(k)                 = max(boundary_resid_field);
          supG_min_raw(k)                 = min(boundary_resid_field);
          
          % Estimated boundary
          observed_boundary_bootstrap       = observed_resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          observed_boundary_resid_field     = sum(observed_boundary_bootstrap, 2)/sqrt(nSubj); 
          % Re-standardizing by bootstrap standard deviation
          observed_boot_std                 = std(observed_boundary_bootstrap, 0, 2);
          observed_boundary_resid_field     = observed_boundary_resid_field./observed_boot_std;
          supG_observed(k)                  = max(abs(observed_boundary_resid_field));
          supG_max_observed(k)            = max(observed_boundary_resid_field);
          supG_min_observed(k)            = min(observed_boundary_resid_field);
      end 
    
    observed_cohen_d     = reshape(observed_cohen_d, dim);
    observed_cohen_d_std = reshape(observed_cohen_d_std, dim);
      
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
    
    supGa_max_raw_80                     = prctile(supG_max_raw, 90);
    supGa_max_raw_90                     = prctile(supG_max_raw, 95);
    supGa_max_raw_95                     = prctile(supG_max_raw, 97.5);
       
    supGa_min_raw_80                     = prctile(supG_min_raw, 10);
    supGa_min_raw_90                     = prctile(supG_min_raw, 5);
    supGa_min_raw_95                     = prctile(supG_min_raw, 2.5);   
       
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
    
    max_min_lower_contour_raw_80             = observed_cohen_d >= thr + supGa_min_raw_80*tau*observed_cohen_d_std;
    max_min_upper_contour_raw_80             = observed_cohen_d >= thr + supGa_max_raw_80*tau*observed_cohen_d_std;
    max_min_lower_contour_raw_80_volume_prct = sum(max_min_lower_contour_raw_80(:))/middle_contour_volume;
    max_min_upper_contour_raw_80_volume_prct = sum(max_min_upper_contour_raw_80(:))/middle_contour_volume;
    max_min_mid_on_upper_raw_80              = max_min_upper_contour_raw_80.*middle_contour;
    max_min_lower_on_mid_raw_80              = middle_contour.*max_min_lower_contour_raw_80;
    max_min_upper_subset_mid_raw_80          = max_min_upper_contour_raw_80 - max_min_mid_on_upper_raw_80;
    max_min_mid_subset_lower_raw_80          = middle_contour - max_min_lower_on_mid_raw_80;
    
    max_min_lower_contour_raw_90             = observed_cohen_d >= thr + supGa_min_raw_90*tau*observed_cohen_d_std;
    max_min_upper_contour_raw_90             = observed_cohen_d >= thr + supGa_max_raw_90*tau*observed_cohen_d_std;
    max_min_lower_contour_raw_90_volume_prct = sum(max_min_lower_contour_raw_90(:))/middle_contour_volume;
    max_min_upper_contour_raw_90_volume_prct = sum(max_min_upper_contour_raw_90(:))/middle_contour_volume;
    max_min_mid_on_upper_raw_90              = max_min_upper_contour_raw_90.*middle_contour;
    max_min_lower_on_mid_raw_90              = middle_contour.*max_min_lower_contour_raw_90;
    max_min_upper_subset_mid_raw_90          = max_min_upper_contour_raw_90 - max_min_mid_on_upper_raw_90;
    max_min_mid_subset_lower_raw_90          = middle_contour - max_min_lower_on_mid_raw_90;
    
    max_min_lower_contour_raw_95             = observed_cohen_d >= thr + supGa_min_raw_95*tau*observed_cohen_d_std;
    max_min_upper_contour_raw_95             = observed_cohen_d >= thr + supGa_max_raw_95*tau*observed_cohen_d_std;
    max_min_lower_contour_raw_95_volume_prct = sum(max_min_lower_contour_raw_95(:))/middle_contour_volume;
    max_min_upper_contour_raw_95_volume_prct = sum(max_min_upper_contour_raw_95(:))/middle_contour_volume;
    max_min_mid_on_upper_raw_95              = max_min_upper_contour_raw_95.*middle_contour;
    max_min_lower_on_mid_raw_95              = middle_contour.*max_min_lower_contour_raw_95;
    max_min_upper_subset_mid_raw_95          = max_min_upper_contour_raw_95 - max_min_mid_on_upper_raw_95;
    max_min_mid_subset_lower_raw_95          = middle_contour - max_min_lower_on_mid_raw_95;    
    
    % Observed boundary
    supGa_observed_80                     = prctile(supG_observed, 80);
    supGa_observed_90                     = prctile(supG_observed, 90);
    supGa_observed_95                     = prctile(supG_observed, 95);
    
    supGa_max_observed_80                     = prctile(supG_max_observed, 90);
    supGa_max_observed_90                     = prctile(supG_max_observed, 95);
    supGa_max_observed_95                     = prctile(supG_max_observed, 97.5);
    
    supGa_min_observed_80                     = prctile(supG_min_observed, 10);
    supGa_min_observed_90                     = prctile(supG_min_observed, 5);
    supGa_min_observed_95                     = prctile(supG_min_observed, 2.5);
       
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
    
    max_min_lower_contour_observed_80             = observed_cohen_d >= thr + supGa_min_observed_80*tau*observed_cohen_d_std;
    max_min_upper_contour_observed_80             = observed_cohen_d >= thr + supGa_max_observed_80*tau*observed_cohen_d_std;
    max_min_lower_contour_observed_80_volume_prct = sum(max_min_lower_contour_observed_80(:))/middle_contour_volume;
    max_min_upper_contour_observed_80_volume_prct = sum(max_min_upper_contour_observed_80(:))/middle_contour_volume;
    max_min_mid_on_upper_observed_80              = max_min_upper_contour_observed_80.*middle_contour;
    max_min_lower_on_mid_observed_80              = middle_contour.*max_min_lower_contour_observed_80;
    max_min_upper_subset_mid_observed_80          = max_min_upper_contour_observed_80 - max_min_mid_on_upper_observed_80;
    max_min_mid_subset_lower_observed_80          = middle_contour - max_min_lower_on_mid_observed_80;
    
    max_min_lower_contour_observed_90             = observed_cohen_d >= thr + supGa_min_observed_90*tau*observed_cohen_d_std;
    max_min_upper_contour_observed_90             = observed_cohen_d >= thr + supGa_max_observed_90*tau*observed_cohen_d_std;
    max_min_lower_contour_observed_90_volume_prct = sum(max_min_lower_contour_observed_90(:))/middle_contour_volume;
    max_min_upper_contour_observed_90_volume_prct = sum(max_min_upper_contour_observed_90(:))/middle_contour_volume;
    max_min_mid_on_upper_observed_90              = max_min_upper_contour_observed_90.*middle_contour;
    max_min_lower_on_mid_observed_90              = middle_contour.*max_min_lower_contour_observed_90;
    max_min_upper_subset_mid_observed_90          = max_min_upper_contour_observed_90 - max_min_mid_on_upper_observed_90;
    max_min_mid_subset_lower_observed_90          = middle_contour - max_min_lower_on_mid_observed_90;
    
    max_min_lower_contour_observed_95             = observed_cohen_d >= thr + supGa_min_observed_95*tau*observed_cohen_d_std;
    max_min_upper_contour_observed_95             = observed_cohen_d >= thr + supGa_max_observed_95*tau*observed_cohen_d_std;
    max_min_lower_contour_observed_95_volume_prct = sum(max_min_lower_contour_observed_95(:))/middle_contour_volume;
    max_min_upper_contour_observed_95_volume_prct = sum(max_min_upper_contour_observed_95(:))/middle_contour_volume;
    max_min_mid_on_upper_observed_95              = max_min_upper_contour_observed_95.*middle_contour;
    max_min_lower_on_mid_observed_95              = middle_contour.*max_min_lower_contour_observed_95;
    max_min_upper_subset_mid_observed_95          = max_min_upper_contour_observed_95 - max_min_mid_on_upper_observed_95;
    max_min_mid_subset_lower_observed_95          = middle_contour - max_min_lower_on_mid_observed_95;

    %
    % Storing all variables of interest
    %
    % True boundary variables
    supG_raw_store(:,t)                                    = supG_raw;
    supG_max_raw_store(:,t)                                = supG_max_raw;
    supG_min_raw_store(:,t)                                = supG_min_raw;
    threshold_raw_80_store(t)                              = supGa_raw_80;
    threshold_max_raw_80_store(t)                          = supGa_max_raw_80;
    threshold_min_raw_80_store(t)                          = supGa_min_raw_80;
    lower_contour_raw_80_volume_prct_store(t)              = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)              = upper_contour_raw_80_volume_prct;
    max_min_lower_contour_raw_80_volume_prct_store(t)      = max_min_lower_contour_raw_80_volume_prct;
    max_min_upper_contour_raw_80_volume_prct_store(t)      = max_min_upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                              = supGa_raw_90;
    threshold_max_raw_90_store(t)                          = supGa_max_raw_90;
    threshold_min_raw_90_store(t)                          = supGa_min_raw_90;
    lower_contour_raw_90_volume_prct_store(t)              = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)              = upper_contour_raw_90_volume_prct;
    max_min_lower_contour_raw_90_volume_prct_store(t)      = max_min_lower_contour_raw_90_volume_prct;
    max_min_upper_contour_raw_90_volume_prct_store(t)      = max_min_upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                              = supGa_raw_95;
    threshold_max_raw_95_store(t)                          = supGa_max_raw_95;
    threshold_min_raw_95_store(t)                          = supGa_min_raw_95;
    lower_contour_raw_95_volume_prct_store(t)              = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)              = upper_contour_raw_95_volume_prct;
    max_min_lower_contour_raw_95_volume_prct_store(t)      = max_min_lower_contour_raw_95_volume_prct;
    max_min_upper_contour_raw_95_volume_prct_store(t)      = max_min_upper_contour_raw_95_volume_prct;
    
    % Observed boundary variables
    supG_observed_store(:,t)                                    = supG_observed;
    supG_max_observed_store(:,t)                                = supG_max_observed;
    supG_min_observed_store(:,t)                                = supG_min_observed;
    threshold_observed_80_store(t)                              = supGa_observed_80;
    threshold_max_observed_80_store(t)                          = supGa_max_observed_80;
    threshold_min_observed_80_store(t)                          = supGa_min_observed_80;
    lower_contour_observed_80_volume_prct_store(t)              = lower_contour_observed_80_volume_prct;
    upper_contour_observed_80_volume_prct_store(t)              = upper_contour_observed_80_volume_prct;
    max_min_lower_contour_observed_80_volume_prct_store(t)      = max_min_lower_contour_observed_80_volume_prct;
    max_min_upper_contour_observed_80_volume_prct_store(t)      = max_min_upper_contour_observed_80_volume_prct;
 
    threshold_observed_90_store(t)                              = supGa_observed_90;
    threshold_max_observed_90_store(t)                          = supGa_max_observed_90;
    threshold_min_observed_90_store(t)                          = supGa_min_observed_90;
    lower_contour_observed_90_volume_prct_store(t)              = lower_contour_observed_90_volume_prct;
    upper_contour_observed_90_volume_prct_store(t)              = upper_contour_observed_90_volume_prct;
    max_min_lower_contour_observed_90_volume_prct_store(t)      = max_min_lower_contour_observed_90_volume_prct;
    max_min_upper_contour_observed_90_volume_prct_store(t)      = max_min_upper_contour_observed_90_volume_prct;

    threshold_observed_95_store(t)                              = supGa_observed_95;
    threshold_max_observed_95_store(t)                          = supGa_max_observed_95;
    threshold_min_observed_95_store(t)                          = supGa_min_observed_95;
    lower_contour_observed_95_volume_prct_store(t)              = lower_contour_observed_95_volume_prct;
    upper_contour_observed_95_volume_prct_store(t)              = upper_contour_observed_95_volume_prct;
    max_min_lower_contour_observed_95_volume_prct_store(t)      = max_min_lower_contour_observed_95_volume_prct;
    max_min_upper_contour_observed_95_volume_prct_store(t)      = max_min_upper_contour_observed_95_volume_prct;
    
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
    lower_condition_80_violations(t) = sum(lower_condition_80_success);
    upper_condition_80_violations(t) = sum(upper_condition_80_success);
    
    lower_condition_90_success = observed_cohen_d_true_boundary_values < lower_condition_90_boundary_values;
    upper_condition_90_success = observed_cohen_d_true_boundary_values >= upper_condition_90_boundary_values;
    lower_condition_90_violations(t) = sum(lower_condition_90_success);
    upper_condition_90_violations(t) = sum(upper_condition_90_success);
    
    lower_condition_95_success = observed_cohen_d_true_boundary_values < lower_condition_95_boundary_values;
    upper_condition_95_success = observed_cohen_d_true_boundary_values >= upper_condition_95_boundary_values;
    lower_condition_95_violations(t) = sum(lower_condition_95_success);
    upper_condition_95_violations(t) = sum(upper_condition_95_success);
    
    % Calculating the subset condition when we obtain max and min
    % distributions in the bootstrap
    max_min_lower_condition_80 = thr + supGa_min_raw_80*tau*observed_cohen_d_std;
    max_min_upper_condition_80 = thr + supGa_max_raw_80*tau*observed_cohen_d_std;
    max_min_lower_condition_90 = thr + supGa_min_raw_90*tau*observed_cohen_d_std;
    max_min_upper_condition_90 = thr + supGa_max_raw_90*tau*observed_cohen_d_std;
    max_min_lower_condition_95 = thr + supGa_min_raw_95*tau*observed_cohen_d_std;
    max_min_upper_condition_95 = thr + supGa_max_raw_95*tau*observed_cohen_d_std;
    
    max_min_lower_condition_80_boundary_values = getBdryvalues(lower_condition_80, Sig_boundary_edges);
    max_min_upper_condition_80_boundary_values = getBdryvalues(upper_condition_80, Sig_boundary_edges);

    max_min_lower_condition_90_boundary_values = getBdryvalues(lower_condition_90, Sig_boundary_edges);
    max_min_upper_condition_90_boundary_values = getBdryvalues(upper_condition_90, Sig_boundary_edges);
    
    max_min_lower_condition_95_boundary_values = getBdryvalues(lower_condition_95, Sig_boundary_edges);
    max_min_upper_condition_95_boundary_values = getBdryvalues(upper_condition_95, Sig_boundary_edges);
    
    max_min_lower_condition_80_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_80_boundary_values;
    max_min_upper_condition_80_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_80_boundary_values;
    max_min_lower_condition_80_violations(t) = sum(max_min_lower_condition_80_success);
    max_min_upper_condition_80_violations(t) = sum(max_min_upper_condition_80_success);
    
    max_min_lower_condition_90_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_90_boundary_values;
    max_min_upper_condition_90_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_90_boundary_values;
    max_min_lower_condition_90_violations(t) = sum(max_min_lower_condition_90_success);
    max_min_upper_condition_90_violations(t) = sum(max_min_upper_condition_90_success);
    
    max_min_lower_condition_95_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_95_boundary_values;
    max_min_upper_condition_95_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_95_boundary_values;
    max_min_lower_condition_95_violations(t) = sum(max_min_lower_condition_95_success);
    max_min_upper_condition_95_violations(t) = sum(max_min_upper_condition_95_success);
                              
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
    lower_condition_80_violations_observed(t) = sum(lower_condition_80_observed_success);
    upper_condition_80_violations_observed(t) = sum(upper_condition_80_observed_success);
    
    lower_condition_90_observed_success = observed_cohen_d_true_boundary_values < lower_condition_90_observed_boundary_values;
    upper_condition_90_observed_success = observed_cohen_d_true_boundary_values >= upper_condition_90_observed_boundary_values;
    lower_condition_90_violations_observed(t) = sum(lower_condition_90_observed_success);
    upper_condition_90_violations_observed(t) = sum(upper_condition_90_observed_success);
    
    lower_condition_95_observed_success = observed_cohen_d_true_boundary_values < lower_condition_95_observed_boundary_values;
    upper_condition_95_observed_success = observed_cohen_d_true_boundary_values >= upper_condition_95_observed_boundary_values;
    lower_condition_95_violations_observed(t) = sum(lower_condition_95_observed_success);
    upper_condition_95_violations_observed(t) = sum(upper_condition_95_observed_success);
    
    % Calculating the subset condition when we obtain max and min
    % distributions in the bootstrap on the observed boundary
    max_min_lower_condition_80_observed = thr + supGa_min_observed_80*tau*observed_cohen_d_std;
    max_min_upper_condition_80_observed = thr + supGa_max_observed_80*tau*observed_cohen_d_std;
    max_min_lower_condition_90_observed = thr + supGa_min_observed_90*tau*observed_cohen_d_std;
    max_min_upper_condition_90_observed = thr + supGa_max_observed_90*tau*observed_cohen_d_std;
    max_min_lower_condition_95_observed = thr + supGa_min_observed_80*tau*observed_cohen_d_std;
    max_min_upper_condition_95_observed = thr + supGa_max_observed_80*tau*observed_cohen_d_std;
    
    max_min_lower_condition_80_observed_boundary_values = getBdryvalues(max_min_lower_condition_80_observed, Sig_boundary_edges);
    max_min_upper_condition_80_observed_boundary_values = getBdryvalues(max_min_upper_condition_80_observed, Sig_boundary_edges);
    
    max_min_lower_condition_90_observed_boundary_values = getBdryvalues(max_min_lower_condition_90_observed, Sig_boundary_edges);
    max_min_upper_condition_90_observed_boundary_values = getBdryvalues(max_min_upper_condition_90_observed, Sig_boundary_edges);
    
    max_min_lower_condition_95_observed_boundary_values = getBdryvalues(max_min_lower_condition_95_observed, Sig_boundary_edges);
    max_min_upper_condition_95_observed_boundary_values = getBdryvalues(max_min_upper_condition_95_observed, Sig_boundary_edges);
    
    max_min_lower_condition_80_observed_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_80_observed_boundary_values;
    max_min_upper_condition_80_observed_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_80_observed_boundary_values;
    max_min_lower_condition_80_violations_observed(t) = sum(max_min_lower_condition_80_observed_success);
    max_min_upper_condition_80_violations_observed(t) = sum(max_min_upper_condition_80_observed_success);
    
    max_min_lower_condition_90_observed_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_90_observed_boundary_values;
    max_min_upper_condition_90_observed_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_90_observed_boundary_values;
    max_min_lower_condition_90_violations_observed(t) = sum(max_min_lower_condition_90_observed_success);
    max_min_upper_condition_90_violations_observed(t) = sum(max_min_upper_condition_90_observed_success);
    
    max_min_lower_condition_95_observed_success = observed_cohen_d_true_boundary_values < max_min_lower_condition_95_observed_boundary_values;
    max_min_upper_condition_95_observed_success = observed_cohen_d_true_boundary_values >= max_min_upper_condition_95_observed_boundary_values;
    lower_condition_95_violations_observed(t) = sum(max_min_lower_condition_95_observed_success);
    upper_condition_95_violations_observed(t) = sum(max_min_upper_condition_95_observed_success);
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the true boundary
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))+sum(upper_condition_80_success)+sum(lower_condition_80_success)==0
      subset_success_vector_raw_80_alternate(t) = 1;
      fprintf('raw nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_80_alternate(t) = 0; 
      fprintf('raw nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))+sum(upper_condition_90_success)+sum(lower_condition_90_success)==0
      subset_success_vector_raw_90_alternate(t) = 1; 
      fprintf('raw nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_90_alternate(t) = 0; 
      fprintf('raw nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))+sum(upper_condition_95_success)+sum(lower_condition_95_success)==0
      subset_success_vector_raw_95_alternate(t) = 1; 
      fprintf('raw nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_raw_95_alternate(t) = 0; 
      fprintf('raw nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the true boundary with max/min distributions
    if sum(max_min_upper_subset_mid_raw_80(:))+sum(max_min_mid_subset_lower_raw_80(:))+sum(max_min_upper_condition_80_success)+sum(max_min_lower_condition_80_success)==0
      max_min_subset_success_vector_raw_80_alternate(t) = 1;
      fprintf('max min raw nominal 80 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_raw_80_alternate(t) = 0; 
      fprintf('max min raw nominal 80 alternate true boundary failure! \n');
    end 

    if sum(max_min_upper_subset_mid_raw_90(:))+sum(max_min_mid_subset_lower_raw_90(:))+sum(max_min_upper_condition_90_success)+sum(max_min_lower_condition_90_success)==0
      max_min_subset_success_vector_raw_90_alternate(t) = 1; 
      fprintf('max min raw nominal 90 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_raw_90_alternate(t) = 0; 
      fprintf('max min raw nominal 90 alternate true boundary failure! \n');
    end 

    if sum(max_min_upper_subset_mid_raw_95(:))+sum(max_min_mid_subset_lower_raw_95(:))+sum(max_min_upper_condition_95_success)+sum(max_min_lower_condition_95_success)==0
      max_min_subset_success_vector_raw_95_alternate(t) = 1; 
      fprintf('max min raw nominal 95 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_raw_95_alternate(t) = 0; 
      fprintf('max min raw nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))+sum(upper_condition_80_observed_success)+sum(lower_condition_80_observed_success)==0
      subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('observed nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))+sum(upper_condition_90_observed_success)+sum(lower_condition_90_observed_success)==0
      subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('observed nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))+sum(upper_condition_95_observed_success)+sum(lower_condition_95_observed_success)==0
      subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('observed nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('observed nominal 95 alternate true boundary failure! \n');
    end 
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary with max/min
    % distributions
    if sum(max_min_upper_subset_mid_observed_80(:))+sum(max_min_mid_subset_lower_observed_80(:))+sum(max_min_upper_condition_80_observed_success)+sum(max_min_lower_condition_80_observed_success)==0
      max_min_subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('max min observed nominal 80 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('max min observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(max_min_upper_subset_mid_observed_90(:))+sum(max_min_mid_subset_lower_observed_90(:))+sum(max_min_upper_condition_90_observed_success)+sum(max_min_lower_condition_90_observed_success)==0
      max_min_subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('max min observed nominal 90 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('max min observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(max_min_upper_subset_mid_observed_95(:))+sum(max_min_mid_subset_lower_observed_95(:))+sum(max_min_upper_condition_95_observed_success)+sum(max_min_lower_condition_95_observed_success)==0
      max_min_subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('max min observed nominal 95 alternate true boundary success! \n');
    else 
      max_min_subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('max min observed nominal 95 alternate true boundary failure! \n');
    end 
                              
end

percentage_success_vector_raw_80_alternate               = mean(subset_success_vector_raw_80_alternate, 1);
percentage_success_vector_raw_90_alternate               = mean(subset_success_vector_raw_90_alternate, 1);
percentage_success_vector_raw_95_alternate               = mean(subset_success_vector_raw_95_alternate, 1);

max_min_percentage_success_vector_raw_80_alternate               = mean(max_min_subset_success_vector_raw_80_alternate, 1);
max_min_percentage_success_vector_raw_90_alternate               = mean(max_min_subset_success_vector_raw_90_alternate, 1);
max_min_percentage_success_vector_raw_95_alternate               = mean(max_min_subset_success_vector_raw_95_alternate, 1);

percentage_success_vector_observed_80_alternate          = mean(subset_success_vector_observed_80_alternate, 1);
percentage_success_vector_observed_90_alternate          = mean(subset_success_vector_observed_90_alternate, 1);
percentage_success_vector_observed_95_alternate          = mean(subset_success_vector_observed_95_alternate, 1);

max_min_percentage_success_vector_observed_80_alternate          = mean(max_min_subset_success_vector_observed_80_alternate, 1);
max_min_percentage_success_vector_observed_90_alternate          = mean(max_min_subset_success_vector_observed_90_alternate, 1);
max_min_percentage_success_vector_observed_95_alternate          = mean(max_min_subset_success_vector_observed_95_alternate, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_observed_80_store threshold_observed_90_store threshold_observed_95_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_observed_80 subset_success_vector_observed_90 subset_success_vector_observed_95 subset_success_vector_raw_80_alternate subset_success_vector_raw_90_alternate subset_success_vector_raw_95_alternate subset_success_vector_observed_80_alternate subset_success_vector_observed_90_alternate subset_success_vector_observed_95_alternate '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_observed_80 percentage_success_vector_observed_90 percentage_success_vector_observed_95 percentage_success_vector_raw_80_alternate percentage_success_vector_raw_90_alternate percentage_success_vector_raw_95_alternate percentage_success_vector_observed_80_alternate percentage_success_vector_observed_90_alternate percentage_success_vector_observed_95_alternate '...
      'supG_raw_store supG_observed_store '...
      'middle_contour_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_observed_80_volume_prct_store lower_contour_observed_90_volume_prct_store lower_contour_observed_95_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_observed_80_volume_prct_store upper_contour_observed_90_volume_prct_store upper_contour_observed_95_volume_prct_store '...
      'monte_carlo_cohen_d_boundary_values ' ...])
