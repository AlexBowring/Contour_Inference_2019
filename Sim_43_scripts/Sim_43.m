function Sim_43(nSubj,SvNm,nRlz)
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
rimFWHM = 2/sqrt(2*log(2)); 				 
stdblk  = prod(dim([1 2])/2);
thr     = 0.8;

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data  = zeros([dim nSubj]);

% Variables for the transformation
a       = sqrt((nSubj-1)/(nSubj -3 ));
b       = sqrt(nSubj)*sqrt((8*nSubj^2 - 17*nSubj + 11)/((4*nSubj-5)^2*(nSubj-3)));
alpha   = b^-1;
beta    = b/a;

transformed_thr = alpha*asinh((1 - (3/(4*nSubj - 5)))^(-1)*thr*beta);
second_deriv_term = (1/(2*nSubj))*b^2*((1 - (3/(4*nSubj - 5)))^(-1)*thr)*(( ((nSubj-1)/(nSubj-3)) + nSubj*(thr^2)*( (8*nSubj^2 - 17*nSubj + 11)/(16*(nSubj-3)*((nSubj - 2)^2)) ))^(-1/2)); 


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
Sig = create_signal(dim, 'ramp', 1);

cohen_d = Sig;

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = cohen_d >= thr;

cohen_d_boundary_edges   = getBdryparams(cohen_d, thr);
  
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
      transformed_observed_cohen_d = alpha*asinh(observed_cohen_d*beta);
      
      observed_residual_std        = sqrt(1+observed_cohen_d.^2*(beta^2));       
      
      observed_boundary_edges   = getBdryparams(reshape(observed_cohen_d, dim), (1 - (3/(4*nSubj - 5)))^(-1)*thr);      
      
      % Residuals
      cohen_resid          = create_resid(observed_data, observed_mean, observed_std, 2);
      standardized_cohen_resid = spdiags(alpha*beta./reshape(observed_residual_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*cohen_resid;
            
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          cohen_resid_bootstrap                  = standardized_cohen_resid*spdiags(signflips, 0, nSubj, nSubj);
          cohen_resid_bootstrap                  = reshape(cohen_resid_bootstrap, [dim nSubj]);
          cohen_resid_field                      = sum(cohen_resid_bootstrap, 3)/sqrt(nSubj);
          % Re-standardizing by bootstrap standard deviation
          boot_std                               = std(cohen_resid_bootstrap, 0, 3);
          cohen_resid_field                      = cohen_resid_field./boot_std;
          
          % Calculating the maximum over the weighted interpolated true boundary edges
          true_boundary_values = getBdryvalues(cohen_resid_field, cohen_d_boundary_edges);
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
       
    lower_contour_raw_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_80*tau;
    upper_contour_raw_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_80*tau;
    lower_contour_raw_80_volume_prct = sum(lower_contour_raw_80(:))/middle_contour_volume;
    upper_contour_raw_80_volume_prct = sum(upper_contour_raw_80(:))/middle_contour_volume;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_90*tau;
    upper_contour_raw_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_90*tau;
    lower_contour_raw_90_volume_prct = sum(lower_contour_raw_90(:))/middle_contour_volume;
    upper_contour_raw_90_volume_prct = sum(upper_contour_raw_90(:))/middle_contour_volume;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;    
    
    lower_contour_raw_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_raw_95*tau;
    upper_contour_raw_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_raw_95*tau;
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
       
    lower_contour_observed_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_80*tau;
    upper_contour_observed_80             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_80*tau;
    lower_contour_observed_80_volume_prct = sum(lower_contour_observed_80(:))/middle_contour_volume;
    upper_contour_observed_80_volume_prct = sum(upper_contour_observed_80(:))/middle_contour_volume;
    mid_on_upper_observed_80              = upper_contour_observed_80.*middle_contour;
    lower_on_mid_observed_80              = middle_contour.*lower_contour_observed_80;
    upper_subset_mid_observed_80          = upper_contour_observed_80 - mid_on_upper_observed_80;
    mid_subset_lower_observed_80          = middle_contour - lower_on_mid_observed_80;
    
    lower_contour_observed_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_90*tau;
    upper_contour_observed_90             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_90*tau;
    lower_contour_observed_90_volume_prct = sum(lower_contour_observed_90(:))/middle_contour_volume;
    upper_contour_observed_90_volume_prct = sum(upper_contour_observed_90(:))/middle_contour_volume;
    mid_on_upper_observed_90              = upper_contour_observed_90.*middle_contour;
    lower_on_mid_observed_90              = middle_contour.*lower_contour_observed_90;
    upper_subset_mid_observed_90          = upper_contour_observed_90 - mid_on_upper_observed_90;
    mid_subset_lower_observed_90          = middle_contour - lower_on_mid_observed_90;    
    
    lower_contour_observed_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term - supGa_observed_95*tau;
    upper_contour_observed_95             = transformed_observed_cohen_d >= transformed_thr - second_deriv_term + supGa_observed_95*tau;
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
    lower_condition_80 = transformed_thr - second_deriv_term - supGa_raw_80*tau;
    upper_condition_80 = transformed_thr - second_deriv_term + supGa_raw_80*tau;
    lower_condition_90 = transformed_thr - second_deriv_term - supGa_raw_90*tau;
    upper_condition_90 = transformed_thr - second_deriv_term + supGa_raw_90*tau;
    lower_condition_95 = transformed_thr - second_deriv_term - supGa_raw_95*tau;
    upper_condition_95 = transformed_thr - second_deriv_term + supGa_raw_95*tau;
    
    transformed_observed_cohen_d_true_boundary_values = getBdryvalues(transformed_observed_cohen_d, cohen_d_boundary_edges);
    
    lower_condition_80_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_80;
    upper_condition_80_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_80;
    
    lower_condition_90_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_90;
    upper_condition_90_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_90;
    
    lower_condition_95_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_95;
    upper_condition_95_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_95;
                              
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the observed boundary
    lower_condition_80_observed = transformed_thr - second_deriv_term - supGa_observed_80*tau;
    upper_condition_80_observed = transformed_thr - second_deriv_term + supGa_observed_80*tau;
    lower_condition_90_observed = transformed_thr - second_deriv_term - supGa_observed_90*tau;
    upper_condition_90_observed = transformed_thr - second_deriv_term + supGa_observed_90*tau;
    lower_condition_95_observed = transformed_thr - second_deriv_term - supGa_observed_95*tau;
    upper_condition_95_observed = transformed_thr - second_deriv_term + supGa_observed_95*tau;
    
    lower_condition_80_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_80_observed;
    upper_condition_80_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_80_observed;
    
    lower_condition_90_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_90_observed;
    upper_condition_90_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_90_observed;
    
    lower_condition_95_observed_success = transformed_observed_cohen_d_true_boundary_values < lower_condition_95_observed;
    upper_condition_95_observed_success = transformed_observed_cohen_d_true_boundary_values >= upper_condition_95_observed;
    
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

eval(['save ' SvNm ' nSubj nRlz dim smo rimFWHM thr nBoot '... 
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
