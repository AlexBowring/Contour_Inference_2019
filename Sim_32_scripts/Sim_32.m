function Sim_32(nSubj,SvNm,nRlz)
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

nBoot   = 5000;
dim       = [124^2 100]; 
dim_const = [124^2 1];
smo     = 3;
mag     = 0.7;
rimFWHM = 2/sqrt(2*log(2)); 				 
stdblk  = prod(dim([1 2])/2);
stdblk_const = prod(dim_const([1 2])/2);
thr     = 0.7;

%-----------Initialization of Some Variables
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);


wdim_const  = dim_const + 2*ceil(rimFWHM*smo*ones(1,2));
trunc_x_const     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim_const(1))};
trunc_y_const     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim_const(2))};
trnind_const = cat(2, trunc_x_const, trunc_y_const);

observed_data  = zeros([dim nSubj]);
observed_data_const = zeros([dim_const nSubj]);

supG_raw                         = zeros(nBoot, 1);
supG_raw_const                   = zeros(nBoot, 1);
supG_raw_stdize_after            = zeros(nBoot, 1);

% Creating linearly increasing signal across columns
Sig = create_signal(dim, 'ramp', [-20, 20]);
Sig_const = mag*ones(dim_const);
% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar

Sig_boundary_edges   = getBdryparams(Sig, thr);

resid_boundary_values = zeros(Sig_boundary_edges.length,nSubj);
resid_boundary_values_stdize_after = zeros(Sig_boundary_edges.length,nSubj);

for t=1:nRlz
    fprintf('.');
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        Noise = create_noise(wdim, 'homo', 1, smo, trnind);
        tImgs = Sig + Noise; % Creates the true image of smoothed signal + smoothed noise
        observed_data(:,:,i) = tImgs;
        
        Noise_const = create_noise(wdim_const, 'homo', 1, smo, trnind_const);
        tImgs_const = Sig_const + Noise_const;
        observed_data_const(:,:,i) = tImgs_const; 
        
      end %========== Loop i (subjects)

      observed_mean = mean(observed_data,3);
      observed_mean_const = mean(observed_data_const,3);

      observed_std = reshape(...
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim);
      observed_std_const = reshape(...
         biasmystd(reshape(observed_data_const,[prod(dim_const) nSubj]),stdblk_const),...
           dim_const);
       
      observed_cohen_d = observed_mean./observed_std;
      observed_cohen_d_const = observed_mean_const./observed_std_const; 
      observed_cohen_d_std        = sqrt(1+observed_cohen_d.^2/2);
      observed_cohen_d_std        = reshape(observed_cohen_d_std, [prod(dim) 1]);
      observed_cohen_d_std_const  = sqrt(1+observed_cohen_d_const.^2/2); 
      observed_cohen_d_std_const  = reshape(observed_cohen_d_std_const, [prod(dim_const) 1]);
       
      % Residuals
      cohen_resid          = create_resid(observed_data, observed_mean, observed_std, 2);
      cohen_resid_const    = create_resid(observed_data_const, observed_mean_const, observed_std_const, 2);
      
      cohen_resid_std = std(cohen_resid,0,2);
      cohen_resid_std_const = std(cohen_resid_const,0,2);
      standardized_cohen_resid = spdiags(1./observed_cohen_d_std, 0,prod(dim),prod(dim))*cohen_resid;
      standardized_cohen_resid_const = spdiags(1./observed_cohen_d_std_const, 0,prod(dim_const),prod(dim_const))*cohen_resid_const;
      
      for i=1:nSubj
          subject_resid_field                    = reshape(standardized_cohen_resid(:,i), dim);
          subject_resid_field_stdize_after       = reshape(cohen_resid(:,i), dim);
          resid_boundary_values(:,i)             = getBdryvalues(subject_resid_field, Sig_boundary_edges);
          resid_boundary_values_stdize_after(:,i) = getBdryvalues(subject_resid_field_stdize_after, Sig_boundary_edges); 
      end
      
      resid_boundary_values_stdize_after = (1/sqrt(1 + (thr^2/2))).*resid_boundary_values_stdize_after;
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          
          % Applying the bootstrap using Gaussian variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          cohen_resid_bootstrap                  = resid_boundary_values*spdiags(signflips, 0, nSubj, nSubj);
          cohen_resid_bootstrap_const            = standardized_cohen_resid_const*spdiags(signflips, 0, nSubj, nSubj);
          cohen_resid_bootstrap_stdize_after     = resid_boundary_values_stdize_after*spdiags(signflips, 0, nSubj, nSubj);
          cohen_resid_field                      = sum(cohen_resid_bootstrap, 2)/sqrt(nSubj);
          cohen_resid_field_const                = sum(cohen_resid_bootstrap_const, 2)/sqrt(nSubj);
          cohen_resid_field_stdize_after         = sum(cohen_resid_bootstrap_stdize_after, 2)/sqrt(nSubj);
          boot_std                               = std(cohen_resid_bootstrap, 0, 2); 
          boot_std_const                         = std(cohen_resid_bootstrap_const, 0, 2); 
          boot_std_stdize_after                  = std(cohen_resid_field_stdize_after, 0, 2);
          cohen_resid_field                      = cohen_resid_field./boot_std;
          cohen_resid_field_const                = cohen_resid_field_const./boot_std_const;
          cohen_resid_field_stdize_after         = cohen_resid_field_stdize_after./boot_std_stdize_after;

          supG_raw(k)              = max(abs(cohen_resid_field)); 
          supG_raw_const(k)        = max(abs(cohen_resid_field_const));
          supG_raw_stdize_after(k) = max(abs(cohen_resid_field_stdize_after));
 
      end                              
end