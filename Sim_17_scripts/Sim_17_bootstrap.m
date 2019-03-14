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
nBoot     = 10000;

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data_120_subj_124_dim = zeros([dim, nSubj]);

idea_one_signflips_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);

idea_one_signflips_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);

idea_f_signflips_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
          
idea_f_signflips_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);

idea_one_signflips_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);

idea_one_signflips_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);

idea_f_signflips_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
          
idea_f_signflips_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);

idea_one_signflips_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);

idea_one_signflips_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_signflips_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_gauss_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_one_mammen_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);

idea_f_signflips_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
          
idea_f_signflips_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_signflips_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_gauss_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
idea_f_mammen_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);

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
      
      observed_cohen_d_std_120_subj_124_dim = reshape(sqrt(1+observed_cohen_d_120_subj_124_dim.^2/2), [prod(dim) 1]);
      observed_cohen_d_std_60_subj_124_dim = reshape(sqrt(1+observed_cohen_d_60_subj_124_dim.^2/2), [prod(dim) 1]);
      observed_cohen_d_std_30_subj_124_dim = reshape(sqrt(1+observed_cohen_d_30_subj_124_dim.^2/2), [prod(dim) 1]);
      
      snr_resid_120_subj_124_dim    = create_resid(observed_data_120_subj_124_dim, observed_mean_120_subj_124_dim, observed_std_120_subj_124_dim, 2);
      snr_resid_60_subj_124_dim     = create_resid(observed_data_60_subj_124_dim, observed_mean_60_subj_124_dim, observed_std_60_subj_124_dim, 2);
      snr_resid_30_subj_124_dim     = create_resid(observed_data_30_subj_124_dim, observed_mean_30_subj_124_dim, observed_std_30_subj_124_dim, 2);
      
      snr_resid_std_120_subj_124_dim = biasmystd(reshape(snr_resid_120_subj_124_dim,[prod(dim) nSubj]),stdblk);
      snr_resid_std_60_subj_124_dim = biasmystd(reshape(snr_resid_60_subj_124_dim,[prod(dim) nSubj/2]),stdblk);
      snr_resid_std_30_subj_124_dim = biasmystd(reshape(snr_resid_30_subj_124_dim,[prod(dim) nSubj/4]),stdblk);
      
      idea_one_residuals_120_subj_124_dim = snr_resid_120_subj_124_dim./observed_cohen_d_std_120_subj_124_dim;
      idea_one_residuals_60_subj_124_dim = snr_resid_60_subj_124_dim./observed_cohen_d_std_60_subj_124_dim;
      idea_one_residuals_30_subj_124_dim = snr_resid_30_subj_124_dim./observed_cohen_d_std_30_subj_124_dim;
      
      idea_f_residuals_120_subj_124_dim = snr_resid_120_subj_124_dim./snr_resid_std_120_subj_124_dim;
      idea_f_residuals_60_subj_124_dim = snr_resid_60_subj_124_dim./snr_resid_std_60_subj_124_dim;
      idea_f_residuals_30_subj_124_dim = snr_resid_30_subj_124_dim./snr_resid_std_30_subj_124_dim;
     
      mammen_120_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_120_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_60_subj                               = signflips_120_subj(1:nSubj/2);
          signflips_30_subj                               = signflips_120_subj(1:nSubj/4);
          
          % Gauss variables 
          gauss_120_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_60_subj                                   = gauss_120_subj(1:nSubj/2);
          gauss_30_subj                                   = gauss_120_subj(1:nSubj/4);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_120_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_120_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_60_subj                                  = mammen_120_subj(1:nSubj/2);
          mammen_30_subj                                  = mammen_120_subj(1:nSubj/4);
          
          %% idea one bootstrap
          idea_one_signflips_bootstrap_120_subj_124_dim = idea_one_residuals_120_subj_124_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_one_signflips_bootstrap_60_subj_124_dim = idea_one_residuals_60_subj_124_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_signflips_bootstrap_30_subj_124_dim = idea_one_residuals_30_subj_124_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_gauss_bootstrap_120_subj_124_dim = idea_one_residuals_120_subj_124_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_one_gauss_bootstrap_60_subj_124_dim = idea_one_residuals_60_subj_124_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_gauss_bootstrap_30_subj_124_dim = idea_one_residuals_30_subj_124_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_mammen_bootstrap_120_subj_124_dim = idea_one_residuals_120_subj_124_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_one_mammen_bootstrap_60_subj_124_dim = idea_one_residuals_60_subj_124_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_mammen_bootstrap_30_subj_124_dim = idea_one_residuals_30_subj_124_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_one_signflips_resid_field_120_subj_124_dim = sum(idea_one_signflips_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_one_signflips_resid_field_60_subj_124_dim = sum(idea_one_signflips_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_one_signflips_resid_field_30_subj_124_dim = sum(idea_one_signflips_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          idea_one_gauss_resid_field_120_subj_124_dim = sum(idea_one_gauss_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_one_gauss_resid_field_60_subj_124_dim = sum(idea_one_gauss_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_one_gauss_resid_field_30_subj_124_dim = sum(idea_one_gauss_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          idea_one_mammen_resid_field_120_subj_124_dim = sum(idea_one_mammen_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_one_mammen_resid_field_60_subj_124_dim = sum(idea_one_mammen_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_one_mammen_resid_field_30_subj_124_dim = sum(idea_one_mammen_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          
          idea_one_signflips_boot_std_120_subj_124_dim = std(idea_one_signflips_resid_field_120_subj_124_dim, 0, 2);
          idea_one_signflips_boot_std_60_subj_124_dim = std(idea_one_signflips_resid_field_60_subj_124_dim, 0, 2);
          idea_one_signflips_boot_std_30_subj_124_dim = std(idea_one_signflips_resid_field_30_subj_124_dim, 0, 2);
          idea_one_gauss_boot_std_120_subj_124_dim = std(idea_one_gauss_resid_field_120_subj_124_dim, 0, 2);
          idea_one_gauss_boot_std_60_subj_124_dim = std(idea_one_gauss_resid_field_60_subj_124_dim, 0, 2);
          idea_one_gauss_boot_std_30_subj_124_dim = std(idea_one_gauss_resid_field_30_subj_124_dim, 0, 2);
          idea_one_mammen_boot_std_120_subj_124_dim = std(idea_one_mammen_resid_field_120_subj_124_dim, 0, 2);
          idea_one_mammen_boot_std_60_subj_124_dim = std(idea_one_mammen_resid_field_60_subj_124_dim, 0, 2);
          idea_one_mammen_boot_std_30_subj_124_dim = std(idea_one_mammen_resid_field_30_subj_124_dim, 0, 2);
          
          idea_one_signflips_tboot_resid_field_120_subj_124_dim = idea_one_signflips_resid_field_120_subj_124_dim./idea_one_signflips_boot_std_120_subj_124_dim;      
          idea_one_signflips_tboot_resid_field_60_subj_124_dim = idea_one_signflips_resid_field_60_subj_124_dim./idea_one_signflips_boot_std_60_subj_124_dim;  
          idea_one_signflips_tboot_resid_field_30_subj_124_dim = idea_one_signflips_resid_field_30_subj_124_dim./idea_one_signflips_boot_std_30_subj_124_dim;
          idea_one_gauss_tboot_resid_field_120_subj_124_dim = idea_one_gauss_resid_field_120_subj_124_dim./idea_one_gauss_boot_std_120_subj_124_dim;      
          idea_one_gauss_tboot_resid_field_60_subj_124_dim = idea_one_gauss_resid_field_60_subj_124_dim./idea_one_gauss_boot_std_60_subj_124_dim;  
          idea_one_gauss_tboot_resid_field_30_subj_124_dim = idea_one_gauss_resid_field_30_subj_124_dim./idea_one_gauss_boot_std_30_subj_124_dim;
          idea_one_mammen_tboot_resid_field_120_subj_124_dim = idea_one_mammen_resid_field_120_subj_124_dim./idea_one_mammen_boot_std_120_subj_124_dim;      
          idea_one_mammen_tboot_resid_field_60_subj_124_dim = idea_one_mammen_resid_field_60_subj_124_dim./idea_one_mammen_boot_std_60_subj_124_dim;  
          idea_one_mammen_tboot_resid_field_30_subj_124_dim = idea_one_mammen_resid_field_30_subj_124_dim./idea_one_mammen_boot_std_30_subj_124_dim;
          
          idea_one_signflips_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_signflips_resid_field_120_subj_124_dim));
          idea_one_signflips_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_signflips_resid_field_120_subj_124_dim);
          idea_one_signflips_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_signflips_resid_field_120_subj_124_dim);
          idea_one_signflips_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_signflips_resid_field_60_subj_124_dim));
          idea_one_signflips_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_signflips_resid_field_60_subj_124_dim);
          idea_one_signflips_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_signflips_resid_field_60_subj_124_dim);
          idea_one_signflips_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_signflips_resid_field_30_subj_124_dim));
          idea_one_signflips_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_signflips_resid_field_30_subj_124_dim);
          idea_one_signflips_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_signflips_resid_field_30_subj_124_dim);
          idea_one_gauss_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_gauss_resid_field_120_subj_124_dim));
          idea_one_gauss_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_gauss_resid_field_120_subj_124_dim);
          idea_one_gauss_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_gauss_resid_field_120_subj_124_dim);
          idea_one_gauss_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_gauss_resid_field_60_subj_124_dim));
          idea_one_gauss_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_gauss_resid_field_60_subj_124_dim);
          idea_one_gauss_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_gauss_resid_field_60_subj_124_dim);
          idea_one_gauss_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_gauss_resid_field_30_subj_124_dim));
          idea_one_gauss_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_gauss_resid_field_30_subj_124_dim);
          idea_one_gauss_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_gauss_resid_field_30_subj_124_dim);
          idea_one_mammen_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_mammen_resid_field_120_subj_124_dim));
          idea_one_mammen_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_mammen_resid_field_120_subj_124_dim);
          idea_one_mammen_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_mammen_resid_field_120_subj_124_dim);
          idea_one_mammen_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_mammen_resid_field_60_subj_124_dim));
          idea_one_mammen_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_mammen_resid_field_60_subj_124_dim);
          idea_one_mammen_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_mammen_resid_field_60_subj_124_dim);
          idea_one_mammen_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_mammen_resid_field_30_subj_124_dim));
          idea_one_mammen_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_mammen_resid_field_30_subj_124_dim);
          idea_one_mammen_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_mammen_resid_field_30_subj_124_dim);
          
          idea_one_signflips_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_120_subj_124_dim));
          idea_one_signflips_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_signflips_tboot_resid_field_120_subj_124_dim);
          idea_one_signflips_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_signflips_tboot_resid_field_120_subj_124_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_60_subj_124_dim));
          idea_one_signflips_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_signflips_tboot_resid_field_60_subj_124_dim);
          idea_one_signflips_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_signflips_tboot_resid_field_60_subj_124_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_30_subj_124_dim));
          idea_one_signflips_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_signflips_tboot_resid_field_30_subj_124_dim);
          idea_one_signflips_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_signflips_tboot_resid_field_30_subj_124_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_120_subj_124_dim));
          idea_one_gauss_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_gauss_tboot_resid_field_120_subj_124_dim);
          idea_one_gauss_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_gauss_tboot_resid_field_120_subj_124_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_60_subj_124_dim));
          idea_one_gauss_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_gauss_tboot_resid_field_60_subj_124_dim);
          idea_one_gauss_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_gauss_tboot_resid_field_60_subj_124_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_30_subj_124_dim));
          idea_one_gauss_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_gauss_tboot_resid_field_30_subj_124_dim);
          idea_one_gauss_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_gauss_tboot_resid_field_30_subj_124_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_120_subj_124_dim));
          idea_one_mammen_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_one_mammen_tboot_resid_field_120_subj_124_dim);
          idea_one_mammen_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_one_mammen_tboot_resid_field_120_subj_124_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_60_subj_124_dim));
          idea_one_mammen_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_one_mammen_tboot_resid_field_60_subj_124_dim);
          idea_one_mammen_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_one_mammen_tboot_resid_field_60_subj_124_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_30_subj_124_dim));
          idea_one_mammen_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_one_mammen_tboot_resid_field_30_subj_124_dim);
          idea_one_mammen_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_one_mammen_tboot_resid_field_30_subj_124_dim);
          
          %% idea Fabian bootstrap
          idea_f_signflips_bootstrap_120_subj_124_dim = idea_f_residuals_120_subj_124_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_f_signflips_bootstrap_60_subj_124_dim = idea_f_residuals_60_subj_124_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_signflips_bootstrap_30_subj_124_dim = idea_f_residuals_30_subj_124_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_gauss_bootstrap_120_subj_124_dim = idea_f_residuals_120_subj_124_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_f_gauss_bootstrap_60_subj_124_dim = idea_f_residuals_60_subj_124_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_gauss_bootstrap_30_subj_124_dim = idea_f_residuals_30_subj_124_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_mammen_bootstrap_120_subj_124_dim = idea_f_residuals_120_subj_124_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_f_mammen_bootstrap_60_subj_124_dim = idea_f_residuals_60_subj_124_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_mammen_bootstrap_30_subj_124_dim = idea_f_residuals_30_subj_124_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_f_signflips_resid_field_120_subj_124_dim = sum(idea_f_signflips_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_f_signflips_resid_field_60_subj_124_dim = sum(idea_f_signflips_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_f_signflips_resid_field_30_subj_124_dim = sum(idea_f_signflips_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          idea_f_gauss_resid_field_120_subj_124_dim = sum(idea_f_gauss_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_f_gauss_resid_field_60_subj_124_dim = sum(idea_f_gauss_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_f_gauss_resid_field_30_subj_124_dim = sum(idea_f_gauss_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          idea_f_mammen_resid_field_120_subj_124_dim = sum(idea_f_mammen_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj);
          idea_f_mammen_resid_field_60_subj_124_dim = sum(idea_f_mammen_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/2);
          idea_f_mammen_resid_field_30_subj_124_dim = sum(idea_f_mammen_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/4);
          
          idea_f_signflips_boot_std_120_subj_124_dim = std(idea_f_signflips_resid_field_120_subj_124_dim, 0, 2);
          idea_f_signflips_boot_std_60_subj_124_dim = std(idea_f_signflips_resid_field_60_subj_124_dim, 0, 2);
          idea_f_signflips_boot_std_30_subj_124_dim = std(idea_f_signflips_resid_field_30_subj_124_dim, 0, 2);
          idea_f_gauss_boot_std_120_subj_124_dim = std(idea_f_gauss_resid_field_120_subj_124_dim, 0, 2);
          idea_f_gauss_boot_std_60_subj_124_dim = std(idea_f_gauss_resid_field_60_subj_124_dim, 0, 2);
          idea_f_gauss_boot_std_30_subj_124_dim = std(idea_f_gauss_resid_field_30_subj_124_dim, 0, 2);
          idea_f_mammen_boot_std_120_subj_124_dim = std(idea_f_mammen_resid_field_120_subj_124_dim, 0, 2);
          idea_f_mammen_boot_std_60_subj_124_dim = std(idea_f_mammen_resid_field_60_subj_124_dim, 0, 2);
          idea_f_mammen_boot_std_30_subj_124_dim = std(idea_f_mammen_resid_field_30_subj_124_dim, 0, 2);
          
          idea_f_signflips_tboot_resid_field_120_subj_124_dim = idea_f_signflips_resid_field_120_subj_124_dim./idea_f_signflips_boot_std_120_subj_124_dim;      
          idea_f_signflips_tboot_resid_field_60_subj_124_dim = idea_f_signflips_resid_field_60_subj_124_dim./idea_f_signflips_boot_std_60_subj_124_dim;  
          idea_f_signflips_tboot_resid_field_30_subj_124_dim = idea_f_signflips_resid_field_30_subj_124_dim./idea_f_signflips_boot_std_30_subj_124_dim;
          idea_f_gauss_tboot_resid_field_120_subj_124_dim = idea_f_gauss_resid_field_120_subj_124_dim./idea_f_gauss_boot_std_120_subj_124_dim;      
          idea_f_gauss_tboot_resid_field_60_subj_124_dim = idea_f_gauss_resid_field_60_subj_124_dim./idea_f_gauss_boot_std_60_subj_124_dim;  
          idea_f_gauss_tboot_resid_field_30_subj_124_dim = idea_f_gauss_resid_field_30_subj_124_dim./idea_f_gauss_boot_std_30_subj_124_dim;
          idea_f_mammen_tboot_resid_field_120_subj_124_dim = idea_f_mammen_resid_field_120_subj_124_dim./idea_f_mammen_boot_std_120_subj_124_dim;      
          idea_f_mammen_tboot_resid_field_60_subj_124_dim = idea_f_mammen_resid_field_60_subj_124_dim./idea_f_mammen_boot_std_60_subj_124_dim;  
          idea_f_mammen_tboot_resid_field_30_subj_124_dim = idea_f_mammen_resid_field_30_subj_124_dim./idea_f_mammen_boot_std_30_subj_124_dim;
          
          idea_f_signflips_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_signflips_resid_field_120_subj_124_dim));
          idea_f_signflips_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_signflips_resid_field_120_subj_124_dim);
          idea_f_signflips_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_signflips_resid_field_120_subj_124_dim);
          idea_f_signflips_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_signflips_resid_field_60_subj_124_dim));
          idea_f_signflips_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_signflips_resid_field_60_subj_124_dim);
          idea_f_signflips_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_signflips_resid_field_60_subj_124_dim);
          idea_f_signflips_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_signflips_resid_field_30_subj_124_dim));
          idea_f_signflips_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_signflips_resid_field_30_subj_124_dim);
          idea_f_signflips_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_signflips_resid_field_30_subj_124_dim);
          idea_f_gauss_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_gauss_resid_field_120_subj_124_dim));
          idea_f_gauss_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_gauss_resid_field_120_subj_124_dim);
          idea_f_gauss_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_gauss_resid_field_120_subj_124_dim);
          idea_f_gauss_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_gauss_resid_field_60_subj_124_dim));
          idea_f_gauss_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_gauss_resid_field_60_subj_124_dim);
          idea_f_gauss_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_gauss_resid_field_60_subj_124_dim);
          idea_f_gauss_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_gauss_resid_field_30_subj_124_dim));
          idea_f_gauss_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_gauss_resid_field_30_subj_124_dim);
          idea_f_gauss_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_gauss_resid_field_30_subj_124_dim);
          idea_f_mammen_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_mammen_resid_field_120_subj_124_dim));
          idea_f_mammen_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_mammen_resid_field_120_subj_124_dim);
          idea_f_mammen_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_mammen_resid_field_120_subj_124_dim);
          idea_f_mammen_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_mammen_resid_field_60_subj_124_dim));
          idea_f_mammen_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_mammen_resid_field_60_subj_124_dim);
          idea_f_mammen_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_mammen_resid_field_60_subj_124_dim);
          idea_f_mammen_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_mammen_resid_field_30_subj_124_dim));
          idea_f_mammen_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_mammen_resid_field_30_subj_124_dim);
          idea_f_mammen_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_mammen_resid_field_30_subj_124_dim);
          
          idea_f_signflips_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_120_subj_124_dim));
          idea_f_signflips_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_signflips_tboot_resid_field_120_subj_124_dim);
          idea_f_signflips_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_signflips_tboot_resid_field_120_subj_124_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_60_subj_124_dim));
          idea_f_signflips_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_signflips_tboot_resid_field_60_subj_124_dim);
          idea_f_signflips_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_signflips_tboot_resid_field_60_subj_124_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_30_subj_124_dim));
          idea_f_signflips_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_signflips_tboot_resid_field_30_subj_124_dim);
          idea_f_signflips_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_signflips_tboot_resid_field_30_subj_124_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_120_subj_124_dim));
          idea_f_gauss_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_gauss_tboot_resid_field_120_subj_124_dim);
          idea_f_gauss_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_gauss_tboot_resid_field_120_subj_124_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_60_subj_124_dim));
          idea_f_gauss_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_gauss_tboot_resid_field_60_subj_124_dim);
          idea_f_gauss_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_gauss_tboot_resid_field_60_subj_124_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_30_subj_124_dim));
          idea_f_gauss_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_gauss_tboot_resid_field_30_subj_124_dim);
          idea_f_gauss_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_gauss_tboot_resid_field_30_subj_124_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_120_subj_124_dim));
          idea_f_mammen_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(idea_f_mammen_tboot_resid_field_120_subj_124_dim);
          idea_f_mammen_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(idea_f_mammen_tboot_resid_field_120_subj_124_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_60_subj_124_dim));
          idea_f_mammen_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(idea_f_mammen_tboot_resid_field_60_subj_124_dim);
          idea_f_mammen_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(idea_f_mammen_tboot_resid_field_60_subj_124_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_30_subj_124_dim));
          idea_f_mammen_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(idea_f_mammen_tboot_resid_field_30_subj_124_dim);
          idea_f_mammen_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(idea_f_mammen_tboot_resid_field_30_subj_124_dim);
          
      end
      
     
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
      
      observed_cohen_d_std_120_subj_60_dim = reshape(sqrt(1+observed_cohen_d_120_subj_60_dim.^2/2), [prod(dim_60) 1]);
      observed_cohen_d_std_60_subj_60_dim = reshape(sqrt(1+observed_cohen_d_60_subj_60_dim.^2/2), [prod(dim_60) 1]);
      observed_cohen_d_std_30_subj_60_dim = reshape(sqrt(1+observed_cohen_d_30_subj_60_dim.^2/2), [prod(dim_60) 1]);
      
      snr_resid_120_subj_60_dim     = create_resid(observed_data_120_subj_60_dim, observed_mean_120_subj_60_dim, observed_std_120_subj_60_dim, 2);
      snr_resid_60_subj_60_dim     = create_resid(observed_data_60_subj_60_dim, observed_mean_60_subj_60_dim, observed_std_60_subj_60_dim, 2);
      snr_resid_30_subj_60_dim     = create_resid(observed_data_30_subj_60_dim, observed_mean_30_subj_60_dim, observed_std_30_subj_60_dim, 2);
      
      snr_resid_std_120_subj_60_dim = biasmystd(reshape(snr_resid_120_subj_60_dim,[prod(dim_60) nSubj]),stdblk_60);
      snr_resid_std_60_subj_60_dim = biasmystd(reshape(snr_resid_60_subj_60_dim,[prod(dim_60) nSubj/2]),stdblk_60);
      snr_resid_std_30_subj_60_dim = biasmystd(reshape(snr_resid_30_subj_60_dim,[prod(dim_60) nSubj/4]),stdblk_60);
      
      idea_one_residuals_120_subj_60_dim = snr_resid_120_subj_60_dim./observed_cohen_d_std_120_subj_60_dim;
      idea_one_residuals_60_subj_60_dim = snr_resid_60_subj_60_dim./observed_cohen_d_std_60_subj_60_dim;
      idea_one_residuals_30_subj_60_dim = snr_resid_30_subj_60_dim./observed_cohen_d_std_30_subj_60_dim;
      
      idea_f_residuals_120_subj_60_dim = snr_resid_120_subj_60_dim./snr_resid_std_120_subj_60_dim;
      idea_f_residuals_60_subj_60_dim = snr_resid_60_subj_60_dim./snr_resid_std_60_subj_60_dim;
      idea_f_residuals_30_subj_60_dim = snr_resid_30_subj_60_dim./snr_resid_std_30_subj_60_dim;
     
      mammen_120_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_120_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_60_subj                               = signflips_120_subj(1:nSubj/2);
          signflips_30_subj                               = signflips_120_subj(1:nSubj/4);
          
          % Gauss variables 
          gauss_120_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_60_subj                                   = gauss_120_subj(1:nSubj/2);
          gauss_30_subj                                   = gauss_120_subj(1:nSubj/4);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_120_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_120_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_60_subj                                  = mammen_120_subj(1:nSubj/2);
          mammen_30_subj                                  = mammen_120_subj(1:nSubj/4);
          
          %% idea one bootstrap
          idea_one_signflips_bootstrap_120_subj_60_dim = idea_one_residuals_120_subj_60_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_one_signflips_bootstrap_60_subj_60_dim = idea_one_residuals_60_subj_60_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_signflips_bootstrap_30_subj_60_dim = idea_one_residuals_30_subj_60_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_gauss_bootstrap_120_subj_60_dim = idea_one_residuals_120_subj_60_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_one_gauss_bootstrap_60_subj_60_dim = idea_one_residuals_60_subj_60_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_gauss_bootstrap_30_subj_60_dim = idea_one_residuals_30_subj_60_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_mammen_bootstrap_120_subj_60_dim = idea_one_residuals_120_subj_60_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_one_mammen_bootstrap_60_subj_60_dim = idea_one_residuals_60_subj_60_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_mammen_bootstrap_30_subj_60_dim = idea_one_residuals_30_subj_60_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_one_signflips_resid_field_120_subj_60_dim = sum(idea_one_signflips_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_one_signflips_resid_field_60_subj_60_dim = sum(idea_one_signflips_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_one_signflips_resid_field_30_subj_60_dim = sum(idea_one_signflips_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          idea_one_gauss_resid_field_120_subj_60_dim = sum(idea_one_gauss_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_one_gauss_resid_field_60_subj_60_dim = sum(idea_one_gauss_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_one_gauss_resid_field_30_subj_60_dim = sum(idea_one_gauss_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          idea_one_mammen_resid_field_120_subj_60_dim = sum(idea_one_mammen_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_one_mammen_resid_field_60_subj_60_dim = sum(idea_one_mammen_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_one_mammen_resid_field_30_subj_60_dim = sum(idea_one_mammen_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          
          idea_one_signflips_boot_std_120_subj_60_dim = std(idea_one_signflips_resid_field_120_subj_60_dim, 0, 2);
          idea_one_signflips_boot_std_60_subj_60_dim = std(idea_one_signflips_resid_field_60_subj_60_dim, 0, 2);
          idea_one_signflips_boot_std_30_subj_60_dim = std(idea_one_signflips_resid_field_30_subj_60_dim, 0, 2);
          idea_one_gauss_boot_std_120_subj_60_dim = std(idea_one_gauss_resid_field_120_subj_60_dim, 0, 2);
          idea_one_gauss_boot_std_60_subj_60_dim = std(idea_one_gauss_resid_field_60_subj_60_dim, 0, 2);
          idea_one_gauss_boot_std_30_subj_60_dim = std(idea_one_gauss_resid_field_30_subj_60_dim, 0, 2);
          idea_one_mammen_boot_std_120_subj_60_dim = std(idea_one_mammen_resid_field_120_subj_60_dim, 0, 2);
          idea_one_mammen_boot_std_60_subj_60_dim = std(idea_one_mammen_resid_field_60_subj_60_dim, 0, 2);
          idea_one_mammen_boot_std_30_subj_60_dim = std(idea_one_mammen_resid_field_30_subj_60_dim, 0, 2);
          
          idea_one_signflips_tboot_resid_field_120_subj_60_dim = idea_one_signflips_resid_field_120_subj_60_dim./idea_one_signflips_boot_std_120_subj_60_dim;      
          idea_one_signflips_tboot_resid_field_60_subj_60_dim = idea_one_signflips_resid_field_60_subj_60_dim./idea_one_signflips_boot_std_60_subj_60_dim;  
          idea_one_signflips_tboot_resid_field_30_subj_60_dim = idea_one_signflips_resid_field_30_subj_60_dim./idea_one_signflips_boot_std_30_subj_60_dim;
          idea_one_gauss_tboot_resid_field_120_subj_60_dim = idea_one_gauss_resid_field_120_subj_60_dim./idea_one_gauss_boot_std_120_subj_60_dim;      
          idea_one_gauss_tboot_resid_field_60_subj_60_dim = idea_one_gauss_resid_field_60_subj_60_dim./idea_one_gauss_boot_std_60_subj_60_dim;  
          idea_one_gauss_tboot_resid_field_30_subj_60_dim = idea_one_gauss_resid_field_30_subj_60_dim./idea_one_gauss_boot_std_30_subj_60_dim;
          idea_one_mammen_tboot_resid_field_120_subj_60_dim = idea_one_mammen_resid_field_120_subj_60_dim./idea_one_mammen_boot_std_120_subj_60_dim;      
          idea_one_mammen_tboot_resid_field_60_subj_60_dim = idea_one_mammen_resid_field_60_subj_60_dim./idea_one_mammen_boot_std_60_subj_60_dim;  
          idea_one_mammen_tboot_resid_field_30_subj_60_dim = idea_one_mammen_resid_field_30_subj_60_dim./idea_one_mammen_boot_std_30_subj_60_dim;
          
          idea_one_signflips_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_signflips_resid_field_120_subj_60_dim));
          idea_one_signflips_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_signflips_resid_field_120_subj_60_dim);
          idea_one_signflips_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_signflips_resid_field_120_subj_60_dim);
          idea_one_signflips_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_signflips_resid_field_60_subj_60_dim));
          idea_one_signflips_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_signflips_resid_field_60_subj_60_dim);
          idea_one_signflips_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_signflips_resid_field_60_subj_60_dim);
          idea_one_signflips_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_signflips_resid_field_30_subj_60_dim));
          idea_one_signflips_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_signflips_resid_field_30_subj_60_dim);
          idea_one_signflips_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_signflips_resid_field_30_subj_60_dim);
          idea_one_gauss_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_gauss_resid_field_120_subj_60_dim));
          idea_one_gauss_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_gauss_resid_field_120_subj_60_dim);
          idea_one_gauss_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_gauss_resid_field_120_subj_60_dim);
          idea_one_gauss_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_gauss_resid_field_60_subj_60_dim));
          idea_one_gauss_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_gauss_resid_field_60_subj_60_dim);
          idea_one_gauss_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_gauss_resid_field_60_subj_60_dim);
          idea_one_gauss_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_gauss_resid_field_30_subj_60_dim));
          idea_one_gauss_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_gauss_resid_field_30_subj_60_dim);
          idea_one_gauss_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_gauss_resid_field_30_subj_60_dim);
          idea_one_mammen_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_mammen_resid_field_120_subj_60_dim));
          idea_one_mammen_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_mammen_resid_field_120_subj_60_dim);
          idea_one_mammen_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_mammen_resid_field_120_subj_60_dim);
          idea_one_mammen_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_mammen_resid_field_60_subj_60_dim));
          idea_one_mammen_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_mammen_resid_field_60_subj_60_dim);
          idea_one_mammen_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_mammen_resid_field_60_subj_60_dim);
          idea_one_mammen_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_mammen_resid_field_30_subj_60_dim));
          idea_one_mammen_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_mammen_resid_field_30_subj_60_dim);
          idea_one_mammen_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_mammen_resid_field_30_subj_60_dim);
          
          idea_one_signflips_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_120_subj_60_dim));
          idea_one_signflips_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_signflips_tboot_resid_field_120_subj_60_dim);
          idea_one_signflips_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_signflips_tboot_resid_field_120_subj_60_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_60_subj_60_dim));
          idea_one_signflips_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_signflips_tboot_resid_field_60_subj_60_dim);
          idea_one_signflips_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_signflips_tboot_resid_field_60_subj_60_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_30_subj_60_dim));
          idea_one_signflips_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_signflips_tboot_resid_field_30_subj_60_dim);
          idea_one_signflips_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_signflips_tboot_resid_field_30_subj_60_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_120_subj_60_dim));
          idea_one_gauss_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_gauss_tboot_resid_field_120_subj_60_dim);
          idea_one_gauss_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_gauss_tboot_resid_field_120_subj_60_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_60_subj_60_dim));
          idea_one_gauss_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_gauss_tboot_resid_field_60_subj_60_dim);
          idea_one_gauss_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_gauss_tboot_resid_field_60_subj_60_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_30_subj_60_dim));
          idea_one_gauss_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_gauss_tboot_resid_field_30_subj_60_dim);
          idea_one_gauss_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_gauss_tboot_resid_field_30_subj_60_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_120_subj_60_dim));
          idea_one_mammen_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_one_mammen_tboot_resid_field_120_subj_60_dim);
          idea_one_mammen_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_one_mammen_tboot_resid_field_120_subj_60_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_60_subj_60_dim));
          idea_one_mammen_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_one_mammen_tboot_resid_field_60_subj_60_dim);
          idea_one_mammen_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_one_mammen_tboot_resid_field_60_subj_60_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_30_subj_60_dim));
          idea_one_mammen_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_one_mammen_tboot_resid_field_30_subj_60_dim);
          idea_one_mammen_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_one_mammen_tboot_resid_field_30_subj_60_dim);
          
          %% idea Fabian bootstrap
          idea_f_signflips_bootstrap_120_subj_60_dim = idea_f_residuals_120_subj_60_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_f_signflips_bootstrap_60_subj_60_dim = idea_f_residuals_60_subj_60_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_signflips_bootstrap_30_subj_60_dim = idea_f_residuals_30_subj_60_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_gauss_bootstrap_120_subj_60_dim = idea_f_residuals_120_subj_60_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_f_gauss_bootstrap_60_subj_60_dim = idea_f_residuals_60_subj_60_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_gauss_bootstrap_30_subj_60_dim = idea_f_residuals_30_subj_60_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_mammen_bootstrap_120_subj_60_dim = idea_f_residuals_120_subj_60_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_f_mammen_bootstrap_60_subj_60_dim = idea_f_residuals_60_subj_60_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_mammen_bootstrap_30_subj_60_dim = idea_f_residuals_30_subj_60_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_f_signflips_resid_field_120_subj_60_dim = sum(idea_f_signflips_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_f_signflips_resid_field_60_subj_60_dim = sum(idea_f_signflips_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_f_signflips_resid_field_30_subj_60_dim = sum(idea_f_signflips_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          idea_f_gauss_resid_field_120_subj_60_dim = sum(idea_f_gauss_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_f_gauss_resid_field_60_subj_60_dim = sum(idea_f_gauss_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_f_gauss_resid_field_30_subj_60_dim = sum(idea_f_gauss_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          idea_f_mammen_resid_field_120_subj_60_dim = sum(idea_f_mammen_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj);
          idea_f_mammen_resid_field_60_subj_60_dim = sum(idea_f_mammen_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/2);
          idea_f_mammen_resid_field_30_subj_60_dim = sum(idea_f_mammen_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/4);
          
          idea_f_signflips_boot_std_120_subj_60_dim = std(idea_f_signflips_resid_field_120_subj_60_dim, 0, 2);
          idea_f_signflips_boot_std_60_subj_60_dim = std(idea_f_signflips_resid_field_60_subj_60_dim, 0, 2);
          idea_f_signflips_boot_std_30_subj_60_dim = std(idea_f_signflips_resid_field_30_subj_60_dim, 0, 2);
          idea_f_gauss_boot_std_120_subj_60_dim = std(idea_f_gauss_resid_field_120_subj_60_dim, 0, 2);
          idea_f_gauss_boot_std_60_subj_60_dim = std(idea_f_gauss_resid_field_60_subj_60_dim, 0, 2);
          idea_f_gauss_boot_std_30_subj_60_dim = std(idea_f_gauss_resid_field_30_subj_60_dim, 0, 2);
          idea_f_mammen_boot_std_120_subj_60_dim = std(idea_f_mammen_resid_field_120_subj_60_dim, 0, 2);
          idea_f_mammen_boot_std_60_subj_60_dim = std(idea_f_mammen_resid_field_60_subj_60_dim, 0, 2);
          idea_f_mammen_boot_std_30_subj_60_dim = std(idea_f_mammen_resid_field_30_subj_60_dim, 0, 2);
          
          idea_f_signflips_tboot_resid_field_120_subj_60_dim = idea_f_signflips_resid_field_120_subj_60_dim./idea_f_signflips_boot_std_120_subj_60_dim;      
          idea_f_signflips_tboot_resid_field_60_subj_60_dim = idea_f_signflips_resid_field_60_subj_60_dim./idea_f_signflips_boot_std_60_subj_60_dim;  
          idea_f_signflips_tboot_resid_field_30_subj_60_dim = idea_f_signflips_resid_field_30_subj_60_dim./idea_f_signflips_boot_std_30_subj_60_dim;
          idea_f_gauss_tboot_resid_field_120_subj_60_dim = idea_f_gauss_resid_field_120_subj_60_dim./idea_f_gauss_boot_std_120_subj_60_dim;      
          idea_f_gauss_tboot_resid_field_60_subj_60_dim = idea_f_gauss_resid_field_60_subj_60_dim./idea_f_gauss_boot_std_60_subj_60_dim;  
          idea_f_gauss_tboot_resid_field_30_subj_60_dim = idea_f_gauss_resid_field_30_subj_60_dim./idea_f_gauss_boot_std_30_subj_60_dim;
          idea_f_mammen_tboot_resid_field_120_subj_60_dim = idea_f_mammen_resid_field_120_subj_60_dim./idea_f_mammen_boot_std_120_subj_60_dim;      
          idea_f_mammen_tboot_resid_field_60_subj_60_dim = idea_f_mammen_resid_field_60_subj_60_dim./idea_f_mammen_boot_std_60_subj_60_dim;  
          idea_f_mammen_tboot_resid_field_30_subj_60_dim = idea_f_mammen_resid_field_30_subj_60_dim./idea_f_mammen_boot_std_30_subj_60_dim;
          
          idea_f_signflips_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_signflips_resid_field_120_subj_60_dim));
          idea_f_signflips_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_signflips_resid_field_120_subj_60_dim);
          idea_f_signflips_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_signflips_resid_field_120_subj_60_dim);
          idea_f_signflips_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_signflips_resid_field_60_subj_60_dim));
          idea_f_signflips_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_signflips_resid_field_60_subj_60_dim);
          idea_f_signflips_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_signflips_resid_field_60_subj_60_dim);
          idea_f_signflips_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_signflips_resid_field_30_subj_60_dim));
          idea_f_signflips_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_signflips_resid_field_30_subj_60_dim);
          idea_f_signflips_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_signflips_resid_field_30_subj_60_dim);
          idea_f_gauss_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_gauss_resid_field_120_subj_60_dim));
          idea_f_gauss_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_gauss_resid_field_120_subj_60_dim);
          idea_f_gauss_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_gauss_resid_field_120_subj_60_dim);
          idea_f_gauss_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_gauss_resid_field_60_subj_60_dim));
          idea_f_gauss_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_gauss_resid_field_60_subj_60_dim);
          idea_f_gauss_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_gauss_resid_field_60_subj_60_dim);
          idea_f_gauss_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_gauss_resid_field_30_subj_60_dim));
          idea_f_gauss_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_gauss_resid_field_30_subj_60_dim);
          idea_f_gauss_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_gauss_resid_field_30_subj_60_dim);
          idea_f_mammen_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_mammen_resid_field_120_subj_60_dim));
          idea_f_mammen_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_mammen_resid_field_120_subj_60_dim);
          idea_f_mammen_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_mammen_resid_field_120_subj_60_dim);
          idea_f_mammen_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_mammen_resid_field_60_subj_60_dim));
          idea_f_mammen_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_mammen_resid_field_60_subj_60_dim);
          idea_f_mammen_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_mammen_resid_field_60_subj_60_dim);
          idea_f_mammen_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_mammen_resid_field_30_subj_60_dim));
          idea_f_mammen_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_mammen_resid_field_30_subj_60_dim);
          idea_f_mammen_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_mammen_resid_field_30_subj_60_dim);
          
          idea_f_signflips_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_120_subj_60_dim));
          idea_f_signflips_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_signflips_tboot_resid_field_120_subj_60_dim);
          idea_f_signflips_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_signflips_tboot_resid_field_120_subj_60_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_60_subj_60_dim));
          idea_f_signflips_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_signflips_tboot_resid_field_60_subj_60_dim);
          idea_f_signflips_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_signflips_tboot_resid_field_60_subj_60_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_30_subj_60_dim));
          idea_f_signflips_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_signflips_tboot_resid_field_30_subj_60_dim);
          idea_f_signflips_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_signflips_tboot_resid_field_30_subj_60_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_120_subj_60_dim));
          idea_f_gauss_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_gauss_tboot_resid_field_120_subj_60_dim);
          idea_f_gauss_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_gauss_tboot_resid_field_120_subj_60_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_60_subj_60_dim));
          idea_f_gauss_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_gauss_tboot_resid_field_60_subj_60_dim);
          idea_f_gauss_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_gauss_tboot_resid_field_60_subj_60_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_30_subj_60_dim));
          idea_f_gauss_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_gauss_tboot_resid_field_30_subj_60_dim);
          idea_f_gauss_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_gauss_tboot_resid_field_30_subj_60_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_120_subj_60_dim));
          idea_f_mammen_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(idea_f_mammen_tboot_resid_field_120_subj_60_dim);
          idea_f_mammen_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(idea_f_mammen_tboot_resid_field_120_subj_60_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_60_subj_60_dim));
          idea_f_mammen_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(idea_f_mammen_tboot_resid_field_60_subj_60_dim);
          idea_f_mammen_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(idea_f_mammen_tboot_resid_field_60_subj_60_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_30_subj_60_dim));
          idea_f_mammen_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(idea_f_mammen_tboot_resid_field_30_subj_60_dim);
          idea_f_mammen_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(idea_f_mammen_tboot_resid_field_30_subj_60_dim);
          
      end   
      
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
      
      observed_cohen_d_std_120_subj_10_dim = reshape(sqrt(1+observed_cohen_d_120_subj_10_dim.^2/2), [prod(dim_10) 1]);
      observed_cohen_d_std_60_subj_10_dim = reshape(sqrt(1+observed_cohen_d_60_subj_10_dim.^2/2), [prod(dim_10) 1]);
      observed_cohen_d_std_30_subj_10_dim = reshape(sqrt(1+observed_cohen_d_30_subj_10_dim.^2/2), [prod(dim_10) 1]);
      
      snr_resid_120_subj_10_dim     = create_resid(observed_data_120_subj_10_dim, observed_mean_120_subj_10_dim, observed_std_120_subj_10_dim, 2);
      snr_resid_60_subj_10_dim     = create_resid(observed_data_60_subj_10_dim, observed_mean_60_subj_10_dim, observed_std_60_subj_10_dim, 2);
      snr_resid_30_subj_10_dim     = create_resid(observed_data_30_subj_10_dim, observed_mean_30_subj_10_dim, observed_std_30_subj_10_dim, 2);
      
      snr_resid_std_120_subj_10_dim = biasmystd(reshape(snr_resid_120_subj_10_dim,[prod(dim_10) nSubj]),stdblk_10);
      snr_resid_std_60_subj_10_dim = biasmystd(reshape(snr_resid_60_subj_10_dim,[prod(dim_10) nSubj/2]),stdblk_10);
      snr_resid_std_30_subj_10_dim = biasmystd(reshape(snr_resid_30_subj_10_dim,[prod(dim_10) nSubj/4]),stdblk_10);
      
      idea_one_residuals_120_subj_10_dim = snr_resid_120_subj_10_dim./observed_cohen_d_std_120_subj_10_dim;
      idea_one_residuals_60_subj_10_dim = snr_resid_60_subj_10_dim./observed_cohen_d_std_60_subj_10_dim;
      idea_one_residuals_30_subj_10_dim = snr_resid_30_subj_10_dim./observed_cohen_d_std_30_subj_10_dim;
      
      idea_f_residuals_120_subj_10_dim = snr_resid_120_subj_10_dim./snr_resid_std_120_subj_10_dim;
      idea_f_residuals_60_subj_10_dim = snr_resid_60_subj_10_dim./snr_resid_std_60_subj_10_dim;
      idea_f_residuals_30_subj_10_dim = snr_resid_30_subj_10_dim./snr_resid_std_30_subj_10_dim;
     
      mammen_120_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_120_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_60_subj                               = signflips_120_subj(1:nSubj/2);
          signflips_30_subj                               = signflips_120_subj(1:nSubj/4);
          
          % Gauss variables 
          gauss_120_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_60_subj                                   = gauss_120_subj(1:nSubj/2);
          gauss_30_subj                                   = gauss_120_subj(1:nSubj/4);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_120_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_120_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_60_subj                                  = mammen_120_subj(1:nSubj/2);
          mammen_30_subj                                  = mammen_120_subj(1:nSubj/4);
          
          %% idea one bootstrap
          idea_one_signflips_bootstrap_120_subj_10_dim = idea_one_residuals_120_subj_10_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_one_signflips_bootstrap_60_subj_10_dim = idea_one_residuals_60_subj_10_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_signflips_bootstrap_30_subj_10_dim = idea_one_residuals_30_subj_10_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_gauss_bootstrap_120_subj_10_dim = idea_one_residuals_120_subj_10_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_one_gauss_bootstrap_60_subj_10_dim = idea_one_residuals_60_subj_10_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_gauss_bootstrap_30_subj_10_dim = idea_one_residuals_30_subj_10_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_one_mammen_bootstrap_120_subj_10_dim = idea_one_residuals_120_subj_10_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_one_mammen_bootstrap_60_subj_10_dim = idea_one_residuals_60_subj_10_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_one_mammen_bootstrap_30_subj_10_dim = idea_one_residuals_30_subj_10_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_one_signflips_resid_field_120_subj_10_dim = sum(idea_one_signflips_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_one_signflips_resid_field_60_subj_10_dim = sum(idea_one_signflips_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_one_signflips_resid_field_30_subj_10_dim = sum(idea_one_signflips_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          idea_one_gauss_resid_field_120_subj_10_dim = sum(idea_one_gauss_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_one_gauss_resid_field_60_subj_10_dim = sum(idea_one_gauss_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_one_gauss_resid_field_30_subj_10_dim = sum(idea_one_gauss_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          idea_one_mammen_resid_field_120_subj_10_dim = sum(idea_one_mammen_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_one_mammen_resid_field_60_subj_10_dim = sum(idea_one_mammen_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_one_mammen_resid_field_30_subj_10_dim = sum(idea_one_mammen_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          
          idea_one_signflips_boot_std_120_subj_10_dim = std(idea_one_signflips_resid_field_120_subj_10_dim, 0, 2);
          idea_one_signflips_boot_std_60_subj_10_dim = std(idea_one_signflips_resid_field_60_subj_10_dim, 0, 2);
          idea_one_signflips_boot_std_30_subj_10_dim = std(idea_one_signflips_resid_field_30_subj_10_dim, 0, 2);
          idea_one_gauss_boot_std_120_subj_10_dim = std(idea_one_gauss_resid_field_120_subj_10_dim, 0, 2);
          idea_one_gauss_boot_std_60_subj_10_dim = std(idea_one_gauss_resid_field_60_subj_10_dim, 0, 2);
          idea_one_gauss_boot_std_30_subj_10_dim = std(idea_one_gauss_resid_field_30_subj_10_dim, 0, 2);
          idea_one_mammen_boot_std_120_subj_10_dim = std(idea_one_mammen_resid_field_120_subj_10_dim, 0, 2);
          idea_one_mammen_boot_std_60_subj_10_dim = std(idea_one_mammen_resid_field_60_subj_10_dim, 0, 2);
          idea_one_mammen_boot_std_30_subj_10_dim = std(idea_one_mammen_resid_field_30_subj_10_dim, 0, 2);
          
          idea_one_signflips_tboot_resid_field_120_subj_10_dim = idea_one_signflips_resid_field_120_subj_10_dim./idea_one_signflips_boot_std_120_subj_10_dim;      
          idea_one_signflips_tboot_resid_field_60_subj_10_dim = idea_one_signflips_resid_field_60_subj_10_dim./idea_one_signflips_boot_std_60_subj_10_dim;  
          idea_one_signflips_tboot_resid_field_30_subj_10_dim = idea_one_signflips_resid_field_30_subj_10_dim./idea_one_signflips_boot_std_30_subj_10_dim;
          idea_one_gauss_tboot_resid_field_120_subj_10_dim = idea_one_gauss_resid_field_120_subj_10_dim./idea_one_gauss_boot_std_120_subj_10_dim;      
          idea_one_gauss_tboot_resid_field_60_subj_10_dim = idea_one_gauss_resid_field_60_subj_10_dim./idea_one_gauss_boot_std_60_subj_10_dim;  
          idea_one_gauss_tboot_resid_field_30_subj_10_dim = idea_one_gauss_resid_field_30_subj_10_dim./idea_one_gauss_boot_std_30_subj_10_dim;
          idea_one_mammen_tboot_resid_field_120_subj_10_dim = idea_one_mammen_resid_field_120_subj_10_dim./idea_one_mammen_boot_std_120_subj_10_dim;      
          idea_one_mammen_tboot_resid_field_60_subj_10_dim = idea_one_mammen_resid_field_60_subj_10_dim./idea_one_mammen_boot_std_60_subj_10_dim;  
          idea_one_mammen_tboot_resid_field_30_subj_10_dim = idea_one_mammen_resid_field_30_subj_10_dim./idea_one_mammen_boot_std_30_subj_10_dim;
          
          idea_one_signflips_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_signflips_resid_field_120_subj_10_dim));
          idea_one_signflips_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_signflips_resid_field_120_subj_10_dim);
          idea_one_signflips_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_signflips_resid_field_120_subj_10_dim);
          idea_one_signflips_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_signflips_resid_field_60_subj_10_dim));
          idea_one_signflips_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_signflips_resid_field_60_subj_10_dim);
          idea_one_signflips_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_signflips_resid_field_60_subj_10_dim);
          idea_one_signflips_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_signflips_resid_field_30_subj_10_dim));
          idea_one_signflips_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_signflips_resid_field_30_subj_10_dim);
          idea_one_signflips_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_signflips_resid_field_30_subj_10_dim);
          idea_one_gauss_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_gauss_resid_field_120_subj_10_dim));
          idea_one_gauss_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_gauss_resid_field_120_subj_10_dim);
          idea_one_gauss_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_gauss_resid_field_120_subj_10_dim);
          idea_one_gauss_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_gauss_resid_field_60_subj_10_dim));
          idea_one_gauss_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_gauss_resid_field_60_subj_10_dim);
          idea_one_gauss_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_gauss_resid_field_60_subj_10_dim);
          idea_one_gauss_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_gauss_resid_field_30_subj_10_dim));
          idea_one_gauss_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_gauss_resid_field_30_subj_10_dim);
          idea_one_gauss_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_gauss_resid_field_30_subj_10_dim);
          idea_one_mammen_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_mammen_resid_field_120_subj_10_dim));
          idea_one_mammen_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_mammen_resid_field_120_subj_10_dim);
          idea_one_mammen_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_mammen_resid_field_120_subj_10_dim);
          idea_one_mammen_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_mammen_resid_field_60_subj_10_dim));
          idea_one_mammen_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_mammen_resid_field_60_subj_10_dim);
          idea_one_mammen_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_mammen_resid_field_60_subj_10_dim);
          idea_one_mammen_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_mammen_resid_field_30_subj_10_dim));
          idea_one_mammen_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_mammen_resid_field_30_subj_10_dim);
          idea_one_mammen_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_mammen_resid_field_30_subj_10_dim);
          
          idea_one_signflips_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_120_subj_10_dim));
          idea_one_signflips_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_signflips_tboot_resid_field_120_subj_10_dim);
          idea_one_signflips_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_signflips_tboot_resid_field_120_subj_10_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_60_subj_10_dim));
          idea_one_signflips_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_signflips_tboot_resid_field_60_subj_10_dim);
          idea_one_signflips_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_signflips_tboot_resid_field_60_subj_10_dim);
          idea_one_signflips_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_signflips_tboot_resid_field_30_subj_10_dim));
          idea_one_signflips_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_signflips_tboot_resid_field_30_subj_10_dim);
          idea_one_signflips_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_signflips_tboot_resid_field_30_subj_10_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_120_subj_10_dim));
          idea_one_gauss_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_gauss_tboot_resid_field_120_subj_10_dim);
          idea_one_gauss_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_gauss_tboot_resid_field_120_subj_10_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_60_subj_10_dim));
          idea_one_gauss_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_gauss_tboot_resid_field_60_subj_10_dim);
          idea_one_gauss_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_gauss_tboot_resid_field_60_subj_10_dim);
          idea_one_gauss_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_gauss_tboot_resid_field_30_subj_10_dim));
          idea_one_gauss_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_gauss_tboot_resid_field_30_subj_10_dim);
          idea_one_gauss_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_gauss_tboot_resid_field_30_subj_10_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_120_subj_10_dim));
          idea_one_mammen_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_one_mammen_tboot_resid_field_120_subj_10_dim);
          idea_one_mammen_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_one_mammen_tboot_resid_field_120_subj_10_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_60_subj_10_dim));
          idea_one_mammen_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_one_mammen_tboot_resid_field_60_subj_10_dim);
          idea_one_mammen_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_one_mammen_tboot_resid_field_60_subj_10_dim);
          idea_one_mammen_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_one_mammen_tboot_resid_field_30_subj_10_dim));
          idea_one_mammen_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_one_mammen_tboot_resid_field_30_subj_10_dim);
          idea_one_mammen_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_one_mammen_tboot_resid_field_30_subj_10_dim);
          
          %% idea Fabian bootstrap
          idea_f_signflips_bootstrap_120_subj_10_dim = idea_f_residuals_120_subj_10_dim*spdiags(signflips_120_subj, 0, nSubj, nSubj);
          idea_f_signflips_bootstrap_60_subj_10_dim = idea_f_residuals_60_subj_10_dim*spdiags(signflips_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_signflips_bootstrap_30_subj_10_dim = idea_f_residuals_30_subj_10_dim*spdiags(signflips_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_gauss_bootstrap_120_subj_10_dim = idea_f_residuals_120_subj_10_dim*spdiags(gauss_120_subj, 0, nSubj, nSubj);
          idea_f_gauss_bootstrap_60_subj_10_dim = idea_f_residuals_60_subj_10_dim*spdiags(gauss_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_gauss_bootstrap_30_subj_10_dim = idea_f_residuals_30_subj_10_dim*spdiags(gauss_30_subj, 0, nSubj/4, nSubj/4);
          idea_f_mammen_bootstrap_120_subj_10_dim = idea_f_residuals_120_subj_10_dim*spdiags(mammen_120_subj, 0, nSubj, nSubj);
          idea_f_mammen_bootstrap_60_subj_10_dim = idea_f_residuals_60_subj_10_dim*spdiags(mammen_60_subj, 0, nSubj/2, nSubj/2);
          idea_f_mammen_bootstrap_30_subj_10_dim = idea_f_residuals_30_subj_10_dim*spdiags(mammen_30_subj, 0, nSubj/4, nSubj/4);
          
          idea_f_signflips_resid_field_120_subj_10_dim = sum(idea_f_signflips_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_f_signflips_resid_field_60_subj_10_dim = sum(idea_f_signflips_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_f_signflips_resid_field_30_subj_10_dim = sum(idea_f_signflips_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          idea_f_gauss_resid_field_120_subj_10_dim = sum(idea_f_gauss_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_f_gauss_resid_field_60_subj_10_dim = sum(idea_f_gauss_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_f_gauss_resid_field_30_subj_10_dim = sum(idea_f_gauss_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          idea_f_mammen_resid_field_120_subj_10_dim = sum(idea_f_mammen_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj);
          idea_f_mammen_resid_field_60_subj_10_dim = sum(idea_f_mammen_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/2);
          idea_f_mammen_resid_field_30_subj_10_dim = sum(idea_f_mammen_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/4);
          
          idea_f_signflips_boot_std_120_subj_10_dim = std(idea_f_signflips_resid_field_120_subj_10_dim, 0, 2);
          idea_f_signflips_boot_std_60_subj_10_dim = std(idea_f_signflips_resid_field_60_subj_10_dim, 0, 2);
          idea_f_signflips_boot_std_30_subj_10_dim = std(idea_f_signflips_resid_field_30_subj_10_dim, 0, 2);
          idea_f_gauss_boot_std_120_subj_10_dim = std(idea_f_gauss_resid_field_120_subj_10_dim, 0, 2);
          idea_f_gauss_boot_std_60_subj_10_dim = std(idea_f_gauss_resid_field_60_subj_10_dim, 0, 2);
          idea_f_gauss_boot_std_30_subj_10_dim = std(idea_f_gauss_resid_field_30_subj_10_dim, 0, 2);
          idea_f_mammen_boot_std_120_subj_10_dim = std(idea_f_mammen_resid_field_120_subj_10_dim, 0, 2);
          idea_f_mammen_boot_std_60_subj_10_dim = std(idea_f_mammen_resid_field_60_subj_10_dim, 0, 2);
          idea_f_mammen_boot_std_30_subj_10_dim = std(idea_f_mammen_resid_field_30_subj_10_dim, 0, 2);
          
          idea_f_signflips_tboot_resid_field_120_subj_10_dim = idea_f_signflips_resid_field_120_subj_10_dim./idea_f_signflips_boot_std_120_subj_10_dim;      
          idea_f_signflips_tboot_resid_field_60_subj_10_dim = idea_f_signflips_resid_field_60_subj_10_dim./idea_f_signflips_boot_std_60_subj_10_dim;  
          idea_f_signflips_tboot_resid_field_30_subj_10_dim = idea_f_signflips_resid_field_30_subj_10_dim./idea_f_signflips_boot_std_30_subj_10_dim;
          idea_f_gauss_tboot_resid_field_120_subj_10_dim = idea_f_gauss_resid_field_120_subj_10_dim./idea_f_gauss_boot_std_120_subj_10_dim;      
          idea_f_gauss_tboot_resid_field_60_subj_10_dim = idea_f_gauss_resid_field_60_subj_10_dim./idea_f_gauss_boot_std_60_subj_10_dim;  
          idea_f_gauss_tboot_resid_field_30_subj_10_dim = idea_f_gauss_resid_field_30_subj_10_dim./idea_f_gauss_boot_std_30_subj_10_dim;
          idea_f_mammen_tboot_resid_field_120_subj_10_dim = idea_f_mammen_resid_field_120_subj_10_dim./idea_f_mammen_boot_std_120_subj_10_dim;      
          idea_f_mammen_tboot_resid_field_60_subj_10_dim = idea_f_mammen_resid_field_60_subj_10_dim./idea_f_mammen_boot_std_60_subj_10_dim;  
          idea_f_mammen_tboot_resid_field_30_subj_10_dim = idea_f_mammen_resid_field_30_subj_10_dim./idea_f_mammen_boot_std_30_subj_10_dim;
          
          idea_f_signflips_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_signflips_resid_field_120_subj_10_dim));
          idea_f_signflips_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_signflips_resid_field_120_subj_10_dim);
          idea_f_signflips_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_signflips_resid_field_120_subj_10_dim);
          idea_f_signflips_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_signflips_resid_field_60_subj_10_dim));
          idea_f_signflips_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_signflips_resid_field_60_subj_10_dim);
          idea_f_signflips_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_signflips_resid_field_60_subj_10_dim);
          idea_f_signflips_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_signflips_resid_field_30_subj_10_dim));
          idea_f_signflips_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_signflips_resid_field_30_subj_10_dim);
          idea_f_signflips_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_signflips_resid_field_30_subj_10_dim);
          idea_f_gauss_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_gauss_resid_field_120_subj_10_dim));
          idea_f_gauss_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_gauss_resid_field_120_subj_10_dim);
          idea_f_gauss_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_gauss_resid_field_120_subj_10_dim);
          idea_f_gauss_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_gauss_resid_field_60_subj_10_dim));
          idea_f_gauss_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_gauss_resid_field_60_subj_10_dim);
          idea_f_gauss_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_gauss_resid_field_60_subj_10_dim);
          idea_f_gauss_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_gauss_resid_field_30_subj_10_dim));
          idea_f_gauss_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_gauss_resid_field_30_subj_10_dim);
          idea_f_gauss_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_gauss_resid_field_30_subj_10_dim);
          idea_f_mammen_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_mammen_resid_field_120_subj_10_dim));
          idea_f_mammen_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_mammen_resid_field_120_subj_10_dim);
          idea_f_mammen_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_mammen_resid_field_120_subj_10_dim);
          idea_f_mammen_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_mammen_resid_field_60_subj_10_dim));
          idea_f_mammen_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_mammen_resid_field_60_subj_10_dim);
          idea_f_mammen_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_mammen_resid_field_60_subj_10_dim);
          idea_f_mammen_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_mammen_resid_field_30_subj_10_dim));
          idea_f_mammen_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_mammen_resid_field_30_subj_10_dim);
          idea_f_mammen_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_mammen_resid_field_30_subj_10_dim);
          
          idea_f_signflips_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_120_subj_10_dim));
          idea_f_signflips_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_signflips_tboot_resid_field_120_subj_10_dim);
          idea_f_signflips_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_signflips_tboot_resid_field_120_subj_10_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_60_subj_10_dim));
          idea_f_signflips_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_signflips_tboot_resid_field_60_subj_10_dim);
          idea_f_signflips_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_signflips_tboot_resid_field_60_subj_10_dim);
          idea_f_signflips_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_signflips_tboot_resid_field_30_subj_10_dim));
          idea_f_signflips_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_signflips_tboot_resid_field_30_subj_10_dim);
          idea_f_signflips_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_signflips_tboot_resid_field_30_subj_10_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_120_subj_10_dim));
          idea_f_gauss_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_gauss_tboot_resid_field_120_subj_10_dim);
          idea_f_gauss_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_gauss_tboot_resid_field_120_subj_10_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_60_subj_10_dim));
          idea_f_gauss_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_gauss_tboot_resid_field_60_subj_10_dim);
          idea_f_gauss_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_gauss_tboot_resid_field_60_subj_10_dim);
          idea_f_gauss_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_gauss_tboot_resid_field_30_subj_10_dim));
          idea_f_gauss_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_gauss_tboot_resid_field_30_subj_10_dim);
          idea_f_gauss_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_gauss_tboot_resid_field_30_subj_10_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_120_subj_10_dim));
          idea_f_mammen_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(idea_f_mammen_tboot_resid_field_120_subj_10_dim);
          idea_f_mammen_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(idea_f_mammen_tboot_resid_field_120_subj_10_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_60_subj_10_dim));
          idea_f_mammen_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(idea_f_mammen_tboot_resid_field_60_subj_10_dim);
          idea_f_mammen_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(idea_f_mammen_tboot_resid_field_60_subj_10_dim);
          idea_f_mammen_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(idea_f_mammen_tboot_resid_field_30_subj_10_dim));
          idea_f_mammen_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(idea_f_mammen_tboot_resid_field_30_subj_10_dim);
          idea_f_mammen_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(idea_f_mammen_tboot_resid_field_30_subj_10_dim);
          
      end
end

%% Saving variables

eval(['save ' SvNm ' nRlz dim smo mag '... 
      'idea_one_signflips_abs_supG_resid_field_120_subj_124_dim idea_one_signflips_max_supG_resid_field_120_subj_124_dim idea_one_signflips_min_supG_resid_field_120_subj_124_dim idea_one_signflips_abs_supG_resid_field_60_subj_124_dim idea_one_signflips_max_supG_resid_field_60_subj_124_dim idea_one_signflips_min_supG_resid_field_60_subj_124_dim idea_one_signflips_abs_supG_resid_field_30_subj_124_dim idea_one_signflips_max_supG_resid_field_30_subj_124_dim idea_one_signflips_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_one_gauss_abs_supG_resid_field_120_subj_124_dim idea_one_gauss_max_supG_resid_field_120_subj_124_dim idea_one_gauss_min_supG_resid_field_120_subj_124_dim idea_one_gauss_abs_supG_resid_field_60_subj_124_dim idea_one_gauss_max_supG_resid_field_60_subj_124_dim idea_one_gauss_min_supG_resid_field_60_subj_124_dim idea_one_gauss_abs_supG_resid_field_30_subj_124_dim idea_one_gauss_max_supG_resid_field_30_subj_124_dim idea_one_gauss_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_one_mammen_abs_supG_resid_field_120_subj_124_dim idea_one_mammen_max_supG_resid_field_120_subj_124_dim idea_one_mammen_min_supG_resid_field_120_subj_124_dim idea_one_mammen_abs_supG_resid_field_60_subj_124_dim idea_one_mammen_max_supG_resid_field_60_subj_124_dim idea_one_mammen_min_supG_resid_field_60_subj_124_dim idea_one_mammen_abs_supG_resid_field_30_subj_124_dim idea_one_mammen_max_supG_resid_field_30_subj_124_dim idea_one_mammen_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_one_signflips_tboot_abs_supG_resid_field_120_subj_124_dim idea_one_signflips_tboot_max_supG_resid_field_120_subj_124_dim idea_one_signflips_tboot_min_supG_resid_field_120_subj_124_dim idea_one_signflips_tboot_abs_supG_resid_field_60_subj_124_dim idea_one_signflips_tboot_max_supG_resid_field_60_subj_124_dim idea_one_signflips_tboot_min_supG_resid_field_60_subj_124_dim idea_one_signflips_tboot_abs_supG_resid_field_30_subj_124_dim idea_one_signflips_tboot_max_supG_resid_field_30_subj_124_dim idea_one_signflips_tboot_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_one_gauss_tboot_abs_supG_resid_field_120_subj_124_dim idea_one_gauss_tboot_max_supG_resid_field_120_subj_124_dim idea_one_gauss_tboot_min_supG_resid_field_120_subj_124_dim idea_one_gauss_tboot_abs_supG_resid_field_60_subj_124_dim idea_one_gauss_tboot_max_supG_resid_field_60_subj_124_dim idea_one_gauss_tboot_min_supG_resid_field_60_subj_124_dim idea_one_gauss_tboot_abs_supG_resid_field_30_subj_124_dim idea_one_gauss_tboot_max_supG_resid_field_30_subj_124_dim idea_one_gauss_tboot_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_one_mammen_tboot_abs_supG_resid_field_120_subj_124_dim idea_one_mammen_tboot_max_supG_resid_field_120_subj_124_dim idea_one_mammen_tboot_min_supG_resid_field_120_subj_124_dim idea_one_mammen_tboot_abs_supG_resid_field_60_subj_124_dim idea_one_mammen_tboot_max_supG_resid_field_60_subj_124_dim idea_one_mammen_tboot_min_supG_resid_field_60_subj_124_dim idea_one_mammen_tboot_abs_supG_resid_field_30_subj_124_dim idea_one_mammen_tboot_max_supG_resid_field_30_subj_124_dim idea_one_mammen_tboot_min_supG_resid_field_30_subj_124_dim '...
      'idea_f_signflips_abs_supG_resid_field_120_subj_124_dim idea_f_signflips_max_supG_resid_field_120_subj_124_dim idea_f_signflips_min_supG_resid_field_120_subj_124_dim idea_f_signflips_abs_supG_resid_field_60_subj_124_dim idea_f_signflips_max_supG_resid_field_60_subj_124_dim idea_f_signflips_min_supG_resid_field_60_subj_124_dim idea_f_signflips_abs_supG_resid_field_30_subj_124_dim idea_f_signflips_max_supG_resid_field_30_subj_124_dim idea_f_signflips_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_f_gauss_abs_supG_resid_field_120_subj_124_dim idea_f_gauss_max_supG_resid_field_120_subj_124_dim idea_f_gauss_min_supG_resid_field_120_subj_124_dim idea_f_gauss_abs_supG_resid_field_60_subj_124_dim idea_f_gauss_max_supG_resid_field_60_subj_124_dim idea_f_gauss_min_supG_resid_field_60_subj_124_dim idea_f_gauss_abs_supG_resid_field_30_subj_124_dim idea_f_gauss_max_supG_resid_field_30_subj_124_dim idea_f_gauss_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_f_mammen_abs_supG_resid_field_120_subj_124_dim idea_f_mammen_max_supG_resid_field_120_subj_124_dim idea_f_mammen_min_supG_resid_field_120_subj_124_dim idea_f_mammen_abs_supG_resid_field_60_subj_124_dim idea_f_mammen_max_supG_resid_field_60_subj_124_dim idea_f_mammen_min_supG_resid_field_60_subj_124_dim idea_f_mammen_abs_supG_resid_field_30_subj_124_dim idea_f_mammen_max_supG_resid_field_30_subj_124_dim idea_f_mammen_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_f_signflips_tboot_abs_supG_resid_field_120_subj_124_dim idea_f_signflips_tboot_max_supG_resid_field_120_subj_124_dim idea_f_signflips_tboot_min_supG_resid_field_120_subj_124_dim idea_f_signflips_tboot_abs_supG_resid_field_60_subj_124_dim idea_f_signflips_tboot_max_supG_resid_field_60_subj_124_dim idea_f_signflips_tboot_min_supG_resid_field_60_subj_124_dim idea_f_signflips_tboot_abs_supG_resid_field_30_subj_124_dim idea_f_signflips_tboot_max_supG_resid_field_30_subj_124_dim idea_f_signflips_tboot_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_f_gauss_tboot_abs_supG_resid_field_120_subj_124_dim idea_f_gauss_tboot_max_supG_resid_field_120_subj_124_dim idea_f_gauss_tboot_min_supG_resid_field_120_subj_124_dim idea_f_gauss_tboot_abs_supG_resid_field_60_subj_124_dim idea_f_gauss_tboot_max_supG_resid_field_60_subj_124_dim idea_f_gauss_tboot_min_supG_resid_field_60_subj_124_dim idea_f_gauss_tboot_abs_supG_resid_field_30_subj_124_dim idea_f_gauss_tboot_max_supG_resid_field_30_subj_124_dim idea_f_gauss_tboot_min_supG_resid_field_30_subj_124_dim ' ...
      'idea_f_mammen_tboot_abs_supG_resid_field_120_subj_124_dim idea_f_mammen_tboot_max_supG_resid_field_120_subj_124_dim idea_f_mammen_tboot_min_supG_resid_field_120_subj_124_dim idea_f_mammen_tboot_abs_supG_resid_field_60_subj_124_dim idea_f_mammen_tboot_max_supG_resid_field_60_subj_124_dim idea_f_mammen_tboot_min_supG_resid_field_60_subj_124_dim idea_f_mammen_tboot_abs_supG_resid_field_30_subj_124_dim idea_f_mammen_tboot_max_supG_resid_field_30_subj_124_dim idea_f_mammen_tboot_min_supG_resid_field_30_subj_124_dim '...
      'idea_one_signflips_abs_supG_resid_field_120_subj_60_dim idea_one_signflips_max_supG_resid_field_120_subj_60_dim idea_one_signflips_min_supG_resid_field_120_subj_60_dim idea_one_signflips_abs_supG_resid_field_60_subj_60_dim idea_one_signflips_max_supG_resid_field_60_subj_60_dim idea_one_signflips_min_supG_resid_field_60_subj_60_dim idea_one_signflips_abs_supG_resid_field_30_subj_60_dim idea_one_signflips_max_supG_resid_field_30_subj_60_dim idea_one_signflips_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_one_gauss_abs_supG_resid_field_120_subj_60_dim idea_one_gauss_max_supG_resid_field_120_subj_60_dim idea_one_gauss_min_supG_resid_field_120_subj_60_dim idea_one_gauss_abs_supG_resid_field_60_subj_60_dim idea_one_gauss_max_supG_resid_field_60_subj_60_dim idea_one_gauss_min_supG_resid_field_60_subj_60_dim idea_one_gauss_abs_supG_resid_field_30_subj_60_dim idea_one_gauss_max_supG_resid_field_30_subj_60_dim idea_one_gauss_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_one_mammen_abs_supG_resid_field_120_subj_60_dim idea_one_mammen_max_supG_resid_field_120_subj_60_dim idea_one_mammen_min_supG_resid_field_120_subj_60_dim idea_one_mammen_abs_supG_resid_field_60_subj_60_dim idea_one_mammen_max_supG_resid_field_60_subj_60_dim idea_one_mammen_min_supG_resid_field_60_subj_60_dim idea_one_mammen_abs_supG_resid_field_30_subj_60_dim idea_one_mammen_max_supG_resid_field_30_subj_60_dim idea_one_mammen_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_one_signflips_tboot_abs_supG_resid_field_120_subj_60_dim idea_one_signflips_tboot_max_supG_resid_field_120_subj_60_dim idea_one_signflips_tboot_min_supG_resid_field_120_subj_60_dim idea_one_signflips_tboot_abs_supG_resid_field_60_subj_60_dim idea_one_signflips_tboot_max_supG_resid_field_60_subj_60_dim idea_one_signflips_tboot_min_supG_resid_field_60_subj_60_dim idea_one_signflips_tboot_abs_supG_resid_field_30_subj_60_dim idea_one_signflips_tboot_max_supG_resid_field_30_subj_60_dim idea_one_signflips_tboot_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_one_gauss_tboot_abs_supG_resid_field_120_subj_60_dim idea_one_gauss_tboot_max_supG_resid_field_120_subj_60_dim idea_one_gauss_tboot_min_supG_resid_field_120_subj_60_dim idea_one_gauss_tboot_abs_supG_resid_field_60_subj_60_dim idea_one_gauss_tboot_max_supG_resid_field_60_subj_60_dim idea_one_gauss_tboot_min_supG_resid_field_60_subj_60_dim idea_one_gauss_tboot_abs_supG_resid_field_30_subj_60_dim idea_one_gauss_tboot_max_supG_resid_field_30_subj_60_dim idea_one_gauss_tboot_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_one_mammen_tboot_abs_supG_resid_field_120_subj_60_dim idea_one_mammen_tboot_max_supG_resid_field_120_subj_60_dim idea_one_mammen_tboot_min_supG_resid_field_120_subj_60_dim idea_one_mammen_tboot_abs_supG_resid_field_60_subj_60_dim idea_one_mammen_tboot_max_supG_resid_field_60_subj_60_dim idea_one_mammen_tboot_min_supG_resid_field_60_subj_60_dim idea_one_mammen_tboot_abs_supG_resid_field_30_subj_60_dim idea_one_mammen_tboot_max_supG_resid_field_30_subj_60_dim idea_one_mammen_tboot_min_supG_resid_field_30_subj_60_dim '...
      'idea_f_signflips_abs_supG_resid_field_120_subj_60_dim idea_f_signflips_max_supG_resid_field_120_subj_60_dim idea_f_signflips_min_supG_resid_field_120_subj_60_dim idea_f_signflips_abs_supG_resid_field_60_subj_60_dim idea_f_signflips_max_supG_resid_field_60_subj_60_dim idea_f_signflips_min_supG_resid_field_60_subj_60_dim idea_f_signflips_abs_supG_resid_field_30_subj_60_dim idea_f_signflips_max_supG_resid_field_30_subj_60_dim idea_f_signflips_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_f_gauss_abs_supG_resid_field_120_subj_60_dim idea_f_gauss_max_supG_resid_field_120_subj_60_dim idea_f_gauss_min_supG_resid_field_120_subj_60_dim idea_f_gauss_abs_supG_resid_field_60_subj_60_dim idea_f_gauss_max_supG_resid_field_60_subj_60_dim idea_f_gauss_min_supG_resid_field_60_subj_60_dim idea_f_gauss_abs_supG_resid_field_30_subj_60_dim idea_f_gauss_max_supG_resid_field_30_subj_60_dim idea_f_gauss_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_f_mammen_abs_supG_resid_field_120_subj_60_dim idea_f_mammen_max_supG_resid_field_120_subj_60_dim idea_f_mammen_min_supG_resid_field_120_subj_60_dim idea_f_mammen_abs_supG_resid_field_60_subj_60_dim idea_f_mammen_max_supG_resid_field_60_subj_60_dim idea_f_mammen_min_supG_resid_field_60_subj_60_dim idea_f_mammen_abs_supG_resid_field_30_subj_60_dim idea_f_mammen_max_supG_resid_field_30_subj_60_dim idea_f_mammen_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_f_signflips_tboot_abs_supG_resid_field_120_subj_60_dim idea_f_signflips_tboot_max_supG_resid_field_120_subj_60_dim idea_f_signflips_tboot_min_supG_resid_field_120_subj_60_dim idea_f_signflips_tboot_abs_supG_resid_field_60_subj_60_dim idea_f_signflips_tboot_max_supG_resid_field_60_subj_60_dim idea_f_signflips_tboot_min_supG_resid_field_60_subj_60_dim idea_f_signflips_tboot_abs_supG_resid_field_30_subj_60_dim idea_f_signflips_tboot_max_supG_resid_field_30_subj_60_dim idea_f_signflips_tboot_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_f_gauss_tboot_abs_supG_resid_field_120_subj_60_dim idea_f_gauss_tboot_max_supG_resid_field_120_subj_60_dim idea_f_gauss_tboot_min_supG_resid_field_120_subj_60_dim idea_f_gauss_tboot_abs_supG_resid_field_60_subj_60_dim idea_f_gauss_tboot_max_supG_resid_field_60_subj_60_dim idea_f_gauss_tboot_min_supG_resid_field_60_subj_60_dim idea_f_gauss_tboot_abs_supG_resid_field_30_subj_60_dim idea_f_gauss_tboot_max_supG_resid_field_30_subj_60_dim idea_f_gauss_tboot_min_supG_resid_field_30_subj_60_dim ' ...
      'idea_f_mammen_tboot_abs_supG_resid_field_120_subj_60_dim idea_f_mammen_tboot_max_supG_resid_field_120_subj_60_dim idea_f_mammen_tboot_min_supG_resid_field_120_subj_60_dim idea_f_mammen_tboot_abs_supG_resid_field_60_subj_60_dim idea_f_mammen_tboot_max_supG_resid_field_60_subj_60_dim idea_f_mammen_tboot_min_supG_resid_field_60_subj_60_dim idea_f_mammen_tboot_abs_supG_resid_field_30_subj_60_dim idea_f_mammen_tboot_max_supG_resid_field_30_subj_60_dim idea_f_mammen_tboot_min_supG_resid_field_30_subj_60_dim '...
      'idea_one_signflips_abs_supG_resid_field_120_subj_10_dim idea_one_signflips_max_supG_resid_field_120_subj_10_dim idea_one_signflips_min_supG_resid_field_120_subj_10_dim idea_one_signflips_abs_supG_resid_field_60_subj_10_dim idea_one_signflips_max_supG_resid_field_60_subj_10_dim idea_one_signflips_min_supG_resid_field_60_subj_10_dim idea_one_signflips_abs_supG_resid_field_30_subj_10_dim idea_one_signflips_max_supG_resid_field_30_subj_10_dim idea_one_signflips_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_one_gauss_abs_supG_resid_field_120_subj_10_dim idea_one_gauss_max_supG_resid_field_120_subj_10_dim idea_one_gauss_min_supG_resid_field_120_subj_10_dim idea_one_gauss_abs_supG_resid_field_60_subj_10_dim idea_one_gauss_max_supG_resid_field_60_subj_10_dim idea_one_gauss_min_supG_resid_field_60_subj_10_dim idea_one_gauss_abs_supG_resid_field_30_subj_10_dim idea_one_gauss_max_supG_resid_field_30_subj_10_dim idea_one_gauss_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_one_mammen_abs_supG_resid_field_120_subj_10_dim idea_one_mammen_max_supG_resid_field_120_subj_10_dim idea_one_mammen_min_supG_resid_field_120_subj_10_dim idea_one_mammen_abs_supG_resid_field_60_subj_10_dim idea_one_mammen_max_supG_resid_field_60_subj_10_dim idea_one_mammen_min_supG_resid_field_60_subj_10_dim idea_one_mammen_abs_supG_resid_field_30_subj_10_dim idea_one_mammen_max_supG_resid_field_30_subj_10_dim idea_one_mammen_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_one_signflips_tboot_abs_supG_resid_field_120_subj_10_dim idea_one_signflips_tboot_max_supG_resid_field_120_subj_10_dim idea_one_signflips_tboot_min_supG_resid_field_120_subj_10_dim idea_one_signflips_tboot_abs_supG_resid_field_60_subj_10_dim idea_one_signflips_tboot_max_supG_resid_field_60_subj_10_dim idea_one_signflips_tboot_min_supG_resid_field_60_subj_10_dim idea_one_signflips_tboot_abs_supG_resid_field_30_subj_10_dim idea_one_signflips_tboot_max_supG_resid_field_30_subj_10_dim idea_one_signflips_tboot_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_one_gauss_tboot_abs_supG_resid_field_120_subj_10_dim idea_one_gauss_tboot_max_supG_resid_field_120_subj_10_dim idea_one_gauss_tboot_min_supG_resid_field_120_subj_10_dim idea_one_gauss_tboot_abs_supG_resid_field_60_subj_10_dim idea_one_gauss_tboot_max_supG_resid_field_60_subj_10_dim idea_one_gauss_tboot_min_supG_resid_field_60_subj_10_dim idea_one_gauss_tboot_abs_supG_resid_field_30_subj_10_dim idea_one_gauss_tboot_max_supG_resid_field_30_subj_10_dim idea_one_gauss_tboot_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_one_mammen_tboot_abs_supG_resid_field_120_subj_10_dim idea_one_mammen_tboot_max_supG_resid_field_120_subj_10_dim idea_one_mammen_tboot_min_supG_resid_field_120_subj_10_dim idea_one_mammen_tboot_abs_supG_resid_field_60_subj_10_dim idea_one_mammen_tboot_max_supG_resid_field_60_subj_10_dim idea_one_mammen_tboot_min_supG_resid_field_60_subj_10_dim idea_one_mammen_tboot_abs_supG_resid_field_30_subj_10_dim idea_one_mammen_tboot_max_supG_resid_field_30_subj_10_dim idea_one_mammen_tboot_min_supG_resid_field_30_subj_10_dim '...
      'idea_f_signflips_abs_supG_resid_field_120_subj_10_dim idea_f_signflips_max_supG_resid_field_120_subj_10_dim idea_f_signflips_min_supG_resid_field_120_subj_10_dim idea_f_signflips_abs_supG_resid_field_60_subj_10_dim idea_f_signflips_max_supG_resid_field_60_subj_10_dim idea_f_signflips_min_supG_resid_field_60_subj_10_dim idea_f_signflips_abs_supG_resid_field_30_subj_10_dim idea_f_signflips_max_supG_resid_field_30_subj_10_dim idea_f_signflips_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_f_gauss_abs_supG_resid_field_120_subj_10_dim idea_f_gauss_max_supG_resid_field_120_subj_10_dim idea_f_gauss_min_supG_resid_field_120_subj_10_dim idea_f_gauss_abs_supG_resid_field_60_subj_10_dim idea_f_gauss_max_supG_resid_field_60_subj_10_dim idea_f_gauss_min_supG_resid_field_60_subj_10_dim idea_f_gauss_abs_supG_resid_field_30_subj_10_dim idea_f_gauss_max_supG_resid_field_30_subj_10_dim idea_f_gauss_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_f_mammen_abs_supG_resid_field_120_subj_10_dim idea_f_mammen_max_supG_resid_field_120_subj_10_dim idea_f_mammen_min_supG_resid_field_120_subj_10_dim idea_f_mammen_abs_supG_resid_field_60_subj_10_dim idea_f_mammen_max_supG_resid_field_60_subj_10_dim idea_f_mammen_min_supG_resid_field_60_subj_10_dim idea_f_mammen_abs_supG_resid_field_30_subj_10_dim idea_f_mammen_max_supG_resid_field_30_subj_10_dim idea_f_mammen_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_f_signflips_tboot_abs_supG_resid_field_120_subj_10_dim idea_f_signflips_tboot_max_supG_resid_field_120_subj_10_dim idea_f_signflips_tboot_min_supG_resid_field_120_subj_10_dim idea_f_signflips_tboot_abs_supG_resid_field_60_subj_10_dim idea_f_signflips_tboot_max_supG_resid_field_60_subj_10_dim idea_f_signflips_tboot_min_supG_resid_field_60_subj_10_dim idea_f_signflips_tboot_abs_supG_resid_field_30_subj_10_dim idea_f_signflips_tboot_max_supG_resid_field_30_subj_10_dim idea_f_signflips_tboot_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_f_gauss_tboot_abs_supG_resid_field_120_subj_10_dim idea_f_gauss_tboot_max_supG_resid_field_120_subj_10_dim idea_f_gauss_tboot_min_supG_resid_field_120_subj_10_dim idea_f_gauss_tboot_abs_supG_resid_field_60_subj_10_dim idea_f_gauss_tboot_max_supG_resid_field_60_subj_10_dim idea_f_gauss_tboot_min_supG_resid_field_60_subj_10_dim idea_f_gauss_tboot_abs_supG_resid_field_30_subj_10_dim idea_f_gauss_tboot_max_supG_resid_field_30_subj_10_dim idea_f_gauss_tboot_min_supG_resid_field_30_subj_10_dim ' ...
      'idea_f_mammen_tboot_abs_supG_resid_field_120_subj_10_dim idea_f_mammen_tboot_max_supG_resid_field_120_subj_10_dim idea_f_mammen_tboot_min_supG_resid_field_120_subj_10_dim idea_f_mammen_tboot_abs_supG_resid_field_60_subj_10_dim idea_f_mammen_tboot_max_supG_resid_field_60_subj_10_dim idea_f_mammen_tboot_min_supG_resid_field_60_subj_10_dim idea_f_mammen_tboot_abs_supG_resid_field_30_subj_10_dim idea_f_mammen_tboot_max_supG_resid_field_30_subj_10_dim idea_f_mammen_tboot_min_supG_resid_field_30_subj_10_dim'
      ])
