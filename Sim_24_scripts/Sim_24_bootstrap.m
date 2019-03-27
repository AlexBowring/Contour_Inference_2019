function Sim_24(nSubj, SvNm, nRlz, mag, smo)
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
% nSubj  = 240_subj_124_dim;
% nRlz = 60_subj_124_dim0;
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

% Variables for the transformation
a_240_subj       = sqrt((nSubj-1)/(nSubj -3 ));
b_240_subj       = sqrt(nSubj)*sqrt((8*nSubj^2 - 17*nSubj + 11)/((5-4*nSubj)^2*(nSubj-3)));
alpha_240_subj   = b_240_subj^-1;
beta_240_subj    = b_240_subj/a_240_subj;

a_120_subj       = sqrt(((nSubj/2)-1)/((nSubj/2) -3 ));
b_120_subj       = sqrt((nSubj/2))*sqrt((8*(nSubj/2)^2 - 17*(nSubj/2) + 11)/((5-4*(nSubj/2))^2*((nSubj/2)-3)));
alpha_120_subj   = b_120_subj^-1;
beta_120_subj    = b_120_subj/a_120_subj;

a_60_subj       = sqrt(((nSubj/4)-1)/((nSubj/4) -3 ));
b_60_subj       = sqrt((nSubj/4))*sqrt((8*(nSubj/4)^2 - 17*(nSubj/4) + 11)/((5-4*(nSubj/4))^2*((nSubj/4)-3)));
alpha_60_subj   = b_60_subj^-1;
beta_60_subj    = b_60_subj/a_60_subj;

a_30_subj       = sqrt(((nSubj/8)-1)/((nSubj/8) -3 ));
b_30_subj       = sqrt((nSubj/8))*sqrt((8*(nSubj/8)^2 - 17*(nSubj/8) + 11)/((5-4*(nSubj/8))^2*((nSubj/8)-3)));
alpha_30_subj   = b_30_subj^-1;
beta_30_subj    = b_30_subj/a_30_subj;

transformed_thr_240_subj = alpha_240_subj*asinh((1 - (3/(4*nSubj - 1)))^(-1)*thr*beta_240_subj);
transformed_thr_120_subj = alpha_120_subj*asinh((1 - (3/(4*(nSubj/2) - 1)))^(-1)*thr*beta_120_subj);
transformed_thr_60_subj  = alpha_60_subj*asinh((1 - (3/(4*(nSubj/4) - 1)))^(-1)*thr*beta_60_subj);
transformed_thr_30_subj  = alpha_30_subj*asinh((1 - (3/(4*(nSubj/8) - 1)))^(-1)*thr*beta_30_subj);

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data_240_subj_124_dim = zeros([dim, nSubj]);

signflips_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);

signflips_tboot_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_240_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_120_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_60_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_30_subj_124_dim = zeros(nBoot, 1);

signflips_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);

signflips_tboot_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_240_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_120_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_60_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_30_subj_60_dim = zeros(nBoot, 1);

signflips_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
signflips_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
signflips_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);

signflips_tboot_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
signflips_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
gauss_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_240_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_120_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_60_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_abs_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_max_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);
mammen_tboot_min_supG_resid_field_30_subj_10_dim = zeros(nBoot, 1);


% Creating linearly increasing signal across columns
Sig = ones(dim).*mag;
  
for t=1:nRlz
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        Noise = create_noise(wdim, 'homo', 1, smo, trnind);
        tImgs = Sig + Noise; % Creates the true image of smoothed signal + smoothed noise
        observed_data_240_subj_124_dim(:,:,i) = tImgs;
        
      end %========== Loop i (subjects)
      
      %% Getting everything for 124 x 124 dimension
      
      observed_data_120_subj_124_dim = observed_data_240_subj_124_dim(:,:,1:(nSubj/2));
      observed_data_60_subj_124_dim = observed_data_240_subj_124_dim(:,:,1:(nSubj/4));
      observed_data_30_subj_124_dim = observed_data_240_subj_124_dim(:,:,1:(nSubj/8));
      
      observed_mean_240_subj_124_dim = mean(observed_data_240_subj_124_dim,3);
      observed_mean_120_subj_124_dim = mean(observed_data_120_subj_124_dim,3);
      observed_mean_60_subj_124_dim = mean(observed_data_60_subj_124_dim,3);
      observed_mean_30_subj_124_dim = mean(observed_data_30_subj_124_dim,3);
      
      observed_std_240_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_240_subj_124_dim,[prod(dim) nSubj]),stdblk),...
           dim);
      observed_std_120_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_124_dim,[prod(dim) nSubj/2]),stdblk),...
           dim);
      observed_std_60_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_124_dim,[prod(dim) nSubj/4]),stdblk),...
           dim);
      observed_std_30_subj_124_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_124_dim,[prod(dim) nSubj/8]),stdblk),...
           dim);
       
      observed_cohen_d_240_subj_124_dim = observed_mean_240_subj_124_dim./observed_std_240_subj_124_dim;
      observed_cohen_d_120_subj_124_dim = observed_mean_120_subj_124_dim./observed_std_120_subj_124_dim;
      observed_cohen_d_60_subj_124_dim = observed_mean_60_subj_124_dim./observed_std_60_subj_124_dim;
      observed_cohen_d_30_subj_124_dim = observed_mean_30_subj_124_dim./observed_std_30_subj_124_dim;
      
      snr_resid_240_subj_124_dim    = create_resid(observed_data_240_subj_124_dim, observed_mean_240_subj_124_dim, observed_std_240_subj_124_dim, 2);
      snr_resid_120_subj_124_dim     = create_resid(observed_data_120_subj_124_dim, observed_mean_120_subj_124_dim, observed_std_120_subj_124_dim, 2);
      snr_resid_60_subj_124_dim     = create_resid(observed_data_60_subj_124_dim, observed_mean_60_subj_124_dim, observed_std_60_subj_124_dim, 2);
      snr_resid_30_subj_124_dim     = create_resid(observed_data_30_subj_124_dim, observed_mean_30_subj_124_dim, observed_std_30_subj_124_dim, 2);
      
      observed_residual_std_240_subj_124_dim = reshape(sqrt(1+observed_cohen_d_240_subj_124_dim.^2*(beta_240_subj^2))./(alpha_240_subj*beta_240_subj), [prod(dim) 1]);
      observed_residual_std_120_subj_124_dim = reshape(sqrt(1+observed_cohen_d_120_subj_124_dim.^2*(beta_120_subj^2))./(alpha_120_subj*beta_120_subj), [prod(dim) 1]);
      observed_residual_std_60_subj_124_dim = reshape(sqrt(1+observed_cohen_d_60_subj_124_dim.^2*(beta_60_subj^2))./(alpha_60_subj*beta_60_subj), [prod(dim) 1]);
      observed_residual_std_30_subj_124_dim = reshape(sqrt(1+observed_cohen_d_30_subj_124_dim.^2*(beta_30_subj^2))./(alpha_30_subj*beta_30_subj), [prod(dim) 1]);
            
      residuals_240_subj_124_dim = snr_resid_240_subj_124_dim./observed_residual_std_240_subj_124_dim;
      residuals_120_subj_124_dim = snr_resid_120_subj_124_dim./observed_residual_std_120_subj_124_dim;
      residuals_60_subj_124_dim = snr_resid_60_subj_124_dim./observed_residual_std_60_subj_124_dim;
      residuals_30_subj_124_dim = snr_resid_30_subj_124_dim./observed_residual_std_30_subj_124_dim;
     
      mammen_240_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_240_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_120_subj                               = signflips_240_subj(1:nSubj/2);
          signflips_60_subj                               = signflips_240_subj(1:nSubj/4);
          signflips_30_subj                               = signflips_240_subj(1:nSubj/8);
          
          % Gauss variables 
          gauss_240_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_120_subj                                   = gauss_240_subj(1:nSubj/2);
          gauss_60_subj                                   = gauss_240_subj(1:nSubj/4);
          gauss_30_subj                                   = gauss_240_subj(1:nSubj/8);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_240_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_240_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_120_subj                                  = mammen_240_subj(1:nSubj/2);
          mammen_60_subj                                  = mammen_240_subj(1:nSubj/4);
          mammen_30_subj                                  = mammen_240_subj(1:nSubj/8);
          
          %% bootstrap
          signflips_bootstrap_240_subj_124_dim = residuals_240_subj_124_dim*spdiags(signflips_240_subj, 0, nSubj, nSubj);
          signflips_bootstrap_120_subj_124_dim = residuals_120_subj_124_dim*spdiags(signflips_120_subj, 0, nSubj/2, nSubj/2);
          signflips_bootstrap_60_subj_124_dim = residuals_60_subj_124_dim*spdiags(signflips_60_subj, 0, nSubj/4, nSubj/4);
          signflips_bootstrap_30_subj_124_dim = residuals_30_subj_124_dim*spdiags(signflips_30_subj, 0, nSubj/8, nSubj/8);
          gauss_bootstrap_240_subj_124_dim = residuals_240_subj_124_dim*spdiags(gauss_240_subj, 0, nSubj, nSubj);
          gauss_bootstrap_120_subj_124_dim = residuals_120_subj_124_dim*spdiags(gauss_120_subj, 0, nSubj/2, nSubj/2);
          gauss_bootstrap_60_subj_124_dim = residuals_60_subj_124_dim*spdiags(gauss_60_subj, 0, nSubj/4, nSubj/4);
          gauss_bootstrap_30_subj_124_dim = residuals_30_subj_124_dim*spdiags(gauss_30_subj, 0, nSubj/8, nSubj/8);
          mammen_bootstrap_240_subj_124_dim = residuals_240_subj_124_dim*spdiags(mammen_240_subj, 0, nSubj, nSubj);
          mammen_bootstrap_120_subj_124_dim = residuals_120_subj_124_dim*spdiags(mammen_120_subj, 0, nSubj/2, nSubj/2);
          mammen_bootstrap_60_subj_124_dim = residuals_60_subj_124_dim*spdiags(mammen_60_subj, 0, nSubj/4, nSubj/4);
          mammen_bootstrap_30_subj_124_dim = residuals_30_subj_124_dim*spdiags(mammen_30_subj, 0, nSubj/8, nSubj/8);
          
          signflips_resid_field_240_subj_124_dim = sum(signflips_bootstrap_240_subj_124_dim, 2)/sqrt(nSubj);
          signflips_resid_field_120_subj_124_dim = sum(signflips_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj/2);
          signflips_resid_field_60_subj_124_dim = sum(signflips_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/4);
          signflips_resid_field_30_subj_124_dim = sum(signflips_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/8);
          gauss_resid_field_240_subj_124_dim = sum(gauss_bootstrap_240_subj_124_dim, 2)/sqrt(nSubj);
          gauss_resid_field_120_subj_124_dim = sum(gauss_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj/2);
          gauss_resid_field_60_subj_124_dim = sum(gauss_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/4);
          gauss_resid_field_30_subj_124_dim = sum(gauss_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/8);
          mammen_resid_field_240_subj_124_dim = sum(mammen_bootstrap_240_subj_124_dim, 2)/sqrt(nSubj);
          mammen_resid_field_120_subj_124_dim = sum(mammen_bootstrap_120_subj_124_dim, 2)/sqrt(nSubj/2);
          mammen_resid_field_60_subj_124_dim = sum(mammen_bootstrap_60_subj_124_dim, 2)/sqrt(nSubj/4);
          mammen_resid_field_30_subj_124_dim = sum(mammen_bootstrap_30_subj_124_dim, 2)/sqrt(nSubj/8);
          
          signflips_boot_std_240_subj_124_dim = std(signflips_bootstrap_240_subj_124_dim, 0, 2);
          signflips_boot_std_120_subj_124_dim = std(signflips_bootstrap_120_subj_124_dim, 0, 2);
          signflips_boot_std_60_subj_124_dim = std(signflips_bootstrap_60_subj_124_dim, 0, 2);
          signflips_boot_std_30_subj_124_dim = std(signflips_bootstrap_30_subj_124_dim, 0, 2);
          gauss_boot_std_240_subj_124_dim = std(gauss_bootstrap_240_subj_124_dim, 0, 2);
          gauss_boot_std_120_subj_124_dim = std(gauss_bootstrap_120_subj_124_dim, 0, 2);
          gauss_boot_std_60_subj_124_dim = std(gauss_bootstrap_60_subj_124_dim, 0, 2);
          gauss_boot_std_30_subj_124_dim = std(gauss_bootstrap_30_subj_124_dim, 0, 2);
          mammen_boot_std_240_subj_124_dim = std(mammen_bootstrap_240_subj_124_dim, 0, 2);
          mammen_boot_std_120_subj_124_dim = std(mammen_bootstrap_120_subj_124_dim, 0, 2);
          mammen_boot_std_60_subj_124_dim = std(mammen_bootstrap_60_subj_124_dim, 0, 2);
          mammen_boot_std_30_subj_124_dim = std(mammen_bootstrap_30_subj_124_dim, 0, 2);
          
          signflips_tboot_resid_field_240_subj_124_dim = signflips_resid_field_240_subj_124_dim./signflips_boot_std_240_subj_124_dim;      
          signflips_tboot_resid_field_120_subj_124_dim = signflips_resid_field_120_subj_124_dim./signflips_boot_std_120_subj_124_dim;  
          signflips_tboot_resid_field_60_subj_124_dim = signflips_resid_field_60_subj_124_dim./signflips_boot_std_60_subj_124_dim;
          signflips_tboot_resid_field_30_subj_124_dim = signflips_resid_field_30_subj_124_dim./signflips_boot_std_30_subj_124_dim;
          gauss_tboot_resid_field_240_subj_124_dim = gauss_resid_field_240_subj_124_dim./gauss_boot_std_240_subj_124_dim;      
          gauss_tboot_resid_field_120_subj_124_dim = gauss_resid_field_120_subj_124_dim./gauss_boot_std_120_subj_124_dim;  
          gauss_tboot_resid_field_60_subj_124_dim = gauss_resid_field_60_subj_124_dim./gauss_boot_std_60_subj_124_dim;
          gauss_tboot_resid_field_30_subj_124_dim = gauss_resid_field_30_subj_124_dim./gauss_boot_std_30_subj_124_dim;
          mammen_tboot_resid_field_240_subj_124_dim = mammen_resid_field_240_subj_124_dim./mammen_boot_std_240_subj_124_dim;      
          mammen_tboot_resid_field_120_subj_124_dim = mammen_resid_field_120_subj_124_dim./mammen_boot_std_120_subj_124_dim;  
          mammen_tboot_resid_field_60_subj_124_dim = mammen_resid_field_60_subj_124_dim./mammen_boot_std_60_subj_124_dim;
          mammen_tboot_resid_field_30_subj_124_dim = mammen_resid_field_30_subj_124_dim./mammen_boot_std_30_subj_124_dim;
          
          signflips_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(signflips_resid_field_240_subj_124_dim));
          signflips_max_supG_resid_field_240_subj_124_dim(k) = max(signflips_resid_field_240_subj_124_dim);
          signflips_min_supG_resid_field_240_subj_124_dim(k) = min(signflips_resid_field_240_subj_124_dim);
          signflips_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(signflips_resid_field_120_subj_124_dim));
          signflips_max_supG_resid_field_120_subj_124_dim(k) = max(signflips_resid_field_120_subj_124_dim);
          signflips_min_supG_resid_field_120_subj_124_dim(k) = min(signflips_resid_field_120_subj_124_dim);
          signflips_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(signflips_resid_field_60_subj_124_dim));
          signflips_max_supG_resid_field_60_subj_124_dim(k) = max(signflips_resid_field_60_subj_124_dim);
          signflips_min_supG_resid_field_60_subj_124_dim(k) = min(signflips_resid_field_60_subj_124_dim);
          signflips_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(signflips_resid_field_30_subj_124_dim));
          signflips_max_supG_resid_field_30_subj_124_dim(k) = max(signflips_resid_field_30_subj_124_dim);
          signflips_min_supG_resid_field_30_subj_124_dim(k) = min(signflips_resid_field_30_subj_124_dim);
          gauss_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(gauss_resid_field_240_subj_124_dim));
          gauss_max_supG_resid_field_240_subj_124_dim(k) = max(gauss_resid_field_240_subj_124_dim);
          gauss_min_supG_resid_field_240_subj_124_dim(k) = min(gauss_resid_field_240_subj_124_dim);
          gauss_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(gauss_resid_field_120_subj_124_dim));
          gauss_max_supG_resid_field_120_subj_124_dim(k) = max(gauss_resid_field_120_subj_124_dim);
          gauss_min_supG_resid_field_120_subj_124_dim(k) = min(gauss_resid_field_120_subj_124_dim);
          gauss_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(gauss_resid_field_60_subj_124_dim));
          gauss_max_supG_resid_field_60_subj_124_dim(k) = max(gauss_resid_field_60_subj_124_dim);
          gauss_min_supG_resid_field_60_subj_124_dim(k) = min(gauss_resid_field_60_subj_124_dim);
          gauss_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(gauss_resid_field_30_subj_124_dim));
          gauss_max_supG_resid_field_30_subj_124_dim(k) = max(gauss_resid_field_30_subj_124_dim);
          gauss_min_supG_resid_field_30_subj_124_dim(k) = min(gauss_resid_field_30_subj_124_dim);
          mammen_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(mammen_resid_field_240_subj_124_dim));
          mammen_max_supG_resid_field_240_subj_124_dim(k) = max(mammen_resid_field_240_subj_124_dim);
          mammen_min_supG_resid_field_240_subj_124_dim(k) = min(mammen_resid_field_240_subj_124_dim);
          mammen_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(mammen_resid_field_120_subj_124_dim));
          mammen_max_supG_resid_field_120_subj_124_dim(k) = max(mammen_resid_field_120_subj_124_dim);
          mammen_min_supG_resid_field_120_subj_124_dim(k) = min(mammen_resid_field_120_subj_124_dim);
          mammen_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(mammen_resid_field_60_subj_124_dim));
          mammen_max_supG_resid_field_60_subj_124_dim(k) = max(mammen_resid_field_60_subj_124_dim);
          mammen_min_supG_resid_field_60_subj_124_dim(k) = min(mammen_resid_field_60_subj_124_dim);
          mammen_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(mammen_resid_field_30_subj_124_dim));
          mammen_max_supG_resid_field_30_subj_124_dim(k) = max(mammen_resid_field_30_subj_124_dim);
          mammen_min_supG_resid_field_30_subj_124_dim(k) = min(mammen_resid_field_30_subj_124_dim);
          
          signflips_tboot_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(signflips_tboot_resid_field_240_subj_124_dim));
          signflips_tboot_max_supG_resid_field_240_subj_124_dim(k) = max(signflips_tboot_resid_field_240_subj_124_dim);
          signflips_tboot_min_supG_resid_field_240_subj_124_dim(k) = min(signflips_tboot_resid_field_240_subj_124_dim);
          signflips_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(signflips_tboot_resid_field_120_subj_124_dim));
          signflips_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(signflips_tboot_resid_field_120_subj_124_dim);
          signflips_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(signflips_tboot_resid_field_120_subj_124_dim);
          signflips_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(signflips_tboot_resid_field_60_subj_124_dim));
          signflips_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(signflips_tboot_resid_field_60_subj_124_dim);
          signflips_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(signflips_tboot_resid_field_60_subj_124_dim);
          signflips_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(signflips_tboot_resid_field_30_subj_124_dim));
          signflips_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(signflips_tboot_resid_field_30_subj_124_dim);
          signflips_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(signflips_tboot_resid_field_30_subj_124_dim);
          gauss_tboot_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(gauss_tboot_resid_field_240_subj_124_dim));
          gauss_tboot_max_supG_resid_field_240_subj_124_dim(k) = max(gauss_tboot_resid_field_240_subj_124_dim);
          gauss_tboot_min_supG_resid_field_240_subj_124_dim(k) = min(gauss_tboot_resid_field_240_subj_124_dim);
          gauss_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(gauss_tboot_resid_field_120_subj_124_dim));
          gauss_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(gauss_tboot_resid_field_120_subj_124_dim);
          gauss_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(gauss_tboot_resid_field_120_subj_124_dim);
          gauss_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(gauss_tboot_resid_field_60_subj_124_dim));
          gauss_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(gauss_tboot_resid_field_60_subj_124_dim);
          gauss_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(gauss_tboot_resid_field_60_subj_124_dim);
          gauss_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(gauss_tboot_resid_field_30_subj_124_dim));
          gauss_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(gauss_tboot_resid_field_30_subj_124_dim);
          gauss_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(gauss_tboot_resid_field_30_subj_124_dim);
          mammen_tboot_abs_supG_resid_field_240_subj_124_dim(k) = max(abs(mammen_tboot_resid_field_240_subj_124_dim));
          mammen_tboot_max_supG_resid_field_240_subj_124_dim(k) = max(mammen_tboot_resid_field_240_subj_124_dim);
          mammen_tboot_min_supG_resid_field_240_subj_124_dim(k) = min(mammen_tboot_resid_field_240_subj_124_dim);
          mammen_tboot_abs_supG_resid_field_120_subj_124_dim(k) = max(abs(mammen_tboot_resid_field_120_subj_124_dim));
          mammen_tboot_max_supG_resid_field_120_subj_124_dim(k) = max(mammen_tboot_resid_field_120_subj_124_dim);
          mammen_tboot_min_supG_resid_field_120_subj_124_dim(k) = min(mammen_tboot_resid_field_120_subj_124_dim);
          mammen_tboot_abs_supG_resid_field_60_subj_124_dim(k) = max(abs(mammen_tboot_resid_field_60_subj_124_dim));
          mammen_tboot_max_supG_resid_field_60_subj_124_dim(k) = max(mammen_tboot_resid_field_60_subj_124_dim);
          mammen_tboot_min_supG_resid_field_60_subj_124_dim(k) = min(mammen_tboot_resid_field_60_subj_124_dim);
          mammen_tboot_abs_supG_resid_field_30_subj_124_dim(k) = max(abs(mammen_tboot_resid_field_30_subj_124_dim));
          mammen_tboot_max_supG_resid_field_30_subj_124_dim(k) = max(mammen_tboot_resid_field_30_subj_124_dim);
          mammen_tboot_min_supG_resid_field_30_subj_124_dim(k) = min(mammen_tboot_resid_field_30_subj_124_dim);
          
      end
       
      %% Getting everything for 60 x 60 dimension
      
      observed_data_240_subj_60_dim = observed_data_240_subj_124_dim(1:60,1:60,1:nSubj);
      observed_data_120_subj_60_dim = observed_data_240_subj_124_dim(1:60,1:60,1:(nSubj/2));
      observed_data_60_subj_60_dim = observed_data_240_subj_124_dim(1:60,1:60,1:(nSubj/4));
      observed_data_30_subj_60_dim = observed_data_240_subj_124_dim(1:60,1:60,1:(nSubj/8));
      
      observed_mean_240_subj_60_dim = mean(observed_data_240_subj_60_dim,3);
      observed_mean_120_subj_60_dim = mean(observed_data_120_subj_60_dim,3);
      observed_mean_60_subj_60_dim = mean(observed_data_60_subj_60_dim,3);
      observed_mean_30_subj_60_dim = mean(observed_data_30_subj_60_dim,3);
      
      observed_std_240_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_240_subj_60_dim,[prod(dim_60) nSubj]),stdblk_60),...
           dim_60);
      observed_std_120_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_60_dim,[prod(dim_60) nSubj/2]),stdblk_60),...
           dim_60);
      observed_std_60_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_60_dim,[prod(dim_60) nSubj/4]),stdblk_60),...
           dim_60);
      observed_std_30_subj_60_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_60_dim,[prod(dim_60) nSubj/8]),stdblk_60),...
           dim_60);
       
      observed_cohen_d_240_subj_60_dim = observed_mean_240_subj_60_dim./observed_std_240_subj_60_dim;
      observed_cohen_d_120_subj_60_dim = observed_mean_120_subj_60_dim./observed_std_120_subj_60_dim;
      observed_cohen_d_60_subj_60_dim = observed_mean_60_subj_60_dim./observed_std_60_subj_60_dim;
      observed_cohen_d_30_subj_60_dim = observed_mean_30_subj_60_dim./observed_std_30_subj_60_dim;
      
      snr_resid_240_subj_60_dim    = create_resid(observed_data_240_subj_60_dim, observed_mean_240_subj_60_dim, observed_std_240_subj_60_dim, 2);
      snr_resid_120_subj_60_dim     = create_resid(observed_data_120_subj_60_dim, observed_mean_120_subj_60_dim, observed_std_120_subj_60_dim, 2);
      snr_resid_60_subj_60_dim     = create_resid(observed_data_60_subj_60_dim, observed_mean_60_subj_60_dim, observed_std_60_subj_60_dim, 2);
      snr_resid_30_subj_60_dim     = create_resid(observed_data_30_subj_60_dim, observed_mean_30_subj_60_dim, observed_std_30_subj_60_dim, 2);
      
      observed_residual_std_240_subj_60_dim = reshape(sqrt(1+observed_cohen_d_240_subj_60_dim.^2*(beta_240_subj^2))./(alpha_240_subj*beta_240_subj), [prod(dim_60) 1]);
      observed_residual_std_120_subj_60_dim = reshape(sqrt(1+observed_cohen_d_120_subj_60_dim.^2*(beta_120_subj^2))./(alpha_120_subj*beta_120_subj), [prod(dim_60) 1]);
      observed_residual_std_60_subj_60_dim = reshape(sqrt(1+observed_cohen_d_60_subj_60_dim.^2*(beta_60_subj^2))./(alpha_60_subj*beta_60_subj), [prod(dim_60) 1]);
      observed_residual_std_30_subj_60_dim = reshape(sqrt(1+observed_cohen_d_30_subj_60_dim.^2*(beta_30_subj^2))./(alpha_30_subj*beta_30_subj), [prod(dim_60) 1]);
            
      residuals_240_subj_60_dim = snr_resid_240_subj_60_dim./observed_residual_std_240_subj_60_dim;
      residuals_120_subj_60_dim = snr_resid_120_subj_60_dim./observed_residual_std_120_subj_60_dim;
      residuals_60_subj_60_dim = snr_resid_60_subj_60_dim./observed_residual_std_60_subj_60_dim;
      residuals_30_subj_60_dim = snr_resid_30_subj_60_dim./observed_residual_std_30_subj_60_dim;
     
      mammen_240_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_240_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_120_subj                              = signflips_240_subj(1:nSubj/2);
          signflips_60_subj                               = signflips_240_subj(1:nSubj/4);
          signflips_30_subj                               = signflips_240_subj(1:nSubj/8);
          
          % Gauss variables 
          gauss_240_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_120_subj                                   = gauss_240_subj(1:nSubj/2);
          gauss_60_subj                                   = gauss_240_subj(1:nSubj/4);
          gauss_30_subj                                   = gauss_240_subj(1:nSubj/8);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_240_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_240_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_120_subj                                  = mammen_240_subj(1:nSubj/2);
          mammen_60_subj                                  = mammen_240_subj(1:nSubj/4);
          mammen_30_subj                                  = mammen_240_subj(1:nSubj/8);
          
          %% bootstrap
          signflips_bootstrap_240_subj_60_dim = residuals_240_subj_60_dim*spdiags(signflips_240_subj, 0, nSubj, nSubj);
          signflips_bootstrap_120_subj_60_dim = residuals_120_subj_60_dim*spdiags(signflips_120_subj, 0, nSubj/2, nSubj/2);
          signflips_bootstrap_60_subj_60_dim = residuals_60_subj_60_dim*spdiags(signflips_60_subj, 0, nSubj/4, nSubj/4);
          signflips_bootstrap_30_subj_60_dim = residuals_30_subj_60_dim*spdiags(signflips_30_subj, 0, nSubj/8, nSubj/8);
          gauss_bootstrap_240_subj_60_dim = residuals_240_subj_60_dim*spdiags(gauss_240_subj, 0, nSubj, nSubj);
          gauss_bootstrap_120_subj_60_dim = residuals_120_subj_60_dim*spdiags(gauss_120_subj, 0, nSubj/2, nSubj/2);
          gauss_bootstrap_60_subj_60_dim = residuals_60_subj_60_dim*spdiags(gauss_60_subj, 0, nSubj/4, nSubj/4);
          gauss_bootstrap_30_subj_60_dim = residuals_30_subj_60_dim*spdiags(gauss_30_subj, 0, nSubj/8, nSubj/8);
          mammen_bootstrap_240_subj_60_dim = residuals_240_subj_60_dim*spdiags(mammen_240_subj, 0, nSubj, nSubj);
          mammen_bootstrap_120_subj_60_dim = residuals_120_subj_60_dim*spdiags(mammen_120_subj, 0, nSubj/2, nSubj/2);
          mammen_bootstrap_60_subj_60_dim = residuals_60_subj_60_dim*spdiags(mammen_60_subj, 0, nSubj/4, nSubj/4);
          mammen_bootstrap_30_subj_60_dim = residuals_30_subj_60_dim*spdiags(mammen_30_subj, 0, nSubj/8, nSubj/8);
          
          signflips_resid_field_240_subj_60_dim = sum(signflips_bootstrap_240_subj_60_dim, 2)/sqrt(nSubj);
          signflips_resid_field_120_subj_60_dim = sum(signflips_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj/2);
          signflips_resid_field_60_subj_60_dim = sum(signflips_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/4);
          signflips_resid_field_30_subj_60_dim = sum(signflips_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/8);
          gauss_resid_field_240_subj_60_dim = sum(gauss_bootstrap_240_subj_60_dim, 2)/sqrt(nSubj);
          gauss_resid_field_120_subj_60_dim = sum(gauss_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj/2);
          gauss_resid_field_60_subj_60_dim = sum(gauss_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/4);
          gauss_resid_field_30_subj_60_dim = sum(gauss_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/8);
          mammen_resid_field_240_subj_60_dim = sum(mammen_bootstrap_240_subj_60_dim, 2)/sqrt(nSubj);
          mammen_resid_field_120_subj_60_dim = sum(mammen_bootstrap_120_subj_60_dim, 2)/sqrt(nSubj/2);
          mammen_resid_field_60_subj_60_dim = sum(mammen_bootstrap_60_subj_60_dim, 2)/sqrt(nSubj/4);
          mammen_resid_field_30_subj_60_dim = sum(mammen_bootstrap_30_subj_60_dim, 2)/sqrt(nSubj/8);
          
          signflips_boot_std_240_subj_60_dim = std(signflips_bootstrap_240_subj_60_dim, 0, 2);
          signflips_boot_std_120_subj_60_dim = std(signflips_bootstrap_120_subj_60_dim, 0, 2);
          signflips_boot_std_60_subj_60_dim = std(signflips_bootstrap_60_subj_60_dim, 0, 2);
          signflips_boot_std_30_subj_60_dim = std(signflips_bootstrap_30_subj_60_dim, 0, 2);
          gauss_boot_std_240_subj_60_dim = std(gauss_bootstrap_240_subj_60_dim, 0, 2);
          gauss_boot_std_120_subj_60_dim = std(gauss_bootstrap_120_subj_60_dim, 0, 2);
          gauss_boot_std_60_subj_60_dim = std(gauss_bootstrap_60_subj_60_dim, 0, 2);
          gauss_boot_std_30_subj_60_dim = std(gauss_bootstrap_30_subj_60_dim, 0, 2);
          mammen_boot_std_240_subj_60_dim = std(mammen_bootstrap_240_subj_60_dim, 0, 2);
          mammen_boot_std_120_subj_60_dim = std(mammen_bootstrap_120_subj_60_dim, 0, 2);
          mammen_boot_std_60_subj_60_dim = std(mammen_bootstrap_60_subj_60_dim, 0, 2);
          mammen_boot_std_30_subj_60_dim = std(mammen_bootstrap_30_subj_60_dim, 0, 2);
          
          signflips_tboot_resid_field_240_subj_60_dim = signflips_resid_field_240_subj_60_dim./signflips_boot_std_240_subj_60_dim;      
          signflips_tboot_resid_field_120_subj_60_dim = signflips_resid_field_120_subj_60_dim./signflips_boot_std_120_subj_60_dim;  
          signflips_tboot_resid_field_60_subj_60_dim = signflips_resid_field_60_subj_60_dim./signflips_boot_std_60_subj_60_dim;
          signflips_tboot_resid_field_30_subj_60_dim = signflips_resid_field_30_subj_60_dim./signflips_boot_std_30_subj_60_dim;
          gauss_tboot_resid_field_240_subj_60_dim = gauss_resid_field_240_subj_60_dim./gauss_boot_std_240_subj_60_dim;      
          gauss_tboot_resid_field_120_subj_60_dim = gauss_resid_field_120_subj_60_dim./gauss_boot_std_120_subj_60_dim;  
          gauss_tboot_resid_field_60_subj_60_dim = gauss_resid_field_60_subj_60_dim./gauss_boot_std_60_subj_60_dim;
          gauss_tboot_resid_field_30_subj_60_dim = gauss_resid_field_30_subj_60_dim./gauss_boot_std_30_subj_60_dim;
          mammen_tboot_resid_field_240_subj_60_dim = mammen_resid_field_240_subj_60_dim./mammen_boot_std_240_subj_60_dim;      
          mammen_tboot_resid_field_120_subj_60_dim = mammen_resid_field_120_subj_60_dim./mammen_boot_std_120_subj_60_dim;  
          mammen_tboot_resid_field_60_subj_60_dim = mammen_resid_field_60_subj_60_dim./mammen_boot_std_60_subj_60_dim;
          mammen_tboot_resid_field_30_subj_60_dim = mammen_resid_field_30_subj_60_dim./mammen_boot_std_30_subj_60_dim;
          
          signflips_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(signflips_resid_field_240_subj_60_dim));
          signflips_max_supG_resid_field_240_subj_60_dim(k) = max(signflips_resid_field_240_subj_60_dim);
          signflips_min_supG_resid_field_240_subj_60_dim(k) = min(signflips_resid_field_240_subj_60_dim);
          signflips_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(signflips_resid_field_120_subj_60_dim));
          signflips_max_supG_resid_field_120_subj_60_dim(k) = max(signflips_resid_field_120_subj_60_dim);
          signflips_min_supG_resid_field_120_subj_60_dim(k) = min(signflips_resid_field_120_subj_60_dim);
          signflips_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(signflips_resid_field_60_subj_60_dim));
          signflips_max_supG_resid_field_60_subj_60_dim(k) = max(signflips_resid_field_60_subj_60_dim);
          signflips_min_supG_resid_field_60_subj_60_dim(k) = min(signflips_resid_field_60_subj_60_dim);
          signflips_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(signflips_resid_field_30_subj_60_dim));
          signflips_max_supG_resid_field_30_subj_60_dim(k) = max(signflips_resid_field_30_subj_60_dim);
          signflips_min_supG_resid_field_30_subj_60_dim(k) = min(signflips_resid_field_30_subj_60_dim);
          gauss_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(gauss_resid_field_240_subj_60_dim));
          gauss_max_supG_resid_field_240_subj_60_dim(k) = max(gauss_resid_field_240_subj_60_dim);
          gauss_min_supG_resid_field_240_subj_60_dim(k) = min(gauss_resid_field_240_subj_60_dim);
          gauss_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(gauss_resid_field_120_subj_60_dim));
          gauss_max_supG_resid_field_120_subj_60_dim(k) = max(gauss_resid_field_120_subj_60_dim);
          gauss_min_supG_resid_field_120_subj_60_dim(k) = min(gauss_resid_field_120_subj_60_dim);
          gauss_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(gauss_resid_field_60_subj_60_dim));
          gauss_max_supG_resid_field_60_subj_60_dim(k) = max(gauss_resid_field_60_subj_60_dim);
          gauss_min_supG_resid_field_60_subj_60_dim(k) = min(gauss_resid_field_60_subj_60_dim);
          gauss_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(gauss_resid_field_30_subj_60_dim));
          gauss_max_supG_resid_field_30_subj_60_dim(k) = max(gauss_resid_field_30_subj_60_dim);
          gauss_min_supG_resid_field_30_subj_60_dim(k) = min(gauss_resid_field_30_subj_60_dim);
          mammen_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(mammen_resid_field_240_subj_60_dim));
          mammen_max_supG_resid_field_240_subj_60_dim(k) = max(mammen_resid_field_240_subj_60_dim);
          mammen_min_supG_resid_field_240_subj_60_dim(k) = min(mammen_resid_field_240_subj_60_dim);
          mammen_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(mammen_resid_field_120_subj_60_dim));
          mammen_max_supG_resid_field_120_subj_60_dim(k) = max(mammen_resid_field_120_subj_60_dim);
          mammen_min_supG_resid_field_120_subj_60_dim(k) = min(mammen_resid_field_120_subj_60_dim);
          mammen_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(mammen_resid_field_60_subj_60_dim));
          mammen_max_supG_resid_field_60_subj_60_dim(k) = max(mammen_resid_field_60_subj_60_dim);
          mammen_min_supG_resid_field_60_subj_60_dim(k) = min(mammen_resid_field_60_subj_60_dim);
          mammen_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(mammen_resid_field_30_subj_60_dim));
          mammen_max_supG_resid_field_30_subj_60_dim(k) = max(mammen_resid_field_30_subj_60_dim);
          mammen_min_supG_resid_field_30_subj_60_dim(k) = min(mammen_resid_field_30_subj_60_dim);
          
          signflips_tboot_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(signflips_tboot_resid_field_240_subj_60_dim));
          signflips_tboot_max_supG_resid_field_240_subj_60_dim(k) = max(signflips_tboot_resid_field_240_subj_60_dim);
          signflips_tboot_min_supG_resid_field_240_subj_60_dim(k) = min(signflips_tboot_resid_field_240_subj_60_dim);
          signflips_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(signflips_tboot_resid_field_120_subj_60_dim));
          signflips_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(signflips_tboot_resid_field_120_subj_60_dim);
          signflips_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(signflips_tboot_resid_field_120_subj_60_dim);
          signflips_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(signflips_tboot_resid_field_60_subj_60_dim));
          signflips_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(signflips_tboot_resid_field_60_subj_60_dim);
          signflips_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(signflips_tboot_resid_field_60_subj_60_dim);
          signflips_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(signflips_tboot_resid_field_30_subj_60_dim));
          signflips_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(signflips_tboot_resid_field_30_subj_60_dim);
          signflips_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(signflips_tboot_resid_field_30_subj_60_dim);
          gauss_tboot_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(gauss_tboot_resid_field_240_subj_60_dim));
          gauss_tboot_max_supG_resid_field_240_subj_60_dim(k) = max(gauss_tboot_resid_field_240_subj_60_dim);
          gauss_tboot_min_supG_resid_field_240_subj_60_dim(k) = min(gauss_tboot_resid_field_240_subj_60_dim);
          gauss_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(gauss_tboot_resid_field_120_subj_60_dim));
          gauss_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(gauss_tboot_resid_field_120_subj_60_dim);
          gauss_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(gauss_tboot_resid_field_120_subj_60_dim);
          gauss_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(gauss_tboot_resid_field_60_subj_60_dim));
          gauss_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(gauss_tboot_resid_field_60_subj_60_dim);
          gauss_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(gauss_tboot_resid_field_60_subj_60_dim);
          gauss_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(gauss_tboot_resid_field_30_subj_60_dim));
          gauss_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(gauss_tboot_resid_field_30_subj_60_dim);
          gauss_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(gauss_tboot_resid_field_30_subj_60_dim);
          mammen_tboot_abs_supG_resid_field_240_subj_60_dim(k) = max(abs(mammen_tboot_resid_field_240_subj_60_dim));
          mammen_tboot_max_supG_resid_field_240_subj_60_dim(k) = max(mammen_tboot_resid_field_240_subj_60_dim);
          mammen_tboot_min_supG_resid_field_240_subj_60_dim(k) = min(mammen_tboot_resid_field_240_subj_60_dim);
          mammen_tboot_abs_supG_resid_field_120_subj_60_dim(k) = max(abs(mammen_tboot_resid_field_120_subj_60_dim));
          mammen_tboot_max_supG_resid_field_120_subj_60_dim(k) = max(mammen_tboot_resid_field_120_subj_60_dim);
          mammen_tboot_min_supG_resid_field_120_subj_60_dim(k) = min(mammen_tboot_resid_field_120_subj_60_dim);
          mammen_tboot_abs_supG_resid_field_60_subj_60_dim(k) = max(abs(mammen_tboot_resid_field_60_subj_60_dim));
          mammen_tboot_max_supG_resid_field_60_subj_60_dim(k) = max(mammen_tboot_resid_field_60_subj_60_dim);
          mammen_tboot_min_supG_resid_field_60_subj_60_dim(k) = min(mammen_tboot_resid_field_60_subj_60_dim);
          mammen_tboot_abs_supG_resid_field_30_subj_60_dim(k) = max(abs(mammen_tboot_resid_field_30_subj_60_dim));
          mammen_tboot_max_supG_resid_field_30_subj_60_dim(k) = max(mammen_tboot_resid_field_30_subj_60_dim);
          mammen_tboot_min_supG_resid_field_30_subj_60_dim(k) = min(mammen_tboot_resid_field_30_subj_60_dim);
          
      end
      
      %% Getting everything for 10 x 10 dimension
      
      observed_data_240_subj_10_dim = observed_data_240_subj_124_dim(1:10,1:10,1:nSubj);
      observed_data_120_subj_10_dim = observed_data_240_subj_124_dim(1:10,1:10,1:(nSubj/2));
      observed_data_60_subj_10_dim = observed_data_240_subj_124_dim(1:10,1:10,1:(nSubj/4));
      observed_data_30_subj_10_dim = observed_data_240_subj_124_dim(1:10,1:10,1:(nSubj/8));
      
      observed_mean_240_subj_10_dim = mean(observed_data_240_subj_10_dim,3);
      observed_mean_120_subj_10_dim = mean(observed_data_120_subj_10_dim,3);
      observed_mean_60_subj_10_dim = mean(observed_data_60_subj_10_dim,3);
      observed_mean_30_subj_10_dim = mean(observed_data_30_subj_10_dim,3);
      
      observed_std_240_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_240_subj_10_dim,[prod(dim_10) nSubj]),stdblk_10),...
           dim_10);
      observed_std_120_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_120_subj_10_dim,[prod(dim_10) nSubj/2]),stdblk_10),...
           dim_10);
      observed_std_60_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_60_subj_10_dim,[prod(dim_10) nSubj/4]),stdblk_10),...
           dim_10);
      observed_std_30_subj_10_dim = reshape(...
         biasmystd(reshape(observed_data_30_subj_10_dim,[prod(dim_10) nSubj/8]),stdblk_10),...
           dim_10);
       
      observed_cohen_d_240_subj_10_dim = observed_mean_240_subj_10_dim./observed_std_240_subj_10_dim;
      observed_cohen_d_120_subj_10_dim = observed_mean_120_subj_10_dim./observed_std_120_subj_10_dim;
      observed_cohen_d_60_subj_10_dim = observed_mean_60_subj_10_dim./observed_std_60_subj_10_dim;
      observed_cohen_d_30_subj_10_dim = observed_mean_30_subj_10_dim./observed_std_30_subj_10_dim;
      
      snr_resid_240_subj_10_dim    = create_resid(observed_data_240_subj_10_dim, observed_mean_240_subj_10_dim, observed_std_240_subj_10_dim, 2);
      snr_resid_120_subj_10_dim     = create_resid(observed_data_120_subj_10_dim, observed_mean_120_subj_10_dim, observed_std_120_subj_10_dim, 2);
      snr_resid_60_subj_10_dim     = create_resid(observed_data_60_subj_10_dim, observed_mean_60_subj_10_dim, observed_std_60_subj_10_dim, 2);
      snr_resid_30_subj_10_dim     = create_resid(observed_data_30_subj_10_dim, observed_mean_30_subj_10_dim, observed_std_30_subj_10_dim, 2);
      
      observed_residual_std_240_subj_10_dim = reshape(sqrt(1+observed_cohen_d_240_subj_10_dim.^2*(beta_240_subj^2))./(alpha_240_subj*beta_240_subj), [prod(dim_10) 1]);
      observed_residual_std_120_subj_10_dim = reshape(sqrt(1+observed_cohen_d_120_subj_10_dim.^2*(beta_120_subj^2))./(alpha_120_subj*beta_120_subj), [prod(dim_10) 1]);
      observed_residual_std_60_subj_10_dim = reshape(sqrt(1+observed_cohen_d_60_subj_10_dim.^2*(beta_60_subj^2))./(alpha_60_subj*beta_60_subj), [prod(dim_10) 1]);
      observed_residual_std_30_subj_10_dim = reshape(sqrt(1+observed_cohen_d_30_subj_10_dim.^2*(beta_30_subj^2))./(alpha_30_subj*beta_30_subj), [prod(dim_10) 1]);
            
      residuals_240_subj_10_dim = snr_resid_240_subj_10_dim./observed_residual_std_240_subj_10_dim;
      residuals_120_subj_10_dim = snr_resid_120_subj_10_dim./observed_residual_std_120_subj_10_dim;
      residuals_60_subj_10_dim = snr_resid_60_subj_10_dim./observed_residual_std_60_subj_10_dim;
      residuals_30_subj_10_dim = snr_resid_30_subj_10_dim./observed_residual_std_30_subj_10_dim;
     
      mammen_240_subj = zeros(nSubj,1);
      for k=1:nBoot
          % Rademacher variables 
          signflips_240_subj                              = randi(2,[nSubj,1])*2-3;
          signflips_120_subj                              = signflips_240_subj(1:nSubj/2);
          signflips_60_subj                               = signflips_240_subj(1:nSubj/4);
          signflips_30_subj                               = signflips_240_subj(1:nSubj/8);
          
          % Gauss variables 
          gauss_240_subj                                  = normrnd(0,1,[nSubj,1]);
          gauss_120_subj                                   = gauss_240_subj(1:nSubj/2);
          gauss_60_subj                                   = gauss_240_subj(1:nSubj/4);
          gauss_30_subj                                   = gauss_240_subj(1:nSubj/8);
          
          % Mammen variables 
          uniform                            = rand(nSubj,1);
          mammen_240_subj(uniform < (sqrt(5)+1)/(2*sqrt(5))) = -(sqrt(5) - 1)/2;
          mammen_240_subj(uniform > (sqrt(5)+1)/(2*sqrt(5))) = (sqrt(5) + 1)/2;
          mammen_120_subj                                  = mammen_240_subj(1:nSubj/2);
          mammen_60_subj                                  = mammen_240_subj(1:nSubj/4);
          mammen_30_subj                                  = mammen_240_subj(1:nSubj/8);
          
          %% bootstrap
          signflips_bootstrap_240_subj_10_dim = residuals_240_subj_10_dim*spdiags(signflips_240_subj, 0, nSubj, nSubj);
          signflips_bootstrap_120_subj_10_dim = residuals_120_subj_10_dim*spdiags(signflips_120_subj, 0, nSubj/2, nSubj/2);
          signflips_bootstrap_60_subj_10_dim = residuals_60_subj_10_dim*spdiags(signflips_60_subj, 0, nSubj/4, nSubj/4);
          signflips_bootstrap_30_subj_10_dim = residuals_30_subj_10_dim*spdiags(signflips_30_subj, 0, nSubj/8, nSubj/8);
          gauss_bootstrap_240_subj_10_dim = residuals_240_subj_10_dim*spdiags(gauss_240_subj, 0, nSubj, nSubj);
          gauss_bootstrap_120_subj_10_dim = residuals_120_subj_10_dim*spdiags(gauss_120_subj, 0, nSubj/2, nSubj/2);
          gauss_bootstrap_60_subj_10_dim = residuals_60_subj_10_dim*spdiags(gauss_60_subj, 0, nSubj/4, nSubj/4);
          gauss_bootstrap_30_subj_10_dim = residuals_30_subj_10_dim*spdiags(gauss_30_subj, 0, nSubj/8, nSubj/8);
          mammen_bootstrap_240_subj_10_dim = residuals_240_subj_10_dim*spdiags(mammen_240_subj, 0, nSubj, nSubj);
          mammen_bootstrap_120_subj_10_dim = residuals_120_subj_10_dim*spdiags(mammen_120_subj, 0, nSubj/2, nSubj/2);
          mammen_bootstrap_60_subj_10_dim = residuals_60_subj_10_dim*spdiags(mammen_60_subj, 0, nSubj/4, nSubj/4);
          mammen_bootstrap_30_subj_10_dim = residuals_30_subj_10_dim*spdiags(mammen_30_subj, 0, nSubj/8, nSubj/8);
          
          signflips_resid_field_240_subj_10_dim = sum(signflips_bootstrap_240_subj_10_dim, 2)/sqrt(nSubj);
          signflips_resid_field_120_subj_10_dim = sum(signflips_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj/2);
          signflips_resid_field_60_subj_10_dim = sum(signflips_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/4);
          signflips_resid_field_30_subj_10_dim = sum(signflips_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/8);
          gauss_resid_field_240_subj_10_dim = sum(gauss_bootstrap_240_subj_10_dim, 2)/sqrt(nSubj);
          gauss_resid_field_120_subj_10_dim = sum(gauss_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj/2);
          gauss_resid_field_60_subj_10_dim = sum(gauss_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/4);
          gauss_resid_field_30_subj_10_dim = sum(gauss_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/8);
          mammen_resid_field_240_subj_10_dim = sum(mammen_bootstrap_240_subj_10_dim, 2)/sqrt(nSubj);
          mammen_resid_field_120_subj_10_dim = sum(mammen_bootstrap_120_subj_10_dim, 2)/sqrt(nSubj/2);
          mammen_resid_field_60_subj_10_dim = sum(mammen_bootstrap_60_subj_10_dim, 2)/sqrt(nSubj/4);
          mammen_resid_field_30_subj_10_dim = sum(mammen_bootstrap_30_subj_10_dim, 2)/sqrt(nSubj/8);
          
          signflips_boot_std_240_subj_10_dim = std(signflips_bootstrap_240_subj_10_dim, 0, 2);
          signflips_boot_std_120_subj_10_dim = std(signflips_bootstrap_120_subj_10_dim, 0, 2);
          signflips_boot_std_60_subj_10_dim = std(signflips_bootstrap_60_subj_10_dim, 0, 2);
          signflips_boot_std_30_subj_10_dim = std(signflips_bootstrap_30_subj_10_dim, 0, 2);
          gauss_boot_std_240_subj_10_dim = std(gauss_bootstrap_240_subj_10_dim, 0, 2);
          gauss_boot_std_120_subj_10_dim = std(gauss_bootstrap_120_subj_10_dim, 0, 2);
          gauss_boot_std_60_subj_10_dim = std(gauss_bootstrap_60_subj_10_dim, 0, 2);
          gauss_boot_std_30_subj_10_dim = std(gauss_bootstrap_30_subj_10_dim, 0, 2);
          mammen_boot_std_240_subj_10_dim = std(mammen_bootstrap_240_subj_10_dim, 0, 2);
          mammen_boot_std_120_subj_10_dim = std(mammen_bootstrap_120_subj_10_dim, 0, 2);
          mammen_boot_std_60_subj_10_dim = std(mammen_bootstrap_60_subj_10_dim, 0, 2);
          mammen_boot_std_30_subj_10_dim = std(mammen_bootstrap_30_subj_10_dim, 0, 2);
          
          signflips_tboot_resid_field_240_subj_10_dim = signflips_resid_field_240_subj_10_dim./signflips_boot_std_240_subj_10_dim;      
          signflips_tboot_resid_field_120_subj_10_dim = signflips_resid_field_120_subj_10_dim./signflips_boot_std_120_subj_10_dim;  
          signflips_tboot_resid_field_60_subj_10_dim = signflips_resid_field_60_subj_10_dim./signflips_boot_std_60_subj_10_dim;
          signflips_tboot_resid_field_30_subj_10_dim = signflips_resid_field_30_subj_10_dim./signflips_boot_std_30_subj_10_dim;
          gauss_tboot_resid_field_240_subj_10_dim = gauss_resid_field_240_subj_10_dim./gauss_boot_std_240_subj_10_dim;      
          gauss_tboot_resid_field_120_subj_10_dim = gauss_resid_field_120_subj_10_dim./gauss_boot_std_120_subj_10_dim;  
          gauss_tboot_resid_field_60_subj_10_dim = gauss_resid_field_60_subj_10_dim./gauss_boot_std_60_subj_10_dim;
          gauss_tboot_resid_field_30_subj_10_dim = gauss_resid_field_30_subj_10_dim./gauss_boot_std_30_subj_10_dim;
          mammen_tboot_resid_field_240_subj_10_dim = mammen_resid_field_240_subj_10_dim./mammen_boot_std_240_subj_10_dim;      
          mammen_tboot_resid_field_120_subj_10_dim = mammen_resid_field_120_subj_10_dim./mammen_boot_std_120_subj_10_dim;  
          mammen_tboot_resid_field_60_subj_10_dim = mammen_resid_field_60_subj_10_dim./mammen_boot_std_60_subj_10_dim;
          mammen_tboot_resid_field_30_subj_10_dim = mammen_resid_field_30_subj_10_dim./mammen_boot_std_30_subj_10_dim;
          
          signflips_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(signflips_resid_field_240_subj_10_dim));
          signflips_max_supG_resid_field_240_subj_10_dim(k) = max(signflips_resid_field_240_subj_10_dim);
          signflips_min_supG_resid_field_240_subj_10_dim(k) = min(signflips_resid_field_240_subj_10_dim);
          signflips_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(signflips_resid_field_120_subj_10_dim));
          signflips_max_supG_resid_field_120_subj_10_dim(k) = max(signflips_resid_field_120_subj_10_dim);
          signflips_min_supG_resid_field_120_subj_10_dim(k) = min(signflips_resid_field_120_subj_10_dim);
          signflips_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(signflips_resid_field_60_subj_10_dim));
          signflips_max_supG_resid_field_60_subj_10_dim(k) = max(signflips_resid_field_60_subj_10_dim);
          signflips_min_supG_resid_field_60_subj_10_dim(k) = min(signflips_resid_field_60_subj_10_dim);
          signflips_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(signflips_resid_field_30_subj_10_dim));
          signflips_max_supG_resid_field_30_subj_10_dim(k) = max(signflips_resid_field_30_subj_10_dim);
          signflips_min_supG_resid_field_30_subj_10_dim(k) = min(signflips_resid_field_30_subj_10_dim);
          gauss_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(gauss_resid_field_240_subj_10_dim));
          gauss_max_supG_resid_field_240_subj_10_dim(k) = max(gauss_resid_field_240_subj_10_dim);
          gauss_min_supG_resid_field_240_subj_10_dim(k) = min(gauss_resid_field_240_subj_10_dim);
          gauss_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(gauss_resid_field_120_subj_10_dim));
          gauss_max_supG_resid_field_120_subj_10_dim(k) = max(gauss_resid_field_120_subj_10_dim);
          gauss_min_supG_resid_field_120_subj_10_dim(k) = min(gauss_resid_field_120_subj_10_dim);
          gauss_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(gauss_resid_field_60_subj_10_dim));
          gauss_max_supG_resid_field_60_subj_10_dim(k) = max(gauss_resid_field_60_subj_10_dim);
          gauss_min_supG_resid_field_60_subj_10_dim(k) = min(gauss_resid_field_60_subj_10_dim);
          gauss_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(gauss_resid_field_30_subj_10_dim));
          gauss_max_supG_resid_field_30_subj_10_dim(k) = max(gauss_resid_field_30_subj_10_dim);
          gauss_min_supG_resid_field_30_subj_10_dim(k) = min(gauss_resid_field_30_subj_10_dim);
          mammen_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(mammen_resid_field_240_subj_10_dim));
          mammen_max_supG_resid_field_240_subj_10_dim(k) = max(mammen_resid_field_240_subj_10_dim);
          mammen_min_supG_resid_field_240_subj_10_dim(k) = min(mammen_resid_field_240_subj_10_dim);
          mammen_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(mammen_resid_field_120_subj_10_dim));
          mammen_max_supG_resid_field_120_subj_10_dim(k) = max(mammen_resid_field_120_subj_10_dim);
          mammen_min_supG_resid_field_120_subj_10_dim(k) = min(mammen_resid_field_120_subj_10_dim);
          mammen_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(mammen_resid_field_60_subj_10_dim));
          mammen_max_supG_resid_field_60_subj_10_dim(k) = max(mammen_resid_field_60_subj_10_dim);
          mammen_min_supG_resid_field_60_subj_10_dim(k) = min(mammen_resid_field_60_subj_10_dim);
          mammen_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(mammen_resid_field_30_subj_10_dim));
          mammen_max_supG_resid_field_30_subj_10_dim(k) = max(mammen_resid_field_30_subj_10_dim);
          mammen_min_supG_resid_field_30_subj_10_dim(k) = min(mammen_resid_field_30_subj_10_dim);
          
          signflips_tboot_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(signflips_tboot_resid_field_240_subj_10_dim));
          signflips_tboot_max_supG_resid_field_240_subj_10_dim(k) = max(signflips_tboot_resid_field_240_subj_10_dim);
          signflips_tboot_min_supG_resid_field_240_subj_10_dim(k) = min(signflips_tboot_resid_field_240_subj_10_dim);
          signflips_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(signflips_tboot_resid_field_120_subj_10_dim));
          signflips_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(signflips_tboot_resid_field_120_subj_10_dim);
          signflips_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(signflips_tboot_resid_field_120_subj_10_dim);
          signflips_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(signflips_tboot_resid_field_60_subj_10_dim));
          signflips_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(signflips_tboot_resid_field_60_subj_10_dim);
          signflips_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(signflips_tboot_resid_field_60_subj_10_dim);
          signflips_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(signflips_tboot_resid_field_30_subj_10_dim));
          signflips_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(signflips_tboot_resid_field_30_subj_10_dim);
          signflips_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(signflips_tboot_resid_field_30_subj_10_dim);
          gauss_tboot_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(gauss_tboot_resid_field_240_subj_10_dim));
          gauss_tboot_max_supG_resid_field_240_subj_10_dim(k) = max(gauss_tboot_resid_field_240_subj_10_dim);
          gauss_tboot_min_supG_resid_field_240_subj_10_dim(k) = min(gauss_tboot_resid_field_240_subj_10_dim);
          gauss_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(gauss_tboot_resid_field_120_subj_10_dim));
          gauss_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(gauss_tboot_resid_field_120_subj_10_dim);
          gauss_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(gauss_tboot_resid_field_120_subj_10_dim);
          gauss_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(gauss_tboot_resid_field_60_subj_10_dim));
          gauss_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(gauss_tboot_resid_field_60_subj_10_dim);
          gauss_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(gauss_tboot_resid_field_60_subj_10_dim);
          gauss_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(gauss_tboot_resid_field_30_subj_10_dim));
          gauss_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(gauss_tboot_resid_field_30_subj_10_dim);
          gauss_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(gauss_tboot_resid_field_30_subj_10_dim);
          mammen_tboot_abs_supG_resid_field_240_subj_10_dim(k) = max(abs(mammen_tboot_resid_field_240_subj_10_dim));
          mammen_tboot_max_supG_resid_field_240_subj_10_dim(k) = max(mammen_tboot_resid_field_240_subj_10_dim);
          mammen_tboot_min_supG_resid_field_240_subj_10_dim(k) = min(mammen_tboot_resid_field_240_subj_10_dim);
          mammen_tboot_abs_supG_resid_field_120_subj_10_dim(k) = max(abs(mammen_tboot_resid_field_120_subj_10_dim));
          mammen_tboot_max_supG_resid_field_120_subj_10_dim(k) = max(mammen_tboot_resid_field_120_subj_10_dim);
          mammen_tboot_min_supG_resid_field_120_subj_10_dim(k) = min(mammen_tboot_resid_field_120_subj_10_dim);
          mammen_tboot_abs_supG_resid_field_60_subj_10_dim(k) = max(abs(mammen_tboot_resid_field_60_subj_10_dim));
          mammen_tboot_max_supG_resid_field_60_subj_10_dim(k) = max(mammen_tboot_resid_field_60_subj_10_dim);
          mammen_tboot_min_supG_resid_field_60_subj_10_dim(k) = min(mammen_tboot_resid_field_60_subj_10_dim);
          mammen_tboot_abs_supG_resid_field_30_subj_10_dim(k) = max(abs(mammen_tboot_resid_field_30_subj_10_dim));
          mammen_tboot_max_supG_resid_field_30_subj_10_dim(k) = max(mammen_tboot_resid_field_30_subj_10_dim);
          mammen_tboot_min_supG_resid_field_30_subj_10_dim(k) = min(mammen_tboot_resid_field_30_subj_10_dim);
          
      end
end

%% Saving variables

eval(['save ' SvNm ' nRlz dim smo mag '... 
      'signflips_abs_supG_resid_field_240_subj_124_dim signflips_max_supG_resid_field_240_subj_124_dim signflips_min_supG_resid_field_240_subj_124_dim signflips_abs_supG_resid_field_120_subj_124_dim signflips_max_supG_resid_field_120_subj_124_dim signflips_min_supG_resid_field_120_subj_124_dim signflips_abs_supG_resid_field_60_subj_124_dim signflips_max_supG_resid_field_60_subj_124_dim signflips_min_supG_resid_field_60_subj_124_dim ' ...
      'gauss_abs_supG_resid_field_240_subj_124_dim gauss_max_supG_resid_field_240_subj_124_dim gauss_min_supG_resid_field_240_subj_124_dim gauss_abs_supG_resid_field_120_subj_124_dim gauss_max_supG_resid_field_120_subj_124_dim gauss_min_supG_resid_field_120_subj_124_dim gauss_abs_supG_resid_field_60_subj_124_dim gauss_max_supG_resid_field_60_subj_124_dim gauss_min_supG_resid_field_60_subj_124_dim ' ...
      'mammen_abs_supG_resid_field_240_subj_124_dim mammen_max_supG_resid_field_240_subj_124_dim mammen_min_supG_resid_field_240_subj_124_dim mammen_abs_supG_resid_field_120_subj_124_dim mammen_max_supG_resid_field_120_subj_124_dim mammen_min_supG_resid_field_120_subj_124_dim mammen_abs_supG_resid_field_60_subj_124_dim mammen_max_supG_resid_field_60_subj_124_dim mammen_min_supG_resid_field_60_subj_124_dim ' ...
      'signflips_tboot_abs_supG_resid_field_240_subj_124_dim signflips_tboot_max_supG_resid_field_240_subj_124_dim signflips_tboot_min_supG_resid_field_240_subj_124_dim signflips_tboot_abs_supG_resid_field_120_subj_124_dim signflips_tboot_max_supG_resid_field_120_subj_124_dim signflips_tboot_min_supG_resid_field_120_subj_124_dim signflips_tboot_abs_supG_resid_field_60_subj_124_dim signflips_tboot_max_supG_resid_field_60_subj_124_dim signflips_tboot_min_supG_resid_field_60_subj_124_dim ' ...
      'gauss_tboot_abs_supG_resid_field_240_subj_124_dim gauss_tboot_max_supG_resid_field_240_subj_124_dim gauss_tboot_min_supG_resid_field_240_subj_124_dim gauss_tboot_abs_supG_resid_field_120_subj_124_dim gauss_tboot_max_supG_resid_field_120_subj_124_dim gauss_tboot_min_supG_resid_field_120_subj_124_dim gauss_tboot_abs_supG_resid_field_60_subj_124_dim gauss_tboot_max_supG_resid_field_60_subj_124_dim gauss_tboot_min_supG_resid_field_60_subj_124_dim ' ...
      'mammen_tboot_abs_supG_resid_field_240_subj_124_dim mammen_tboot_max_supG_resid_field_240_subj_124_dim mammen_tboot_min_supG_resid_field_240_subj_124_dim mammen_tboot_abs_supG_resid_field_120_subj_124_dim mammen_tboot_max_supG_resid_field_120_subj_124_dim mammen_tboot_min_supG_resid_field_120_subj_124_dim mammen_tboot_abs_supG_resid_field_60_subj_124_dim mammen_tboot_max_supG_resid_field_60_subj_124_dim mammen_tboot_min_supG_resid_field_60_subj_124_dim '...
      'signflips_abs_supG_resid_field_240_subj_60_dim signflips_max_supG_resid_field_240_subj_60_dim signflips_min_supG_resid_field_240_subj_60_dim signflips_abs_supG_resid_field_120_subj_60_dim signflips_max_supG_resid_field_120_subj_60_dim signflips_min_supG_resid_field_120_subj_60_dim signflips_abs_supG_resid_field_60_subj_60_dim signflips_max_supG_resid_field_60_subj_60_dim signflips_min_supG_resid_field_60_subj_60_dim ' ...
      'gauss_abs_supG_resid_field_240_subj_60_dim gauss_max_supG_resid_field_240_subj_60_dim gauss_min_supG_resid_field_240_subj_60_dim gauss_abs_supG_resid_field_120_subj_60_dim gauss_max_supG_resid_field_120_subj_60_dim gauss_min_supG_resid_field_120_subj_60_dim gauss_abs_supG_resid_field_60_subj_60_dim gauss_max_supG_resid_field_60_subj_60_dim gauss_min_supG_resid_field_60_subj_60_dim ' ...
      'mammen_abs_supG_resid_field_240_subj_60_dim mammen_max_supG_resid_field_240_subj_60_dim mammen_min_supG_resid_field_240_subj_60_dim mammen_abs_supG_resid_field_120_subj_60_dim mammen_max_supG_resid_field_120_subj_60_dim mammen_min_supG_resid_field_120_subj_60_dim mammen_abs_supG_resid_field_60_subj_60_dim mammen_max_supG_resid_field_60_subj_60_dim mammen_min_supG_resid_field_60_subj_60_dim ' ...
      'signflips_tboot_abs_supG_resid_field_240_subj_60_dim signflips_tboot_max_supG_resid_field_240_subj_60_dim signflips_tboot_min_supG_resid_field_240_subj_60_dim signflips_tboot_abs_supG_resid_field_120_subj_60_dim signflips_tboot_max_supG_resid_field_120_subj_60_dim signflips_tboot_min_supG_resid_field_120_subj_60_dim signflips_tboot_abs_supG_resid_field_60_subj_60_dim signflips_tboot_max_supG_resid_field_60_subj_60_dim signflips_tboot_min_supG_resid_field_60_subj_60_dim ' ...
      'gauss_tboot_abs_supG_resid_field_240_subj_60_dim gauss_tboot_max_supG_resid_field_240_subj_60_dim gauss_tboot_min_supG_resid_field_240_subj_60_dim gauss_tboot_abs_supG_resid_field_120_subj_60_dim gauss_tboot_max_supG_resid_field_120_subj_60_dim gauss_tboot_min_supG_resid_field_120_subj_60_dim gauss_tboot_abs_supG_resid_field_60_subj_60_dim gauss_tboot_max_supG_resid_field_60_subj_60_dim gauss_tboot_min_supG_resid_field_60_subj_60_dim ' ...
      'mammen_tboot_abs_supG_resid_field_240_subj_60_dim mammen_tboot_max_supG_resid_field_240_subj_60_dim mammen_tboot_min_supG_resid_field_240_subj_60_dim mammen_tboot_abs_supG_resid_field_120_subj_60_dim mammen_tboot_max_supG_resid_field_120_subj_60_dim mammen_tboot_min_supG_resid_field_120_subj_60_dim mammen_tboot_abs_supG_resid_field_60_subj_60_dim mammen_tboot_max_supG_resid_field_60_subj_60_dim mammen_tboot_min_supG_resid_field_60_subj_60_dim '...
      'signflips_abs_supG_resid_field_240_subj_10_dim signflips_max_supG_resid_field_240_subj_10_dim signflips_min_supG_resid_field_240_subj_10_dim signflips_abs_supG_resid_field_120_subj_10_dim signflips_max_supG_resid_field_120_subj_10_dim signflips_min_supG_resid_field_120_subj_10_dim signflips_abs_supG_resid_field_60_subj_10_dim signflips_max_supG_resid_field_60_subj_10_dim signflips_min_supG_resid_field_60_subj_10_dim ' ...
      'gauss_abs_supG_resid_field_240_subj_10_dim gauss_max_supG_resid_field_240_subj_10_dim gauss_min_supG_resid_field_240_subj_10_dim gauss_abs_supG_resid_field_120_subj_10_dim gauss_max_supG_resid_field_120_subj_10_dim gauss_min_supG_resid_field_120_subj_10_dim gauss_abs_supG_resid_field_60_subj_10_dim gauss_max_supG_resid_field_60_subj_10_dim gauss_min_supG_resid_field_60_subj_10_dim ' ...
      'mammen_abs_supG_resid_field_240_subj_10_dim mammen_max_supG_resid_field_240_subj_10_dim mammen_min_supG_resid_field_240_subj_10_dim mammen_abs_supG_resid_field_120_subj_10_dim mammen_max_supG_resid_field_120_subj_10_dim mammen_min_supG_resid_field_120_subj_10_dim mammen_abs_supG_resid_field_60_subj_10_dim mammen_max_supG_resid_field_60_subj_10_dim mammen_min_supG_resid_field_60_subj_10_dim ' ...
      'signflips_tboot_abs_supG_resid_field_240_subj_10_dim signflips_tboot_max_supG_resid_field_240_subj_10_dim signflips_tboot_min_supG_resid_field_240_subj_10_dim signflips_tboot_abs_supG_resid_field_120_subj_10_dim signflips_tboot_max_supG_resid_field_120_subj_10_dim signflips_tboot_min_supG_resid_field_120_subj_10_dim signflips_tboot_abs_supG_resid_field_60_subj_10_dim signflips_tboot_max_supG_resid_field_60_subj_10_dim signflips_tboot_min_supG_resid_field_60_subj_10_dim ' ...
      'gauss_tboot_abs_supG_resid_field_240_subj_10_dim gauss_tboot_max_supG_resid_field_240_subj_10_dim gauss_tboot_min_supG_resid_field_240_subj_10_dim gauss_tboot_abs_supG_resid_field_120_subj_10_dim gauss_tboot_max_supG_resid_field_120_subj_10_dim gauss_tboot_min_supG_resid_field_120_subj_10_dim gauss_tboot_abs_supG_resid_field_60_subj_10_dim gauss_tboot_max_supG_resid_field_60_subj_10_dim gauss_tboot_min_supG_resid_field_60_subj_10_dim ' ...
      'mammen_tboot_abs_supG_resid_field_240_subj_10_dim mammen_tboot_max_supG_resid_field_240_subj_10_dim mammen_tboot_min_supG_resid_field_240_subj_10_dim mammen_tboot_abs_supG_resid_field_120_subj_10_dim mammen_tboot_max_supG_resid_field_120_subj_10_dim mammen_tboot_min_supG_resid_field_120_subj_10_dim mammen_tboot_abs_supG_resid_field_60_subj_10_dim mammen_tboot_max_supG_resid_field_60_subj_10_dim mammen_tboot_min_supG_resid_field_60_subj_10_dim'...
     ])
