function Y = create_resid(observed_data, observed_mean, observed_std, snr_transform)

% Generates Signal
% 	Input:
%		observed_data:		Array of all subjects observed data.
%		observed_mean:		Array of the mean of all subjects data.
%		observed_std:		Array of the standard devations of all subjects data.
%		snr_transform:      A binary value specifying whether the SNR transformation is applied to the residuals for the cohens d case. 
%
%	Output
%		Y:					Vector of standardized residuals, of size (prod(dim), nSubj). 
%

resid = bsxfun(@minus, observed_data, observed_mean);

% Reshaping to save memory
resid = reshape(resid, [prod(dim) nSubj]);
observed_mean = reshape(observed_mean, [prod(dim) nSubj]);
observed_std = reshape(observed_std, [prod(dim) nSubj]);


% Cohens d case
if snr_transform == 1
	cohen_d = observed_mean./observed_std;
	Y = (resid./observed_std - cohen_d/2.*(resid./observed_std))./sqrt(1+cohen_d.^2/2);
else
	Y = spdiags(1./observed_std, 0, prod(dim), prod(dim))*resid;
end 