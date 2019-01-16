function Y = create_signal(dim, signal_type, param)

% Generates Signal
% 	Input:
%		dim:				Dimension.
%		signal_type:		Shape of the signal. 2D options are 'ramp' or 'circle'. 3D options are 'sphere'.
%		param:				Parameters of the signal. 
$							If 'signal_type' is 'ramp', then param is either a number X, in which case the signal
%							is a ramp from 0 to X, or a vector [X,Y], in which case the signal a ramp from X to Y in the x direction.
%							If 'signal_type' is 'circle', then param is a vector [M, R], where M is the magnitude of the signal and R is the radius.
%							If 'signal_type' is 'sphere', then param is a vector [M, R], where M is the magnitude of the signal and R is the radius.  
%	Output
%		Y:					Array of size dim containing the signal
%

switch(signal_type)
	case 'ramp'
		if size(param,2) == 1
			Y = repmat(linspace(0, param), dim(2), 1);
		elseif size(param,2) == 2
			Y = repmat(linspace(param(1), param(2)), 1);
		else
			error('param must be a scalar X or vector [X,Y]')
		end
	case 'circle'
		Y = CircularSignal(dim, rag, mag, 0);  
end

