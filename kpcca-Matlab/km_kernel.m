function K = km_kernel(X1,X2,ktype,kpar, w,perm)
%% NOTE: this method has been modified by Laleh Soltan Ghoraie from
%% the original version to pass two extra parameters: weights (w) 
%% and permutation vector (perm) from km_kernel_icd method
%% the 'von' kernel has also been added

% KM_KERNEL calculates the kernel matrix between two data sets.
% Input:	- X1, X2: data matrices in row format (data as rows)
%			- ktype: string representing kernel type
%			- kpar: vector containing the kernel parameters
% Output:	- K: kernel matrix
% USAGE: K = km_kernel(X1,X2,ktype,kpar, w, perm)
%
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2010.
% Id: km_kernel.m v1.0
% This file is part of the Kernel Methods Toolbox (KMBOX) for MATLAB.
% http://sourceforge.net/p/kmbox
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, version 3 (as included and available at
% http://www.gnu.org/licenses).

switch ktype
	case 'gauss'	% Gaussian kernel
		% this calculation introduces a small numerical error w.r.t. to the
		% direct calculation, but it is much faster.
		sgm = kpar;	% kernel width
		
		dim1 = size(X1,1);
		dim2 = size(X2,1);
		
% 		if (dim1*dim2 >= 1E6)	% warning for large computations
% 			fprintf('Warning: computation of %d kernel elements might be slow. Press key to continue?\n',dim1*dim2);
% 			pause
% 		end
		
		norms1 = sum(X1.^2,2);
		norms2 = sum(X2.^2,2);
		
		mat1 = repmat(norms1,1,dim2);
		mat2 = repmat(norms2',dim1,1);
		
		distmat = mat1 + mat2 - 2*X1*X2';	% full distance matrix
		K = exp(-distmat/(2*sgm^2));
		
	case 'gauss-diag'	% only diagonal of Gaussian kernel
		sgm = kpar;	% kernel width
		K = exp(-sum((X1-X2).^2,2)/(2*sgm^2));
		
	case 'poly'	% polynomial kernel
		p = kpar(1);	% polynome order
		c = kpar(2);	% additive constant
		
		K = (X1*X2' + c).^p;

 
     case 'von'  
		dim1 = size(X1,1);
		dim2 = size(X2,1);
		
		mat1 = repmat(X1,dim2,1);
		mat2 = repmat(X2,dim1,1);
		
		xd = mat1-mat2;                
        xdrad = degtorad(xd); 
       %% set right dimension for the kappa vector
        kappa = kpar(1:size(xd,2));

        if (~isempty(perm))
          i = size(w,1) - size(xdrad,1);
          w1w2 = bsxfun(@times,w(perm(i)),w(perm(i+1:end)));
          logw1w2 = log(w1w2);
          a = zeros(size(w1w2,1),1);
          x = cos(xdrad)';          
          for l=1:size(w1w2,1)
              a(l) = (kappa*x(:,l)) +logw1w2(l); 
          end          
          K = exp(a)'; 

         else
          distmat = kappa*cos(xdrad)';
          K = exp(distmat);
         end

        
	case 'von-diag'	   
        kappa = kpar(1:size(X1,2));

       %% note: compute log of weights 
       %% note: we need to compute i (from the km_kernel_icd.m: i=1:m)
        if ((~isempty(w)) && (~isempty(perm)))
          i = size(w,1) - size(X1,1);      
          ww = bsxfun(@times,w(perm(i+1:end)),w(perm(i+1:end))); %% weighted cosine
          logww = log(ww); 
          x = cos(zeros(size(X1,1),size(X1,2)));
          a = zeros(size(ww,1),1);
          for l=1:size(ww,1)
              a(l) = (x(l,:)*kappa') +logww(l); 
          end          
          K = exp(a); 
        else
          K = exp((cos(zeros(size(X1,1),size(X1,2)))*kappa'));             
        end

        
	otherwise	% default case
		error ('unknown kernel type')
end









