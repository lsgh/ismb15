% kernelized partial canonical correlation analysis
% Author: Laleh Soltan Ghoraie (lsoltang *at* uwaterloo.ca), 2015.
%
% This file uses methods implemented in the Kernel Methods Toolbox
% (KMBOX) for MATLAB by:
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2012.
%==============================================================================================
%% NOTE: 
%% INPUT: this script requires the following files to work properly:
%% list.txt : contains list of PDB Ids. One Id per row
%% dim_'id'.txt: contains number of dihedral angles for side chains of a protein...
%% For side chains of ALA and GLY, 0 is stored in this file  
%% dih_'id'.txt: contains dihedral angle infomration of each structure in the 
%% generated ensemble. One structure per row. Note that the file does not contain 
%% the no-rotamer reidues (ALA and GLY); but for all others it contains 4 dihedral angles
%% i.e., for residue type ILE with two dihedral angles: the first two stored columns are non-zero
%% following with two zero columns ... 
%% 'id'-normWeights.txt: contains normalized weights for structures in the 
%% generated ensemble. The weights are inverse of calculated energies for structures

%% OUTPUT: the residues' pairwise kpcca scores are stored in the following file: 
%% kpcca_'id'.txt
%==============================================================================================


close all; clear all;


%% read a list of protein (PDB) ids
%% read id from the list.txt
ids{1}='';

fid = fopen(strcat('list.txt'));
i = 0;
while(~feof(fid))
    s = fgetl(fid);
    if s
         parts = textscan(s, '%5s');
         i = i+1;
         ids{i} = parts{1};
         disp(parts{1});
    end      
end
fclose(fid);
no = i;


%% loop over all proteins ... 
for i=1:no
    
    id = cell2mat(ids{i});
    %% open the dimensions file (no. of dihedral angles for each residue) and read them in
    dims = load(strcat('dim_',id,'.txt')); 
    resNum = length(dims);    
    
    %% open and read in the dihedral angles/ delete zeros
    dihs = load(strcat('dih_',id,'.txt'));
    dihs(:,all(~any(dihs),1)) = [];    

    %% important checks: do the size-check 
    %% total dihedral angles must be equal to the sum of dimensions from the dim file
    [smplno,dihangls] = size(dihs); 
    assert(dihangls==sum(dims)); 
    %% total number of the structures generated in the ensemble
    assert(smplno==((resNum*(resNum-1))/2));
             
    %% load the structure weights ....  
    weights = load(strcat(id,'-normWeights.txt'));
    weights = weights';
    %% size-check: make sure we have one weight per structure
    assert(size(weights,1)==smplno);

    %% matrix to store the kpcca scores in ...
    canCorrMat = zeros(resNum,resNum);

    %% kernel settings ...    
    reg = 0.0000001; 
    k1{1} = 'von'; 
    k1{2} = [8,8,4,2];
    k2 = k1;
        
    %% settings for parallel computing in matlab
    out = findResource();
    if (matlabpool('size')~=0) 
            matlabpool close;
    end;
    matlabpool(out.ClusterSize);     
    
    %% loop over the residues to compute regressions (for partial correlation) 
    for k=1:resNum
       if dims(k)~=0
           
          %% pick r_i (residue i)
           ri = zeros(smplno,dims(k));
           ri(:,1:dims(k)) = dihs(:,sum(dims(1:k-1))+1:sum(dims(1:k-1))+dims(k));          
           assert(size(ri,2)==dims(k)); 
           
           dim = dims(k);          
           parfor j=k+1:resNum
              if dims(j)~=0                     
                    
                  %% pick r_j (residue j)
                    rj = zeros(smplno,dims(j));
                    rj(:,1:dims(j)) = dihs(:,sum(dims(1:j-1))+1:sum(dims(1:j-1))+dims(j));
                    assert(size(rj,2)==dims(j));

                 
                  %% Generate the R_0
                    R0 = dihs;
                  %% exclude (r_i) and (r_j) from X -> get all others
                    R0(:, [sum(dims(1:k-1))+1:sum(dims(1:k)),sum(dims(1:j-1))+1:sum(dims(1:j))]) = [];                    
                  %% add the ones column
                    R0 = [ones(smplno,1) R0];
                 
                  %% regression r_i to R_0
                    bOLS = R0\ri;
                    E1 = (ri-R0*bOLS);
                    assert(size(E1,2)==dims(k));
                    
                  %% regression x_j to R_0
                    bOLS = R0\rj;
                    E2 = (rj-R0*bOLS);
                    assert(size(E2,2)==dims(j));
                         
                  %% KERNEL CCA
                    [~,~,beta] = km_kcca(E1,E2,k1,k2, reg, 'ICD', 0, weights);

                    canCorrMat(k,j) = beta;

              end
          end
       end
    end
    
     matlabpool close;

    
    
     %% we can make a symmetric canCorrMat ...
     %% but it's not necessary

    
    %% save the matrix in a file
     corf = fopen(strcat('kpcca_', id, '.txt'), 'wt');
     for j=1:resNum
         fprintf(corf, '%d ', canCorrMat(j,:));
         fprintf(corf, '\n');
     end
     fclose(corf);

    
end

display('done')



