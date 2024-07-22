% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function log_conformation = calc_log_conformation_par(conformationTensor,nproc)

delete(gcp('nocreate'))
parpool('Processes',nproc)


    [n1 n2 n3]=size(  conformationTensor.Cxx)
	parfor k=1:n3
	end

    Cxxl = 0*conformationTensor.Cxx;  % remove ghost points in x, z
    Cyyl = 0*conformationTensor.Cyy;
    Czzl = 0*conformationTensor.Czz;
    Cxyl = 0*conformationTensor.Cxy;
    Cxzl = 0*conformationTensor.Cxz;
    Cyzl = 0*conformationTensor.Cyz;
    %c=zeros(3,3);
    
    parfor k =1:n3
    %	tic
    	c=zeros(3,3);

	for j =1:n2
		for i =1:n1

			    %coordi=[i j k]
			    c(1,1)=conformationTensor.Cxx(i,j,k);
			    c(1,2)=conformationTensor.Cxy(i,j,k);
			    c(1,3)=conformationTensor.Cxz(i,j,k);
			    c(2,1)=c(1,2);
                            c(2,2)=conformationTensor.Cyy(i,j,k);
                            c(2,3)=conformationTensor.Cyz(i,j,k);
			    c(3,1)=c(1,3);
			    c(3,2)=c(2,3);
			    c(3,3)=conformationTensor.Czz(i,j,k);
			    lc=logm(c);
			    Cxxl(i,j,k) =lc(1,1); 
    			    Cyyl(i,j,k) =lc(2,2); 
    			    Czzl(i,j,k) =lc(3,3); 
    			    Cxyl(i,j,k) =lc(1,2); 
    			    Cxzl(i,j,k) =lc(1,3); 
    			    Cyzl(i,j,k) =lc(2,3); 
		end
	end
%	k
%	toc
    end
    

    log_conformation.Cxx = Cxxl;  % remove ghost points in x, z
    log_conformation.Cyy = Cyyl;
    log_conformation.Czz = Czzl;
    log_conformation.Cxy = Cxyl;
    log_conformation.Cxz = Cxzl;
    log_conformation.Cyz = Cyzl;

delete(gcp('nocreate'))
end
