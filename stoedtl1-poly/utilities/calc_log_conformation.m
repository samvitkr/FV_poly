% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function log_conformation = calc_log_conformation(conformationTensor)
    [n1 n2 n3]=size(  conformationTensor.Cxx)

    Cxx = conformationTensor.Cxx;  
    Cyy = conformationTensor.Cyy;
    Czz = conformationTensor.Czz;
    Cxy = conformationTensor.Cxy;
    Cxz = conformationTensor.Cxz;
    Cyz = conformationTensor.Cyz;
    
    Cxxl = 0*Cxx;
    Cyyl = 0*Cyy;
    Czzl = 0*Czz;
    Cxyl = 0*Cxy;
    Cxzl = 0*Cxz;
    Cyzl = 0*Cyz;

    c=zeros(3,3);
    for k =1:n3
%	    tic
	    for j =1:n2
		    %j
		    for i =1:n1
			    %coordi=[i j k]
			    c(1,1)=Cxx(i,j,k);
			    c(1,2)=Cxy(i,j,k);
			    c(1,3)=Cxz(i,j,k);
			    c(2,1)=c(1,2);
                            c(2,2)=Cyy(i,j,k);
                            c(2,3)=Cyz(i,j,k);
			    c(3,1)=c(1,3);
			    c(3,2)=c(2,3);
			    c(3,3)=Czz(i,j,k);
			    lc=real(logm(c));
			    Cxxl(i,j,k) =lc(1,1); 
    			    Cyyl(i,j,k) =lc(2,2); 
    			    Czzl(i,j,k) =lc(3,3); 
    			    Cxyl(i,j,k) =lc(1,2); 
    			    Cxzl(i,j,k) =lc(1,3); 
    			    Cyzl(i,j,k) =lc(2,3); 
		    end
	    end
%	    k 
%	    toc
    end
	log_conformation.Cxx = Cxxl;
	log_conformation.Cyy = Cyyl;
	log_conformation.Czz = Czzl;
	log_conformation.Cxy = Cxyl;
	log_conformation.Cxz = Cxzl;
	log_conformation.Cyz = Cyzl;
end
