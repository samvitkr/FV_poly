tic
n = 200;
A = 500;
a = zeros(n);
delete(gcp('nocreate'))
parpool('Processes',48)
parfor i = 1:n
    
end
toc
parfor i = 1:n
    i
    for j =1:10
    a(i) = min(abs(eig(rand(A))));
    end
end
toc
delete(gcp('nocreate'))
