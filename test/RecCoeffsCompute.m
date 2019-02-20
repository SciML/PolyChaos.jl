low = 1;
step = 20;
high = 200;

path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/config.txt';
fileID = fopen(path,'w');
fprintf(fileID,'%i\n',[low;step;high]);

%% hermite coefficients

for n =low:step:high
    for mu = 0.0:0.1:1AME may also start wit
        path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/hermite',num2str(n),'mu',num2str(mu,'%1.1f'),'.txt'];
        fileID = fopen(path,'w');
        ab=r_hermite(n,mu);
        fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
    end
end

%% logistic coefficients

for n =low:step:high
    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/log',num2str(n),'.txt'];
    fileID = fopen(path,'w');
    ab=r_logistic(n);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
end
%% jacobi coefficients

for n =low:step:high
    for al = 0.0:0.2:2
        for be = 0.0:0.2:2
            path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/jac',num2str(n),'al',num2str(al,'%1.1f'),'be',num2str(be,'%1.1f'),'.txt'];
            fileID = fopen(path,'w');
            ab=r_jacobi(n,al,be);
            fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
        end
    end
end

%% jacobi01 coefficients

for n =low:step:high
    for al = 0.0:0.2:2
        for be = 0.0:0.2:2
            path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/jac01',num2str(n),'al',num2str(al,'%1.1f'),'be',num2str(be,'%1.1f'),'.txt'];
            fileID = fopen(path,'w');
            ab=r_jacobi01(n,al,be);
            fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
        end
    end
end

%% laguerre coefficients

for n =low:step:high
    for mu = 0.0:0.1:1
        path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/laguerre',num2str(n),'a',num2str(mu,'%1.1f'),'.txt'];
        fileID = fopen(path,'w');
        ab=r_laguerre(n,mu);
        fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
    end
end

%% meixner pollaczek coefficients

for n =low:step:high
    for al = 0.1:0.2:2
        for be = 0.1:0.2:2
            path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataRecCoeffs/meixpol',num2str(n),'la',num2str(al,'%1.1f'),'phi',num2str(be,'%1.1f'),'.txt'];
            fileID = fopen(path,'w');
            ab=r_meixner_pollaczek(n,al,be);
            fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
        end
    end
end