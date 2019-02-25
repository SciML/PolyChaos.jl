clearvars

low = 2;
step = 25;
high = 100;

path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/config.txt';
fileID = fopen(path,'w');
fprintf(fileID,'%i\n',[low;step;high]);

%% Fejer rule

for n =low:step:high
    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/fejer',num2str(n),'.txt'];
    fileID = fopen(path,'w');
    uv=flip(fejer(n),1);
    fprintf(fileID,'%12.10f\n',reshape(uv,numel(uv),1));
end

%% Fejer 2 rule

for n =low:step:high
    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/fejer2_',num2str(n),'.txt'];
    fileID = fopen(path,'w');
    uv=fejer2(n);
    fprintf(fileID,'%12.10f\n',reshape(uv,numel(uv),1));
end

%% Clenshaw-Curtis rule

for n =low:step:high
    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/cc',num2str(n),'.txt'];
    fileID = fopen(path,'w');
    uv=clenshaw_curtis(n);
    fprintf(fileID,'%12.10f\n',reshape(uv,numel(uv),1));
end

%% create recurrence coefficients test sets

noCoeffs = high+5;
% hermite
mu = 2;
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/hermiteHigh.txt'];
fileID = fopen(path,'w');
ab=r_hermite(noCoeffs,mu);
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

% logistic
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/logHigh.txt'];
fileID = fopen(path,'w');
    ab=r_logistic(noCoeffs);
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

% jacobi
abLow = 0.0;
abStep = 0.2;
abHigh = 2;
    for al = abLow:abStep:abHigh
        for be = abLow:abStep:abHigh
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/jacHighal',num2str(al,'%1.1f'),'be',num2str(be,'%1.1f'),'.txt'];
fileID = fopen(path,'w');
            ab=r_jacobi(noCoeffs,al,be);
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
        end
    end


% Laguerre
muLag = 0.5;
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/laguerre',num2str(high),'.txt'];
fileID = fopen(path,'w');
        ab=r_laguerre(high,muLag);
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

%endpoints
endPoints = [2;5;10];
path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/endPts.txt';
fileID = fopen(path,'w');
fprintf(fileID,'%f \n',endPoints);

%names
path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/namesConfig.txt';
fileID = fopen(path,'w');
names =["hermite";"log";"jac";"laguerre"];
fprintf(fileID,'%s \n',names);

%% gauss quadrature rule

path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussLog',num2str(low),'.txt'];
fileID = fopen(path,'w');
ab=gauss(low,r_logistic(noCoeffs));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
 
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussLog',num2str(step),'.txt'];
fileID = fopen(path,'w');
ab=gauss(step,r_logistic(noCoeffs));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussLog',num2str(high),'.txt'];
fileID = fopen(path,'w');
ab=gauss(high,r_logistic(noCoeffs));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussHerm',num2str(low),'.txt'];
fileID = fopen(path,'w');
ab=gauss(low,r_hermite(noCoeffs,mu));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));
 
path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussHerm',num2str(step),'.txt'];
fileID = fopen(path,'w');
ab=gauss(step,r_hermite(noCoeffs,mu));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussHerm',num2str(high),'.txt'];
fileID = fopen(path,'w');
ab=gauss(high,r_hermite(noCoeffs,mu));
fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussConfig.txt';
fileID = fopen(path,'w');
names =["log";"hermite"];
fprintf(fileID,'%s \n',names);

%% radau quadrature rule
for i = 1:length(endPoints)
    endpt=endPoints(i);

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(low),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(low,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(step),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(step,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(high),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(high,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(low),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(low,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(step),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(step,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(high),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
      ab=radau(high,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussConfig.txt';
    fileID = fopen(path,'w');
    names =["log";"hermite"];
    fprintf(fileID,'%s \n',names);
    
        endpt=-endPoints(i);

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(low),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(low,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(step),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(step,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauLog',num2str(high),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(high,r_logistic(noCoeffs),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(low),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(low,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(step),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
    ab=radau(step,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = ['/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/radauHerm',num2str(high),'pt',num2str(endpt,'%1.1f'),'.txt'];
    fileID = fopen(path,'w');
      ab=radau(high,r_hermite(noCoeffs,mu),endpt);
    fprintf(fileID,'%12.10f\n',reshape(ab,numel(ab),1));

    path = '/home/ws/uxsng/JuliaDev/PolyChaos.jl/test/dataQuadratureRules/gaussConfig.txt';
    fileID = fopen(path,'w');
    names =["log";"hermite"];
    fprintf(fileID,'%s \n',names);
end

%% radau jacobi