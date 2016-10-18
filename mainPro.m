function [] = mainPro(ID)
% main function has been written in C. This is function is now the function
% for generating simulated data and unit sphere
% main function to call the tests
% run for job ID = 1:(plen*deltalen)
global alpha useparfor M Xs Ys N

alpha = 0.05;
useparfor = 0; % need parallel toolbox purchased
verbose = 1;

deltavec = [0, 0.05:0.3:1.90]; 
deltavec0 = [0, linspace(0.01,  1/2, 7)]; %(1-1/sqrt(2))
deltalen = numel(deltavec); 
pvec = [20]; %2:2:10; 
plen = length(pvec);
B = 1e2; %2e3; % replicates of simulated data for calculating power
M = 1e2; %1e3; % bootstrap replicates

if useparfor == 1
    poolsize = matlabpool('size');
    if poolsize == 0
        matlabpool local
    end
end

% translation region: for high performance runing
ch = str2double(num2str(ID)); ch0 = ch;
nChain = 1; % no MCMC !
design = ceil(ch/(nChain*B));
ch = ch - (design-1)*(nChain*B);
rep = ceil(ch/nChain);
ch = ch - (rep-1)*nChain;
fprintf('design = %d, replication = %d, chain = %d:\n', [design, rep, ch])

% set random seeds for reprodicibility of synthetic data
% use randn('state',seeds) and rand('twister',seeds) for older MATLAB version
rng('default')
rng(20)

% save unit sphere
nets.eps = 0.01; % nets.len = 2e3; % nets.eps makes no difference when nets.len <= 1e3
for np = pvec
    nam = strcat('unitsphere',num2str(np),'.mat');
    if ~exist(nam, 'file')
        % N = net(np, nets.eps, 1);
        N = net(np, 6);
        %         if np == 3 % check if the approximation is okay
        %             plot3(N(:,1), N(:,2), N(:,3), '.')
        %             grid on
        %         elseif np == 2
        %             plot(N(:,1), N(:,2),'b.')
        %         end
        save(nam, 'N')
    end
end


error('No need to go futher')


% simulate data
% this part is now moved to R as well
ndesign = 4;
for mydesign = 3:ndesign
    nam = strcat('dat_simu',num2str(mydesign),'.mat');
    n = 50; m = 20; myflag = 0;
    if ~exist(nam, 'file') % need to create the data once
        myflag = 1;
        Xs = cell(1, plen); Ys = cell(1, plen);
        for k = 1:plen
            Xs{k} = nan(n, pvec(k), B, deltalen);
            Ys{k} = nan(m, pvec(k), B, deltalen);
            for j = 1:deltalen
                for i = 1:B
                    switch mydesign
                        case 1 % location alternatives
                            Xs{k}(:, :, i, j) = randn(n, pvec(k)); % std normal
                            tmpY = randn(m, pvec(k)); tmpY(:, 1) = tmpY(:, 1) + deltavec(j);
                        case 2 % dispersion alternatives
                            Xs{k}(:, :, i, j) = randn(n, pvec(k)); % std normal
                            if j == 1
                                tmpY = randn(m, pvec(k));
                            else
                                tmpY = deltavec(j) * randn(m, pvec(k)); % no need to sqrt(delta)
                            end
                        case 3 % dependence alternatives
                            Xs{k}(:, :, i, j) = randn(n, pvec(k)); % std normal
                            A = eye(pvec(k)); 
                            %                             for ii = 1:2 %(pvec(k)-1) %-1 for p=2
                            %                                 tmpU = randn(1, pvec(k));
                            %                                 A(ii, :) = tmpU/norm(tmpU);
                            %                             end
                            A(1, 1:2) = sqrt([1-deltavec0(j), deltavec0(j)]); 
                            A(2, 1:2) = sqrt([deltavec0(j), 1-deltavec0(j)]); 
                            tmpY = A*randn(pvec(k), m); tmpY = tmpY';
                        case 4 % dependence alternatives for multivariate T
                            Xs{k}(:, :, i, j) = mvtrnd(eye(pvec(k)), 5, n); % std normal
                            A = eye(pvec(k));
                            %                             for ii = 1:(pvec(k)-1) %-1 for p=2
                            %                                 tmpU = randn(1, pvec(k));
                            %                                 A(ii, :) = tmpU/norm(tmpU);
                            %                             end
                            A(1, 1:2) = sqrt([1-deltavec0(j), deltavec0(j)]); 
                            A(2, 1:2) = sqrt([deltavec0(j), 1-deltavec0(j)]); 
                            tmpY = A*mvtrnd(eye(pvec(k)), 5, m)'; tmpY = tmpY';
                    end
                    Ys{k}(:, :, i, j) = tmpY;
                end
            end
        end
        save(strcat('dat_simu',num2str(mydesign),'.mat'), 'Xs', 'Ys','pvec','deltavec')
    end
    if myflag == 1 && mydesign == ndesign
        error('Stopped intendedly as the synthetic data were created')
    end
end


load(strcat('dat_simu',num2str(design),'.mat'))

% run the full simulation with B replicated data sets
rng(ch0*11); % set random seeds for reproducibility (no effect if using parfor)
sigs = nan(deltalen, plen);

for np = 1:plen
    load(strcat('unitsphere',num2str(pvec(np)),'.mat'))
    for ndelta = 1:deltalen
        if verbose == 1
            fprintf('(%1d,%2d) ', [np, ndelta])
            if ndelta == deltalen
                fprintf('\n')
            end
        end
        
        
        X_d = Xs{np}(:, :, rep, ndelta);
        Y_d = Ys{np}(:, :, rep, ndelta);
        
        rng('default')
        rng(20); Es = randn(n, M);save('Es.mat','Es')
        
        tic
        [out] = proposed(X_d, Y_d, Es);
        disp([out.stats, mean(out.T1)])
        toc
        error('stop here')
        sigs(ndelta, np) = out.sig; % record the significance
    end
end
toc

save(strcat('out',num2str(design),'_',num2str(rep),'.mat'), 'sigs')

if useparfor == 1
    matlabpool close
end

end

function [out] = proposed(X_d, Y_d, Es)
% function for the proposed smoothing methods
global alpha useparfor M N
d = 4;

% % use randn('state',seeds) and rand('twister',seeds) for older MATLAB version
% rng('default')
% rng(20) % set random seeds

[n, p] = size(X_d);
m = size(Y_d, 1);
len = size(N, 1);

% tic
for i = 1:len
    u = N(i,:)';
    X = X_d*u;
    Y = Y_d*u;
    % Z = ksdensity(X, Y, 'kernel','epanechnikov','function','cdf');
    Z = mean(repmat(X', [m,1]) <= repmat(Y, [1, n]), 2)';
    Wavelet_Z = zeros([m, d]);
    for j = 1:d
        Wavelet_Z(:,j) = double(Fourier(Z,j));
    end
    Phi1 = sum(Wavelet_Z, 1);
    T = max(abs(Phi1))/sqrt(m);
    if i == 1
        T0 = T; u0 = u;
    else
        if T > T0
            T0 = T; u0 = u;
        end
    end
end
% toc
out.stats = T0*sqrt(n/(m+n));

% tic
T1 = nan(1, M);
if useparfor == 1
    parfor k = 1:M
        e = Es(:,k); %e = randn(n, 1);
        for i = 1:len
            u = N(i,:)';
            X = X_d*u;
            Z = mean(repmat(X', [n,1]) <= repmat(X, [1, n]), 2)';
            Wavelet_Z = zeros([n, d]);
            for j = 1:d
                Wavelet_Z(:,j) = double(Fourier(Z,j));
            end
            Phi1 = sum( repmat(e, [1,d]) .* Wavelet_Z, 1); % multiplier
            parT = max(abs(Phi1))/sqrt(n);
            if i == 1
                parT0 = parT;
            else
                if parT > parT0
                    parT0 = parT;
                end
            end
        end
        T1(k) = parT0;
    end
else
    for k = 1:M
        e = Es(:,k); %randn(n, 1);
        for i = 1:len
            u = N(i,:)';
            X = X_d*u;
            Z = mean(repmat(X', [n,1]) <= repmat(X, [1, n]), 2)';
            Wavelet_Z = zeros([n, d]);
            for j = 1:d
                Wavelet_Z(:,j) = double(Fourier(Z,j));
            end
            Phi1 = sum( repmat(e, [1,d]) .* Wavelet_Z, 1); % multiplier
            parT = max(abs(Phi1))/sqrt(n);
            if i == 1
                parT0 = parT;
            else
                if parT > parT0
                    parT0 = parT;
                end
            end
        end
        T1(k) = parT0;
    end
end
% toc

% % if you want to check how empirical distribution of Phi_MB compare to
% % Phi_MS, set a break point below and uncomment the following
% hist(T1)
% line([out.stats out.stats],get(gca,'YLim'),'Color',[1 0 0])

out.T1 = T1; 
% out.CV = quantile(T1, 1-alpha);
% out.sig = (out.stats >= out.CV); % =1: --> H_0: F=G rejected

end

% function [ y ] = Fourier(x,k)
% %Normalized Fourier Base function,
% %input value x and order k, output the result y
% 
% y=sqrt(2)*cos(x*pi*k);
% end

function N = net_old(p,eps,type)
% used for first submission

% net   gernate an eps-net of the unit sphere with cardinality T in random
% N = net(p,eps,T) all points on the unit sphere are collected as rows in N

% % randn('state',sum(100*clock))
% % set the random seeds in the main function for reproducibility
% if T > 1000
%     T = round((1/eps)^p); %round((1+2/eps)^(p));
% end

T = min(round(1/eps)^(p), 2e4); % maximally consider 20,000


if type==1
    N = randn(T,p);
else
    N = laprnd(0,1,T,p);       % Double exponential
end

for t = 1:T
    N(t,:) = N(t,:)/norm(N(t,:));
end

N = [N; eye(p)]; % add the unit direction

if p == 3;
    x = N(:,1);
    y = N(:,2);
    z = N(:,3);
    plot3(x,y,z,'.')
    grid on
end

end


