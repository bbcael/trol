clear all; close all; clc; load trol.mat; % load time-series

%% 4x

for k = 1:26; % for all models

f = ones(150,1).*F2x(k)';
t = T4x(:,k)';

c10 = 5; % initial parameter values
c20 = 50;
l10 = 1;
d0 = 0.5;
c100 = 5;
c200 = 50;
l100 = 1;
d00 = .5;

nb = 1000; % parameter search iterations
snr = .1; % characteristic jump size in parameter space

r0 = Inf;

for i = 1:nb;

    c1b = abs(c10+snr.*randn(1).*c100); % randomly perturb parameters
    c2b = abs(c20+snr.*randn(1).*c200);
    l1b = abs(l10+snr.*randn(1).*l100);
    db = mod(d0+snr.*randn(1).*d00,1);

    T1 = 0; % initialise at zero
    T2 = 0;

    % put any model equations in here to forward-step
    for j = 2:151;
        T1(j) = T1(j-1) + f(j-1)./c1b - l1b.*T1(j-1)./c1b;
        T2(j) = T2(j-1) + f(j-1)./c2b;
    end
    T = (1-db./2).*T1 + (db./2).*T2;
    T = T(2:end);

    rb = sqrt(sum((T-t).^2)./length(T)); % calculate error

    if rb<=r0; % update parameters if better fit
        c10 = c1b;
        c20 = c2b;
        l10 = l1b;
        d0 = db;
        r0 = rb;
    end

end
R4C1(k) = c10;
R4C2(k) = c20;
R4L1(k) = l10;
R4D(k) = d0;
R4R(k) = r0;
[k 1]
end

%% 1%

for k = 1:26;

f = linspace(0,150/70,150).*F2x(k)';
t = T1pct(:,k)';

c10 = 5; % initial parameter values
c20 = 50;
l10 = 1;
d0 = 0.5;
c100 = 5;
c200 = 50;
l100 = 1;
d00 = .5;

nb = 1000; % parameter search iterations
snr = .1; % characteristic jump size in parameter space

r0 = Inf;

for i = 1:nb;

    c1b = abs(c10+snr.*randn(1).*c100); % randomly perturb parameters
    c2b = abs(c20+snr.*randn(1).*c200);
    l1b = abs(l10+snr.*randn(1).*l100);
    db = mod(d0+snr.*randn(1).*d00,1);

    T1 = 0; % initialise at zero
    T2 = 0;

    % put any model equations in here to forward-step
    for j = 2:151;
        T1(j) = T1(j-1) + f(j-1)./c1b - l1b.*T1(j-1)./c1b;
        T2(j) = T2(j-1) + f(j-1)./c2b;
    end
    T = (1-db./2).*T1 + (db./2).*T2;
    T = T(2:end);

    rb = sqrt(sum((T-t).^2)./length(T)); % calculate error

    if rb<=r0; % update parameters if better fit
        c10 = c1b;
        c20 = c2b;
        l10 = l1b;
        d0 = db;
        r0 = rb;
    end

end

R1C1(k) = c10;
R1C2(k) = c20;
R1L1(k) = l10;
R1A(k) = d0;
R1R(k) = r0;
[k 2]
end

%% historical

for k = 1:2237;

f = Fh(:,k)';
t = Th(:,mod(k,200)+1)';
f = f(101:end)-mean(f(1:100));
t = t(1:end-2) - mean(t(1:50));

c10 = 5; % initial parameter values
c20 = 50;
l10 = 1;
d0 = 0.5;
c100 = 5;
c200 = 50;
l100 = 1;
d00 = .5;

nb = 1000; % parameter search iterations
snr = .1; % characteristic jump size in parameter space

r0 = Inf;

for i = 1:nb;

    c1b = abs(c10+snr.*randn(1).*c100); % randomly perturb parameters
    c2b = abs(c20+snr.*randn(1).*c200);
    l1b = abs(l10+snr.*randn(1).*l100);
    db = mod(d0+snr.*randn(1).*d00,1);

    T1 = 0; % initialise at zero
    T2 = 0;

    % put any model equations in here to forward-step
    for j = 2:171;
        T1(j) = T1(j-1) + f(j-1)./c1b - l1b.*T1(j-1)./c1b;
        T2(j) = T2(j-1) + f(j-1)./c2b;
    end
    T = (1-db./2).*T1 + (db./2).*T2;
    T = T(2:end);

    rb = sqrt(sum((T-t).^2)./length(T)); % calculate error

    if rb<=r0; % update parameters if better fit
        c10 = c1b;
        c20 = c2b;
        l10 = l1b;
        d0 = db;
        r0 = rb;
    end

end

RhC1(k) = c10;
RhC2(k) = c20;
RhL1(k) = l10;
DhR(k) = d0;
RhR(k) = r0;
[k 3]
end

clearvars -EXCEPT R*