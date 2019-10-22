% Reflective Index
n1 = 1.5;    n2 = 1.45;
c = 3e8;

% a planar dielectric waveguide of thickness d = 2a
a = 10e-6; % Condition 1
% a = 20e-6; % Condition 2

vp = c/n1;
cutoff = @(m) (pi*vp*m)./(2*a*sqrt(n1^2-n2^2));
% beta = @(w,m) w/vp.*sqrt(1-(cutoff(m)./w).^2);

step = 2e11;
endpoint = 3e14;

figure(1)
for mode = 0:2
    ww = cutoff(mode):step:endpoint;
    th = 90:-0.001:0;
    RH = real(sqrt((sind(th)).^2-(n2/n1)^2)./cosd(th));
    idx = find(RH == 0, 1);
    th = th(1:idx);
    RH = RH(1:idx);
    
    theta = zeros(size(ww));
    beta = zeros(size(ww));
    omega = zeros(size(ww));
    th_idx = 1;
    for w = ww
        temp = a*w/vp.*cosd(th)-mode*pi/2;
        idx = temp>=0 & temp<=pi/2;
        if isempty(idx)
            continue
        end
        temp = temp(idx);
        X = th(idx);
        RHX = RH(idx);
        LF = tan(temp);
        idx = find(RHX-LF<=0, 1);
        if isempty(idx)
            continue
        end
        theta(th_idx) = X(idx);   
        beta(th_idx) = w/vp*sind(theta(th_idx));
        omega(th_idx) = w;
        th_idx = th_idx+1;
    end
    if th_idx > 1
        beta = beta(1:th_idx-1);
        omega = omega(1:th_idx-1);
        plot(beta, omega);
        hold on
        disp(mode);
    end
end
ww = 0:step:endpoint;
slope1 = @(w) n2/c*w;
slope2 = @(w) n1/c*w;
plot(slope1(ww), ww, 'k:');
plot(slope2(ww), ww, 'k:');
hold off
xlabel('\beta_m');
ylabel('\omega');
%%


figure(1)
hold on
plot(slope1(w), w, slope2(w), w);
% hold on
hold off
