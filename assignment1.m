% Reflective Index
n1 = 1.5;    n2 = 1.45;

% Wavelength
lambda = 1550e-9;

% wave number
k = 2*pi*n1/lambda;

% a planar dielectric waveguide of thickness d = 2a
% a = 10e-6; % Condition 1
a = 20e-6; % Condition 2

% Combining
% The expression for phi in TIR for the TE mode.
% The waveguide condition for phi.
LF = @(th, mode) a*k*cosd(th)-mode*pi/2;    % shoud func tanget.
RH = @(th) real(sqrt((sind(th)).^2-(n2/n1)^2)./cosd(th));

th = 0:1e-4:90; th = th(1:end-1);
fRH = RH(th);
idx = find(fRH == 0, 1, 'last');
th = th(idx-1:end);
fRH = fRH(idx-1:end);

figure(1)
plt0 = plot(th, fRH);
hold on

mode = 30;
theta = zeros(1,mode+1);
cLF = cell(2,mode+1);
for m = 0:mode
    fLF = LF(th, m);
    idx = find(fLF>=0 & fLF<pi/2);
    if isempty(idx)
        fprintf("break at %d\n", m);
        cLF = cLF(:,1:m);
        theta = theta(:,1:m);
        break
    end
    xLF = th(idx);
    fLF = tan(fLF(idx));
    
    plt = plot(xLF, fLF, 'k');
    if mod(m,2)==1
        plt.LineStyle = '--';
    end
    if m == 0
        even = plt;
    elseif m == 1
        odd = plt;
    end
    cLF(:,m+1) = {xLF; fLF};
    
    % ±³Á¡
    fRH0 = fRH(idx);
    for x = 1:length(idx)
        if fLF(x)-fRH0(x) <= 0
            theta(m+1) = xLF(x);
            plot(xLF(x), fLF(x), 'ro');
            break
        end
    end
end
hold off

beta = k*sind(theta);

legend([plt0, even, odd], {'f(\theta_m)', 'Even Mode', 'Odd Mode'}, 'Location','northwest');
xlim([th(1) 90]); ylim([0 30]);
% title('a = 10\mum')
title('a = 20\mum')
xlabel('\theta_m (deg)')
ylabel('tan(ak_1cos\theta_m-m\pi/2)')

for k = 1:3
    fprintf("mode: %d / ¥è: %.4f / ¥â: %.4f \n", k-1, theta(k), beta(k));
%     gtext(sprintf("%.4f¨¬", theta(k)))
end