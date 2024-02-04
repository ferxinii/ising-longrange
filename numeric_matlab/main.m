% Code to explore numerically the equations arising in the super long-range
% ising model. The plots produced are the ones I used in the project.

global T H
% Attention!! variables T and H are global!!!

addpath(genpath('Newton'))

Tvec = linspace(0.1, 2.5, 1000);


%% m, cV, xM, E for zero field

H=0;

Mvec = []; Cvec=[]; Xvec=[]; Evec=[];
ii = 0;
for T = Tvec
    m = abs(m_general()); %Note that T and H are GLOBALS!!
    cv = cv_general(T, H);
    xm = xm_general(T, H);
    E = energy(m, T, H);
    
    Mvec = [Mvec m];
    Cvec = [Cvec cv];
    Xvec = [Xvec xm];
    Evec = [Evec E];
    ii = ii + 1 %This is an indicator
end

fig = figure();
tit = "Mean magnetization per spin |m|, as a function of T/T_C";
nice_plot(Tvec, Mvec, tit, "0", 1);
xlabel("T / T_C"); ylabel("|m|");
saveas(gcf,'m.png');

fig = figure;
tit = "Specific heat capacity (reduced) as a function of T/T_C";
color = "#0072BD";
nice_plot(Tvec(1:375), Cvec(1:375), tit, color, 0); hold on;
nice_plot(Tvec(376:end), Cvec(376:end), tit, color, 1);
xlabel("T / T_C"); ylabel("c_V / k_B"); 
saveas(gcf,'cv.png');

fig = figure;
nice_plot(Tvec, Xvec, tit, "0", 0);
grid on; xlabel("T / T_C"); ylabel("X_M * (k_B T_C)"); ylim([0,10]);
nice_plot([],[],tit,"0",1);
saveas(gcf,'xm.png');

data = [Tvec; Mvec; Evec; Cvec; Xvec];
save("data_numerical_H0.mat", "data");

%% Helmholtz free energy

mvec = linspace(-1, 1, 100);

% Zero field
H=0;
fig = figure;
for T = [0.6, 0.8, 1.2]
    Avec = [];
    for m = mvec
        A = helmholtz(m, T, H);
        Avec = [Avec A];
    end
    nice_plot(mvec, Avec - Avec(round(length(Avec)/2)), "", "0", 0); hold on;
end
tit = "Helmholtz free energy A(m), zero field";
ylim([-0.3,0.5]);
nice_plot([],[],tit,"0",1)
xlabel("m"); ylabel("(A - A_0)/(N J)");
legend("T(k_B/J)=0.6", "T(k_B/J)=0.8", "T(k_B/J)=1.2");
saveas(gcf,'A_H0.png');

% Non-zero field
H=0.1;
fig = figure;
for T = [0.6, 0.8, 1.2]
    Avec = [];
    for m = mvec
        A = helmholtz(m, T, H);
        Avec = [Avec A];
    end
    nice_plot(mvec, Avec - Avec(round(length(Avec)/2)), "", "0", 0); hold on;
end
tit = "Helmholtz free energy A(m), H/J=0.1";
ylim([-0.3,0.5]);
nice_plot([],[],tit,"0",1)
xlabel("m"); ylabel("(A - A_0)/(N J)");
legend("T(k_B/J)=0.6", "T(k_B/J)=0.8", "T(k_B/J)=1.2");
saveas(gcf,'A_H01.png');



%% Calculate surfaces as a function of (T, H)

Hvec = linspace(-0.5, 0.5, 101);
Tvec = linspace(0, 3, 100);

Msurf=[]; Esurf=[]; Csurf=[]; Xsurf=[];
for H = Hvec
    Mvec=[]; Evec=[]; Cvec=[]; Xvec=[];
    for T = Tvec
        m = m_general();
        Mvec = [Mvec m];
        
        E = energy(m,T,H);
        Evec = [Evec E];

        cv = cv_general(T, H);
        Cvec = [Cvec, cv];
        
        xm = xm_general(T,H);
        Xvec = [Xvec xm];
    end
    Msurf = [Msurf, Mvec'];     Esurf = [Esurf, Evec'];
    Csurf = [Csurf, Cvec'];     Xsurf = [Xsurf, Xvec'];
    H
end


%% 3D PLOTS

[Hsurf, Tsurf] = meshgrid(Hvec, Tvec);

figure;
tit = "surface m(T,H)";
nice_plot_3D(Hsurf, Tsurf, Msurf, tit, 10, 6);
xlabel("H/(k_B T_C)"); ylabel("T/T_C"); zlabel("m"); 

figure;
tit = "surface cv(T,H)";
nice_plot_3D(Hsurf, Tsurf, Csurf, tit, 10, 6);
xlabel("H/(k_B T_C)"); ylabel("T/T_C"); zlabel("cv/k_B"); 

figure;
tit = "surface E(T,H)";
nice_plot_3D(Hsurf, Tsurf, Esurf, tit, 10, 6);
xlabel("H/(k_B T_C)"); ylabel("T/T_C"); zlabel("E/(N J)");

figure;
Xsurf(Xsurf > 10) = 11;
tit = "surface \chi_M(T,H)";
nice_plot_3D(Hsurf, Tsurf, Xsurf, tit, 10, 6);
xlabel("H/(k_B T_C)"); ylabel("T/T_C"); zlabel("\chi_M * (k_B T_C)"); %
% zlim([-1,10]);


%% Functions

function E = energy(m, T, H)
    E = -m*H - 1/2*m^2;
end

function A = helmholtz(m, T, H)
    aux = -1/2*(m^2) - H*m; 
    aux2 = T/2*log( 1/2 * (1-m)^(1-m) .* (1+m).^(1+m) );

    A = aux + aux2;
end

function xm = xm_general(T, H)
    m = m_general();
    b = 1/T;

    num = 1-m^2;
    den = m^2 - 1 + T;
    xm = num/den;
end

function cv = cv_general(T, H)
    m = m_general();
    b = 1/T;
    
    num = 1/T*(1-m^2)*(m^2+H*m);
    den = m^2 - 1 + T;
    cv = num/den;
end

function m = m_general()
    global T H
    %Solves the transcendental equation and returns stable value
    tol = 1e-5;
    itmax = 1000;
    precision = 4; %Number of digits after the comma for m

    N = 100;
    mvec = zeros(1, N);
    ii = 1;
    for x0 = linspace(-1, 1, N)
        [XK, ~, it] = newtonn(x0,tol,itmax,@implicit_m);
        if it ~= itmax 
            % Only save those that have converged
            mvec(ii) = round(XK(end), precision);
        end
        ii = ii+1;
    end
    
    %When there are multiple solutions, the stable is the one with A min
    munique = unique(mvec);
    A = [];
    for mm = munique
        A = [A helmholtz(mm, T, H)];
    end

    if length(munique) > 1
        [~, id] = min(A);
        m = munique(id);
    else
        m = munique;
    end
end

function res = implicit_m(m)
%Implicit function m(T,H), implemented to use Newton. T, H are GLOBAL
    global T H
    res = m - tanh(1/T*(m+H));
end

function nice_plot(x, y, tit, color, indicator_title)
%color == "0" : default, indicator_title is 0 or 1
    if color ~= "0"
        plot(x, y, 'LineWidth', 3, 'Color', color); grid on; 
    else
        plot(x, y, 'LineWidth', 3); grid on; 
    end

    if indicator_title == 1
        h_title = title(tit, 'FontSize', 15);
        
        title_pos = get(gca, 'Title').Position; % Get the current position
        title_pos(2) = title_pos(2) + 0.01; % Increase the vertical position by 0.05
        set(h_title, 'Position', title_pos, 'VerticalAlignment', 'bottom'); % Set the new position
        set(gca,'FontWeight','bold', 'GridLineWidth', 1.4, 'GridAlpha', 0.3)
        fontsize(gca, scale=1.2)
    end
end


function nice_plot_3D(xsurf, ysurf, zsurf, tit, nx, ny)

s = surf(xsurf, ysurf, zsurf);
shading interp;

x=s.XData; y=s.YData; z=s.ZData;
x=x(1,:); y=y(:,1);
xspacing = (round(length(x)/nx)); yspacing = round(length(y)/ny);
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end

h_title = title(tit, 'FontSize', 15);

title_pos = get(gca, 'Title').Position; % Get the current position
title_pos(2) = title_pos(2) + 0.01; % Increase the vertical position by 0.05
set(h_title, 'Position', title_pos, 'VerticalAlignment', 'bottom'); % Set the new position
set(gca,'FontWeight','bold', 'GridLineWidth', 1.4, 'GridAlpha', 0.3)
fontsize(gca, scale=1.2)

end