% Script to read the results and produce the good plots

colors = ["#0072BD"; 
          "#0072BD"; 
          "#D95319";
          "#D95319";
          "#EDB120";
          "#EDB120"];
figure;
colororder(colors)

N = 500;
data = readmatrix("data_N"+string(N));
plot_data(data, N)

N = 100;
data = readmatrix("data_N"+string(N));
plot_data(data, N)
hold on

N = 50;
data = readmatrix("data_N"+string(N));
plot_data(data, N)

saveas(gca, "summary.png");

%% Functions

function plot_data(data, N)
    Tvec = data(1,:); Mvec = data(2,:); Evec = data(3,:) / N; % E/N to compare
    Cvec = data(4,:); Xvec = data(5,:);
    leg = ["N=500", "", "N=100", "", "N=50", "Analytic"];
    
    subplot(2,2,1); tit = "Magnetization";
    nice_plot(Tvec, Mvec, tit, "0", 1); hold on;
    plot_numerical(1);
    ylabel("<|m|>");
    h = legend(leg);
    set(h,'FontSize',11);

    subplot(2,2,2); tit = "Energy";
    nice_plot(Tvec, Evec, "", "0", 0); hold on;
    plot_numerical(2);
    nice_plot([], [], tit, "0", 1); hold on;
    ylabel("<E>/N");
    
    subplot(2,2,3); tit = "Specific heat capacity";
    nice_plot(Tvec, Cvec, "", "0", 0); ylim([0,1.5]); hold on;
    plot_numerical(3);
    nice_plot([], [], tit, "0", 1); hold on;
    ylabel("c_V / k_B"); xlabel("T/T_C");
    
    subplot(2,2,4); tit = "Magnetic susceptibility";
    nice_plot(Tvec, Xvec, "", "0", 0); hold on;
    plot_numerical(4);
    nice_plot([], [], tit, "0", 1); hold on; 
    ylabel("\chi_M * (k_B T_C)"); xlabel("T/T_C");

end


function plot_numerical(id)
    data = load("data_numerical_H0.mat");
    data = data.data;
    
    Tvec = data(1,:); Mvec = data(2,:); Evec = data(3,:); %Already normalised
    Cvec = data(4,:); Xvec = data(5,:);
    
    if id == 1
        plot(Tvec, Mvec, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--');
    elseif id == 2
        plot(Tvec, Evec, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--');
    elseif id == 3
        plot(Tvec, Cvec, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--');
    elseif id == 4
        plot(Tvec, Xvec, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--');
        ylim([0,10]);
    end
end

function nice_plot(x, y, tit, color, indicator_title)
%color == "0" : default, indicator_title is 0 or 1
    if color ~= "0"
        plot(x, y, 'LineWidth', 3, 'Color', color); grid on; 
    else
        plot(x, y, 'LineWidth', 3); grid on; 
    end

    xlim([0,2.5]);

    if indicator_title == 1
        h_title = title(tit, 'FontSize', 12);
        
        title_pos = get(gca, 'Title').Position; % Get the current position
        title_pos(2) = title_pos(2) + 0.01; % Increase the vertical position by 0.05
        set(h_title, 'Position', title_pos, 'VerticalAlignment', 'bottom'); % Set the new position
        set(gca,'FontWeight','bold', 'GridLineWidth', 1.4, 'GridAlpha', 0.3)
        fontsize(gca, scale=1)
    end
end