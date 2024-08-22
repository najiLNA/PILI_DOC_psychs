function optimalG = displayKSDistance(KSdist,G_range,options)

% Plot the KS distances (Check out this guy's functions on MATLAB central
% https://www.mathworks.com/matlabcentral/fileexchange/?term=profileid:5962292)
% options.x_axis = 0:0.01:1;
% options.handle     = figure(1);
% % options.color_area = [243 169 114]./255;    % Orange theme
% % options.color_line = [236 112  22]./255;
% options.color_area = [128 193 219]./255;    % Blue theme
% options.color_line = [ 52 148 186]./255;
% options.alpha      = 0.5;
% options.line_width = 2;
% options.error      = 'std';

plot_areaerrorbar(KSdist.',options)
ylabel('KS distance')
xlabel('Global Coupling Strength')
set(gca,'fontsize', 14)
ylim([0 0.8])

% Get the optimal G value that minimizes KS distance
KSdist_average = mean(KSdist,2); 
[~,index] = min(KSdist_average);
optimalG = G_range(index);