function f = custom_plot(y_lower,y_upper,qoi_HDM,qoi_ROM,...
    mean_qoi_SROM,t,indicate_dof_type)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
th = 1.1; % thickness of the line
width_plot = 595;
height_plot = width_plot/2.1;
f = figure('Color',[1 1 1],'units','points','position',[0,0,width_plot,height_plot]);
h = fill([t, flip(t)], [y_lower, flip(y_upper)],'c');  % plot filled area
h.FaceColor = '#a6cce3';
h.EdgeColor = "none";
hold on
p1 = plot(t,qoi_HDM,'k','LineWidth',th);
p2 = plot(t,qoi_ROM,'LineWidth',th);
p2.Color = '#ee3a2b';
p3 = plot(t,mean_qoi_SROM,'LineWidth',th);
p3.Color = "#1e78b3";
xlim([0 77])
xlabel('Time (ms)')
if indicate_dof_type ==1
    legend([p1,p2,p3,h],{'HDM','ROM','Mean SROM','SROM 95$\%$ PI'},'location','southwest',Box='off',Interpreter='latex')
    ylim([-9 7.2])
    ylabel('Velocity in X (in/s)')
    % create a new pair of axes inside current figure
    axes('position',[.643 .185 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 40) & (t >= 30); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    xticks(30:5:40)
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.36 .185 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 10) & (t >= 0); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor ='#a6cce3';
    xticks(0:5:10)
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.643 .665 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <=76.8 ) & (t >= 66.8); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color = "#1e78b3";
    axis tight
elseif indicate_dof_type ==2
    ylabel('Acceleration in X ($10^4$ in/s$^2$)', Interpreter='latex')
    legend([p1,p2,p3,h],{'HDM','ROM','Mean SROM','SROM 95$\%$ PI'},'Position',[0.45 0.715 0.1 0.2],Box='off')
    % create a new pair of axes inside current figure
    axes('position',[.643 .185 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 40) & (t >= 30); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.36 .185 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 10) & (t >= 0); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor ='#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.63 .665 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <=76.8 ) & (t >= 66.8); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    % ylim([-3e4 3e4])
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color = "#1e78b3";
    axis tight
elseif indicate_dof_type ==3
    ylabel('Displacement in X ($10^{-3}$ in)',Interpreter='latex')
    legend([p1,p2,p3,h],{'HDM','ROM','Mean SROM','SROM 95$\%$ PI'},'Position',[0.2 0.4 0.1 0.2],Box='off',Interpreter='latex')
    % create a new pair of axes inside current figure
    axes('position',[.643 .5 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 40) & (t >= 30); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.2 .655 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <= 10) & (t >= 0); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor ='#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color =  "#1e78b3";
    axis tight
    % create a new pair of axes inside current figure
    axes('position',[.643 .185 .25 .25])
    box on % put box around new pair of axes
    indexOfInterest = (t <=76.8 ) & (t >= 66.8); % range of t near perturbation
    h = fill([t(indexOfInterest), flip(t(indexOfInterest))], [y_lower(indexOfInterest), flip(y_upper(indexOfInterest))],'y','FaceAlpha',1);  % plot filled area
    h.EdgeColor = "none";
    h.FaceColor = '#a6cce3';
    hold on
    plot(t(indexOfInterest),qoi_HDM(indexOfInterest),'k','LineWidth',th) % plot on new axes
    p = plot(t(indexOfInterest),qoi_ROM(indexOfInterest),'LineWidth',th);
    p.Color = '#ee3a2b';
    p = plot(t(indexOfInterest),mean_qoi_SROM(indexOfInterest),'LineWidth',th);
    p.Color = "#1e78b3";
    axis tight
end
end