function [] = Plotter(saveFig,t,idx_r,idx_m,nu,comparison,par_HBC, emaxSelect)

nr = length(idx_r); nm = length(idx_m); n = nr+nm;

% =======================================================
% -- fig1: plot trust and sets for agent i
close all
col = turbo(nr);
figure('color','w'); 
hold on;
inuc = [1,10,length(nu)]; 
t1 = 1;
for inu = 1:length(inuc)
    sp0 = subplot(2,length(inuc),inu); 
    hold on;
    epsilon = par_HBC.epsilon;
    xm = comparison.x{inuc(inu)}(idx_m,:);
    plot(t(t1:length(t)),xm(:,t1:length(t)),'-','color',[0.9 0.9 0.9],'Marker','.','MarkerSize',2,'LineWidth',.8);
    ylim([-6,6]); 
    xlim([t(1),t(end)]);
    for idx_i = idx_r
        xi = comparison.x{inuc(inu)}(idx_i,:);
        c = fill([t,flip(t)],[xi-epsilon,flip(xi+epsilon)],'r');
        c.EdgeColor = 'none'; 
        c.FaceColor = col(idx_i,:); 
        c.FaceAlpha = .1;
    end
    for idx_i = idx_r
        xi = comparison.x{inuc(inu)}(idx_i,:);
        plot(t(t1:length(t)),xi(t1:length(t)),'-','color',col(idx_i,:),'Marker','.','MarkerSize',5,'LineWidth',.8);
    end
    set(gca,'FontSize',15)
    if inu == 1 
        ylabel(sprintf('states $x(t)$'),'FontSize',15)
    end
    title(['$\nu = ',sprintf('%.2f$',nu(inuc(inu)))],'FontSize',15)
    xlabel('t','FontSize',15); 
end 

% Create a color map as a subplot
sp_bar = subplot(2,length(inuc),[4,5,6]); 
hold on; 
set(gca,'FontSize',15); 

% Width of each color map
w = 1/2;

% Color map of agents
for in = 1:n
    
    % Specify the (x,y) position of the agent's color map
    X = in + w*[-1,1,1,-1];
    Y = w*[-1,-1,1,1];
    % Color it
    c = patch(X,Y,'r');
    
    if in <= nr
        % Cooperative agents
        c.EdgeColor = 'none'; 
        c.FaceColor = col(in,:); 
    else
        % Non-cooperative agents
        c.EdgeColor = [1 1 1]; 
        c.FaceColor = [0.9 0.9 0.9]; 
    end
end

% Set the figure specifications
xlim([1-w,n+w]); 
ylim([-w,w]); 
sp_bar.XTick = 1:n; 
sp_bar.YTickLabel = [];
sp_bar.Position(2) = 0.40; % 0.12; Specify the y position of colormap
sp_bar.Position(4) = 0.05; % 0.05; Specify the height of colormap
xlabel('agents color map','FontSize',15); 

% Save the figure if required
if saveFig
    if emaxSelect == 1
        export_fig('x_plot_eps_050.png','-r600');  
    elseif emaxSelect == 2
        export_fig('x_plot_eps_100.png','-r600');
    elseif emaxSelect == 3
        export_fig('x_plot_eps_150.png','-r600');
    else
        export_fig('x_plot_T_5.png','-r600');
    end 
end

% =========================================================================
% -- fig2: End Behaviour
% =========================================================================
col = turbo(nr);
figure('color','w'); 
hold on;
sp1 = subplot(2,1,1); 
hold on;
for inu = 1:length(nu)
    for idx_i = idx_r
        plot(nu(inu),comparison.x_end(idx_i,inu),'.','color',col(idx_i,:),'MarkerSize',16,'LineWidth',1.1);
    end
    set(sp1,'FontSize',15);    
    xlabel('$\nu$','FontSize',15); 
    ylabel(sprintf('x($t=%d$)',t(end)),'FontSize',15)    
end
sp1.YLim = sp0.YLim;

% Create a color map as a subplot
sp_bar = subplot(2,1,2); 
hold on; 
set(gca,'FontSize',15); 

% Width of each color map
w = 1/2;

% Color map of agents
for in = 1:n
    
    % Specify the (x,y) position of the agent's color map
    X = in + w*[-1,1,1,-1];
    Y = w*[-1,-1,1,1];
    % Color it
    c = patch(X,Y,'r');
    
    if in <= nr
        % Cooperative agents
        c.EdgeColor = 'none'; 
        c.FaceColor = col(in,:); 
    else
        % Non-cooperative agents
        c.EdgeColor = [1 1 1]; 
        c.FaceColor = [0.9 0.9 0.9]; 
    end
end

% Set the figure specifications
xlim([1-w,n+w]); 
ylim([-w,w]); 
sp_bar.XTick = 1:n; 
sp_bar.YTickLabel = [];
sp_bar.Position(2) = 0.40; % 0.12; Specify the y position of colormap
sp_bar.Position(4) = 0.05; % 0.05; Specify the height of colormap
xlabel('agents color map','FontSize',15); 

% Save the figure if required
if saveFig    
    if emaxSelect == 1
        export_fig('x_endplot_eps_050.png','-r600');
    elseif emaxSelect == 2
        export_fig('x_endplot_eps_100.png','-r600');
    elseif emaxSelect == 3
        export_fig('x_endplot_eps_150.png','-r600');
    else
        export_fig('x_endplot_T_5.png','-r600'); 
    end 
end

% =========================================================================
% -- fig3: Trust
% =========================================================================

clear sp
idx_i = 1; 
idx_j = [2,3]; 
comp = [2,idx_m];
col = turbo(nr);
figure('color','w'); 
hold on; 
for idx_c = 1:length(comp)
    Wend = zeros(nr,length(nu));
    sp(idx_c) = subplot(3,2,idx_c); 
    hold on; 
    set(gca,'FontSize',15);
    for inu = 1:length(nu)
        Wend(:,inu)= comparison.W{inu}{length(t)-1}(idx_r,comp(idx_c));        
    end
    for ia = 1:nr
        plot(nu,Wend(ia,:),'-','Marker','.','MarkerSize',10,'LineWidth',.8,'color',col(ia,:));
    end
    if idx_c>2
        xlabel('$\nu$','FontSize',15); 
    end
    ylabel(sprintf('$w_{i%d}(t=200)$',comp(idx_c)),'FontSize',15);
end

% Create a color map
sp_bar=subplot(3,2,[5,6]); 
hold on; 
set(gca,'FontSize',15); 
w = 1/2;
for in = 1:nr
    X = in+w*[-1,1,1,-1];
    Y = w*[-1,-1,1,1];
    c = patch(X,Y,'r'); 
    c.EdgeColor = 'none'; 
    c.FaceColor = col(in,:); 
end
for in = (nr+1):n
    X = in+w*[-1,1,1,-1];
    Y = w*[-1,-1,1,1];
    c = patch(X,Y,'r'); 
    c.EdgeColor = [1 1 1]; 
    c.FaceColor = [0.9 0.9 0.9]; 
end
xlim([1-w,n+w]); 
sp_bar.XTick = 1:n; 
ylim([-w,w]); 
sp_bar.YTickLabel = [];
xlabel('agents','FontSize',15); 
sp_bar.Position(4) = 0.05; 
sp_bar.Position(2) = 0.12; 
%sp2.Position(3)=.95;
for idx_c = length(comp):-1:(length(comp)-1)
    pos = sp(idx_c).Position;
    sp(idx_c).Position(2) = sp_bar.Position(2) + sp_bar.Position(4) + 0.15;
    sp(idx_c).Position(4) = sp(idx_c).Position(4) + (pos(2) - sp(idx_c).Position(2))/2;
    sp(idx_c).Position(1) = (-1)^mod(idx_c,2)*0.01 + sp(idx_c).Position(1);
end
for idx_c = 1:2
    pos = sp(idx_c+2).Position;
    sp(idx_c).Position(2) = pos(2) + pos(4) + 0.1;
    sp(idx_c).Position(4) = pos(4);
    sp(idx_c).Position(1) = (-1)^mod(idx_c,2)*0.01+sp(idx_c).Position(1);
end

% Save the figure if needed
if saveFig
    if emaxSelect == 1
        export_fig('w_endplot_eps_050.png','-r600');
    elseif emaxSelect == 2
        export_fig('w_endplot_eps_100.png','-r600');
    elseif emaxSelect == 3
        export_fig('w_endplot_eps_150.png','-r600');
    else
        export_fig('w_endplot_T_5.png','-r600'); 
    end     
    disp('saved')
end


end