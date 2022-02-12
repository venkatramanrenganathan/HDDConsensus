function [] = f_fig_nu(saveFig,t,idx_r,idx_m,nu,comparison,par_HBC)

nr = length(idx_r); nm = length(idx_m); n = nr+nm;

% =======================================================
% -- fig1: plot trust and sets for agent i
close all
col = turbo(nr);
figure('color','w'); 
hold on;
inuc = [1,10,length(nu)]; 
t1 = 1;
%find(t==0);
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
sp1 = subplot(2,length(inuc),[(length(inuc)+1):(length(inuc)*2)]); 
hold on;
for inu = 1:length(nu)
    for idx_i = idx_r
        plot(nu(inu),comparison.x_end(idx_i,inu),'.','color',col(idx_i,:),'MarkerSize',16,'LineWidth',1.1);
    end
    set(sp1,'FontSize',15); 
    %xlim([t(1),t(end)])
    xlabel('$\nu$','FontSize',15); 
    ylabel(sprintf('x($t=%d$)',t(end)),'FontSize',15)    
end
sp1.YLim = sp0.YLim;
if saveFig
    export_fig('ex-nu-config-eq.png','-r600')    
end

% =======================================================
% -- fig2: trust

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
if saveFig
    export_fig('ex-nu-W.png','-r600'); 
    disp('saved')
end

% =======================================================
% -- other old **** things 
% figure('color','w');hold on; 
% 
% inu_c = [1,8,length(nu)];
% for ix = 1:3
%     sp(ix)=subplot(2,3,ix); hold on; plot(t,comparison.x{inu_c(ix)}(idx_r,:),'k','LineWidth',1.1);
%     set(gca,'FontSize',15); xlabel('$t$','FontSize',15);
%     if ix ==1 
%         ylabel('$x_i(t)$, $i$ cooperative','FontSize',15)
%     end
%     title(['$\nu = ',sprintf('%.2f$',nu(inu_c(ix)))],'FontSize',15)
% end
% subplot(2,3,[4:6]); hold on
% for inu = 1:length(nu)
%     
%     subplot(221); hold on; plot(t,comparison.x{inu}(idx_r,:),'color',col(inm,:),'LineWidth',1.1);
%     subplot(222); hold on; plot(t,comparison.x{inu}(idx_m ,:),'color',col(inm,:),'LineWidth',1.1);
% 
%     subplot(223); hold on; 
%     plot(nu(inu),comparison.x_end(:,inu),'k.','MarkerSize',20,'LineWidth',1.1);
%     set(gca,'FontSize',15)
%     xlabel('$\nu$','FontSize',15); ylabel(sprintf('configuration at $t=%d$',t(end)),'FontSize',15)    
%     plot(t,abs(comparison_nm.x{inm}(idx_r,:)-comparison_nm.x_end(idx_r,end)),'color',col(inm,:),'LineWidth',1.1);
%     bound = (nm+nr-1)*(1-par_HBC.nu.^(t-par_HBC.T+2))/(1-par_HBC.nu)* eps_R .*(eps_rho.^t);
%     y = abs(comparison.x{inm}(idx_r,:)-comparison.x_end(idx_r,end));
%     bound =eps_R.*(eps_rho.^t);
%     plot(t,bound-y,'color',col(inm,:),'LineWidth',1.1);
% end
% ax = gca; ax.YLim = sp(3).YLim;
% if saveFig
%     export_fig('ex-nu-config-eq.png','-r600')
% end
%xlim([0,t(end)]); xticks([0:20:t(end)]);
%subplot(121);plot(t(1:end-1),mu(1:n,:),'-','Marker','.','MarkerSize',6,'LineWidth',.8);
%subplot(122);plot(t(1:end-1),mu(n,:),'-','Marker','.','MarkerSize',6,'LineWidth',.8);
%str=cell(1,length(nu));
%figure; hold on
%for inu = 1:length(nu)
%     for it = find(t==0):(length(t)-1)
%         comparison.trust{inu}.delta{it}{idx_i}(idx_Ni,:)
%         mu(it,:) = comparison.trust{inu}.mu{it}{idx_i}(idx_Ni);
%         plot(it,mu,'o','LineWidth',1.1);
%     end
    %plot(k,nu(inu).^(t-k),'*-','LineWidth',1.1);
    %str{inu}=sprintf('%.2f',nu(inu));
%end
%legend(str)
% 
%{    
    idx_m = nr+(1:nm);
    subplot(221); hold on; plot(t,comparison.x{inm}(idx_r,:),'color',col(inm,:),'LineWidth',1.1);
    subplot(222); hold on; plot(t,comparison.x{inm}(idx_m ,:),'color',col(inm,:),'LineWidth',1.1);

    subplot(223); hold on; plot(inm,comparison_nm.x_end(:,inm),'o','LineWidth',1.1);

    subplot(224); hold on;
    %plot(t,abs(comparison_nm.x{inm}(idx_r,:)-comparison_nm.x_end(idx_r,end)),'color',col(inm,:),'LineWidth',1.1);
    %bound = (nm+nr-1)*(1-par_HBC.nu.^(t-par_HBC.T+2))/(1-par_HBC.nu)* eps_R .*(eps_rho.^t);
    y = abs(comparison_nm.x{inm}(idx_r,:)-comparison_nm.x_end(idx_r,end));
    bound =eps_R.*(eps_rho.^t);
    plot(t,bound-y,'color',col(inm,:),'LineWidth',1.1);
%} 
% end
end