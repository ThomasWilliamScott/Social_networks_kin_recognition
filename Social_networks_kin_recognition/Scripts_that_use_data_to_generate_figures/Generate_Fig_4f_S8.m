% This script uses data from the no tag mutation version of the agent-based
% simulation described in Supp. Info. 3c. It plots the ratio of tag 
% fixation time under selection to tag fixation time under neutrality (y 
% axis, for different population sizes (x axis). It does this for different
% parameter regimes, and generates the three plots in Fig S8 (panel a is 
% also used as Fig. 4f).

close all
clearvars

load('../Saved_data/Balancing_sel_finite_pop_data.mat') % This loads the results matrix ('resfix') that stores fixation times for different runs of the agent based simulation.

nIndR = unique(resfix(:,2)) ; % This generates an array for all the different population sizes that feature in the dataset. This range of population sizes will be the x axis of our plots.

for i=[1:numel(unique(resfix(:,2)))]
      
    NN= nIndR(i); % This gives the population size
    
    neu(i) = mean(resfix(resfix(:,1)==0 & resfix(:,2)==NN & resfix(:,4)==0 & resfix(:,6)==0.005 ,5)); % This gives the mean neutral fixation time for a given population size (NN).
    
    sela1(i) = mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==1& resfix(:,4)==0 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=1, csearch=0, c=0.5, for a given population size (NN).
    sela05(i)= mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==0.5& resfix(:,4)==0 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=0.5, csearch=0, c=0.5, for a given population size (NN).
    sela0(i) = mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==0& resfix(:,4)==0 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=0, csearch=0, c=0.5, for a given population size (NN).
    
    sela09csL(i) = mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==0.9&  resfix(:,4)==0.02 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=0.9, csearch=0.02, c=0.5, for a given population size (NN).
    sela09csS(i) = mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==0.9&  resfix(:,4)==0.01 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=0.9, csearch=0.01, c=0.5, for a given population size (NN).
    sela09(i) = mean(resfix(resfix(:,1)==0.5 & resfix(:,2)==NN & resfix(:,3)==0.9& resfix(:,4)==0& resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=0.9, csearch=0, c=0.5, for a given population size (NN).
        
    selMa1(i) = mean(resfix(resfix(:,1)==0.4 & resfix(:,2)==NN & resfix(:,3)==1& resfix(:,4)==0 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=1, csearch=0, c=0.4 (& b=9c, as in all cases), for a given population size (NN). 
    selSa1(i) = mean(resfix(0.2999<resfix(:,1)<0.3001 & resfix(:,2)==NN & resfix(:,3)==1& resfix(:,4)==0 & resfix(:,6)==0.005,5)); % This gives the mean fixation time, when alpha=1, csearch=0, c=0.3 (& b=9c, as in all cases), for a given population size (NN).
    
end   

% The next 8 lines of code serve to relatavise the tag fixation times, by 
% writing them as a fraction of the tag fixation times under neutrality. 
sela1 = sela1 ./ neu;
sela05= sela05./ neu;
sela0 = sela0 ./ neu;
sela09csL= sela09csL./ neu;
sela09csS= sela09csS./ neu;
sela09= sela09./ neu;
selMa1 = selMa1 ./ neu;
selSa1 = selSa1 ./ neu;


% The next 10 lines of code generate and format Fig. S8a (replicated as 
% Fig. 4f).
figure
plot(nIndR,sela1,'-o','Color','k','LineWidth',2)
hold on
plot(nIndR,sela05,'-o','Color','k','LineWidth',2)
plot(nIndR,sela0,'-o','Color','k','LineWidth',2)
yline(1,'--','LineWidth',3,'Color','k')
hold off
xlim([30*7 50*7])
box off
set(gcf,'color','white')


% The next 10 lines of code generate and format Fig. S8b
figure
plot(nIndR,sela09,'-o','Color','k','LineWidth',2)
hold on
plot(nIndR,sela09csL,'-o','Color','k','LineWidth',2)
plot(nIndR,sela09csS,'-o','Color','k','LineWidth',2)
yline(1,'--','LineWidth',3,'Color','k')
hold off
xlim([30*7 50*7])
box off
set(gcf,'color','white')

% The next 10 lines of code generate and format Fig. S8c
figure
plot(nIndR,sela1,'-o','Color','k','LineWidth',2)
hold on
plot(nIndR,selMa1,'-o','Color','k','LineWidth',2)
plot(nIndR,selSa1,'-o','Color','k','LineWidth',2)
yline(1,'--','LineWidth',3,'Color','k')
hold off
xlim([30*7 50*7])
box off
set(gcf,'color','white')