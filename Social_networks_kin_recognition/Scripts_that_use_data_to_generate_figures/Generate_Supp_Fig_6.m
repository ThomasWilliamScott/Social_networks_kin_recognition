% This uses data saved from iterating our 'Simple model' of kin
% discrimination (Supp. Info. 3). It uses the data to generate 5 figures,
% corresponding to Supp. Fig. 6a, b, c, d & e.


for part = [1:5] % This loop is to load data for parts a, b, c, d & e of Supp. Fig. 6 in turn.
    
    if part == 1
        load ../Saved_data/Data_Supp_Fig_6a.mat
    
    elseif part==2
        load ../Saved_data/Data_Supp_Fig_6b.mat
        
    elseif part==3
        load ../Saved_data/Data_Supp_Fig_6c.mat
        
    elseif part==4
        load ../Saved_data/Data_Supp_Fig_6d.mat
        
    elseif part==5
        load ../Saved_data/Data_Supp_Fig_6e.mat
        
    end

figure
set(gcf,'color','white')


% Panel i: equilibrium number of tags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
subplot(5,1,1) % specify subplot position
imagesc(thetaR,rR,1./avgtagfreq) % plot as heatmap

% The following 12 lines are for formatting the heatmap.
set(gca,'YDir','normal')
axis([min(thetaR) max(thetaR) min(rR) max(rR)])
yspace=min(rR):max(rR)/5:max(rR);
yticks([yspace]);
yticklabels({yspace});
xspace = min(thetaR):max(thetaR)/5:max(thetaR);
xticks([xspace])
xticklabels({xspace})
colormap(gca, summer)
colorbar
caxis([1 tag]);
colorbar('Ticks',[1,1+(tag-1)/2,tag])

hold on
% The following two lines demarcate the regions of parameter space where, 
% respectively: (i) indiscriminate cooperation is favoured over defection;
% (ii) conditional cooperation (kin discrimination) is favoured over
% defection.
line([c/b,c/b],[0,1],'color','k','Linestyle',':','LineWidth',2) % invasion of indiscrim altruist
line([c/(c-c*tag+tag*b),c/(c-c*tag+tag*(b))],[0,1],'color','b','Linestyle',':','LineWidth',2) % invasion of discrim altruist
hold off

% the following line annotates the heatmap to show the r & m combination 
% that the single trial in panel v is based on.
rectangle('Position',[thetas(1) - (( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1))/2  , rs(1) - (( rR(end)-rR(1) ) / (numel(rR)-1))/2 , ( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1) , ( rR(end)-rR(1) ) / (numel(rR)-1)],'EdgeColor','r','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Panel ii: equilibrium tag polymorphism %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,2) % specify subplot position
imagesc(thetaR,rR,poly) % plot as heatmap

% The following 12 lines are for formatting the heatmap.
set(gca,'YDir','normal')
axis([min(thetaR) max(thetaR) min(rR) max(rR)])
yspace=min(rR):max(rR)/5:max(rR);
yticks([yspace]);
yticklabels({yspace});
xspace = min(thetaR):max(thetaR)/5:max(thetaR);
xticks([xspace])
xticklabels({xspace})
colormap(gca, summer)
colorbar
caxis([0 1]);
colorbar('Ticks',[0,0.5,1]);

hold on
% The following two lines demarcate the regions of parameter space where, 
% respectively: (i) indiscriminate cooperation is favoured over defection;
% (ii) conditional cooperation (kin discrimination) is favoured over
% defection.
line([c/b,c/b],[0,1],'color','k','Linestyle',':','LineWidth',2) % invasion of indiscrim altruist
line([c/(c-c*tag+tag*b),c/(c-c*tag+tag*b)],[0,1],'color','b','Linestyle',':','LineWidth',2) % invasion of discrim altruist
hold off

% the following line annotates the heatmap to show the r & m combination 
% that the single trial in panel v is based on.
rectangle('Position',[thetas(1) - (( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1))/2  , rs(1) - (( rR(end)-rR(1) ) / (numel(rR)-1))/2 , ( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1) , ( rR(end)-rR(1) ) / (numel(rR)-1)],'EdgeColor','r','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Panel iii: equilibrium tag fluctuation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,3)  % specify subplot position
imagesc(thetaR,rR,fluc) % plot as heatmap

% The following 12 lines are for formatting the heatmap.
set(gca,'YDir','normal')
axis([min(thetaR) max(thetaR) min(rR) max(rR)])
yspace=min(rR):max(rR)/5:max(rR);
yticks([yspace]);
yticklabels({yspace});
xspace = min(thetaR):max(thetaR)/5:max(thetaR);
xticks([xspace])
xticklabels({xspace})
colormap(gca, summer)
colorbar
caxis([0 1]);
colorbar('Ticks',[0,0.5,1]);

hold on
% The following two lines demarcate the regions of parameter space where, 
% respectively: (i) indiscriminate cooperation is favoured over defection;
% (ii) conditional cooperation (kin discrimination) is favoured over
% defection.
line([c/b,c/b],[0,1],'color','k','Linestyle',':','LineWidth',2) % invasion of indiscrim altruist
line([c/(c-c*tag+tag*b),c/(c-c*tag+tag*(b))],[0,1],'color','b','Linestyle',':','LineWidth',2) % invasion of discrim altruist
hold off

% the following line annotates the heatmap to show the r & m combination 
% that the single trial in panel v is based on.
rectangle('Position',[thetas(1) - (( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1))/2  , rs(1) - (( rR(end)-rR(1) ) / (numel(rR)-1))/2 , ( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1) , ( rR(end)-rR(1) ) / (numel(rR)-1)],'EdgeColor','r','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Panel iv: equilibrium cooperator frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,4) % specify subplot position
imagesc(thetaR,rR,resalt) % plot as heatmap

% The following 12 lines are for formatting the heatmap.
set(gca,'YDir','normal')
axis([min(thetaR) max(thetaR) min(rR) max(rR)])
yspace=min(rR):max(rR)/5:max(rR);
yticks([yspace]);
yticklabels({yspace});
xspace = min(thetaR):max(thetaR)/5:max(thetaR);
xticks([xspace])
xticklabels({xspace})
colormap(gca, summer)
colorbar
caxis([0 1]);
colorbar('Ticks',[0,0.5,1])

hold on
% The following two lines demarcate the regions of parameter space where, 
% respectively: (i) indiscriminate cooperation is favoured over defection;
% (ii) conditional cooperation (kin discrimination) is favoured over
% defection.
line([c/b,c/b],[0,1],'color','k','Linestyle',':','LineWidth',2) % invasion of indiscrim altruist
line([c/(c-c*tag+tag*b),c/(c-c*tag+tag*(b))],[0,1],'color','b','Linestyle',':','LineWidth',2) % invasion of discrim altruist
hold off

% the following line annotates the heatmap to show the r & m combination 
% that the single trial in panel v is based on.
rectangle('Position',[thetas(1) - (( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1))/2  , rs(1) - (( rR(end)-rR(1) ) / (numel(rR)-1))/2 , ( thetaR(end)-thetaR(1) ) / (numel(thetaR)-1) , ( rR(end)-rR(1) ) / (numel(rR)-1)],'EdgeColor','r','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Panel v: plot a single trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before plotting panel v, we first need to generate the data it is based 
% on. We therefore undertake a single run of our 'simple' population 
% genetic model of kin discrimination, and track tag frequencies over time.
% Panel v then plots how the tag frequencies change over time in a single
% run.

T= Tsingle ; % Tsingle is the number of generations displayed in the single trial plot (panel v).
theta = thetas; % thetas is the theta value (population viscosity) used to generate the single trial plot.
r = rs; % rs is the r value (recombination) used to generate the single trial plot.
tagfreqs=zeros(tag,T); % this empty matrix will be populated during the trial run, to track the freqeuncy of each tag over generations.

% The population is tracked using a 4 dimensional matrix 'pop'. The 4th
% dimension denotes the generation (ranging from 1 to T+1). The 3rd 
% dimension denotes the tag (ranging from 1 to tag, where tag is called 
% L_{max} in the main text). The 2nd dimention (rows) denotes the trait 
% (first row for cooperators; second row for defectors). The 1st dimension
% (columns) denotes: (column 1) genotype frequency at start of generation;
% (column 2) trait identity (1 if cooperator; 0 if defector; note that
% column 2 is technically superfluous given that trait identity is also
% reflected by the row); (column 3) genotype frequency after selection;
% (column 4) genotype frequency after recombination.
pop = zeros(2,4,tag,T+1);
pop(:,:,:,1) = popIni; % set inital genotype frequencies
pop(2,2,:,:) = 1;  % fill in column 2 as 1 for cooperators (bottom row) and 0 for defectors (top row).

for t=1:T % iterate over generations

tagfreqs(:,t) =   reshape(sum(pop(:,1,:,t)),[tag,1]); % populate 'tagfreqs' matrix with tag frequencies at the start of the generation.
    
% The below entry populates the 3rd column of 'pop' with genotype
% frequencies after selection.
pop(:,3,:,t) = pop(:,1,:,t) .* (1 - csearch .* (1-theta) .* alpha .* (1-sum(pop(:,1,:,t)))  + ...  % this line applies the cost of partner search
                  b .*  ((theta .* pop(:,2,:,t) + (1-theta) .* (sum(pop(:,1,:,t) .* pop(:,2,:,t)) )) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) - ...   % this line applies the benefit received from being the recipient of cooperation
                  c .*  pop(:,2,:,t) .* ((theta + (1-theta) .* sum(pop(:,1,:,t)) ) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) );  % this line applies the cost of giving cooperation
             

% This next line ensures that genotype frequencies sum to 1. It
% means that (global) competition is incorporated in the fitness function,
% and it stops the proliferation of rounding errors.
pop(:,3,:,t) = pop(:,3,:,t) ./ sum(sum(pop(:,3,:,t))) ; 
  
% The below entry populates the 4th column of 'pop' with genotype
% frequencies after recombination.
pop(:,4,:,t) = pop(:,3,:,t) .* (pop(:,3,:,t) + (sum(pop(:,3,:,t))-pop(:,3,:,t)) + (sum(pop(:,3,:,t),3) - pop(:,3,:,t)) + ...  % individual with both tag and trait mating with individual with one or both of tag and trait 
    (1-r) .* (1-sum(pop(:,3,:,t))-sum(pop(:,3,:,t),3) + pop(:,3,:,t))) + ... % individual with tag and trait mating but not recombining with individual lacking both tag and trait
    r .* (sum(pop(:,3,:,t))-pop(:,3,:,t)) .* (sum(pop(:,3,:,t),3) - pop(:,3,:,t)); % individual with one of tag and trait mating and recombining with individual with the other of the tag / trait

% This next line is just to stop the proliferation of rounding errors.
pop(:,4,:,t) = pop(:,4,:,t) ./ sum(sum(pop(:,4,:,t))) ;

% The below entry populates the 1st column of 'pop' in the next generation 
% (t+1) with genotype frequencies, accounting for trait mutation at the end
% of the previous generation.
pop(1,1,:,t+1) = pop(2,4,:,t) .* (mu) + pop(1,4,:,t) .* (1-mu);  
pop(2,1,:,t+1) = pop(1,4,:,t) .* (mu) + pop(2,4,:,t) .* (1-mu);  

% This next line is just to stop the proliferation of rounding errors.
pop(:,1,:,t+1) = pop(:,1,:,t+1) ./ sum(sum(pop(:,1,:,t+1))) ;

end

% We can now plot 'tagfreqs' in panel v to display how tag frequencies 
% change over generations.
subplot(5,1,5) % specify subplot position
hold on
for l=1:tag
plot(1:T,tagfreqs(l,:))
end
xlim([1 T])
ylim([0 1])
hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%