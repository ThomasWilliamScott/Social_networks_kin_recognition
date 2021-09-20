% This script generates data for a single run of the 'simple' model of
% genetic kin discrimination described in Supp. Info. 3. It then plots how
% tag frequencies change over time, alongside the average-over-time tag
% frequencies, and the average-over-time-and-tags tag frequency. This plot 
% is used to illustrate our summary statistics (Supp. Fig. 4).

clearvars
close all
clc

alpha = 0; % Partner recycling.
mu = 0.0005; % Trait mutation rate.
b = 3; % Benefit of cooperation.
c = 1; % Cost of cooperation.
helpini = 0.1; % Initial cooperator frequency.
csearch=0; % Cost of partner recycling.
r=0.26; % Recombination rate.
theta=0.22; % Population viscosity.
tag=2; % Maximum number of tags (this is L_{max} in the text).
T=4000; % Generations
tagSkewIni=0.999; % Initial tag frequency distribution. When 0, tags are approximately equal frequency. When 1, one tag is at fixation. 

% The next 7 lines are to set the initial genotype frequencies.
maj = (1/tag) * (1 - tagSkewIni) + tagSkewIni ;
popIni = zeros(2,4,tag);
popIni(1,1,:) = rand(1,tag);
popIni(1,1,:) = popIni(1,1,:) ./ ( sum(sum(popIni(1,1,:)))-popIni(1,1,1) ) .* (1-maj)  ;
popIni(1,1,1) = maj ;
popIni(1,1,:) = popIni(1,1,:) .* (1-helpini);
popIni(2,1,:) = (popIni(1,1,:) ./ (1-helpini) ) .* helpini;

% The population is tracked using a 4 dimensional matrix 'pop'. The 4th
% dimension denotes the generation (ranging from 1 to T+1). The 3rd 
% dimension denotes the tag (ranging from 1 to tag, where tag is called 
% L_{max} in the text). The 2nd dimention (rows) denotes the trait 
% (first row for cooperators; second row for defectors). The 1st dimension
% (columns) denotes: (column 1) genotype frequency at start of generation;
% (column 2) trait identity (1 if cooperator; 0 if defector; note that
% column 2 is technically superfluous given that trait identity is also
% reflected by the row); (column 3) genotype frequency after selection;
% (column 4) genotype frequency after recombination.

pop = zeros(2,4,tag,T+1); % This empty matrix will be filled up to track the evolutionary process.
pop(:,:,:,1) = popIni; % This enters initial genotype frequencies.  
pop(2,2,:,:) = 1;  % This enters trait ID.

for t=1:T % Iterate over generations.

tagfreqstrack(:,t) =   reshape(sum(pop(:,1,:,t)),[tag,1]); % This matrix is updated each generation to track how tag frequencies change over generations.
    
% The below entry populates the 3rd column of 'pop' with genotype
% frequencies after selection.
pop(:,3,:,t) = pop(:,1,:,t) .* (1 - csearch .* (1-theta) .* alpha .* (1-sum(pop(:,1,:,t)))  + ...   % this line applies the cost of partner search
                  b .*  ((theta .* pop(:,2,:,t) + (1-theta) .* (sum(pop(:,1,:,t) .* pop(:,2,:,t)) )) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) - ... % this line applies the benefit received from being the recipient of cooperation
                  c .*  pop(:,2,:,t) .* ((theta + (1-theta) .* sum(pop(:,1,:,t)) ) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) ); % this line applies the cost of giving cooperation
              
% This next line ensures that genotype frequencies sum to 1. It
% means that (global) competition is incorporated in the fitness function,
% and it stops the proliferation of rounding errors.
pop(:,3,:,t) = pop(:,3,:,t) ./ sum(sum(pop(:,3,:,t))) ; % this is just to stop rounding errors from proliferating over iterations.
  
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

% We define the following three variables because we will make use of them
% when calculating our summary statistics.
tagfreqs = sum(pop(:,1,:,round(T/2): T)); % This matrix gives the frequency of each tag, in each generation between T/2 and T (the first T/2 generations are discounted, because equilibrium may not have been reached yet).
propgens = 1/numel(round(T/2): T); % This is the number of generations between T/2 and T.
tagavgtime = mean(tagfreqs,4); % This array gives the average frequency of each tag (averages are taken between generations T/2 and T).

% We then calculate our 4 summary statistics for this run.
avgtagfreq = sum(sum((tagfreqs.^2) .* propgens)); % avg tag frequency.
poly = sqrt(sum(tagavgtime.*((tagavgtime - avgtagfreq).^2))); %  avg deviation of 'avg tag freq over time' from mean (measures polymorpshim - tag divergence)
fluc = sqrt(sum(sum(propgens.*tagfreqs.*((tagfreqs- tagavgtime).^2)))); % avg distance of tag from its 'over time mean freq' (measures oscillation)
resalt = mean(sum(sum ( pop(:,1,:,round(T/2) : T) .* pop(:,2,:,round(T/2) : T) ))) ; % equilibrium cooperator frequency

% We then generate a figure. This figure will show how tag frequencies
% change over time. It will also plot multiple (coloured) horizontal lines,
% showing the average frequency (over time) of each tag. It will also plot
% a further (grey) line,  showing the average tag frequency (averaged 
% across all individuals, over time).

figure
hold on
for l=1:tag % This for loop plots all tag frequencies over time.
plot(1:T,tagfreqstrack(l,:))
end
hold off

% The following three lines of code plot three horizontal lines. The first
% gives average tag frequency over all individuals. The second gives
% average frequency of tag 1. The third gives average frequency of tag 2.
yline(avgtagfreq,'LineStyle','--','LineWidth',2);
yline(tagavgtime(1),'LineStyle','--','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
yline(tagavgtime(2),'LineStyle','--','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);

% The following 5 lines of code are for formatting the graph.
xlim([1 T])
ylim([0 1])
set(gcf,'color','white')
box off
set(gca,'fontsize', 18)

% Finally, we print the numerical values of our summary statistics for this
% run, which are cited in the legend of Supp. Fig. 4.

avgtagfreq % Equilibrium tag frequency.
poly % Equilibrium polymorphism.
fluc % Equilibrium fluctuation
resalt % Equilibrium cooperator frequency.
