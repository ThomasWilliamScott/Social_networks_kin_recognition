% This script runs the 'simple' model of kin discrimination described in
% Supp. Info. 3. Its outputs are the following summary statistics:
% equilibrium tag frequency (inverse is the number of tags maintained at 
% equilibrium); equilibrium polymorphism; equilibrium fluctuation;
% equilibrium cooperator frequency. Outputs of this script are saved in the 
% as 'Saved_Data' folder as: 'Data_Supp_Fig_6a', 'Data_Supp_Fig_6b', 
% 'Data_Supp_Fig_6c', 'Data_Supp_Fig_6d' & 'Data_Supp_Fig_6e'. Outputs of
% this script are plotted in Supp. Info. Fig. 6.  

clearvars
close all
clc

% Parameter Specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tagSkewIni = 0.9; % Initial tag frequency distribution. Ranges from zero to one. When 0, all tags have roughly equal frequency. When 1, one tag starts at fixation. When intermediate, the dominant tag has some intermediate initial bias.
helpini = 0.1; % Initial cooperator frequency. 
alpha = 1; % Partner recycling.
b=0.3; % Benefit of cooperation.
c=0.1; % Cost of cooperation.
csearch= 0.01; % Cost of partner recycling.
T=100000; % Generations.
thetaR=[0:0.01:0.5]; % Range of theta values (population viscosity).
rR=[0:0.025:0.5]; % Range of r values (recombination).
tag=100; % Maximum number of segregating tags (this is L_{max} in the text).
mu=0.01; % Trait mutation rate.
% The following three parameters are not used in the following iterations.
% They are saved here because they are used in some other scripts where
% single trials are plotted.
thetas = 0.3; % 'Single trial' theta (population viscosity) value. This value is saved, and used in illustrative plots of single trials in Supp. Fig. 6.
rs = 0.3; % 'Single trial' r (recombination) value. This value is saved, and used in illustrative plots of single trials in Supp. Fig. 6.
Tsingle = 20000; % 'Single trial' run time (generations). This value is saved, and used in illustrative plots of single trials in Supp. Fig. 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following 4 empty arrays are populated to give our 4 summary 
% statistics. 
resalt  = zeros(length(rR),length(thetaR)); % equilibrium cooperator frequency
fluc = zeros(length(rR),length(thetaR)); % equilibrium fluctuation
poly  = zeros(length(rR),length(thetaR)); % equilibrium polymorphism
avgtagfreq  = zeros(length(rR),length(thetaR)); % equilibrium tag frequency (the inverse of this gives the number of tags maintained at equilibrium).

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

% The following 7 lines establish the initial population genotype
% frequencies.
maj = (1/tag) * (1 - tagSkewIni) + tagSkewIni ;
popIni = zeros(2,4,tag);
popIni(1,1,:) = rand(1,tag);
popIni(1,1,:) = popIni(1,1,:) ./ ( sum(sum(popIni(1,1,:)))-popIni(1,1,1) ) .* (1-maj)  ;
popIni(1,1,1) = maj ;
popIni(1,1,:) = popIni(1,1,:) .* (1-helpini);
popIni(2,1,:) = (popIni(1,1,:) ./ (1-helpini) ) .* helpini;

for cur_r = 1:length(rR) % obtain results across the range of r values
    
    r = rR(cur_r);
    
for cur_theta = 1:length(thetaR) % obtain results across the range of theta values
    
    theta = thetaR(cur_theta);
    
    pop = zeros(2,4,tag,T+1);
    pop(:,:,:,1) = popIni; % This enters initial genotype frequencies  
    pop(2,2,:,:) = 1;  % This enters trait ID
    
for t=1:T % iterate over generations

% The below entry populates the 3rd column of 'pop' with genotype
% frequencies after selection.
pop(:,3,:,t) = pop(:,1,:,t) .* (1 - csearch .* (1-theta) .* alpha .* (1-sum(pop(:,1,:,t)))  + ... % this line applies the cost of partner search
                  b .*  ((theta .* pop(:,2,:,t) + (1-theta) .* (sum(pop(:,1,:,t) .* pop(:,2,:,t)) )) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) - ...    % this line applies the benefit received from being the recipient of cooperation
                  c .*  pop(:,2,:,t) .* ((theta + (1-theta) .* sum(pop(:,1,:,t)) ) ./ (1-alpha .*(1-sum(pop(:,1,:,t))).*(1-theta))) ); % this line applies the cost of giving cooperation
            
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

% We define the following three variables because we will make use of them
% when calculating our summary statistics.
tagfreqs = sum(pop(:,1,:,round(T/2): T)); % This matrix gives the frequency of each tag, in each generation between T/2 and T (the first T/2 generations are discounted, because equilibrium may not have been reached yet).
propgens = 1/numel(round(T/2): T); % This is the number of generations between T/2 and T.
tagavgtime = mean(tagfreqs,4); % This array gives the average frequency of each tag (averages are taken between generations T/2 and T).

% We then calculate our 4 summary statistics, for each r & theta
% combination.
avgtagfreq(cur_r,cur_theta) = sum(sum((tagfreqs.^2) .* propgens)); % avg tag frequency.
poly(cur_r,cur_theta) = sqrt(sum(tagavgtime.*((tagavgtime - avgtagfreq(cur_r,cur_theta)).^2))); %  avg deviation of 'avg tag freq over time' from mean (measures polymorpshim - tag divergence)
fluc(cur_r,cur_theta) = sqrt(sum(sum(propgens.*tagfreqs.*((tagfreqs- tagavgtime).^2)))); % avg distance of tag from its 'over time mean freq' (measures oscillation)
resalt(cur_r,cur_theta) = mean(sum(sum ( pop(:,1,:,round(T/2) : T) .* pop(:,2,:,round(T/2) : T) ))) ; % equilibrium cooperator frequency

clear tagfreqs propgens tagavgtime

end
end

clear list pop 