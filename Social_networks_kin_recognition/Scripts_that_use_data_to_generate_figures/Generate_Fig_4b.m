% This script generates Fig. 4b (Main text). It generates data for the case
% where there are only a maximum of 2 segregating tags (tag=2). It then
% uses this generated data to plot how tag cooperator proportion deviates
% from the average cooperator proportion, and how this changes over time,
% alongside changes in tag frequency.

clearvars
close all
clc

tagSkewIni = 0.99; % This parameter captures the initial tag frequency distribution. It ranges from zero to one. When 0, all tag groups are roughly equal frequency initially. When 1, one tag starts at fixation. When intermediate, the dominant tag has some intermediate initial skew.
helpini = 0.01; % proportion altruists within tag groups at start of run.
alpha = 1; % partner recycling probability
b=0.015; % benefit of cooperation
c=0.005; % cost of cooperation
T = 10000000; % generations
m=0.3; % migration rate
r=0.1; % recombination rate
tag = 2; % maximum number of tags (this is labelled as L_{max} in the text).
mu =  0.0001 ; % trait mutation rate
N = 30 ; % deme size

% the following 7 lines of code set up the initial genotype frequencies;
maj = (1/tag) * (1 - tagSkewIni) + tagSkewIni ;
popIni = zeros(2,4,tag);
popIni(1,1,:) = rand(1,tag);
popIni(1,1,:) = popIni(1,1,:) ./ ( sum(sum(popIni(1,1,:)))-popIni(1,1,1) ) .* (1-maj)  ;
popIni(1,1,1) = maj ;
popIni(1,1,:) = popIni(1,1,:) .* (1-helpini);
popIni(2,1,:) = (popIni(1,1,:) ./ (1-helpini) ) .* helpini;

% 'pop' (below) is a 4 dimensional matrix that is initially empty, but 
% will be populated to track the evolutionary process. The 2 x 4 x tag x 
% T+1 matrix can be understood as follows. The fourth dimension denotes the
% generation (ranging from generation 1 to generation T+1). In a given
% generation, the third dimension of the matrix denotes tag ID. This ranges
% from 1 to 'tag' (where 'tag' is labelled as Lmax in the text). The first
% dimension of the matrix (rows) denotes trait ID. This ranges from 1 
% (defectors) to 2 (helpers). The second dimension of the matrix (columns)
% is populated as follows. The first column denotes genotype frequency at
% the start of the generation. The second column reflects trait ID, and is
% 1 for helpers (i.e. second row of the second column of the matrix is 
% always 1) and zero for defectors (i.e. first row the second column of 
% the matrix is always 0). The third column denotes genotype frequency
% after selection. The forth column denotes genotype frequency after
% recombination. The matrix therefore tracks all genotype frequencies for
% each generation (1 to T+1) and at each lifecycle stage within each
% generation.
    
pop = zeros(2,4,tag,T+1); 
pop(:,:,:,1) = popIni; % This populates initial genotype frequencies.
pop(2,2,:,:) = 1;  % As mentioned above, the second row of the second column is set to one (indicating helping).

% gamma (below) is the probability that 3 individuals coalesce at 2 loci
% (see Supp. Info. 4b).
gamma = ((-1+m)^3*(-12*(-1+r)^3-16*m^7*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*m^8*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^6*(-2+N)*(-1+N)^2*(-1+r)*(168*(-1+r)^2+N*(-221+(448-225*r)*r)+N^2*(55+56*(-2+r)*r))-m^5*(-1+N)^2*(-1+r)*(-672*(-1+r)^2+N^3*(106+r*(-225+113*r))+4*N*(299+r*(-616+311*r))-2*N^2*(321+r*(-671+338*r)))+m^4*(-1+N)*(-1+r)*(840*(-1+r)^2+N^3*(-883+(1951-989*r)*r)+N^4*(122+r*(-281+143*r))-10*N*(229+r*(-476+241*r))+N^2*(2209+r*(-4750+2417*r)))+m^2*(-1+N)*(336*(-1+r)^3-N^3*(-1+r)*(260+r*(-756+397*r))-2*N*(-1+r)*(431+r*(-952+491*r))+N^2*(-1+r)*(751+r*(-1884+989*r))+N^4*(-24+r*(122+r*(-158+57*r))))-m^3*(-1+N)*(672*(-1+r)^3-4*N^3*(-1+r)*(159+2*r*(-193+99*r))-8*N*(-1+r)*(223+r*(-476+243*r))+4*N^2*(-1+r)*(415+r*(-947+488*r))+N^4*(-80+r*(298+r*(-334+115*r))))+m*(96*(-1+r)^3-4*N^5*r*(5+2*r*(-5+2*r))+22*N^2*(-1+r)*(19+r*(-49+26*r))-4*N*(-1+r)*(83+r*(-184+95*r))-N^3*(-1+r)*(230+r*(-743+403*r))+N^4*(-48+r*(278+r*(-364+131*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r))))));
% phi (below) is the probability that 2 individuals coalesce at 2 loci
% (see Supp. Info. 4b).
phi = -(((-1+m)^2*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+35*N^3*(-1+r)^2+N^2*(-211-212*(-2+r)*r)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^3*(-1+r)*(-76+77*r)+N^2*(-466+(945-472*r)*r)+14*N*(63+r*(-129+65*r)))+m^5*(-2+N)*(-1+N)*(-1+r)*(756*(-1+r)^2+2*N^4*(-1+r)*(-49+52*r)+N^3*(-731-753*(-2+r)*r)-42*N*(48+r*(-99+50*r))+N^2*(1893+2*r*(-1977+998*r)))-m^4*(-1+N)*(-1512*(-1+r)^3+2*N^5*(-1+r)^2*(-36+43*r)+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+2*r*(-1439+732*r))+N^3*(-1+r)*(2932+r*(-6321+3194*r))+N^4*(754+r*(-2371+(2441-823*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+8*N^5*(-1+r)^2*(-3+5*r)+140*N*(-1+r)*(21+r*(-45+23*r))-2*N^2*(-1+r)*(1592+3*r*(-1203+619*r))+N^3*(-1+r)*(1583+2*r*(-1837+938*r))+N^4*(348+r*(-1168+(1253-429*r)*r)))+m^2*(-432*(-1+r)^3-8*N^6*(-1+r)^2*r+N^3*(-1+r)*(1642+15*r*(-279+146*r))+12*N*(-1+r)*(135+r*(-291+149*r))-2*N^2*(-1+r)*(1173+r*(-2748+1429*r))+N^4*(556+r*(-2097+(2362-815*r)*r))+N^5*(-72+r*(304+r*(-374+135*r))))+m*(108*(-1+r)^3-N^3*(-1+r)*(293-924*r+500*r^2)-2*N^5*r*(12+5*r*(-5+2*r))-4*N*(-1+r)*(96+r*(-213+110*r))+N^2*(-1+r)*(503+r*(-1292+685*r))+N^4*(-66+r*(360+r*(-465+167*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r)))))));
% phi (below) is the probability that 2 individuals coalesce at 1 locus
% (see Supp. Info. 4b).
F = (1-m)^2 / ( (1-m)^2 + N*(1-(1-m)^2)) ;

tagfreq = zeros(2,T); % This matrix is used to track the frequency of each of the two tags over generations.
tagfreq(:,1) = reshape(sum(pop(:,1,:,1) ),[2,1]) ; % This gives the frequencies of the two tags initially.

deviationsWeight = zeros(2,T); % This matrix is used to track the deviations of each tag's cooperator proportion from the population average cooperator proportion. The deviation is weighted by tag frequency.
deviationsWeight(:,1) = reshape((pop(2,1,:,1) ./ sum(pop(:,1,:,1)) - sum(pop(2,1,:,1)) ).*sum(pop(:,1,:,1)),[2,1]); % This gives the initial deviations.

for t=1:T % This for loop iterates the population genetic model given in the main text and Supp. Info. 4b, over T generations. 

fi = sum(pop(:,1,:,t) ); % This array gives the population frequency of each tag
pi = sum(pop(:,1,:,t) .* pop(:,2,:,t)) ./ fi; % This array gives the cooperator proportion for each tag (i.e. it gives, for each tag, the proportion of individuals who are cooperators).
denom = (1 - alpha .* (1-fi) .* (1-F)); % This variable (short for denominator) is defined as it is used in multiple calculations below.
M1 = ((F+(1-F).* fi) ./ denom); % This array gives the per-generation probability of obtaining a social interaction.
M2 = ((pop(:,2,:,t) .* (phi + (F-phi).* fi) + (F-phi) .* pi + (1-2*F+phi) .* pi .* fi) ./ denom); % This array gives the per-generation probability of receiving help.

% The below entry (into column 3) gives genotype frequencies after 
% selection (but before recombination and mutation).
pop(:,3,:,t) = pop(:,1,:,t) .* ( 1  + ... 
                  b .*  M2   ...
                - c .* M1 .* pop(:,2,:,t)   ... 
                - (1-m) .* (b-c)  .* ( (M1.*pop(:,2,:,t))/N  + M2/N + (N-2)/N .* ... 
               (gamma .* pop(:,2,:,t) + (F-gamma) .* pi + ...
               + (F-gamma) .* pop(:,2,:,t) .* sum((fi .* (fi + (1-fi).*alpha.*M1))) ...
               + (1-2*F+gamma) .* sum((fi .* pi .* (fi + (1-fi).*alpha.*M1))))) ...
               - m .* (b-c) .*  sum((fi .* pi .* (fi + (1-fi).*alpha.*M1))));
 
% The below calculation (dividing all genotype frequencies by the sum of 
% genotype frequencies) serves to ensure that genotype frequencies sum to 
% 1, which prevents the proliferation of rounding errors over successive 
% generations.          
pop(:,3,:,t) = pop(:,3,:,t) ./ sum(sum(pop(:,3,:,t))) ;
  
% The below entry (into column 4) gives genotype frequencies after 
% recombination (but before mutation).
pop(:,4,:,t) =  pop(:,3,:,t) .* (2*F-phi) .* (1-m)^2 + (1-(2*F-phi).*(1-m)^2) .* ...
    (pop(:,3,:,t) .* (pop(:,3,:,t) + (sum(pop(:,3,:,t))-pop(:,3,:,t)) + (sum(pop(:,3,:,t),3) - pop(:,3,:,t)) + ...  % individual with both tag and trait mating with individual with one or both of tag and trait 
    (1-r) .* (1-sum(pop(:,3,:,t))-sum(pop(:,3,:,t),3) + pop(:,3,:,t))) + ... % individual with tag and trait mating but not recombining with individual lacking both tag and trait
    r .* (sum(pop(:,3,:,t))-pop(:,3,:,t)) .* (sum(pop(:,3,:,t),3) - pop(:,3,:,t))); % individual with one of tag and trait mating and recombining with individual with the other of the tag / trait

% The below calculation (dividing all genotype frequencies by the sum of 
% genotype frequencies) serves to prevent rounding errors (as above). 
pop(:,4,:,t) = pop(:,4,:,t) ./ sum(sum(pop(:,4,:,t))) ;

% The below entry gives genotype frequencies after trait mutation. This is
% entered into column 1 for the next generation (t+1).
pop(1,1,:,t+1) = pop(2,4,:,t) .* (mu) + pop(1,4,:,t) .* (1-mu);  
pop(2,1,:,t+1) = pop(1,4,:,t) .* (mu) + pop(2,4,:,t) .* (1-mu);  

% The below calculation (dividing all genotype frequencies by the sum of 
% genotype frequencies) serves to prevent rounding errors (as above).
pop(:,1,:,t+1) = pop(:,1,:,t+1) ./ sum(sum(pop(:,1,:,t+1))) ;

% The following two lines of code serve to record our summary statistics
% for each generation: tag frequencies, and weighted tag cooperator
% frequency deviations.
tagfreq(:,t+1) = reshape(sum(pop(:,1,:,t) ),[2,1]) ; 
deviationsWeight(:,t+1) = reshape((pop(2,1,:,t+1) ./ sum(pop(:,1,:,t+1)) - sum(pop(2,1,:,t+1))).*sum(pop(:,1,:,t+1)),[2,1]);

end

figure
% The following line plots tag 1 deviation over generations.
plot(1:T+1,deviationsWeight(1,1:T+1) ,'LineWidth',2,'LineStyle',':','Color',[0, 0.4470, 0.7410]	)
hold on
% The following line plots tag 2 deviation over generations.
plot(1:T+1,deviationsWeight(2,1:T+1) ,'LineWidth',2,'LineStyle',':','Color',[1.0000    0.6235    0.4588]	)

% The following 7 lines of code are for formatting the plot.
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off
box off
set(gcf,'color','white')
xlim([1,7000000])

yyaxis right
% The following line plots tag 1 frequency over generations.
plot(1:T+1,tagfreq(1,1:T+1),'LineWidth',2,'LineStyle','-','Color',[0, 0.4470, 0.7410])
hold on
% The following line plots tag 2 frequency over generations.
plot(1:T+1,tagfreq(2,1:T+1),'LineWidth',2,'LineStyle','-','Color',[1.0000    0.6235    0.4588]	)

% The following 9 lines of code are for formatting the plot.
hold off
box off
set(gcf,'color','white')
xlim([1,7000000])
ylim([-1,1])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca,'Layer','top')