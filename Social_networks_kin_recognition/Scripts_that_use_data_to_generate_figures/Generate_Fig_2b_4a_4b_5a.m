% This script generates Fig. 2b, 4a, 4b & 5a (the "igloo" plots). It generates
% these plots by loading data that has been collected by iterating our
% weak-selection mathematical model described in the main text and Supp. 
% Info. 3, for the following parameter values: tag = 100 (maximum tags), 
% N=30 (deme size), T = 100000 (generations), mutC = 0.001 (trait mutation 
% rate), tagSkewIni = 0.9 (initial frequeny of one dominant tag), b = 0.3 
% (benefit), c = 0.1 (cost), helpini = 0.1 (initial helper frequency).

clearvars
close all
clc

for aaaa = 1:4
    
    if aaaa == 1
        load ../Saved_data/weak_sel_alpha=1.csearch=0.mat % This loads data for alpha=1, used to generate Fig. 2b.

    elseif aaaa == 2        
        load ../Saved_data/weak_sel_alpha=0.975.csearch=0.mat % This loads data for alpha=0.975, used to generate Fig. 4a.
        
    elseif aaaa == 3        
        load ../Saved_data/weak_sel_alpha=0.csearch=0.mat % This loads data for alpha=0, used to generate Fig. 4b.
       
    elseif aaaa == 4  
        load ../Saved_data/weak_sel_alpha=0.999.csearch=0.0009.mat % This loads data for alpha=0.999 & csearch=0.0009, used to generate Fig. 5a.
               
    end      

figure

set(gcf,'color','white')

rR2=0:0.001:0.5;
mR2= 0:0.002:1;

for cur_r = 1:length(rR2)
    
    r = rR2(cur_r);
    
for cur_m = 1:length(mR2)
    
    m = mR2(cur_m);
    
gamma = ((-1+m)^3*(-12*(-1+r)^3-16*m^7*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*m^8*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^6*(-2+N)*(-1+N)^2*(-1+r)*(168*(-1+r)^2+N*(-221+(448-225*r)*r)+N^2*(55+56*(-2+r)*r))-m^5*(-1+N)^2*(-1+r)*(-672*(-1+r)^2+N^3*(106+r*(-225+113*r))+4*N*(299+r*(-616+311*r))-2*N^2*(321+r*(-671+338*r)))+m^4*(-1+N)*(-1+r)*(840*(-1+r)^2+N^3*(-883+(1951-989*r)*r)+N^4*(122+r*(-281+143*r))-10*N*(229+r*(-476+241*r))+N^2*(2209+r*(-4750+2417*r)))+m^2*(-1+N)*(336*(-1+r)^3-N^3*(-1+r)*(260+r*(-756+397*r))-2*N*(-1+r)*(431+r*(-952+491*r))+N^2*(-1+r)*(751+r*(-1884+989*r))+N^4*(-24+r*(122+r*(-158+57*r))))-m^3*(-1+N)*(672*(-1+r)^3-4*N^3*(-1+r)*(159+2*r*(-193+99*r))-8*N*(-1+r)*(223+r*(-476+243*r))+4*N^2*(-1+r)*(415+r*(-947+488*r))+N^4*(-80+r*(298+r*(-334+115*r))))+m*(96*(-1+r)^3-4*N^5*r*(5+2*r*(-5+2*r))+22*N^2*(-1+r)*(19+r*(-49+26*r))-4*N*(-1+r)*(83+r*(-184+95*r))-N^3*(-1+r)*(230+r*(-743+403*r))+N^4*(-48+r*(278+r*(-364+131*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r))))));
phi = -(((-1+m)^2*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+35*N^3*(-1+r)^2+N^2*(-211-212*(-2+r)*r)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^3*(-1+r)*(-76+77*r)+N^2*(-466+(945-472*r)*r)+14*N*(63+r*(-129+65*r)))+m^5*(-2+N)*(-1+N)*(-1+r)*(756*(-1+r)^2+2*N^4*(-1+r)*(-49+52*r)+N^3*(-731-753*(-2+r)*r)-42*N*(48+r*(-99+50*r))+N^2*(1893+2*r*(-1977+998*r)))-m^4*(-1+N)*(-1512*(-1+r)^3+2*N^5*(-1+r)^2*(-36+43*r)+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+2*r*(-1439+732*r))+N^3*(-1+r)*(2932+r*(-6321+3194*r))+N^4*(754+r*(-2371+(2441-823*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+8*N^5*(-1+r)^2*(-3+5*r)+140*N*(-1+r)*(21+r*(-45+23*r))-2*N^2*(-1+r)*(1592+3*r*(-1203+619*r))+N^3*(-1+r)*(1583+2*r*(-1837+938*r))+N^4*(348+r*(-1168+(1253-429*r)*r)))+m^2*(-432*(-1+r)^3-8*N^6*(-1+r)^2*r+N^3*(-1+r)*(1642+15*r*(-279+146*r))+12*N*(-1+r)*(135+r*(-291+149*r))-2*N^2*(-1+r)*(1173+r*(-2748+1429*r))+N^4*(556+r*(-2097+(2362-815*r)*r))+N^5*(-72+r*(304+r*(-374+135*r))))+m*(108*(-1+r)^3-N^3*(-1+r)*(293-924*r+500*r^2)-2*N^5*r*(12+5*r*(-5+2*r))-4*N*(-1+r)*(96+r*(-213+110*r))+N^2*(-1+r)*(503+r*(-1292+685*r))+N^4*(-66+r*(360+r*(-465+167*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r)))))));
F = (1-m)^2 / ( (1-m)^2 + N*(1-(1-m)^2)) ;

denom = (1 - alpha .* (1-1/tag) .* (1-F)) ;
M1 = ((F+(1-F).* 1/tag) ./ denom); % prob give help (requires tag match and that you are a helper)
M2 = ( phi + (F-phi).* 1/tag)   ./ denom; % prob of receiving help when you are in the role of a recipeint

Rtag = M2/M1;
Rdisp = Rtag/N + 1/N + (1-2/N) * ((gamma  + (F-gamma).* (1/tag .* (1/tag + (1-1/tag).*alpha.*M1))) / M1) ;

if b .*  (M2/M1)  - c    - (1-m) .* (b-c)  .* ( 1/N  + (M2/M1)/N + ((N-2)/N) .* (gamma  + (F-gamma).* ( (1/tag + (1-1/tag).*alpha.*M1)) )./M1) > 0            
resrel(cur_r,cur_m) = 1 ;
else
resrel(cur_r,cur_m) = 0 ; 
end
% 'resrel' indicates when kin selection (Hamilton's Rule) suggests kin
% recognition is favoured. It will be used to generate the solid 
% lines for the "igloo" plots.
    
end
end

res = 1./avgtagfreq > 10 & resalt>0.4; % This criterion indicates when (for given m and r values) kin discrimination based on genetic cues evolves. It is assumed (arbitrarily) that kin discrimination based on genetic cues corresponds to an excess of 10 tags maintained at equilibrium alongside a population cooperator frequency of over 0.4.
[row, column] = find(res == 1);
scatter( mR(column), rR(row),[],[0.5 0.5 0.5],'filled') % This generates a scatter plot, where each point corresponds to parameter (m & r) combinations where kin discrimination based on genetic cues evolves.
box on

axis([min(mR) max(mR) min(rR) max(rR)])
yspace=min(rR):max(rR)/5:max(rR);
yticks([yspace]);
yticklabels({yspace});
xspace = min(mR):max(mR)/5:max(mR);
xticks([xspace])
xticklabels({xspace})
hold on

% The following lines of code are used to generate the solid lines of the
% "igloo" plots, which show when kin selection (Hamilton's Rule) predicts
% that kin discrimination is favoured. These solid lines are superimposed
% onto the scatter plots, which show when kin discrimination based on
% genetic cues is stable.
X = mR2;
YYY = resrel - [resrel(2:length(resrel(:,1)),:); zeros(1,length(resrel(1,:)))];
Y=zeros(1,length(resrel(1,:)));
for ii = 1: length(resrel(1,:))
    if find(YYY(:,ii))>0
Y(ii) =  rR2(find(YYY(:,ii)));
    else
        Y(ii) =  0;
    end
end
plot(X,Y,'color','k','Linestyle','-','LineWidth',2)

hold off

end
