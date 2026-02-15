function calculateAndPlotPhytoplanktonSizeClasses(HPLC)

% CALCULATEANDPLOTPHYTOPLANKTONSIZECLASSES Using pigment composition from 
% HPLC data, calculate fractional contribution of three main phytoplankton 
% size classes using published fit parameters as well as our own, and plot 
% the results using various visuals.
%
%   INPUT: 
%       HPLC - HPLC data filtered
%            
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 7 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

fprintf('\nCalculating phytoplankton size classes...\n')

%% Weights

nWeightStudies = 7; % 7 studies featured in our collection
W = zeros(nWeightStudies,7);

% Weights from published studies
W(1,:) = [1.41,  1.41,  1.27,  0.35,  0.60,  1.01,  0.86];  % Uitz et al. (2006) - global (N=2419)
W(2,:) = [1.65,  1.30,  0.83,  0.78,  0.73,  0.77,  1.29];  % Uitz et al. (2008)
W(3,:) = [1.72,  1.27,  0.68,  1.42,  4.96,  0.81,  1.28];  % Brewin et al. (2014) - North Atl., South Atl. (N=466)
W(4,:) = [1.51,  1.35,  0.95,  0.85,  2.71,  1.27,  0.93];  % Brewin et al. (2015) - global (N=5841)
W(5,:) = [1.65,  1.04,  0.78,  1.19,  3.14,  1.38,  1.02];  % Brewin et al. (2017) - North Atl. (N=2791)
W(6,:) = [1.554, 0.413, 0.855, 1.174, 2.387, 1.062, 2.037]; % Soppa et al. (2014) - global (N=3988)

% Weights estimated from the CEFAS HPLC dataset, follow Chase et al. (2020)
% Help: https://uk.mathworks.com/matlabcentral/answers/588217-how-can-i-solve-the-set-of-liner-equations-using-lsqlin
corrCoeff = [HPLC.Fuc_ug_L(:),HPLC.Per_ug_L(:),HPLC.Hex_ug_L(:),HPLC.But_ug_L(:),HPLC.Allo_ug_L(:),HPLC.Chlb_ug_L(:),HPLC.Zea_ug_L(:)];
R = HPLC.TChlA_ug_L(:);
W(7,:) = lsqlin(corrCoeff,R,[],[],[],[],0.1*ones(7,1),5*ones(7,1));

%% Estimate phytoplankton size classes (% fraction and their chla concentration)

nMethods = 4; % 4 methods to estimate chla and phytoplankton size classes
nDataPoints = height(HPLC);

estim.TChla = zeros(nDataPoints,nWeightStudies);

% percentual fraction of phytoplankton size class
estim.fmicro = zeros(nDataPoints,nWeightStudies,nMethods);
estim.fnano = zeros(nDataPoints,nWeightStudies,nMethods);
estim.fpico = zeros(nDataPoints,nWeightStudies,nMethods);

% Chla concentration by phytoplankton size class
estim.microTChla = zeros(nDataPoints,nWeightStudies,nMethods);
estim.nanoTChla = zeros(nDataPoints,nWeightStudies,nMethods);
estim.picoTChla = zeros(nDataPoints,nWeightStudies,nMethods);

for i = 1:nDataPoints
    
    for j = 1:nWeightStudies
                   
        estim.TChla(i,j) = W(j,1)*HPLC.Fuc_ug_L(i) + W(j,2)*HPLC.Per_ug_L(i) +...
                           W(j,3)*HPLC.Hex_ug_L(i) + W(j,4)*HPLC.But_ug_L(i) +... 
                           W(j,5)*HPLC.Allo_ug_L(i) + W(j,6)*HPLC.Chlb_ug_L(i) +... 
                           W(j,7)*HPLC.Zea_ug_L(i);
        
        for k = 1:nMethods

            if (k == 1) % Method 1: Uitz et al. (2006) (as in Brewin et al. 2014, Method A)

                estim.fmicro(i,j,k) = (W(j,1)*HPLC.Fuc_ug_L(i) + W(j,2)*HPLC.Per_ug_L(i))/estim.TChla(i,j);
                estim.fnano(i,j,k) = (W(j,3)*HPLC.Hex_ug_L(i) + W(j,4)*HPLC.But_ug_L(i) + W(j,5)*HPLC.Allo_ug_L(i))/estim.TChla(i,j);
                estim.fpico(i,j,k) = (W(j,6)*HPLC.Chlb_ug_L(i) + W(j,7)*HPLC.Zea_ug_L(i))/estim.TChla(i,j);

            elseif (k == 2) % Method 2: Brewin et al. (2010) (as in Brewin et al. 2014, Method B) - derived for oligotrophic conditions (Chla < 0.08 ug/L)
    
                estim.fmicro(i,j,k) = estim.fmicro(i,j,1); % as in Method 1

                if (HPLC.TChlA_ug_L(i) <= 0.08)

                    estim.fnano(i,j,k) = ((12.5 * TChlA_ug_L(i) * W(j,3)*HPLC.Hex_ug_L(i))/estim.TChla(i,j)) +...
                        ((W(j,4)*HPLC.But_ug_L(i) + W(j,5)*HPLC.Allo_ug_L(i))/estim.TChla(i,j));

                     estim.fpico(i,j,k) = (((-12.5 * TChlA_ug_L(i) + 1) * W(j,3)*HPLC.Hex_ug_L(i))/estim.TChla(i,j)) +...
                        ((W(j,6)*HPLC.Chlb_ug_L(i) + W(j,7)*HPLC.Zea_ug_L(i))/estim.TChla(i,j));

                else
                    estim.fnano(i,j,k) = estim.fnano(i,j,1); % as in Method 1
                    estim.fpico(i,j,k) = estim.fpico(i,j,1); % as in Method 1
                end
                
            elseif (k == 3) % Method 3: Devred et al. (2011), fucoxanthin adjustment - derived for oligotrophic conditions (Chla < 0.08 ug/L)
    
                % Fucoxanthin adjustment
                nanoFucoxanthin = 10^(0.531*log10(HPLC.Hex_ug_L(i)) + 0.708*log10(HPLC.But_ug_L(i)));
                microFucoxanthin = HPLC.Fuc_ug_L(i) - nanoFucoxanthin;
                if (microFucoxanthin < 0)
                    microFucoxanthin = 0;
                    nanoFucoxanthin = HPLC.Fuc_ug_L(i);
                end

                estim.fmicro(i,j,k) = ((W(j,1)*HPLC.Fuc_ug_L(i) + W(j,2)*HPLC.Per_ug_L(i)) - W(j,1)*nanoFucoxanthin)/estim.TChla(i,j);
                estim.fnano(i,j,k) = ((W(j,3)*HPLC.Hex_ug_L(i) + W(j,4)*HPLC.But_ug_L(i) + W(j,5)*HPLC.Allo_ug_L(i)) + W(j,1)*nanoFucoxanthin)/estim.TChla(i,j);
                estim.fpico(i,j,k) = estim.fpico(i,j,1); % as in Method 1
    
            elseif (k == 4) % Method 4: Hirata et al. (2008)
                
                estim.fmicro(i,j,k) = estim.fmicro(i,j,1); % as in Method 1

                if (HPLC.TChlA_ug_L(i) >= 0.25)
                    estim.fnano(i,j,k) = ((W(j,3)*HPLC.Hex_ug_L(i) + W(j,4)*HPLC.But_ug_L(i) + W(j,5)*HPLC.Allo_ug_L(i)) + W(j,6)*HPLC.Chlb_ug_L(i))/estim.TChla(i,j);
                    estim.fpico(i,j,k) = (W(j,7)*HPLC.Zea_ug_L(i))/estim.TChla(i,j);
                elseif (HPLC.TChlA_ug_L(i) < 0.25)
                    estim.fnano(i,j,k) = estim.fnano(i,j,k); % as in Method 1
                    estim.fpico(i,j,k) = estim.fpico(i,j,k); % as in Method 1
                end
                
            end
    
            estim.microTChla(i,j,k) = estim.fmicro(i,j,k) * HPLC.TChlA_ug_L(i);
            estim.nanoTChla(i,j,k) = estim.fnano(i,j,k) * HPLC.TChlA_ug_L(i);
            estim.picoTChla(i,j,k) = estim.fpico(i,j,k) * HPLC.TChlA_ug_L(i);
            
        end % k (nMethods)
        
    end % j (nWeightStudies)
    
end % i (nDataPoints)

fprintf('\n...done.\n')

%% Statistical tests

stats.phi = zeros(nWeightStudies,1); % RMSE
stats.delta = zeros(nWeightStudies,1); % bias
stats.corrCoeff = zeros(nWeightStudies,1);
stats.pvalue = zeros(nWeightStudies,1);

for j = 1:nWeightStudies                        
    [stats.phi(j),~] = calcRootMeanSquaredError(HPLC.TChlA_ug_L(:),estim.TChla(:,j)); % centre-patterned (or unbiased) root-mean-square error
    [stats.delta(j),~] = calcMeanError(HPLC.TChlA_ug_L(:),estim.TChla(:,j));
    [stats.corrCoeff(j),pvalue(j)] = corr(estim.TChla(:,j),HPLC.TChlA_ug_L(:),'Type','Pearson','Rows','complete'); % Pearson linear correlation coefficient using only the rows that contain no missing values.
end

% What does RMSE indicate?:
%   It indicates the absolute fit of the model to the data.
%   Provides average model prediction error in units of the variable of interest.
%   They are negatively-oriented scores, which means lower values are better.

%% Plots

plotTChlaObservedVsEstimated(HPLC,nWeightStudies,estim,stats)
plotTernaryPlots(nMethods,nWeightStudies,estim)
plotChlaVsPercentagePhytoplanktonSizeClass(estim)

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function plotTChlaObservedVsEstimated(P,nWeightStudies,estim,stats)

isUseLogScale = 1;

labelWeightStudies = {'Uitz et al. (2006)','Uitz et al. (2008)',...
    'Brewin et al. (2014)','Brewin et al. (2015)','Brewin et al. (2017)',...
    'Soppa et al. (2014)','CEFAS'};

%%%%%%%%% Separated plots %%%%%%%%%

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.70],'Color','w') 
haxis = zeros(nWeightStudies,1);

for iSubplot = 1:nWeightStudies

    haxis(iSubplot) = subaxis(4,2,iSubplot,'Spacing',0.03,'Padding',0.03,'Margin', 0.04);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position'); 

    if (iSubplot == 1)
        % CEFAS weights
        scatter(haxis(iSubplot),P.TChlA_ug_L(:),squeeze(estim.TChla(:,end)),...
            15,[0.8 0.8 0.8],'filled'); hold on; 
        % Publication weights
        scatter(haxis(iSubplot),P.TChlA_ug_L(:),squeeze(estim.TChla(:,iSubplot)),...
            10,'+','red'); hold on;
    elseif (iSubplot > 1 && iSubplot < nWeightStudies)
        % CEFAS weights
        scatter(haxis(iSubplot),P.TChlA_ug_L(:),squeeze(estim.TChla(:,end)),...
            15,[0.8 0.8 0.8],'filled','HandleVisibility','off'); hold on; 
        % Publication weights
        scatter(haxis(iSubplot),P.TChlA_ug_L(:),squeeze(estim.TChla(:,iSubplot)),...
            10,'+','red','HandleVisibility','off'); hold on;
    else
        scatter(haxis(iSubplot),P.TChlA_ug_L(:),squeeze(estim.TChla(:,iSubplot)),...
            10,'+','red','HandleVisibility','off'); hold on;
    end

    if (isUseLogScale)
        
        set(gca,'Xscale','log','Yscale','log')

        xlim([0.06 50]); ylim([0.06 50])   
        xticks([1e-1 1 10]); yticks([1e-1 1 10])
        xticklabels({'0.1','1','10'}); yticklabels({'0.1','1','10'})
        
        % Plot 1:1 reference line
        xref = 10.^(-1.8:.1:1.8); yref = xref;
        plot(xref,yref,'k-')
        hold off;
        
        xt = min(xlim)+0.20*min(xlim); 
        yt = max(ylim)-0.20*max(ylim);

    else

        xlim([0 45]); ylim([0 45])   
        xticks([0 20 40]); yticks([0 20 40])
        xticklabels({'0','20','40'}); yticklabels({'0','20','40'})
        
        % Plot 1:1 reference line
        hline = refline(haxis(iSubplot),1,0); 
        hline.Color = 'k';
        hline.LineStyle = '-';
        hold off;
        
        xt = min(xlim)+1;
        yt = max(ylim)-0.05*max(ylim);
    
    end
    
    % Add statistics
    text(xt, yt,... % text position relative to axis
        {['\bf ',labelWeightStudies{iSubplot}, ' \rm'],...
         [strcat('{\it r}=',num2str(stats.corrCoeff(iSubplot,1),'%.2f'))],...
         [strcat('RMSE=',num2str(stats.phi(iSubplot),'%.2f'))],...
         [strcat('\delta=',num2str(stats.delta(iSubplot),'%.2f'))]},...
        'FontSize',8,'Horiz','left','Vert','top');
         
    ylabel('Estimated TChla (\mug/L)')
    xlabel('TChla (\mug/L)')
    
    box on

    if (iSubplot == 1)
        lg = legend(haxis(iSubplot),...
            'CEFAS weights','Publication weights','1:1 line','Location','southeast');
        lg.Position(1) = 0.62; lg.Position(2) = 0.11;
        lg.ItemTokenSize = [15,1];
        lg.FontSize = 12;
        set(lg,'Box','off')
    end
    
end

ax(1).pos(1) = ax(1).pos(1)+0.02; 
ax(2).pos(1) = ax(2).pos(1)+0.06;
ax(3).pos(1) = ax(3).pos(1)+0.02; 
ax(4).pos(1) = ax(4).pos(1)+0.06;
ax(5).pos(1) = ax(5).pos(1)+0.02; 
ax(6).pos(1) = ax(6).pos(1)+0.06;
ax(7).pos(1) = ax(7).pos(1)+0.02; 

for iSubplot = 1:nWeightStudies
    set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
end

exportgraphics(gcf,fullfile('.','figures','HPLC_chla_vs_estim_split_log.png'),'Resolution',600)   
clear ax cb

%%%%%%%%% One plot %%%%%%%%%

coloursStudies = jet(nWeightStudies);
sizeScatter = [56, 48, 40, 32, 24, 16, 8]; % use different dot size so that they don't overlap when I plot them altogether

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.40],'Color','w') 
haxis = zeros(1,1);

haxis(1) = subaxis(1,1,1,'Spacing',0.01,'Padding',0.02,'Margin', 0.17);
ax(1).pos = get(haxis(1),'Position'); 
for j = 1:nWeightStudies 
    scatter(haxis(1),P.TChlA_ug_L(:),squeeze(estim.TChla(:,j)),...
        sizeScatter(j),coloursStudies(j,:),'filled'); hold on;
end
box on;

if (isUseLogScale)

    set(gca,'Xscale','log','Yscale','log')

    xlim([0.06 50]); ylim([0.06 50])   
    xticks([1e-1 1 10]); yticks([1e-1 1 10])
    xticklabels({'0.1','1','10'}); yticklabels({'0.1','1','10'})

    % Plot 1:1 reference line
    xref = 10.^(-1.8:.1:1.8); yref = xref;
    plot(xref,yref,'k-')
    hold off;

    xt = min(xlim)+0.20*min(xlim); 
    yt = max(ylim)-0.20*max(ylim);

else

    xlim([0 45]); ylim([0 45])   
    xticks([0 20 40]); yticks([0 20 40])
    xticklabels({'0','20','40'}); yticklabels({'0','20','40'})

    % Plot 1:1 reference line
    hline = refline(haxis(1),1,0); 
    hline.Color = 'k';
    hline.LineStyle = '-';
    hold off;

    xt = min(xlim)+1;
    yt = max(ylim)-0.05*max(ylim);

end

xlabel('TChla (\mug/L)')
ylabel('Estimated TChla (\mug/L)')

lg = legend(haxis(1),labelWeightStudies);
lg.Position(1) = 0.72; lg.Position(2) = 0.57;
lg.ItemTokenSize = [15,1];
lg.FontSize = 11;
set(lg,'Box','off')

ax(1).pos(1) = ax(1).pos(1)-0.10; 
set(haxis(1),'Position',ax(1).pos) 

exportgraphics(gcf,fullfile('.','figures','HPLC_chla_vs_estim_one_log.png'),'Resolution',600)    
clear ax cb

end % plotTChlaObservedVsEstimated

% *************************************************************************

function plotTernaryPlots(nMethods,nWeightStudies,estim)

labelWeightStudies = {'W: Uitz et al. (2006)','W: Uitz et al. (2008)',...
    'W: Brewin et al. (2014)','W: Brewin et al. (2015)','W: Brewin et al. (2017)',...
    'W: Soppa et al. (2014)','W: CEFAS'};

nSubplots = nMethods*nWeightStudies;

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.58 0.58],'Color','w') 
haxis = zeros(nSubplots,1);

iWeightStudy = 1;
iMethod = 1;

for iSubplot = 1:nSubplots

    if (mod(iSubplot,((nWeightStudies*iMethod)+1)) == 0)
        iWeightStudy = 1;
        iMethod = iMethod + 1;
    end
    
    haxis(iSubplot) = subaxis(nMethods,nWeightStudies,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin', 0.02);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position'); 
        
    hp = ternplot(estim.fpico(:,iWeightStudy,iMethod),estim.fnano(:,iWeightStudy,iMethod),estim.fmicro(:,iWeightStudy,iMethod),...
        'r.','majors',5);
    ternlabel('pico %','nano %','micro %') 

    
    iWeightStudy = iWeightStudy + 1;
    
    if (iSubplot == 1)
        ht = text(-0.23,-0.04,'M: Uitz et al. (2006)','FontSize',10,'FontWeight','bold');
    elseif (iSubplot == 8)
        ht = text(-0.23,-0.12,'M: Brewin et al. (2010)','FontSize',10,'FontWeight','bold');
    elseif (iSubplot == 15)
        ht = text(-0.23,-0.12,'M: Devred et al. (2011)','FontSize',10,'FontWeight','bold');
    elseif (iSubplot == 22)
        ht = text(-0.23,-0.12,'M: Hirata et al. (2008)','FontSize',10,'FontWeight','bold');
    end
    
    set(ht,'Rotation',90)
    
    if (iSubplot >= 1 && iSubplot <= 7)
        title(labelWeightStudies(iSubplot),'FontSize',10)
    end
    
end

% set(gcf,'PaperPositionMode','manual')
orient(gcf,'landscape')
exportgraphics(gcf,fullfile('.','figures','HPLC_ternary_plots.png'),'Resolution',600)    
clear ax cb

end % plotTernaryPlots

% *************************************************************************

function plotChlaVsPercentagePhytoplanktonSizeClass(estim)

iWeightStudy = 7; % CEFAS fit
iMethod = 3; % Devred et al. (2010)
nBins = 7;

% maxChla = max(TChla_estim(:,iWeightStudy),[],'all');
% minChla = min(TChla_estim(:,iWeightStudy),[],'all');
% edges = linspace(log10(minChla),log10(maxChla),nBins);
% Y = discretize(TChla_estim(:,iWeightStudy),10.^edges);

[N,edges,bins] = histcounts(log10(estim.TChla(:,iWeightStudy)),nBins); % histcounts provides better binning
binCenters = edges(2:end) - (edges(2)-edges(1))./2; 

labelSubplots = {'pico %','nano %','micro %'};

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.60],'Color','w') 
haxis = zeros(3,1);

for iSubplot = 1:3

    haxis(iSubplot) = subaxis(3,1,iSubplot,'Spacing',0.02,'Padding',0.02,'Margin', 0.08);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position'); 

    if (iSubplot == 1)
        x = estim.fpico(:,iWeightStudy,iMethod);
    elseif (iSubplot == 2)
        x = estim.fnano(:,iWeightStudy,iMethod);
    elseif (iSubplot == 3)
        x = estim.fmicro(:,iWeightStudy,iMethod);
    end

    for iBin = 1:nBins
        binPhyto = x(bins==iBin);
        scatter(haxis(iSubplot),binCenters(iBin),binPhyto,22,'red'); hold on;
    end
    hold on
    meanBinPhyto = splitapply(@mean,x,bins);
    plot(haxis(iSubplot),binCenters,meanBinPhyto,'b--o','LineWidth',1.5,...
        'MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue'); hold off
        
    xlim([edges(1) edges(end)])
    xticks([-1 0 1])
    xticklabels({'0.1','1','10'})
    ylim([0 1])
    ylabel(labelSubplots(iSubplot))
    box on

    if (iSubplot == 3)
        xlabel('Estimated TChla (\mug/L)')
    end

end

exportgraphics(gcf,fullfile('.','figures','HPLC_frac_phyto_vs_estim_TChla.png'),'Resolution',600)  
clear ax cb

end % plotChlaVsPercentagePhytoplanktonSizeClass

% *************************************************************************

end % calculatePhytoplanktonSizeClasses
