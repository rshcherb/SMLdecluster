function ax = plot_cumul_synth_rates(testCats,BkgrdInds,AshkInds,Model,OptArgs)
%
%   Plots the cumulative event rates for the synthetic clustered and background events
%   vCat - catalogue of events
%   BkgrdInds, AshkInds - indeces for background and aftershock events 
%   Model - structure that defines the seismicity model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 October 2024
%   ...
%   version 1.0.0, 21 October 2024
%
    arguments
        testCats table
        BkgrdInds logical
        AshkInds logical
        Model struct
        OptArgs.DateFormat char = 'day'   % 'date', 'day' format to use for x-axis 
    end

    % vTbgrd = vCat(BkgrdInds,1);
    % vTashk = vCat(AshkInds,1);
    % t0 = vCat(find(BkgrdInds,1),1);
    % t1 = vCat(find(BkgrdInds,1,'last'),1);
    % if strcmp(OptArgs.DateFormat,'date')
    %     vTbgrd = datenum(Model.vDateStart) + vTbgrd;
    %     vTashk = datenum(Model.vDateStart) + vTashk;
    %     t0 = datenum(Model.vDateStart) + t0;
    %     t1 = datenum(Model.vDateStart) + t1;
    % end

    catalog_list_synth = unique(testCats.Catalog);
    
    % Set plotting visuals
    synthBkgrdAlpha = 0.15;
    synthAshkAlpha = 0.2;
    synthBkgrdLineWidth = 0.1;
    synthAshkLineWidth = 0.1;
    
    % Set plotting visuals
    clrs = lines(5);
    colors = [clrs(1,:);clrs(2,:);[0 0 0]];
    
    for i = 1:length(catalog_list_synth)
        cat_idx = testCats.Catalog == catalog_list_synth(i);
        tot_rate = cumsum(ones(sum(cat_idx),1))/sum(cat_idx);
        ashk_idx = cat_idx & AshkInds;
        ashk_rate = cumsum(ones(sum(ashk_idx),1))/sum(ashk_idx);
        bkgrd_idx = cat_idx & BkgrdInds;
        bkgrd_rate = cumsum(ones(sum(bkgrd_idx),1))/sum(bkgrd_idx);
        p(1) = plot(testCats.Time(bkgrd_idx),bkgrd_rate,'-','LineWidth',synthBkgrdLineWidth,...
            'Color',[colors(1,:) synthBkgrdAlpha],'DisplayName','Synth Bkgrd');
        p(2) = plot(testCats.Time(ashk_idx),ashk_rate,'-','LineWidth',synthAshkLineWidth,...
            'Color',[colors(2,:) synthAshkAlpha],'DisplayName','Synth Ashk');
%         plot(synthCat.Date(cat_idx),tot_rate,'-','LineWidth',0.1,'Color',[colors(2,:) 0.2])
    end
end

