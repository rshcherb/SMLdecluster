function plot_nnd_clusters(DiG,NNDStats,vCat,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 3 June 2022
%   ...
%   version 1.3.0, 1 December 2022
%
    bAll = true;
    bFamilyTrees = false;
    bClustStruct = false;
    bClustMag = false;
    bSortClustMap = false;
    bClustMap = false;
    bSaveFig = true;
    for k = 1:length(varargin)
        if strcmp('FamilyTrees',varargin{k})
            bFamilyTrees = true;
            bAll = false;
        end
        if strcmp('ClusterStructure',varargin{k})
            bClustStruct = true;
            bAll = false;
        end
        if strcmp('ClusterMagnitude',varargin{k})
            bClustMag = true;
            bAll = false;
        end
        if strcmp('SortedClusterMap',varargin{k})
            bSortClustMap = true;
            bAll = false;
        end
        if strcmp('ClusterMap',varargin{k})
            bClustMap = true;
            bAll = false;
        end
    end
    
    sPlateNum = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)',...
                 'm)', 'n)', 'o)', 'p)', 'q)', 'r)', 's)', 't)', 'u)', 'v)', 'w)', 'x)'};

    nComp_indx = NNDStats.nComp_indx; % the index of each connected cluster in the graph DiG
    nComp_size = NNDStats.nComp_size; % the size of each connected cluster
    ClstProp   = NNDStats.ClstProp;   % properties of each cluster
    nClst      = NNDStats.nClst;      % number of clusters
    nSingle    = NNDStats.nSingle;    % number of singleton events

    if bAll || bFamilyTrees
        figure('Name','Plot event family trees','Position',[550 300 900 800]);
        set(gcf,'Color','w'); % background color for the figure
    
        subplot(2,2,1)
        plot(DiG);
        hold on
        title(['Family trees: ',num2str(nClst),'; singletons: ',num2str(nSingle)]);
        text(0.0,0.9,sPlateNum{1},'FontSize',12,'Units','normalized');
        axis off
        hold off
    
        idx = nComp_size(nComp_indx) > 1; % all clusters except singletons
        SG1 = subgraph(DiG,idx);
        subplot(2,2,2)
        plot(SG1,'Layout','layered');
        hold on
        title(['Family trees without singletons: ',num2str(nClst - nSingle)]);
        text(0.0,0.9,sPlateNum{2},'FontSize',12,'Units','normalized');
        axis off
        hold off
    
        idx = nComp_size(nComp_indx) >= Model.nNode_min;
        SGs = subgraph(DiG,idx);
        nClst_thrs = numel(find(nComp_size >= Model.nNode_min));
        subplot(2,2,3)
        plot(SGs);
        %plot(SGs,'Layout','layered');
        hold on
        title(['Family trees with number of events above ',num2str(Model.nNode_min),': ',num2str(nClst_thrs)]);
        text(0.0,0.9,sPlateNum{3},'FontSize',12,'Units','normalized');
        axis off
        hold off
        
        idx = nComp_size(nComp_indx) == max(nComp_size);
        SGmax = subgraph(DiG,idx);
        subplot(2,2,4)
        plot(SGmax,'NodeLabel',SGmax.Nodes.Name,'Layout','layered');
        hold on
        title(['The largest tree with number of events ',num2str(max(nComp_size))]);
        text(0.0,0.9,sPlateNum{4},'FontSize',12,'Units','normalized');
        axis off
        hold off
        sgtitle(Model.sTitle,'FontSize',12);
    
        if bSaveFig
            save_cf(gcf,[Model.sFileName,'_family_tree_struct'],'fig','png','pdf');
        end
    end

    % plot the first Model.nClst_largest clusters
    fwx = 900;
    fwy = 800;
    nClst_plot = Model.nClst_largest;
    if nClst_plot > nClst
        nClst_plot = nClst;
    end
    if nClst_plot <= 6
        nx = 3;
        ny = 2;
        fwy = 600;
    elseif nClst_plot > 6 && nClst_plot <= 9
        nx = 3;
        ny = 3;
    elseif nClst_plot > 9 && nClst_plot <= 12
        nx = 3;
        ny = 4;
    elseif nClst_plot > 12 && nClst_plot <= 16
        nx = 4;
        ny = 4;
    elseif nClst_plot > 16 && nClst_plot <= 20
        nx = 4;
        ny = 5;
    end

    [nComp_sort, idx_sort] = sort(nComp_size,'descend'); % sort the components from the largest to smallest

    if bAll || bClustStruct
        figure('Name','Sorted cluster structure','Position',[550 300 fwx fwy]);
        set(gcf,'Color','w'); % background color for the figure
        for nc = 1:nClst_plot
            subplot(ny,nx,nc);
            idx = nComp_indx == idx_sort(nc); % mark nodes belonging to a selected cluster
            Td = ClstProp.Td(idx_sort(nc));
            dm = ClstProp.dm(idx_sort(nc));
            avdepth = ClstProp.avdepth(idx_sort(nc));
            delta = ClstProp.delta(idx_sort(nc));
            BI = ClstProp.BI(idx_sort(nc));
            SG_comp = subgraph(DiG,idx);
            plot(SG_comp,'NodeLabel',SG_comp.Nodes.Magnitude,'Layout','layered');
            str = ['N = ',num2str(nComp_sort(nc)),'; T_d = ',num2str(Td,'%.2f'),' d; \Deltam = ',num2str(dm,'%.2f')];
            text(0.05,0.95,str,'Units','normalized','FontSize',8);
            str = ['<d> = ',num2str(avdepth,'%.3f'),'; \delta = ',num2str(delta,'%.3f'),'; B_I = ',num2str(BI,'%.3f')];
            text(0.05,-0.08,str,'Units','normalized','FontSize',8);
            text(0.0,1.1,sPlateNum{nc},'FontSize',12,'Units','normalized');
            title(['Cluster #',num2str(idx_sort(nc))]);
            axis off
        end
        sgtitle(Model.sTitle,'FontSize',12);
    
        if bSaveFig
            save_cf(gcf,[Model.sFileName,'_sorted_clusters_struct'],'fig','png','pdf');
        end
    end

    if bAll || bClustMag
        figure('Name','Sorted cluster magnitudes','Position',[550 300 fwx fwy]);
        set(gcf,'Color','w'); % background color for the figure
    
        for nc = 1:nClst_plot
            subplot(ny,nx,nc);
            idx = nComp_indx == idx_sort(nc); % mark nodes belonging to a selected cluster
            Td = ClstProp.Td(idx_sort(nc));
            dm = ClstProp.dm(idx_sort(nc));
            vTime = vCat(idx,1);
            vMag = vCat(idx,4);
            SG_comp = subgraph(DiG,idx);
            %plot(SG_comp,'XData',vTime,'YData',vMag,'NodeLabel',{});
            %plot(SG_comp,'XData',vTime,'YData',vMag,'NodeLabel',SG_comp.Nodes.Magnitude);
            plot(SG_comp,'XData',vTime,'YData',vMag,'NodeLabel',{},'Marker','none');
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            hold on
            sz = (2 - Model.fMc + vMag).^4;
            scatter(vTime,vMag,sz,'filled','MarkerFaceColor','#0072BD');
            rootidx = find(strcmp(SG_comp.Nodes.Type,'Root'));  % find the root node
            scatter(vTime(rootidx),vMag(rootidx),sz(rootidx),'filled','MarkerFaceColor','#A2142F'); % plot the root node
            %scatter(vTime(rootidx),vMag(rootidx),sz(rootidx),'MarkerEdgeColor','#D95319'); % plot the root node
            %str = ['N = ',num2str(nComp_sort(nc)),'; T_d = ',num2str(Td),' d; \Deltam = ',num2str(dm)];
            text(0.6,0.95,['N = ',num2str(nComp_sort(nc))],'Units','normalized','FontSize',8);
            text(0.6,0.86, ['T_d = ',num2str(Td,'%.2f'),' d'],'Units','normalized','FontSize',8);
            text(0.6,0.77,['\Deltam = ',num2str(dm,'%.2f')],'Units','normalized','FontSize',8);
            text(0.0,1.1,sPlateNum{nc},'FontSize',12,'Units','normalized');
            %xlabel('days');
            %ylabel('m');
            if ceil(nc/nx) == ny
                xlabel('days');
            end
            if mod(nc-1+nx,nx) == 0
                ylabel('m');
            end
            title(['Cluster #',num2str(idx_sort(nc))]);
            hold off
        end
        sgtitle(Model.sTitle,'FontSize',12);
    
        if bSaveFig
            save_cf(gcf,[Model.sFileName,'_sorted_clusters_magnitude'],'fig','png','pdf');
        end
    end

    if bAll || bSortClustMap
        figure('Name','Sorted cluster maps','Position',[550 300 fwx fwy]);
        set(gcf,'Color','w'); % background color for the figure
    
        for nc = 1:nClst_plot
            subplot(ny,nx,nc);
            idx = nComp_indx == idx_sort(nc); % mark nodes belonging to a selected cluster
            vLat = vCat(idx,2);
            vLon = vCat(idx,3);
            vTime = vCat(idx,1);
            vMag = vCat(idx,4);
            SG_comp = subgraph(DiG,idx);
    %         if hascycles(SG_comp)
    %             disp(['has cycles: ',num2str(nc)])
    %         end
            plot(SG_comp,'XData',vLon,'YData',vLat,'NodeLabel',{},'Marker','none');
            hold on
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            axop = get(ax,'OuterPosition');
            axop(2) = axop(2)*0.95;
            axop(3) = axop(3)*1.1;
            set(ax,'OuterPosition',axop);
            sz = (2 - Model.fMc + vMag).^4;
            scatter(vLon,vLat,sz,vTime,'filled');
            rootidx = find(strcmp(SG_comp.Nodes.Type,'Root'));  % find the root node
            scatter(vLon(rootidx),vLat(rootidx),sz(rootidx),vTime(rootidx),'filled','MarkerFaceColor','#A2142F'); % plot the root node
            %scatter(vLon(rootidx),vLat(rootidx),sz(rootidx),vTime(rootidx),'MarkerEdgeColor','#D95319'); % plot the root node
            c = colorbar;
            xc = get(c,'Position');
            xax = get(ax,'Position');
            xc(3) = xc(3)*0.5; % reduce the width in half
            set(c,'Position',xc);
            set(ax,'Position',xax);
            if mod(nc,nx) == 0
                c.Label.String = 'time (d)';
            end
            if ceil(nc/nx) == ny
                xlabel('longitude');
            end
            if mod(nc-1+nx,nx) == 0
                ylabel('latitude');
            end
            text(0.0,1.1,sPlateNum{nc},'FontSize',12,'Units','normalized');
            title(['Cluster #',num2str(idx_sort(nc))]);
            hold off
        end
        sgtitle(Model.sTitle,'FontSize',12);
    
        if bSaveFig
            save_cf(gcf,[Model.sFileName,'_sorted_clusters_map'],'fig','png','pdf');
        end
    end

    if bAll || bClustMap
        figure('Name','Cluster map','Position',[550 300 fwx fwy]);
        set(gcf,'Color','w'); % background color for the figure
        co = get(0,'DefaultAxesColorOrder');
    
        idx = nComp_size(nComp_indx) == 1; % all singletons
        SG_single = subgraph(DiG,idx);
        vLat = vCat(idx,2);
        vLon = vCat(idx,3);
        vMag = vCat(idx,4);
        vgpl(1) = plot(SG_single,'XData',vLon,'YData',vLat,'NodeLabel',{},'MarkerSize',2,'NodeColor','black');
        hold on
        ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
        sz = (2 - Model.fMc + vMag).^4;
        scatter(vLon,vLat,sz,'black','filled');
        cLegend = {'1'};
    
        idxc = nComp_size > 1; % all clusters except singletons
        vnc = unique(nComp_size(idxc));
        %vcol = lines; % define the colormap
        vcol = [1.071*co; co; 1.2*co/1.4; 1.2*co/1.6; 1.2*co/1.8; 1.2*co/2.0; 1.2*co/2.2; 1.2*co/2.4]; % define the colormap
        ncol = 1;
        for nc = vnc
            idx = nComp_size(nComp_indx) == nc; % a cluster of size nc
            vLat = vCat(idx,2);
            vLon = vCat(idx,3);
            vMag = vCat(idx,4);
            SG_comp = subgraph(DiG,idx);
            vgpl(ncol+1) = plot(SG_comp,'XData',vLon,'YData',vLat,'NodeLabel',{},'MarkerSize',1,'NodeColor',vcol(ncol,:),'EdgeColor',vcol(ncol,:));
            sz = (2 - Model.fMc + vMag).^4;
            scatter(vLon,vLat,sz,vcol(ncol,:),'filled');
            ncol = ncol + 1;
            cLegend = [cLegend, {num2str(nc)}];
        end
        legend(vgpl,cLegend,'Location','northwest','Box','off');
        xlabel('longitude');
        ylabel('latitude');
        title(Model.sTitle,'FontSize',12);
        hold off
    
        if bSaveFig
            save_cf(gcf,[Model.sFileName,'_cluster_map'],'fig','png','pdf');
        end
    end
end

