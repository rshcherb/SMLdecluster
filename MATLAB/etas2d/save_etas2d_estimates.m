function save_etas2d_estimates(vCat,Model,Bckgrnd,varargin)
%
%   Saves the results of the ETAS fitting and stochastic declustering
%     - background probability associated with each event in the catalogue vCat
%     - stochastically declustered events
%     - earthquake catalogue with events marked as background (1) or aftershocks (0)
%     - variable kernel estimate of the spatial background rate
%
%   Model   - structure that defines the model
%   Bckgrnd - structure that defines the background rate 
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 30 December 2021
%   ...
%   version: 1.0.0, 3 January 2022
%
    bSaveResults = false;
    for k = 1:length(varargin)
        if strcmp('SaveResults',varargin{k})
            bSaveResults = true;
            sFileName = varargin{k+1};
        end
    end

    if bSaveResults
        fid = fopen([sFileName,'_background_prob.dat'],'w');
        for n = 1:length(Bckgrnd.vBgProb)
            fprintf(fid,'%10e %8.4f %9.4f %10e\n',vCat(n,1:3),Bckgrnd.vBgProb(n));
        end
        fclose(fid);

        fid = fopen([sFileName,'_declustered_events.dat'],'w');
        for n = 1:length(Bckgrnd.vBgDeclstr(:,1))
            fprintf(fid,'%10e %8.4f %9.4f %.2f\n',Bckgrnd.vBgDeclstr(n,:));
        end
        fclose(fid);

        fid = fopen([sFileName,'_declustered_catalog.dat'],'w');
        for n = 1:length(Bckgrnd.indxBg(:,1))
            fprintf(fid,'%10e %8.4f %9.4f %.2f %d\n',vCat(n,1:4),Bckgrnd.indxBg(n,1));
        end
        fclose(fid);

        fid = fopen([sFileName,'_est_background_rate.dat'],'w');
        for iy = 1:length(Model.vY)
            for ix = 1:length(Model.vX)
                fprintf(fid,'%d %d %f %f %e\n',ix,iy,Model.vX(ix),Model.vY(iy),Bckgrnd.mMu(iy,ix));
            end
        end
        fclose(fid);
    end
end
