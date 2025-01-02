function plot_confusion_chart(true_labels,predicted_labels)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        true_labels
        predicted_labels
    end
    
    correct_0_inds = true_labels == 0 & predicted_labels == 0;
    correct_1_inds = true_labels == 1 & predicted_labels == 1;
    missed_0_inds = true_labels == 0 & predicted_labels == 1;
    missed_1_inds = true_labels == 1 & predicted_labels == 0;
    
    correct_0_perc = 100*(sum(correct_0_inds)/sum(true_labels == 0));
    correct_1_perc = 100*(sum(correct_1_inds)/sum(true_labels == 1));
    missed_0_perc = 100 - correct_0_perc;
    missed_1_perc = 100 - correct_1_perc;
    
    % Set plotting visuals
    clrs = lines(4);
%     clrs = [rgb('gray');rgb('gray')];
    txtSize = 12;
    txtColor = 'w';
    patchFaceAlpha = 0.9;
    imagesc([reshape(clrs(1,:),1,1,3),reshape(clrs(4,:),1,1,3);...
        reshape(clrs(3,:),1,1,3),reshape(clrs(2,:),1,1,3)],...
        'AlphaData',[patchFaceAlpha,patchFaceAlpha;patchFaceAlpha,patchFaceAlpha]);
    
    text(1,1,sprintf('%s\n\\textbf{%.1f\\%%}',addComma(sum(correct_0_inds)),correct_0_perc),...
        'HorizontalAlignment', 'center','FontSize',txtSize,'Color',txtColor,'Interpreter','latex')
    text(2,2,sprintf('%s\n\\textbf{%.1f\\%%}',addComma(sum(correct_1_inds)),correct_1_perc),...
        'HorizontalAlignment', 'center','FontSize',txtSize,'Color',txtColor,'Interpreter','latex')
    text(2,1,sprintf('%s\n\\textbf{%.1f\\%%}',addComma(sum(missed_0_inds)),missed_0_perc),...
        'HorizontalAlignment', 'center','FontSize',txtSize,'Color',txtColor,'Interpreter','latex')
    text(1,2,sprintf('%s\n\\textbf{%.1f\\%%}',addComma(sum(missed_1_inds)),missed_1_perc),...
        'HorizontalAlignment', 'center','FontSize',txtSize,'Color',txtColor,'Interpreter','latex')
    
    set(gca, 'XTick', [1,2], ...                             % Change the axes tick marks
         'XTickLabel', {'Bkgrd', 'Ashk'}, ...  %   and tick labels
         'YTick', [1,2], ...
         'YTickLabel', {'Bkgrd', 'Ashk'}, ...
         'TickLength', [0 0]);
end

