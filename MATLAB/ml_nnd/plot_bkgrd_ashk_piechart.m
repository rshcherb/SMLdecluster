function pax = plot_bkgrd_ashk_piechart(NAshks,NBkgrd,OptArgs)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        NAshks (1,1)
        NBkgrd (1,1)
        OptArgs.ax = [];
        OptArgs.PieSize double = 0.075;
    end
    
    % Set plotting visuals
    clrs = lines(2);
    colors = [clrs(1,:);clrs(2,:)];
    patchFaceAlpha = 0.6;
    textFontSize = 9;
        
    if isempty(OptArgs.ax)
        ax = gca;
    else
        ax = OptArgs.ax;
    end

    pax = axes();
    pax.Position = [ax.Position(1)+0.05*ax.Position(3),ax.Position(2)+0.6*ax.Position(4),OptArgs.PieSize,OptArgs.PieSize];

    p = pie(pax,[NAshks,NBkgrd],[1,1]);
    p(1).FaceColor = colors(2,:);
    p(3).FaceColor = colors(1,:);
    p(1).FaceAlpha = patchFaceAlpha;
    p(3).FaceAlpha = patchFaceAlpha;
    [x1,y1] = centroid(polyshape(p(1).Vertices));
    [x2,y2] = centroid(polyshape(p(3).Vertices));
    p(2).Position = [x1-0.03,y1,0];
    p(4).Position = [x2+0.08,y2,0];
    p(2).FontSize = textFontSize;
    p(4).FontSize = textFontSize;
    p(2).FontWeight = 'b';
    p(4).FontWeight = 'b';
    p(2).Color = 'k';
    p(4).Color = 'k';
    p(2).HorizontalAlignment='center';
    p(4).HorizontalAlignment='center';

    set(pax,'Tag','Pie Chart');
end