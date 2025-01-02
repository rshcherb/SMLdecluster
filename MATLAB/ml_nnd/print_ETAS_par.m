function [ETASout] = print_ETAS_par(ETASpar,OptArgs)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        ETASpar table
        OptArgs.delimiter char = ', '
    end

    ETASparnames = string(ETASpar.Properties.VariableNames);
    greekVars = matches(ETASparnames,["mu","alpha","gamma"]);
    ETAStxt = ETASparnames;
    if any(greekVars)
        ETAStxt(greekVars) = "\" + ETASparnames(greekVars);
    end
    
    ETASout = strings(1,length(ETASparnames));
    if size(ETASpar,1) > 1
        for i = 1:length(ETASparnames)
            ETASout(i) = ETAStxt(i) + "_{av}" + " = " + ...
                    num2str(mean(ETASpar.(ETASparnames(i))));
        end
    else
        for i = 1:length(ETASparnames)
            ETASout(i) = ETAStxt(i) + " = " + ...
                    num2str(ETASpar.(ETASparnames(i)));
        end
    end
                
%     ETASout = strings(size(ETASpar,1),length(ETASparnames));
%     for i = 1:length(ETASparnames)
%         for j = 1:size(ETASpar,1)
%             ETASout(j,i) = ETAStxt(i) + " = " + ...
%                     num2str(ETASpar.(ETASparnames(i))(j));
%         end
%     end
    
    ETASout = join(ETASout,OptArgs.delimiter);

end