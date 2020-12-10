function plabel = pitchlabel(F_q, stdfreq)
    [~, I] = find(F_q == stdfreq);
    
    plabel = {'A  ', 'B  ', 'H  ', 'C  ', 'Des', 'D  ', 'Es ', 'E  ', 'F  ', 'Ges', 'G  ', 'As '};
    plabel = repmat(plabel, [1 8]);
    plabel = [plabel{end-mod(I, 12)+2:end} plabel];
    
    for l=1:length(F_q)
        plabel{l} = strcat(num2str(F_q(l)), {' - '}, plabel(l));
        plabel{l} = plabel{l}{1};
    end
    
    plabel(l+1:end) = [];
end