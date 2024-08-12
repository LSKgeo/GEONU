function tab2 = createLatexTable2(x,h)
% Always put a '&' after so it goes to the next value
for i = 1:size(h.rows,1)

    if i == 1 || i == 6 
        x(i).row = sprintf('%s &',strrep(h.rows{i},'_',' '));
    else
        x(i).row = sprintf('%s &',h.rows{i}); 
    end

    if i == 1 || i == 6 || i == 7
    x(i).sU = latexUncSymAsym(h.U238(i,:),2);
    x(i).sTh = latexUncSymAsym(h.Th232(i,:),2);
    x(i).sTotal = latexUncSymAsym(h.total(i,:),2,1);
    else
        x(i).sU = latexUncSymAsym(h.U238(i,:),1);
        x(i).sTh = latexUncSymAsym(h.Th232(i,:),1);
        x(i).sTotal = latexUncSymAsym(h.total(i,:),1,1);
    end


if i == size(h.rows,1)
    x(i).end = '\\ \hline';
else
    x(i).end = '\\';
end
        
   

end

y = strcat(x.row,x.sU,x.sTh,x.sTotal,x.end); 

% Create Name of file
str.flux = strcat('tflux_',str1,'.txt');

% Ammend Initial lines with detector
z = MASTER.detector;
start{1,1} = sprintf('\\multirow{2}{*}{} & \\multicolumn{3}{c}{%s} \\\\',z.Properties.RowNames{1});
start{2,1} = sprintf('& \\multicolumn{3}{c}{%0.2f\\textdegree, %0.2f\\textdegree} \\\\ \\hline',z{1,2},z{1,1});
start{3,1} = ' & S(U) & S(Th) & S(U+Th) \\'; 

% Ammend final lines with caption of information (so we can identify the table)
final{1,1} = '\end{tabular}}';
z = strrep(str.flux,'_',' '); z = strrep(z,'+','$+$'); 
final{2,1} = sprintf('\\caption{%s}',z);



% Write to text file
fid = fopen(str.flux,'wt');
fprintf(fid,'%s\n',start{:,:}); 
fclose(fid);

fid = fopen(str.flux,'at');
fprintf(fid, '%s\n', y{:,1});
%fclose(fid);

fid = fopen(str.flux,'at');
fprintf(fid,'%s\n',final{:,:});

fclose(fid);    
fclose('all'); % Fully releases the file (otherwise Matlab won't for some reason)
disp('6) Saved Table: Flux Results')
fprintf('     %s \n',str.flux)



disp(' ')
disp(' ')
disp('------------------------------------------------------')
disp('---------------------- Run End -----------------------')
fprintf('End date and time:         %s \n',datetime('now','format','mmmm dd, yyyy HH:MM AM')); 
fprintf('Total Run time (min):      %0.2f \n',toc(tstart)/60)
disp(' ')
disp('------------------------------------------------------')


% Turn Diary off
diary off
