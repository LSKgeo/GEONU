function [table1] = createLatexTable1(x, h)

% Always put a '&' after so it goes to the next value
for i = 1:size(h.rows,1)
    % Row Name
if i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 7 || i == 8 %push it over 1 column
    x(i).row = sprintf(' & %s &',h.rows{i}); 
else 
    x(i).row = ''; % for DM,EM, and BSE
end

    % Density
x(i).rho = sprintf('& %0.2f &',h.rho(i,1));

    % Thickness
    if i == 9  % DM
            x(i).thick = sprintf(' %0.0f $%s$%0.0f &',h.thick(i,1),'\pm',h.thick(i,2));
    elseif i == 10 || i == 11
            x(i).thick = sprintf(' %0.0f &',h.thick(i,1));
    else    
            x(i).thick = sprintf(' %0.2f $%s$%0.1f &',h.thick(i,1),'\pm',h.thick(i,2));
    end
    
    % Mass
x(i).mass.mass = sprintf(' %0.1f $%s$%0.1f &',h.mass.mass(i,1),'\pm',h.mass.mass(i,2));    

    % U abundance % THIS ONE
x(i).abund.U = latexUncSymAsym(h.abund.U(i,:),2);
    
    % Th abundance
x(i).abund.Th = latexUncSymAsym(h.abund.Th(i,:),2);    

    % K abundance
x(i).abund.K = latexUncSymAsym(h.abund.K(i,:),2);
 
    % U Mass
x(i).mass.U = latexUncSymAsym(h.mass.U(i,:),1);

    % Th Mass
x(i).mass.Th = latexUncSymAsym(h.mass.Th(i,:),1);

    % K Mass
x(i).mass.K = latexUncSymAsym(h.mass.K(i,:),1);

    % Heat Production
x(i).hp = latexUncSymAsym(h.hp(i,:),1,1);

% Ammend latex format stuff onto front of strings
if i == 1
   x(i).front = '\multicolumn{1}{|c|}{\multirow{6}{*}{CC}}';
elseif i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 8
   x(i).front = '\multicolumn{1}{|c|}{} ';
elseif i == 7
   x(i).front = '\multicolumn{1}{|c|}{\multirow{2}{*}{OC}}';
elseif i == 9
   x(i).front = '\multicolumn{2}{|c}{DM} & '; 
elseif i == 10
   x(i).front = '\multicolumn{2}{|c}{EM} & ';   
elseif i == 11
   x(i).front = '\multicolumn{2}{|c}{BSE} & ';  
end

% Amend latex format stuff onto end of strings
if i == 1 || i == 2 || i == 3 || i == 4|| i == 5 || i == 7
    x(i).end = '\\ \cline{2-13}'; 
else 
    x(i).end = '\\ \hline';
end


end

y = strcat(x.front,x.row,x.rho,x.thick,x.mass.mass,x.abund.U,...
   x.abund.Th,x.abund.K,x.mass.U,x.mass.Th,x.mass.K,...
   x.hp,x.end); 


% Create Name of file
str.geo = strcat('tgeo_',str1,'.txt');

% Ammend final lines with caption of information (so we can identify the table)
final{1,1} = '\end{tabular}}';
x = strrep(str.geo,'_',' '); x = strrep(x,'+','$+$'); 
final{2,1} = sprintf('\\caption{%s}',x);


% Write to text file
fid = fopen(str.geo,'wt');
fprintf(fid, '%s\n', y{1,1});
fclose(fid);

fid = fopen(str.geo,'at');
fprintf(fid, '%s\n', y{2:end,1});
%fclose(fid);

fid = fopen(str.geo,'at');
fprintf(fid,'%s\n',final{:,:});

fclose(fid);    
fclose('all'); % Fully releases the file (otherwise Matlab won't for some reason)
    

disp('5) Saved Table: Geo-Physical/Chemical Results')
fprintf('     %s \n',str.geo)
end