function gpreproc(pathname,filename,variables,c1,loc)
% general preprocessing routine for ascii and mat files (like load)
% function gpreproc(pathname,filename,variables,c1)
% 
%
% pathname = drive name to target folder, ex. 'G:\awiggins\VCTA\Cases\Exp4\Low RE\0916tests\'
%            will be prompted if an empty array [] is used
%
% filename = filename including extension. ex. 'exp4_100_00.ASC'
%            will be prompted if an empty array [] is used. Directories will be filtered by extension
%            (.txt,.mat.,.asc) when left empty.
% 
% variable = cell array of strings containing variable names
%            i.e. variables = {'lift','drag,'poo','etc'};
%
% c1 = program option
%       c1 = 1
%               uses variables contained in 'variables' for assignments
%       c1 = 2
%               user selected variable names
%       c1 = 3 or anything else non-empty
%               default column names such as col1,col2,col3
%
% Note: If variables is defined, c1 is automatically set to 1 
%--------------obtain filename if not input-----------
% if isempty(pathname) ==  1 | isempty(filename) == 1
%     [filename, pathname] = uigetfile( ...
%         {'*.txt','Labview ASCII (*.txt)';...
%             '*.asc','Dasylab ASCII (*.asc)';...
%             '*.mat','Matlab Data Matrix (*.mat)'; ...   
%             '*.comp','Run particulars (*.comp)'}, ... 
%         'Load Time Series Data');
% end
% %------------check path and filename validity---------
% if isequal(filename,0)|isequal(pathname,0)
%     disp('File not found')
%     return %needed if canceled once routine is started
% else
%     disp(['File ', pathname, filename, ' found'])
% end

% loc = 1, variables will be in base workspace
% loc = 0, variables will be in caller workspace, base is default

if isempty(pathname) ==  1 | isempty(filename) == 1
    [pathname,filename,found] = getfile('*.*');
end

fullfilename = fullfile(pathname,filename);
[dir,fname,ext,ver] = fileparts(filename);

if ~isempty(variables)
    c1 = 1;
end

if isempty(loc) 
    loc = 1;
end

if loc == 0
    space = 'caller';
else
    space = 'base';
end


    

if ~isempty(strmatch(ext,strvcat('.mat','.comp'),'exact')) 
    
    fid = fopen(fullfilename,'r'); % input file
    
    if fid == -1
        disp('File could not be opened')
        return
    end
    
    fclose(fid);
    
    
    load(fullfilename,'-mat');
    S = whos('-file',fullfilename);  
    names = {S.name};
    
    
    if isempty(c1)== 1
        c1 = input(['There are ' num2str(length(names)) ' variables in this file, would you like to... \n','1) Output variables listed in program \n', '2) List file contents and choose variables \n', '3) Spit them all out \n']);
    end
    
    switch c1
        case 1
            if isempty(variables) == 1
                variables = {'cylOD','command','blendedpos','hippo'}; 
            end
            
            for j = 1:length(variables)
                
                if isempty(strmatch(variables{j},names,'exact')) == 1
                    disp([variables{j} ' not found in file, set c1 = 2 and list file contents'])
                    variables = variables(1:j-1);
                    break
                end
            end
            
        case 2
            
            disp('Here are the variables in your file, type in the variables you want to export, press enter if you are done')
            disp(names)
            
            for k = 1:length(names)
                variables{k} = input('Variable to export ','s');
                if max(size(variables{k})) == 0
                    variables = variables(1:k-1);
                    break
                end
            end
            
        otherwise
            
            variables = names; 
            
    end
    
    
    for i = 1:length(variables)
        assignin(space,variables{i},eval(variables{i}));
    end
    
    % evalin('base',) 
    %     S = whos('-file',fullfilename);
    %     names = {S.name};
    %     for i = 1:length(names)
    %         temp(i) = max(S(i).size);
    %     end
    %     
    %     testcase = char(names(find(temp == max(temp))));
else 
    
    fid = fopen(fullfilename,'r'); % input file
    
    if fid == -1
        disp('File could not be opened')
        return
    end
    
    testfirst = []; 
    no_lines = 0;
    max_line = 0;
    line = fgetl(fid);
    [testfirst, ncols, errmsg, nxtindex] = sscanf(line, '%f');
    
    while isempty(testfirst)|(nxtindex==1) %this searches for the first line of data
        max_line = max([max_line, length(line)]);
        eval(['line', num2str(no_lines), '=line;']);
        firstline = fgetl(fid); % obtains one line
        [testfirst, ncols, errmsg, nxtindex] = sscanf(firstline, '%f');
    end  
    colnum = ncols;
    fseek(fid,-nxtindex+1,'cof');
    % helpdlg(warning,'File Scan Complete');
    
    
    
    if isempty(c1)== 1
        c1 = input(['There are ' num2str(colnum) ' columns in this data, would you like to... \n','1) Use names from program \n', '2) Assign here \n', '3) Use Default (i.e. col1, col2 etc) \n']);
    end
    
    
    
    switch c1
        case 1
            if isempty(variables) == 1
                
                variables = {'time','pos','drag1','lift2'}; 
            end
            
            if colnum == length(variables)
                colnames = variables;
            end
            
            if colnum > length(variables)
                disp('More columns than variable names, only columns with variables will be processed')
                colnum = length(variables)
            end 
            
            if colnum < length(variables)
                disp('More variables names than columns')
            end
            
            
            
        case 2
            
            for i = 1:colnum
                colnames{i} = input(['What is column ' num2str(i) ' ? '],'s');
            end
            
            
        otherwise 
            
            for i = 1:colnum
                col  = ['col_' num2str((i))];   
                colnames{i} = sprintf(col,'%s');
            end
            
    end
    
    
    fmt = '';
    
    for i = 1:colnum
        format = '%f ';
        fmt = [fmt format];
    end
    
    
    fmt = [fmt '\n'];
    
    
    bigX = fscanf(fid,fmt,[colnum,inf])';
    
    %     for j = 1:length(colnames) %should be bursts 
    %         %callervariable = ['var', int2str(k)];
    %         assignin('base',colnames{j}, bigX(:,j)); 
    %     end
    S =evalin(space,'who');
    
    for j = 1:colnum
        
        if isempty(strmatch(colnames{j},S))==1
            assignin(space,colnames{j},[] );
        end 
        
        temp = evalin(space,colnames{j});
        
        if isempty(temp) == 0
            endpoint = min([length(temp) length(bigX(:,j))]);
            
            assignin(space,colnames{j}, [temp(1:endpoint) bigX((1:endpoint),j)]); 
        else
            assignin(space,colnames{j}, [temp bigX(:,j)]); 
            
        end
        
    end
    
    fclose(fid);
    
    
end

