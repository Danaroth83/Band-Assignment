function matrix2latex_mod3(matrix, filename, varargin)

% function: matrix2latex(...)
% Author:   M. Koehler
% Contact:  koehler@in.tum.de
% Version:  1.1
% Date:     May 09, 2004

% This software is published under the GNU GPL, by the free software
% foundation. For further reading see: http://www.gnu.org/licenses/licenses.html#GPL

% Usage:
% matrix2latex(matrix, filename, varargs)
% where
%   - matrix is a 2 dimensional numerical or cell array
%   - filename is a valid filename, in which the resulting latex code will
%   be stored
%   - varargs is one ore more of the following (denominator, value) combinations
%      + 'rowLabels', array -> Can be used to label the rows of the
%      resulting latex table
%      + 'columnLabels', array -> Can be used to label the columns of the
%      resulting latex table
%      + 'alignment', 'value' -> Can be used to specify the alginment of
%      the table within the latex document. Valid arguments are: 'l', 'c',
%      and 'r' for left, center, and right, respectively
%      + 'format', 'value' -> Can be used to format the input data. 'value'
%      has to be a valid format string, similar to the ones used in
%      fprintf('format', value);
%      + 'size', 'value' -> One of latex' recognized font-sizes, e.g. tiny,
%      HUGE, Large, large, LARGE, etc.
%
% Example input:
%   matrix = [1.5 1.764; 3.523 0.2];
%   rowLabels = {'row 1', 'row 2'};
%   columnLabels = {'col 1', 'col 2'};
%   matrix2latex(matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%
% The resulting latex file can be included into any latex document by:
% /input{out.tex}
%
% Enjoy life!!!


    
    width = size(matrix, 2);
    height = size(matrix, 1);

    rowLabels = [];
    colLabels = [];
    alignment = 'l';
    format = [];
    textsize = [];
    highlight='n';
    groups_r=ones(1,height);
    groups_c=ones(1,width);
    groups_hr=height;
    groups_hc=width;
    if (rem(nargin,2) == 1 || nargin < 2)
        error(['matrix2latex: ', 'Incorrect number of arguments to %s.'], mfilename);
    end
    
    okargs = {'rowlabels','columnlabels', 'alignment', 'format', 'size','highlight','gRow','gCol','hRow','hCol'};
    for j=1:2:(nargin-2)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,4));
        if isempty(k)
            error(['matrix2latex: ', 'Unknown parameter name: %s.', pname]);
        elseif length(k)>1
            error(['matrix2latex: ', 'Ambiguous parameter name: %s.', pname]);
        else
            switch(k)
                case 1  % rowlabels
                    rowLabels = pval;
                    if isnumeric(rowLabels)
                        rowLabels = cellstr(num2str(rowLabels(:)));
                    end
                case 2  % column labels
                    colLabels = pval;
                    if isnumeric(colLabels)
                        colLabels = cellstr(num2str(colLabels(:)));
                    end
                case 3  % alignment
                    alignment = lower(pval);
                    if strcmp(alignment,'right')
                        alignment = 'r';
                    end
                    if strcmp(alignment,'left')
                        alignment = 'l';
                    end
                    if strcmp(alignment,'center')
                        alignment = 'c';
                    end
                    if alignment ~= 'l' && alignment ~= 'c' && alignment ~= 'r'
                        alignment = 'l';
                        warning(['matrix2latex: ', 'Unknown alignment. (Set it to \''left\''.)']);
                    end
                case 4  % format
                    format = lower(pval);
                case 5  % format
                    textsize = pval;
                case 6  % highlight
                    highlight=pval;
                    if highlight(1)~='c' && highlight(1)~='r' && highlight(1)~='n'
                        error('First character must be ''c'' for highlighting columns, ''r'' for rows, ''n'' for no alignment');
                    else
                        strrep(highlight(2:end),'l','s');
                        strrep(highlight(2:end),'b','g');
                        if ~all(highlight(2:end)=='s' | highlight(2:end)=='g')
                            error('Use ''g'' for highlighting greatest value, ''s'' for smallest');
                        end
                    end
                case 7 % Grouping of rows
                    if sum(pval)~=height, error('Row grouping is wrong'); end
                    groups_r=ismember(1:height,cumsum(pval));
                case 8 % Grouping of columns
                    if sum(pval)~=width, error('Column grouping is wrong'); end
                    groups_c=ismember(1:width,cumsum(pval));
                case 9 % Grouping of highlights for rows
                    if sum(pval)~=height, error('Row highlighting grouping is wrong'); end
                    groups_hr=pval;
                case 10 % Grouping of highlights for columns
                    if sum(pval)~=width, error('Column highlighting grouping is wrong'); end
                    groups_hc=pval;                   
            end
        end
    end
    
    if (highlight(1)=='c' && sum(groups_hc)~=width) || (highlight(1)=='r' && sum(groups_hr)~=height)
        error('Highlight grouping is wrong');
    end
    if (highlight(1)=='c' && length(groups_hc)~=length(highlight)-1) || (highlight(1)=='r' && length(groups_hr)~=length(highlight)-1)
        error('Highlighting dimensions are not correct');
    end

    fid = fopen(filename, 'a');
    
    matrix_orig=matrix;

    if isnumeric(matrix)
        matrix = num2cell(matrix);
        for h=1:height
            for w=1:width
                if(~isempty(format))
                    matrix{h, w} = num2str(matrix{h, w}, format);
                else
                    matrix{h, w} = num2str(matrix{h, w});
                end
            end
        end
    end
    
    if(~isempty(textsize))
        fprintf(fid, '\\begin{%s}', textsize);
    end

    fprintf(fid, '\\begin{tabular}{|');

    if(~isempty(rowLabels))
        fprintf(fid, 'l|');
        fprintf(fid, '%s',repmat('l|',[1,length(strfind(colLabels{1},'&'))]));
    end
    for i=1:width
            fprintf(fid, '%c', alignment);
        if groups_c(i)==1, fprintf(fid, '|'); end
    end
    fprintf(fid, '}\r\n');
    
    fprintf(fid, '\\hline\r\n');
    
    if(~isempty(colLabels))
        if(~isempty(rowLabels))
            fprintf(fid, '&');
        end
        for w=1:width
            title=colLabels{w};
            idx=strfind(title,'&');
            if ~isempty(idx), idx=idx(end); else idx=0; end
            if isempty(regexp(title(idx+1:end),'[_^=]','once'))
                fprintf(fid, [title(1:idx),'\\textbf{%s}'], title(idx+1:end));
            else
                fprintf(fid, [title(1:idx),'$\\mathbf{%s}$'], title(idx+1:end));
            end
            if w==width
                fprintf(fid, '\\\\\\hline\r\n');
            else
                fprintf(fid, '&');
            end
        end
    end
    
    highlight_matrix=zeros(height,width);
    highlight_matrix2=zeros(height,width);
    start_groups_hr=[1,cumsum(groups_hr)+1]; end_groups_hr=start_groups_hr(2:end)-1;
    start_groups_hc=[1,cumsum(groups_hc)+1]; end_groups_hc=start_groups_hc(2:end)-1;
    if ~strcmpi(highlight,'n')
        for i1=1:length(groups_hr)
            idx1=start_groups_hr(i1):end_groups_hr(i1);
            for i2=1:length(groups_hc)
                idx2=start_groups_hc(i2):end_groups_hc(i2);
                matrix_block=matrix_orig(idx1,idx2);
                if highlight(1)=='r'
                    ordering=highlight(i1+1);
                else
                    ordering=highlight(i2+1);
                end
                if ordering=='g'
                    matrix_ordered=sort(matrix_block(:),'descend');
                else
                    matrix_ordered=sort(matrix_block(:),'ascend');
                end
                if numel(matrix_block)>1
                    highlight_matrix(idx1,idx2)=(matrix_block==matrix_ordered(1));
                end
                if numel(matrix_block)>2
                    highlight_matrix2(idx1,idx2)=(matrix_block==matrix_ordered(2));
                end
            end
        end
    end
            
    
    for h=1:height
        if(~isempty(rowLabels))
            title=rowLabels{h};
            idx=strfind(title,'&');
            if ~isempty(idx), idx=idx(end); else idx=0; end
            if isempty(regexp(title(idx+1:end),'[_^=]','once'))
                fprintf(fid, [title(1:idx),'\\textbf{%s}&'], title(idx+1:end));
            else
                fprintf(fid, [title(1:idx),'$\\mathbf{%s}$&'], title(idx+1:end));
            end
        end
        for w=1:width
            if highlight_matrix(h,w)==1
                fprintf(fid, '\\textbf{%s}', matrix{h, w});
            elseif highlight_matrix2(h,w)==1
                fprintf(fid, '\\uline{%s}', matrix{h, w});
            else
            	fprintf(fid, '%s', matrix{h, w});
            end
            if w==width
                if groups_r(h)==1
                    fprintf(fid,'\\\\\\hline\r\n');
                else
                    fprintf(fid,'\\\\\r\n');
                end
            else
                fprintf(fid,'&');
            end
        end
    end

    fprintf(fid, '\\end{tabular}\r\n');
    
    if(~isempty(textsize))
        fprintf(fid, '\\end{%s}', textsize);
    end

    fclose(fid);