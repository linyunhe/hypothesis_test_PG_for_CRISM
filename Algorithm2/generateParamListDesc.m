function [ desc_text ] = generateParamListDesc( params )
%GENERATEPARAMLISTDESC Builds description text including all params
%   Creates a description containing all parameter values. For inclusion in
%   the ENVI-style headers of data produced by the algorithm

newline = sprintf('\r\n');

desc_text = [ 'File Generated in CRISM MLM algorithm.' newline newline];

desc_text = [ desc_text 'Parameters:' newline ];

for field = fieldnames(params)'
    desc_text = [desc_text paramValueReport(params,field{1}) newline];
end

end

function txt = paramValueReport(params,field)
    value = params.(field);
    if isstruct(value)
        % we must go deeper (this will look odd if there is more than one
        % layer of substructures under params-proper)
        valueAsStr = '{';
        for innerfield = fieldnames(value)'
            valueAsStr = [valueAsStr paramValueReport(value,innerfield{1}) ', '];
        end
        valueAsStr = [valueAsStr(1:end-2) '}']; % remove trailing comma before closing brace
    else
        % typical case: just turn value into a string and print the pair
        valueAsStr = mat2str(value);
    end

    txt = sprintf('%s: %s',field,valueAsStr);
    
end
    