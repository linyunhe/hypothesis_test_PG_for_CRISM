function [ desc_text ] = generateParamListDesc( params, version, algo, bounds )
%GENERATEPARAMLISTDESC Builds description text including all params
%   Creates a description containing all parameter values. For inclusion in
%   the ENVI-style headers of data produced by the algorithm

newline = sprintf('\r\n'); % TODO: should this be operating system dependent?

desc_text = [ 'File Generated in CRISM ' algo ' algorithm, version ' version newline newline];

row_min = bounds.rowBounds(1); row_max = bounds.rowBounds(2);
col_min = bounds.colBounds(1); col_max = bounds.colBounds(2);
band_min = bounds.bandBounds(1); band_max = bounds.bandBounds(2);

desc_text = [ desc_text ...
    sprintf('Ran rows %i to %i (%i), cols %i to %i (%i), bands %i to %i (%i) (all 1-indexed).\n', ...
        row_min,row_max,row_max-row_min+1, ...
        col_min,col_max,col_max-col_min+1, ...
        band_min,band_max,band_max-band_min+1 ...
    ),...
    newline ...
];

desc_text = [ desc_text 'Parameters used:' newline ];

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
    