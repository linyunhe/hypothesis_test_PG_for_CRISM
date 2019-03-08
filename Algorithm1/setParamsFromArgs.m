function [ params ] = setParamsFromArgs( args, params )
%SETPARAMSFROMARGS Update structure fields in params using named args pairs
%   Expects pairs of cells in args that go 'FieldName','FieldValue'. The
%   value of the field can be a string or a double.

% SAFTIES THIS DOESN'T HAVE YET:
% 1. enforcing types for parameters

nargs = length(args);

for i=1:2:nargs-1
    % We expect that a parameter name is at i, and the value is at i+1
    
    if ~ischar(args{i}) % if paramter name isn't a string, abort
        error('Expected argument #%i to be a string (command field name): instead it''s a %s',...
            i,class(args{i}));
    end
    
    if isfield(params,args{i})
        params.(args{i}) = args{i+1};
    else
        % Check one layer deep into nested structs in params
        match = regexp(args{i},'^(\w+)\.(\w+)$','tokens'); % TODO: what if it doesn't match this pattern?
        substruct_name = match{1}{1};
        substruct_member = match{1}{2};
        
        if isfield(params.(substruct_name),substruct_member)
            params.(substruct_name).(substruct_member) = args{i+1};
        else
            error('Parameter "%s" does not exist; it cannot be set through function arguments!',args{i});
        end
    end
end

