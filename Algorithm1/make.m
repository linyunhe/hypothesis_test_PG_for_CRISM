function make(varargin)
%MAKE A utility to do the heavy lifting of Matlab MEX compilation
%
%   This function simplifies compiling MEX files for Matlab by automating
%   the process of running several compilation jobs with different files
%   with the same compiler flags. This information is always read from a
%   file called mexMakeFile in the current directory.
%   Different groups of compiler flags are specified by grouping them into
%   'targets.' These are specified in lines in mexMakeFile formatted as:
%   
%       !<targetName>: <flags>
%
%   When calling this function, the flags for this target are applied when
%   the call takes the form:
%   
%       make <targetName> [<possibly other targets>]
%
%   Multiple targets may be used, in which case the flags for all of the
%   targets will be concatenated together and added to each call to mex()
%   which this function makes.
%
%   The special target 'all' can be used in mexMakeFile to specify flags
%   which should always be added to mex() calls regardless of whether or
%   not any targets are specified.
%
%   All lines of mexMakeFile that don't begin with an exclamation point are
%   assumed to be a filename/list of filenames that should be passed to one
%   call of mex().
%
%   This function was inspired by the much-renowned Unix tool make, though
%   it is not nearly as flexible as its predecessor.

fid = fopen('mexMakeFile');
if fid==-1
    error('mexMakeFile not found in current directory');
end

commands = {};
targets = struct();

errnum = 0;
while errnum==0
    line = strtrim(fgetl(fid));
    if isempty(line)
        continue; % nothing here, skip
    end
    
    % Does this line specify a target? If so, this regex will match it
    targetLineMatch = regexp(line,'!(\S+):(.*)','tokens');
    if ~isempty(targetLineMatch)
        targetName = targetLineMatch{1}{1};
        targetOptions = targetLineMatch{1}{2};
        
        targets.(targetName) = targetOptions;
    else
        % must be file arguments for a mex() call
        command = strtrim(line);
        commands{end+1} = command;
    end
    
    [~,errnum] = ferror(fid); % -4 is EOF reached

end

fclose(fid);

% collect a string of all the flags specified by the targets the user wants
opts = '';
if isfield(targets,'all')
    varargin{end+1} = 'all'; % always run the target 'all' if it's supported
end
for selectedTarget = varargin
    selectedTarget = selectedTarget{1};
    if isfield(targets,selectedTarget)
        opts = [opts getfield(targets,selectedTarget)];
    else
        error('Current mexMakeFile does not support target ''%s''',selectedTarget);
    end
end

% run each command
for filename = commands
    filename = filename{1};
    
    command = ['mex ',filename,' ',opts];
    
    fprintf('* Now building %s ...\n',filename);
    fprintf('Running command ''%s'' ...\n',command);
    
    % run the thing
    eval(command);
end

end