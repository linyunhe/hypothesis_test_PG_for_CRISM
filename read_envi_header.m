function [ headerData ] = read_envi_header( filename )
%READ_ENVI_HEADER Extracts simple line-items from an ENVI header
%   Extracts the fields samples,lines,bands,interleave,and data type.
%   We might eventually be interested in the offset as well, but not
%   currently supported.
fid = fopen(filename);

if fid==-1
    error('The file %s does not exist!\n',filename);
end

totalData = '';

% read lines until file is consumed
errnum = 0;
while errnum == 0
    nextLine = fgetl(fid);
    
    if isempty(nextLine) || nextLine(1)~=';' % ignore ENVI full-line comments
        totalData = [totalData,nextLine,sprintf('\r')];
    end
    
    [~,errnum] = ferror(fid); % -4 is EOF reached
end
fclose(fid);

totalData = [totalData sprintf('\rEND')];

headerData.samples = getValueForKey(totalData,'samples','\d+');
headerData.lines = getValueForKey(totalData,'lines','\d+');
headerData.bands = getValueForKey(totalData,'bands','\d+');
headerData.interleave = getValueForKey(totalData,'interleave','\w{3}');

dataTypeNum = getValueForKey(totalData,'data type','\d+');

% get name of data type from the numeric code given in the header
switch dataTypeNum
    case 4
        headerData.dataType = 'float'; % not 'single', because these types meant for multibandread
    case 5
        headerData.dataType = 'double';
    otherwise
        error('Header has unknown data type %i. You might consider modifying this function to support this type.', dataTypeNum);
end
        
end

function val = getValueForKey(lines,key,valuePattern)
    captured = regexp(lines,sprintf('%s\\W*=\\W*(%s)',key,valuePattern),'tokens');
    val = str2double(captured{1}); % try parsing as number (will be NaN if it fails)
    
    if isnan(val) % not number; extract string value from cell
        val = captured{1}{1};
    end
end