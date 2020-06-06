% A more convenient directory listing command.  Only outputs the directories within the passed folder.
% Files and directories starting with '.' are ignored.
% Output is a cell array of directory names.
% Cross Platform.
function fileNames = list_directories(parameterString)

if nargin==0
    parameterString = '.';
end

temp = dir(parameterString);
fileNames = {};
[fileNames{1:length(temp),1}] = deal(temp.name);

% Find hidden files.
hiddenFiles = nan(size(fileNames));
for fileIndex = 1:length(fileNames)
    tempName = fileNames{fileIndex};
    hiddenFiles(fileIndex) = tempName(1) == '.';
end

% Find which entries in the listing are directories.
directoryFlags = {};
[directoryFlags{1:length(temp),1}] = deal(temp.isdir);
directoryFlags = cell2mat(directoryFlags);

tempFlags = (~hiddenFiles) & directoryFlags;
fileNames = fileNames(tempFlags);
