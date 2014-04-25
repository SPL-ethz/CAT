function makeDocs

%% Make Docs
%
% Make documentation for publishing on web/pdf

% Make PDF files
opts.format = 'pdf';

% Put in the docs folder
opts.outputDir = [fileparts(mfilename('fullpath')) filesep 'docs'];

% Don't show or evaluate code
opts.showCode = false;
opts.evalCode = false;

codefile = [fileparts(mfilename('fullpath')) filesep 'Install.m'];

fprintf('Creating %s file for %s in directory %s\n',opts.format,codefile,opts.outputDir);

publish(codefile,opts);