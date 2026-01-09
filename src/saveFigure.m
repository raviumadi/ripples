function saveFigure(figHandle, destDir, baseName, saveFigFile)
    % saveFigure Save figure as .pdf (and optionally .fig) to destination
    %
    %   saveFigure(figHandle, destDir, baseName)
    %   saveFigure(figHandle, destDir, baseName, saveFigFile)
    %
    %   figHandle     - handle to figure (use gcf for current)
    %   destDir       - string, path to destination folder
    %   baseName      - base filename without extension
    %   saveFigFile   - true/false (default: true)

    if nargin < 1 || isempty(figHandle)
        figHandle = gcf;
    end
    if nargin < 2 || isempty(destDir)
        destDir = pwd;
    end
    if nargin < 3 || isempty(baseName)
        baseName = ['figure_' datestr(now, 'yyyymmdd_HHMMSS')];
    end
    if nargin < 4
        saveFigFile = true;
    end

    % Ensure directory exists
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end

    % Construct filenames
    figFile = fullfile(destDir, [baseName '.fig']);
    pdfFile = fullfile(destDir, [baseName '.pdf']);

    % Try saving .fig if requested
    if saveFigFile
        try
            savefig(figHandle, figFile);
            fprintf('Saved: %s\n', figFile);
        catch ME
           warning(ME.identifier, 'Could not save .fig file: %s', ME.message);
        end
    end

    % Save as vector PDF
    try
        exportgraphics(figHandle, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'none', ...
            'Resolution', 300);
        fprintf('Saved: %s\n', pdfFile);
    catch ME
       warning(ME.identifier, 'Could not save .fig file: %s', ME.message);
    end
end