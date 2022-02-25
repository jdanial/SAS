function loadParamAuxFn(app)
% loadParamAuxFn() -
% loads BAS parameters.
%
% Syntax -
% exportFn(app).
%
% Parameters -
% - app: BAS UI class

% saving parameters
[importFile,importPath] = uigetfile('*.sp');

% replacing extension with .mat
newFilePath = strrep(fullfile(importPath,importFile),'sp','.mat');

% copying file
copyfile(fullfile(importPath,importFile),newFilePath,'f');

% loading data
try
    paramTemp = load(newFilePath);
catch
    app.msgBox.Value = sprintf('%s',['Error: cannot load data from (' importFile ').']);
    return;
end

% deleting copied file
delete(newFilePath);

% copying paramTemp
app.param = paramTemp.param;

% listing parameters
app.DetectSwitch.Value = app.param.detection.detect;
app.MaximumSigmaEditField.Value = app.param.detection.maxSigma;
app.LocalizeCheckBox.Value = app.param.detection.localize;
app.ROIRadiusEditField.Value = app.param.detection.roiRadius;
app.AnalyzeSwitch.Value = app.param.analysis.analyze;
app.NumSubUnitsPerCalibComplexEditField.Value = app.param.analysis.numSubUnitsPerCalibComplex;
app.TimeSliceEditField.Value = app.param.analysis.timeSlice;
app.BinSizeEditField.Value = app.param.analysis.binSize;
end