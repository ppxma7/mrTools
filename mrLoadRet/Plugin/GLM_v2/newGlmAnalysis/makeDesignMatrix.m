% makeDesignMatrix.m
%
%        $Id$
%      usage: makeDesignMatrix(d,params,verbose)
%         by: farshad moradi, modified by julien besle
%       date: 06/14/07, 11/02/2010
%       e.g.: makeDesignMatrix(d,params,verbose)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it, as well as a model hrf
%              optional parameters can be passed in the params structure 
%
function d = makeDesignMatrix(d,params,verbose, scanNum)

if ~any(nargin == [1 2 3 4 5])
   help makeDesignMatrix;
   return
end

if ieNotDefined('params')
  params=struct;
end
if fieldIsNotDefined(params,'scanParams')
  params.scanParams{1}=struct;
  scanNum=1;
end
if ieNotDefined('verbose')
  verbose = 1;
end

% check if the hrf starts from zero (except if it is the identity matrix, the deconvolution case)
if verbose && any(d.hrf(1,:)>.05) && (size(d.hrf,1)~=size(d.hrf,2) || ~isempty(find(d.hrf^-1-d.hrf>1e-6, 1)))
   mrWarnDlg(['(makeDesignMatrix) HRF does not start from zero (hrf(0) = ' mat2str(d.hrf(1,:)) ')']);
end

if isfield(params,'nonLinearityCorrection') && params.nonLinearityCorrection && isfield(params.hrfParams,'maxModelHrf')
   saturationThreshold = params.saturationThreshold*params.hrfParams.maxModelHrf;
else
   saturationThreshold = Inf(1,size(d.hrf,2));
end

if ~fieldIsNotDefined(params.scanParams{scanNum},'estimationSupersampling')
  estimationSupersampling = params.scanParams{scanNum}.estimationSupersampling;
else
  estimationSupersampling = 1;
end
if ~fieldIsNotDefined(params.scanParams{scanNum},'acquisitionDelay')
  acquisitionDelay = params.scanParams{scanNum}.acquisitionDelay;
else
  acquisitionDelay = d.tr/2;
end
if ~fieldIsNotDefined(params.scanParams{scanNum},'acquisitionDuration')
  acquisitionDuration = params.scanParams{scanNum}.acquisitionDuration;
else
  acquisitionDuration = d.tr;
end

if isfield(params.scanParams{scanNum},'stimToEVmatrix') && ~isempty(params.scanParams{scanNum}.stimToEVmatrix)
  %match stimNames in params to stimNames in structure d
  [isInMatrix,whichStims] = ismember(d.stimNames,params.scanParams{scanNum}.stimNames);
  stimToEVmatrix = zeros(length(d.stimvol),size(params.scanParams{scanNum}.stimToEVmatrix,2));
  stimToEVmatrix(isInMatrix,:) = params.scanParams{scanNum}.stimToEVmatrix(whichStims(isInMatrix),:);
  if size(stimToEVmatrix,1)~=length(d.stimvol)
    mrWarnDlg('(makeDesignMatrix) EV combination matrix is incompatible with number of event types');
    d.scm = [];
    return;
  end
else
  stimToEVmatrix = eye(length(d.stimvol));
end

% if we have only a single run then we set
% the runTransitions for that single run
if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
  runTransition = [1 d.dim(4)];
else
  runTransition = d.concatInfo.runTransition;
end

numberSamples = diff(runTransition,1,2)+1;
runTransition(:,1) = ((runTransition(:,1)-1)*round(d.designSupersampling)+1);
runTransition(:,2) = runTransition(:,2)*round(d.designSupersampling);

%apply duration and convert to matrix form
stimMatrix = stimCell2Mat(d.stimvol,d.stimDurations,runTransition);
%if design sampling is larger than estimation sampling, we need to correct the amplitude of the hrf 
stimMatrix = stimMatrix*estimationSupersampling/d.designSupersampling; 

% apply EV combination matrix
d.EVmatrix = stimMatrix*stimToEVmatrix;
downsampleFactor = d.designSupersampling/estimationSupersampling;
acquisitionSubSamples = round(rem(acquisitionDelay-acquisitionDuration/2,d.tr/estimationSupersampling)*downsampleFactor/d.tr)+1:...
                        round(rem(acquisitionDelay+acquisitionDuration/2,d.tr/estimationSupersampling)*downsampleFactor/d.tr-eps*2);
if isempty(acquisitionSubSamples) % this means the acquisition period is between two subsamples and
  %is less than 1 subsample long, so taking the subsample closest to the middle of the aquisition period (acquisitionDelay)
  acquisitionSubSamples = floor(rem(acquisitionDelay,d.tr/estimationSupersampling)*downsampleFactor/d.tr)+1;
end
paddingSamples = 1-min(1,min(acquisitionSubSamples)); % if some subsamples are outside the design matrix (less than 1)
% if any(acquisitionSubSamples<1|acquisitionSubSamples>downsampleFactor)
% %   mrErrorDlg('(makeDesignMatrix) Acquisition period is outside frame period');
% %   keyboard
% end
allscm = [];
for iRun = 1:size(runTransition,1)
  scm = [];
  thisEVmatrix = d.EVmatrix(runTransition(iRun,1):runTransition(iRun,2),:);
  samplesToKeep = [];
  for iSubSample = acquisitionSubSamples
    samplesToKeep = [samplesToKeep;(iSubSample:downsampleFactor:size(thisEVmatrix,1)+iSubSample)+paddingSamples];
  end
  samplesToKeep = samplesToKeep(:,1:estimationSupersampling*numberSamples(iRun));
  % make stimulus convolution matrix
  for iEV = 1:size(thisEVmatrix,2)
      m = convn(thisEVmatrix(:,iEV), d.hrf);
      %apply saturation
      m = min(m,repmat(saturationThreshold,size(m,1),1));
      % downsample to estimation sampling rate 
      % (we keep extra samples at the end because downsample can remove samples)
      m = [zeros(paddingSamples,size(m,2));m];  %first pad with zeroes as necessary
      m = sum(m(samplesToKeep'),2);    %sum across acquisition sub-samples
      % Previously, downsampling was done by integrating, downsampling and differentiating, in order to conserve
      % the integral of each EV, which is equivalent to downsampling after a moving sum (low-pass filtering) on the full frame period
      % However this is incorrect in some situations (e.g. sparse imaging or slice-time-corrected timeseries):
      % In these case, the summation should not be done over the entire frame period but over the acquisition period
      % (which can now be specified). Old code:
      % m = mrDownsample(m, downsampleFactor, floor(rem(acquisitionDelay,d.tr/estimationSupersampling)*downsampleFactor/d.tr)+1);
      %only keep acquisition samples
      m = m(floor(acquisitionDelay/d.tr*estimationSupersampling)+1:estimationSupersampling:end,:);
      % remove mean 
      m = m-repmat(mean(m), size(m,1), 1); %DOES IT CHANGE ANYTHING IF I REMOVE THIS ?
      % apply the same filter as original data
      if isfield(d,'concatInfo') 
         % apply hipass filter
         if isfield(d.concatInfo,'hipassfilter') && ~isempty(d.concatInfo.hipassfilter{iRun}) ...
             && isfield(params.scanParams{scanNum},'highpassDesign') && params.scanParams{scanNum}.highpassDesign
           m = real(ifft(fft(m) .* repmat(d.concatInfo.hipassfilter{iRun}', 1, size(m,2)) ));
         end
         % project out the mean vector
         if isfield(d.concatInfo,'projection') && ~isempty(d.concatInfo.projection{iRun})...
              && isfield(params.scanParams{scanNum},'projectOutGlobalComponent') && params.scanParams{scanNum}.projectOutGlobalComponent
           projectionWeight = d.concatInfo.projection{iRun}.sourceMeanVector * m;
           m = m - d.concatInfo.projection{iRun}.sourceMeanVector'*projectionWeight;
         end
      end
      % stack stimmatrices horizontally
      scm =  [scm m];
   end
   % stack this run's stimcmatrix on to the last one
   allscm = [allscm;scm];
end

% set values
d.nhdr = size(stimToEVmatrix,2);
d.scm = allscm;
d.hdrlen = size(d.hrf,1);
d.nHrfComponents = size(d.hrf,2);
d.runTransitions = runTransition;
