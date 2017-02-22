function [er, rrr, intEvidence] = LEI2_perfectMulti(x, pa)
% [er, rrr, intEvidence] = LEI2_perfect_check(x, pa)
% x   - parameter vector
% pa  - data structure (with fields: stis, runs, envSyllType, meanResp)

% weight for each syllable, w(2) is the block syllable with fixed weight of 1.0
w = [x(1) 1 x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9)];
s = x(end-2);      % noise level (std)
thres1 = x(end-1); % thres for pos response
thres2 = x(end);   % thres for neg response

% noise-less cumulative evidence
evidence = mapVal(pa.envSyllType, [unique(pa.envSyllType) w']'); % set weight for each syllable
intEvidence = cumsum(evidence);                 % integrate
intEvidence = intEvidence(1:33,:);              % cut after stimulus

cumNoise = s.*pa.cumNoise;% scale by noise std %

rrr = nan(1,pa.stis); % preallocate response for all stims
allSti = randperm(length(pa.meanResp), ceil(length(pa.meanResp)*pa.batch));% shuffle noise

for sti = 1:length(allSti) % process each stimulus
   
   % add (cumulated) noise to cumulated evidence
   evi = bsxfun(@plus, cumNoise(:,:,sti) , intEvidence(:,allSti(sti)));
   
   % set all superthres values to 1, 0 otherwise
   % append 1 at the end to mark non-crossing rows
   eviPos = evi> thres1;
   eviPos(end+1,:) = 1;
   eviNeg = evi<-thres2;
   eviNeg(end+1,:) = 1;
   
   % get indices of threshold crossings (i.e. index of first max value (1 since its a logical matrix)
   idxPos = argmax(eviPos);
   idxNeg = argmax(eviNeg);
   
   % set response to whatever threshold crossing came first
   resp = sign(idxNeg-idxPos)>0;
   
   % if thresh was simultaneous, there was no thres cross -
   % hence response is sign of evidence at the end
   resp(idxPos==idxNeg) = sign(evi(end,idxPos==idxNeg))>0;
   
   % get response probability by averaging over all noise instantiations
   rrr(allSti(sti)) = mean(resp);
end

% mse over the stimulus set
er = nanmean((rrr(allSti)-pa.meanResp(allSti)').^2);
