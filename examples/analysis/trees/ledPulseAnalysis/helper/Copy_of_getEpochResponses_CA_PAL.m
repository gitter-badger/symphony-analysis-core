function outputStruct = getEpochResponses_CA_PAL(cellData, epochInd, varargin)%adopted getEpochResponses_CA_PAL

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('DeviceName', 'Amplifier_Ch1', @(x)ischar(x));
ip.addParameter('BaselineTime', 0, @(x)isnumeric(x)); %ms - overrides preTime for baseline interval calculation, use if preTime is missing
ip.addParameter('StartTime', 0, @(x)isnumeric(x)); %ms with 0 as stimulus start time
ip.addParameter('EndTime', 0, @(x)isnumeric(x)); %ms
ip.addParameter('LowPassFreq', 30, @(x)isnumeric(x)); %Hz
ip.addParameter('BinWidth', 10, @(x)isnumeric(x)); %ms
ip.addParameter('EndOffset', 0, @(x)isnumeric(x)); %ms
ip.parse(varargin{:});

epochLength = length(epochInd);
NO_OF_FEATURES = 15;

epochFeatures = cell(epochLength, NO_OF_FEATURES);
epochSummaryFeatures = cell(NO_OF_FEATURES);

epoch = cellData.epochs(epochInd(1));
intervalStart = ip.Results.StartTime * 1e-3;
% Default to end at stimEnd
if ip.Results.EndTime == 0
    intervalEnd = (epoch.get('stimTime') + ip.Results.EndOffset) * 1E-3; %s
else
    intervalEnd = ip.Results.EndTime * 1E-3; %s
end

[~, xvals] = epoch.getData(ip.Results.DeviceName);
[~, stimStart] = min(abs(xvals));

stimStart = stimStart(1);
baselineStart = xvals(1);
baselineEnd = ip.Results.BaselineTime * 1e-3;
responseIntervalLen = intervalEnd - intervalStart;
baselineIntervalLen = baselineEnd - baselineStart;
poststimIntervalLen = max(xvals) - intervalEnd;

[spikeAmps_all, ~, averageWaveform] = getSpikeAmplitudesForEpochs(cellData, epochInd, ip.Results.DeviceName);

% First pass through epochs is for getting baseline spikes, ISIs.
% These will be used later for figuring out the response probability

for i = 1 : epochLength
   
    epoch = cellData.epochs(epochInd(i));
    spikeTimes = getSpikeTimes(epoch, ip.Results.DeviceName);
    spikeAmps = spikeAmps_all{i} / max(spikeAmps_all{i});
    
    baselineSpikeTimes = spikeTimes(spikeTimes < baselineEnd);
    baselineIsi = diff(baselineSpikeTimes);
    poststimSpikeTimes = spikeTimes(spikeTimes > intervalEnd);
    baselineIsi2 = baselineIsi(1 : end - 1) + baselineIsi(2 : end);
    
    % Now we go through each response type in its own block
    % count spikes in stimulus interval
    count = length(find(spikeTimes >= intervalStart & spikeTimes < intervalEnd));
    epochFeatures{i, Features.SPIKE_COUNT_STIM_INTERVAL.id} = count;
    
    % Count spikes after stimulus
    count = length(find(spikeTimes >= intervalEnd));
    epochFeatures{i, Features.SPIKE_COUNT_POST_STIM.id} = count;
    epochFeatures{i, Features.SPIKE_COUNT_POST_STIM_BASE_LINE_SUB.id} = count - epochFeatures{i, Features.SPIKE_COUNT_BASE_LINE.Id};
    
    % Count spikes in 400 ms after onset and offset
    if responseIntervalLen >= 0.4
        spikeCount = length(find(spikeTimes >= intervalStart & spikeTimes < intervalStart + 0.4));
        epochFeatures{i, Features.SPIKE_COUNT_ONSET_400_MS.id} = spikeCount;
    end
    if intervalEnd + 0.4 <= xvals(end)
        spikeCount = length(find(spikeTimes >= intervalEnd & spikeTimes < intervalEnd + 0.4));
        epochFeatures{i, Features.SPIKE_COUNT_OFFSET_400_MS.id} = spikeCount;
    end
    
    epochFeatures{i, Features.ISI.id} = cumsum(diff(spikeTimes));
    epochFeatures{i, Features.SPIKE_AMP_DIFF.id} = cumsum(diff(spikeAmps));
    epochFeatures{i, Features.BASE_LINE_RATE.id} = length(baselineSpikeTimes) / baselineIntervalLen; 
    epochFeatures{i, Features.SPIKE_COUNT_BASE_LINE.id} = length(baselineSpikeTimes);                    
    epochFeatures{i, Features.POST_STIMULS_RATE.id} = length(poststimSpikeTimes) / poststimIntervalLen; 
    epochFeatures{i, Features.BASE_LINE_ISI.id} = baselineIsi;                                   
    epochFeatures{i, Features.BASE_LINE_ISI_TWO_SPIKES.id} = baselineIsi2;                                     
    epochFeatures{i, Features.SPIKE_AMP.id} = spikeAmps; 
    epochFeatures{i, Features.SPIKE_TIMES.id} = spikeTimes;        
end

epochSummaryFeatures{Features.STIM_ON_SET.id} = intervalStart;
epochSummaryFeatures{Features.STIM_OFF_SET.id} = intervalEnd;
epochSummaryFeatures{Features.RECORDING_ON_SET.id} = baselineStart;
epochSummaryFeatures{Features.RECORDING_OFF_SET.id} = max(xvals);
epochSummaryFeatures{Features.AVG_WAVE_FORM.id} = averageWaveform;
epochSummaryFeatures{Features.FULL_ISI.id} = flatten(epochFeatures{:, Features.ISI.id});
epochSummaryFeatures{Features.FULL_SPIKE_AMP_DIFF.id} = flatten(epochFeatures{:, Features.SPIKE_AMP_DIFF.id});


% Get response ISI threshold
[blistQshort, baselineISI, blistQ10, blistQ90] = getResponseISIThreshold(blISI, baselineIsi2);

epochSummaryFeatures{Features.BLIST_Q10} = blistQ10;
epochSummaryFeatures{Features.BLIST_Q90} = blistQ90;
epochSummaryFeatures{Features.FULL_BASE_LINE_ISI} = baselineISI;

% Get meanBaseline: to be used throughout
meanBaselineRate = mean(epochFeatures{:, Features.BASE_LINE_RATE.id});


% Loop over epochs to calculate responses
for i = 1 : epochLength
  
    spikeTimes = getSpikeTimes(epoch, ip.Results.DeviceName);
 
    % Subtract baseline
    spikeCount_baselineSubtracted = spikeCount - meanBaselineRate/responseIntervalLen;
    epochFeatures{i, Features.SPIKE_COUNT_STIM_INTERVAL_BASE_LINE_SUB.id} = spikeCount_baselineSubtracted;
    epochFeatures{i, Features.SPIKE_RATE_STIM_INTERVAL_BASE_LINE_SUB.id} = spikeCount_baselineSubtracted / responseIntervalLen;
    
    % Subtract baseline
    epochFeatures{i, Features.SPIKE_COUNT_ONSET_400_MS_BASE_LINE_SUB.id} = epochFeatures{i, Features.SPIKE_COUNT_ONSET_400_MS.id} - meanBaselineRate/0.4;
    epochFeatures{i, Features.SPIKE_COUNT_OFFSET_400_MS_BASE_LINE_SUB.id} = epochFeatures{i, Features.SPIKE_COUNT_OFFSET_400_MS.id} - meanBaselineRate/0.4;
 
    % find response start and end times based on ISIs: This is for a
    % stimulus that starts and ends and a particular time
    % save all values for calculating PSTH stuff at the end
    
    [startOnSet, endOnSet, startOffSet, endOffSet] = enhancedFiringResponse(spikeTimes, intervalStart, intervalEnd, blistQshort);
    
    %onset latency, spike count, duration, fullISI, and mean rate
    if ~ isempty(startOnSet)
        epochFeatures{i, Features.ONSET_LATENCY.id} = startOnSet - intervalStart;
        
        
        if ~ isempty(endOnSet)
            onSetSpikes = sum((spikeTimes >= startOnSet) & (spikeTimes <= endOnSet));
            epochFeatures{i, Features.ONSET_RESP_DURATION.id} = endOnSet - startOnSet;
            epochFeatures{i, Features.ONSET_RESP_RATE.id} =  onSetSpikes / (endOnSet - startOnSet);
            epochFeatures{i, Features.ONSET_RESP_RATE_BASE_LINE_SUB.id} = epochFeatures{i, Features.ONSET_RESP_RATE.id} - meanBaselineRate;
            epochFeatures{i, Features.ONSET_ISI.id} = diff(spikeTimes((spikeTimes >= startOnSet) & (spikeTimes <= endOnSet)));
        end
        
        epochFeatures{i, Features.ON_SET_SPIKES.id} = onSetSpikes;
        endOnsetBurstTime = burstInResponse(spikeTimes, startOnSet, endOnSet);
        
        if ~ isnan(endOnsetBurstTime)
            
            burstSpikesOnset =  sum((spikeTimes >= startOnSet) & (spikeTimes <= endOnsetBurstTime));
            nonBurstSpikesOnset =  sum((spikeTimes > endOnsetBurstTime) & (spikeTimes <= endOnSet));
            
            epochFeatures{i, Features.ONSET_BURST_SPIKES.id} = burstSpikesOnset;
            epochFeatures{i, Features.ONSET_NON_BURST_SPIKES.id} = nonBurstSpikesOnset;
            epochFeatures{i, Features.ONSET_BURST_DURATION.id} = endOnsetBurstTime - startOnSet;
            epochFeatures{i, Features.ONSET_NON_BURST_DURATION.id} = endOnSet - endOnsetBurstTime;
            epochFeatures{i, Features.ONSET_BURST_RATE.id} = burstSpikesOnset / epochFeatures{i, Features.ONSET_BURST_DURATION.id};
            epochFeatures{i, Features.ONSET_NON_BURST_RATE.id} = nonBurstSpikesOnset / epochFeatures{i, Features.ONSET_NON_BURST_DURATION.id};
            epochFeatures{i, Features.ONSET_BURST_NON_BURST_RATIO.id} = burstSpikesOnset ./ nonBurstSpikesOnset;
            epochFeatures{i, Features.ONSET_BURST_NON_BURST_RATIO_DURATION.id} =  epochFeatures{i, Features.ONSET_BURST_DURATION.id} /epochFeatures{i, Features.ONSET_NON_BURST_DURATION.id};
            epochFeatures{i, Features.ONSET_BURST_NON_BURST_RATIO_RATE.id} =  epochFeatures{i, Features.ONSET_BURST_RATE.id} /epochFeatures{i, Features.ONSET_NON_BURST_RATE.id};
            
        end
    end
    
    
    %offset latency, spike count, duration, fullISI, and mean rate
    if ~ isempty(startOffSet)
        
        epochFeatures{i, Features.OFFSET_LATENCY.id} = startOffSet - intervalEnd;
        
        if ~ isempty(OFFSETresponseEndTime)
            offSetSpikes = sum((spikeTimes >= startOffSet) & (spikeTimes <= endOffSet));
            
            epochFeatures{i, Features.OFFSET_RESP_DURATION.id} = endOffSet - startOffSet;
            epochFeatures{i, Features.OFFSET_RESP_RATE.id} = offSetSpikes / (endOffSet - startOffSet);
            epochFeatures{i, Features.OFFSET_RESP_RATE_BASE_LINE_SUB.id} = epochFeatures{i, Features.OFFSET_RESP_RATE.id} - meanBaselineRate;
            epochFeatures{i, Features.OFFSET_ISI.id} = offSetSpikes;
        end
        epochFeatures{i, Features.OFFSET_SPIKES.id} = offSetSpikes;
        endOffsetBurstTime = burstInResponse(spikeTimes, startOnSet, endOnSet);
        
        if ~ isnan(endOffsetBurstTime)
            
            burstSpikesOffset =  sum((spikeTimes >= startOffSet) & (spikeTimes <= endOffsetBurstTime));
            nonBurstSpikesOffset =  sum((spikeTimes > endOffsetBurstTime) & (spikeTimes <= endOffSet));
            
            epochFeatures{i, Features.OFFSET_BURST_SPIKES.id} = burstSpikesOffset;
            epochFeatures{i, Features.OFFSET_NON_BURST_SPIKES.id} = nonBurstSpikesOffset;
            epochFeatures{i, Features.OFFSET_BURST_DURATION.id} = endOffsetBurstTime - startOffSet;
            epochFeatures{i, Features.OFFSET_NON_BURST_DURATION.id} = endOffSet - endOffsetBurstTime;
            epochFeatures{i, Features.OFFSET_BURST_RATE.id} = burstSpikesOffset / epochFeatures{i, Features.OFFSET_BURST_DURATION.id};
            epochFeatures{i, Features.OFFSET_NON_BURST_RATE.id} = nonBurstSpikesOffset / epochFeatures{i, Features.OFFSET_NON_BURST_DURATION.id};
            epochFeatures{i, Features.OFFSET_BURST_NON_BURST_RATIO.id} = burstSpikesOffset ./ nonBurstSpikesOffset;
            epochFeatures{i, Features.OFFSET_BURST_NON_BURST_RATIO_DURATION.id} =  epochFeatures{i, Features.OFFSET_BURST_DURATION.id} /epochFeatures{i, Features.OFFSET_NON_BURST_DURATION.id};
            epochFeatures{i, Features.OFFSET_BURST_NON_BURST_RATIO_RATE.id} =  epochFeatures{i, Features.OFFSET_BURST_RATE.id} /epochFeatures{i, Features.OFFSET_NON_BURST_RATE.id};
        end
    end
    
    epochFeatures{i, Features.ON_OFF_INDEX.id} = onSetSpikes - offSetSpikes / (onSetSpikes + offSetSpikes);
  
    %peak instantaneous firing rate
    ONISI = outputStruct.ONSET_ISI_full.value;
    if  length(ONISI) >=2
        ONISI2 = ( ONISI(1:end-1) + ONISI(2:end) );
        outputStruct.ONSETpeakInstantaneousFR.value(i) = max(2./ONISI2);
    end
    OFFISI = outputStruct.OFFSET_ISI_full.value;
    if  length(OFFISI) >=2
        OFFISI2 = ( OFFISI(1:end-1) + OFFISI(2:end) );
        outputStruct.OFFSETpeakInstantaneousFR.value(i) = max(2./OFFISI2);
    end
    
end %end loop over epochs

%PSTH-based parameters
smoothingWindow = 0;%msec
PSTH_binwidth = 10;%msec
ONSETresponseStartTime_min = min(ONSETresponseStartTime_all);
ONSETresponseEndTime_max = max(ONSETresponseEndTime_all);
OFFSETresponseStartTime_min = min(OFFSETresponseStartTime_all);
OFFSETresponseEndTime_max = max(OFFSETresponseEndTime_all);
%[psth, xvals] = cellData.getPSTH(epochInd, ip.Results.BinWidth, ip.Results.DeviceName, smoothingWindow);
[psth, xvals] = cellData.getPSTH(epochInd, PSTH_binwidth, ip.Results.DeviceName, smoothingWindow);
outputStruct.PSTH.xvalue{1} = xvals;
outputStruct.PSTH.value{1} = psth;
%outputStruct.PSTH.smathvalue{1} = psth_smth;
outputStruct.PSTH.binWidth = PSTH_binwidth;
outputStruct.PSTH.smoothingWindow = smoothingWindow;

%ONSET
if ONSETresponseEndTime_max > ONSETresponseStartTime_min
    xvals_onset = xvals(xvals >= ONSETresponseStartTime_min & xvals < ONSETresponseEndTime_max);
    psth_onset = psth(xvals >= ONSETresponseStartTime_min & xvals < ONSETresponseEndTime_max);
    outputStruct.ONSETpsth.value = psth_onset;
    [outputStruct.ONSET_FRmax.value, maxLoc] = max(psth_onset);
    if ~isempty(maxLoc)
        maxLoc = maxLoc(1); 
        outputStruct.ONSET_FRmaxLatency.value = xvals_onset(maxLoc);
    end
    outputStruct.ONSET_FRrampLatency.value = outputStruct.ONSET_FRmaxLatency.value - nanmedian(outputStruct.ONSETlatency.value); %latency from start to peak
    outputStruct.ONSET_FRrange.value = outputStruct.ONSET_FRmax.value - min(psth_onset(maxLoc:end)); %range from max to end
    outputStruct.ONSET_FRrangeFrac.value = outputStruct.ONSET_FRrange.value / outputStruct.ONSET_FRmax.value;
end
%OFFSET
if OFFSETresponseEndTime_max > OFFSETresponseStartTime_min
    xvals_offset = xvals(xvals >= OFFSETresponseStartTime_min & xvals < OFFSETresponseEndTime_max);
    psth_offset = psth(xvals >= OFFSETresponseStartTime_min & xvals < OFFSETresponseEndTime_max);
    outputStruct.OFFSETpsth.value = psth_offset;
    [outputStruct.OFFSET_FRmax.value, maxLoc] = max(psth_offset);
    if ~isempty(maxLoc)
        maxLoc = maxLoc(1); 
        outputStruct.OFFSET_FRmaxLatency.value = xvals_offset(maxLoc) - OFFSETresponseStartTime;
    end
    outputStruct.OFFSET_FRrampLatency.value = outputStruct.OFFSET_FRmaxLatency.value - nanmedian(outputStruct.OFFSETlatency.value); %latency from start to peak
    outputStruct.OFFSET_FRrange.value = outputStruct.OFFSET_FRmax.value - min(psth_offset(maxLoc:end)); %range from max to end
    outputStruct.OFFSET_FRrangeFrac.value = outputStruct.OFFSET_FRrange.value / outputStruct.OFFSET_FRmax.value;
end

end


function spikeTimes = getSpikeTimes(epoch, deviceName)

% Get spike times (in units of seconds from startTime
spikeTimes = epoch.getSpikes(deviceName);
spikeTimes = spikeTimes - stimStart;
spikeTimes = spikeTimes / sampleRate;

end
