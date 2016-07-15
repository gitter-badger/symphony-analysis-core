classdef Features < handle
   
   enumeration
       ISI(1, 'sec')
       SPIKE_AMP_DIFF(2, '')
       BASE_LINE_RATE(3, 'hz', 'byEpoch')
       SPIKE_COUNT_BASE_LINE(4, 'spikes', 'byEpoch')
       POST_STIMULS_RATE(5, 'hz', 'byEpoch')
       BASE_LINE_ISI(6, '')
       BASE_LINE_ISI_TWO_SPIKES(7, '')
       SPIKE_AMP(8, '')
       SPIKE_TIMES(9, '')
       SPIKE_COUNT_STIM_INTERVAL(10, '', '')
       SPIKE_COUNT_POST_STIM(11, '', '')
       SPIKE_COUNT_POST_STIM_BASE_LINE_SUB(12, '', '')
       SPIKE_COUNT_ONSET_400_MS(13, '', '')
       SPIKE_COUNT_OFFSET_400_MS(14, '', '')
       SPIKE_COUNT_STIM_INTERVAL_BASE_LINE_SUB(15, '', '')
       SPIKE_RATE_STIM_INTERVAL_BASE_LINE_SUB(16, '', '')
       SPIKE_COUNT_ONSET_400_MS_BASE_LINE_SUB(17, '', '')
       SPIKE_COUNT_OFFSET_400_MS_BASE_LINE_SUB(18, '', '')
       ONSET_LATENCY(19, '', '')
       ONSET_RESP_DURATION(20, '', '')
       ONSET_RESP_RATE(21, '', '')
       ONSET_RESP_RATE_BASE_LINE_SUB(22, '', '')
       ON_SET_SPIKES(23, '', '')
       ONSET_ISI(24, '', '')
       
       OFFET_LATENCY(25, '', '')
       OFFSET_RESP_DURATION(26, '', '')
       OFFSET_RESP_RATE(27, '', '')
       OFFSET_RESP_RATE_BASE_LINE_SUB(28, '', '')
       OFFSET_SPIKES(29, '', '')
       OFFSET_ISI(30, '', '')
       
       ON_OFF_INDEX(31, '', '')
       ONSET_BURST_SPIKES(32, '', '')
       ONSET_NON_BURST_SPIKES(33, '', '')
       ONSET_BURST_DURATION(34, '', '')
       ONSET_NON_BURST_DURATION(35, '', '')
       ONSET_BURST_RATE(36, '', '')
       ONSET_NON_BURST_RATE(37, '', '')
       ONSET_BURST_NON_BURST_RATIO(38, '', '')
       ONSET_BURST_NON_BURST_RATIO_DURATION(39, '', '')
       ONSET_BURST_NON_BURST_RATIO_RATE(40, '', '')
   end
   
   enumeration
       STIM_ON_SET(1, 's', 'combinedAcrossEpochs')
       STIM_OFF_SET(2, 's', 'combinedAcrossEpochs')
       STIMULS_START(3, 's', 'combinedAcrossEpochs')
       RECORDING_ON_SET(4, 's', 'combinedAcrossEpochs')
       RECORDING_OFF_SET(5, 's', 'combinedAcrossEpochs')
       AVG_WAVE_FORM(6, '', 'combinedAcrossEpochs')
       FULL_ISI(7, 's', 'combinedAcrossEpochs')
       FULL_SPIKE_AMP_DIFF(8, 'combinedAcrossEpochs')
       FULL_BASE_LINE_ISI(9, 's', 'combinedAcrossEpochs')
       BLIST_Q10(10, 's', 'singleValue')
       BLIST_Q90(11, 's', 'singleValue')
   end
   
   properties
       id
       units
       type
   end
   
   methods
       
       function obj = Features(id, units)
           obj.id = id;
           obj.units = units;
       end
       
   end
    
end

