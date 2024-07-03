function [Clean_TimeSeries,AR1_Stats] = AR1_Clean(TimeSeries)
%Written by David Rothlein
% Correct for n-1 autoregression in timeseries
TimeSeries=TimeSeries(:);
% n-1 timeseries deplicates first timepoint and removes last.
TimeSeriesM1=[TimeSeries(1,1);TimeSeries(1:end-1)];
ConstInt=ones(length(TimeSeries),1);
[~,~,Clean_TimeSeries,~,AR1_Stats]=regress(TimeSeries,[TimeSeriesM1,ConstInt]);
end

