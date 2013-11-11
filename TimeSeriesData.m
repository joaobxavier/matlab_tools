classdef TimeSeriesData
    % defines objects that read a Tecan file and interprete the
    % time series data
    
    % properties
    properties (SetAccess = protected)
        numberOfWavelengths = 0;
        numberOfTimePoints = 0;
        wavelenghtData = [];
        samples = [];
        map = [];
    end
    
    % Method definisions
    methods
        % constructor: creates a new instance of time series
        % reads all data from the matrices
        function tsd = TimeSeriesData(filename)
            fid = fopen(filename);
            
            tline = fgetl(fid);
            while ischar(tline)
                % get the number of time points
                firstCell = TimeSeriesData.textInCell(tline, 1);
                if ~isempty(strfind(firstCell, 'Cycles'))
                    tsd.numberOfTimePoints =...
                        num2str(TimeSeriesData.textInCell(tline, 5));
                end;
                
                % get the data matrices
                if ~isempty(strfind(firstCell, 'Cycle Nr.'))
                    % increment the wavelength count
                    tsd.numberOfWavelengths = tsd.numberOfWavelengths + 1;
                    % get the wavelength name from previous line
                    tsd.wavelenghtData(tsd.numberOfWavelengths).name =...
                        TimeSeriesData.textInCell(previousLine, 1);
                    % extract well map
                    commas = strfind(tline, ','); % find indexes of commas
                    map = [];
                    column = 1;
                    for j = 4:length(commas)
                        map{end+1} = tline(commas(j-1)+1:commas(j)-1);
                        column = column + 1;
                    end;
                    % get the cells in last columns
                    map{column} = tline(commas(j)+1:end);
                    tsd.map = map;
                    % jump to next line
                    tline = fgetl(fid);
                    % start reading and compiling data until finding
                    % empty cell
                    n = 0;
                    while TimeSeriesData.numberInCell(tline, 1) == n+1
                        %check if this line has any data
                        if TimeSeriesData.checkIfLineHasData(tline),
                            % extract time, convert to hour
                            tsd.wavelenghtData(tsd.numberOfWavelengths)...
                                .times(n+1) =...
                                TimeSeriesData.numberInCell(tline, 2) / 3600;
                            % extract data
                            tsd.wavelenghtData(tsd.numberOfWavelengths)...
                                .data(n+1, :) =...
                                TimeSeriesData.numberInLine(tline, 4);
                            % update line number
                            n = TimeSeriesData.numberInCell(tline, 1);
                        end;
                        % proceed to next line in spreadsheet
                        tline = fgetl(fid);
                    end;
                end;
                % proceed to next line
                previousLine = tline;
                tline = fgetl(fid);
            end
            
            fclose(fid);
        end
        
        % Adds a second spreadsheet to this TimeSeriesData.
        % This can be used to join to consecutive time series.
        % Must be added before the samples are set.
        function tsd = appendMoreTimeSeriesData(tsd, filename, gap)
            if nargin == 2
                gap = 0;
            end
            newData = TimeSeriesData(filename);
            for i=1:length(tsd.wavelenghtData)
                endTime = tsd.wavelenghtData(i).times(end);
                % append the time array of part 2adding the final time of
                % part 1
                tsd.wavelenghtData(i).times = [tsd.wavelenghtData(i).times,...
                    (endTime + newData.wavelenghtData(i).times)+gap];
                % append the data of part 2
                tsd.wavelenghtData(i).data = [tsd.wavelenghtData(i).data;...
                    newData.wavelenghtData(i).data]
            end
        end
        
        % set the name of samples in a column
        % 'lines' is an optional argument to set which lines to include
        % if 'lines' is not defined, then assume lines = 1:8
        function tsd = setSampleName(tsd, name, column, lines)
            if nargin == 3
                lines = 1:8;
            end
            for i = lines,
                for j = 1:length(column)
                    well = [char(64 + i) num2str(column(j))];
                    tsd = setWell(tsd, name, well);
                end;
            end;
        end
        
        
        % assign a well to a given sample
        function tsd = setWell(tsd, sampleName, well)
            s = tsd.getSampleNumber(sampleName);
            % if sample doesn't exist yet, add it
            if (s == 0),
                tsd.samples(end+1).name  = sampleName;
                s = length(tsd.samples);
                tsd.samples(end).wells = [];
                tsd.samples(end).columns = [];
            end;
            tsd.samples(s).wells{end+1} = well;
            % convert wells to columns
            tsd.samples(s).columns(end+1) = tsd.findColumnNumber(well);
            % append data to the sample
            for w = 1:length(tsd.wavelenghtData)
                tsd.samples(s).wavelength(w).data    =...
                    tsd.wavelenghtData(w).data(:,...
                    tsd.samples(s).columns);
            end;
        end
        
        
        % get a plate data for a given cycle number
        function data = getPlateMatrix(tsd, wavelength, cycle, iOrders)
            data = tsd.wavelenghtData(wavelength).data(cycle, :);
            if nargin > 3
                data = data(iOrders);
            end
            data = reshape(data, 8, 12);
        end
        
        function data = getPlateMatrixSmall(tsd, wavelength, cycle, iOrders)
            data = tsd.wavelenghtData(wavelength).data(cycle, :);
            if nargin > 3
                data = data(iOrders);
            end
            data = reshape(data, 6, 10);
        end
        
        % get the matrix od data for a sample number
        function data = getData(tsd, wavelength, sampleNumber)
            data = tsd.samples(sampleNumber).wavelength(wavelength).data;
        end
        
        % get the array of times
        function times = getTimes(tsd, wavelength)
            times = tsd.wavelenghtData(wavelength).times;
        end
        
        % get the data for a sample name
        function data = getDataFromName(tsd, wavelength, sampleName)
            % find the number of sample
            for s = 1:length(tsd.samples),
                if strcmp(tsd.samples(s).name, sampleName)
                    break
                end
            end
            % get the data
            data = tsd.getData(wavelength, s);
        end
        
        
        % do blank correction using median over entire blak replicates
        function tsd = performBlankCorrection(tsd, blankName)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelenghtData)
                % get the blankData
                blankData = tsd.getDataFromName(w, blankName);
                % calculate median along 2nd dimension (along replicates)
                blankData = median(blankData, 2);
                % subtract that from all data
                for s = 1:length(tsd.samples),
                    % get the data
                    data =  tsd.samples(s).wavelength(w).data;
                    % subtract blank to each column
                    for i = 1:size(data, 2)
                        data(:, i) = data(:, i) - blankData;
                    end
                    % rewrite the variable
                    tsd.samples(s).wavelength(w).data = data;
                end
            end
            tsd = tsd.updateSampleData;
        end
        
        % do blank correction but take as correction the same line in
        % the blank sample rather than median over entire blak replicates
        % This function has two ways to be used used.
        % If the argument 'samples' is not supplied, then the function
        % preforms blank correction for all the samples.
        % If 'samples' is supplied as a vector of sample numbers, then the
        % correction is carried out only for the samples in the array.
        function tsd = performBlankCorrectionPerLine(tsd, blankName,...
                samples)
            if (nargin == 2)
                samples = 1:length(tsd.samples);
            end;
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelenghtData)
                % get the blankData
                blankData = tsd.getDataFromName(w, blankName);
                % get median accross replicates
                blankData = median(blankData, 2);
                % subtract that from all data
                for s = samples,
                    % get the data
                    data =  tsd.samples(s).wavelength(w).data;
                    for i = 1:size(data, 2)
                        data(:, i) =  data(:, i) - blankData;
                    end;
                    % rewrite the variable
                    tsd.samples(s).wavelength(w).data = data;
                end
            end
            tsd = tsd.updateSampleData;
        end
        
           % correction is carried out only for the samples in the array.
        function tsd = subtractBackgroundValue(tsd, value,...
                samples, w1, w2)
            % subtract that from all data
            for s = samples,
                % get the data
                data =  tsd.samples(s).wavelength(w1).data;
                for i = 1:size(data, 2)
                    data(:, i) =  data(:, i) - value *...
                        tsd.samples(s).wavelength(w2).data(:, i);
                end;
                % rewrite the variable
                tsd.samples(s).wavelength(w1).data = data;
            end
            tsd = tsd.updateSampleData;
        end
        
        % do blank correction taking as background the first 10 point
        function tsd = performBlankCorrection10Points(tsd)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelenghtData)
                % get the data
                data =  tsd.wavelenghtData(w).data;
                % define the blankData
                blankData = median(data(1:10, :), 1);
                blankData = repmat(blankData, size(data, 1), 1);
                % subtract
                data =  data - blankData;
                % rewrite the variable
                tsd.wavelenghtData(w).data = data;
            end
            tsd = tsd.updateSampleData;
        end
        
        % use after blanking to make sure sample data is up to date
        function tsd = updateSampleData(tsd)
            % subtract that from all data
            for s = 1:length(tsd.samples),
                for w = 1:length(tsd.wavelenghtData)
                    tsd.wavelenghtData(w).data(:,...
                        tsd.samples(s).columns)=...
                        tsd.samples(s).wavelength(w).data;
                end
            end
        end
        
        
        % move the curves in the x-axis by a time step tau
        function tsd = moveInTime(tsd, samples, tau)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelenghtData)
                time = tsd.getTimes(w);
                indexValid = find( and(time>tau, time<(time(end)-tau)) );
                % correct all samples in array
                for s = samples,
                    % get the data
                    data =  tsd.samples(s).wavelength(w).data;
                    for i = 1:size(data, 2)
                        newTime  = time(indexValid)-tau;
                        dataMoved = interp1(newTime, data(indexValid, i), time);
                        data(:, i) =  dataMoved;
                    end;
                    % rewrite the variable
                    tsd.samples(s).wavelength(w).data = data;
                end
            end
        end;
        
        % plot the time series of the median for a given sample
        function h = plotMedian(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            h(1) = plot(times, median(data, 2), 'k-');
            set(h(1), 'Linewidth', 3);
        end
        
        % plot the time series of ranges for a given sample
        function h = plotRangesAsLines(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            h(1) = plot(times, max(data, [], 2), 'k-');
            hold on;
            h(2) = plot(times, min(data, [], 2), 'k-');
            hold off;
        end
        
        % plot the time series of ranges for a given sample
        function h = plotRangesAsArea(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            maxData = max(data, [], 2);
            minData = min(data, [], 2);
            range   = maxData - minData;
            % take out NaN's, otherwise area won't work properly
            nans = isnan(minData) | isnan(range) | (minData <= 0);
            minData = minData(~nans);
            range = range(~nans);
            times = times(~nans);
            % plot the area
            h = area(times, [minData, range]);
            set(h, 'EdgeColor', 'none');
            % erase the first area
            set(h(1), 'FaceColor', 'none');
            % return only the handle to the second area
            h = h(2);
        end
        
        
        % plot the time series of the median for a given sample
        % with ranges overlapped
        function h = plotMedianWithRanges(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            h(1) = tsd.plotMedian(wavelength, sampleNumber);
            hold on;
            h2 = tsd.plotRangesAsLines(wavelength, sampleNumber);
            hold off;
            h = [h h2(1) h2(2)];
            set(h(1), 'Linewidth', 3);
            set(gca, 'YScale', 'log');
            title(tsd.samples(sampleNumber).name);
            set(gca, 'XLim', [0 24], 'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelenghtData(wavelength).name);
        end
        
        
        % plots OD and GFP in the same plot, normalized by the entire
        % sample
        function plotNormalizedODAndGFP(tsd, sampleNumber)
            %
            odData = (tsd.getData(1, sampleNumber));
            maxVal = max(odData(~isnan(odData)));
            minVal = min(odData(~isnan(odData)));
            odData = (odData - minVal) ./ (maxVal - minVal);
            % gfp (second wavelength)
            gfpData = tsd.getData(2, sampleNumber);
            maxVal = max(gfpData(~isnan(gfpData)));
            minVal = min(gfpData(~isnan(gfpData)));
            gfpData = (gfpData - minVal) ./ (maxVal - minVal);
            % time
            times = tsd.getTimes(1);
            maxTime = max(times(~isnan(times)));
            minTime = min(times(~isnan(times)));
            plot(times, median(odData, 2), 'k', 'LineWidth', 2);
            hold on;
            plot(tsd.getTimes(2), median(gfpData, 2),...
                'Color', [0 0.8 0], 'LineWidth', 2);
            plot(times, max(odData, [], 2), 'k-');
            plot(times, min(odData, [], 2), 'k-');
            plot(times, max(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(times, min(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            hold off;
            set(gca, 'YLim', [0 1],...
                'XLim', [minTime maxTime], 'YTickLabel', []);
        end
        
        
        % plots three wavelengths in the same plot, normalized by the entire
        % sample. the example here is for od, gfp and pyoverdine
        function plotThreeWavelengths(tsd, sampleNumber)
            % od (first wavelength)
            odData = (tsd.getData(1, sampleNumber));
            maxVal = max(odData(~isnan(odData)));
            minVal = min(odData(~isnan(odData)));
            odData = (odData - minVal) ./ (maxVal - minVal);
            % gfp (second wavelength)
            gfpData = tsd.getData(2, sampleNumber);
            maxVal = max(gfpData(~isnan(gfpData)));
            minVal = min(gfpData(~isnan(gfpData)));
            gfpData = (gfpData - minVal) ./ (maxVal - minVal);
            % pyo (third wavelength)
            pyoData = tsd.getData(3, sampleNumber);
            maxVal = max(pyoData(~isnan(pyoData)));
            minVal = min(pyoData(~isnan(pyoData)));
            pyoData = (pyoData - minVal) ./ (maxVal - minVal);
            % time
            times = tsd.getTimes(1);
            maxTime = max(times(~isnan(times)));
            minTime = min(times(~isnan(times)));
            plot(times, median(odData, 2), 'k', 'LineWidth', 2);
            hold on;
            plot(tsd.getTimes(2), median(gfpData, 2),...
                'Color', [0 0.8 0], 'LineWidth', 2);
            plot(tsd.getTimes(3), median(pyoData, 2),...
                'Color', [0 0 1], 'LineWidth', 2);
            plot(times, max(odData, [], 2), 'k-');
            plot(times, min(odData, [], 2), 'k-');
            plot(times, max(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(times, min(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(times, max(pyoData, [], 2), '-', 'Color', [0 0 1]);
            plot(times, min(pyoData, [], 2), '-', 'Color', [0 0 1]);
            hold off;
            set(gca, 'YLim', [0 1],...
                'XLim', [minTime maxTime], 'YTickLabel', []);
        end
        
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        function plotManySamples(tsd, wavelength, sampleNumbers, rangeFlag)
            % construct a colormap
            labels = [];
            cmap = jet(length(sampleNumbers));
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotMedian(wavelength, sampleNumbers(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(sampleNumbers(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
            if ((nargin == 4) && (rangeFlag ~= false)) || (nargin == 3)
                for i = 1:length(sampleNumbers)
                    hold on;
                    %h = tsd.plotRangesAsArea(wavelength, sampleNumbers(i));
                    %set(h, 'FaceColor', cmap(i, :));
                    h = tsd.plotRangesAsLines(wavelength, sampleNumbers(i));
                    set(h, 'Color', cmap(i, :));
                end;
            end
            set(gca, 'YScale', 'log', 'XLim', [0 24], 'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelenghtData(wavelength).name);
        end;
        
        % plot wavelength w2 as function of wavelength w2
        % for sample s
        function h = plotWavelength1VsWavelength2(tsd, w1, w2, s)
            d1 = tsd.getData(w1, s);
            d2 = tsd.getData(w2, s);
            h = plot(d1(:), d2(:), '.');
            xlabel(tsd.wavelenghtData(w1).name);
            ylabel(tsd.wavelenghtData(w2).name);
        end;
        
        % same as plotWavelength1VsWavelength2 but works with array
        % of samples
        function plotW1VsW2ManySamples(tsd, w1, w2, samples)
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            for i = 1:length(samples)
                hold on;
                h = tsd.plotWavelength1VsWavelength2(w1, w2, samples(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            legend(labels, 'Location', 'Best');
        end
        
        
        % plts the specific growth rate time series
        function h = plotSpecificGrowthW1(tsd, w1, s)
            d1 = tsd.getData(w1, s);
            time = tsd.getTimes(w1);
            d1 = TimeSeriesData.calculateSpecificRate(time, d1);
            h = plot(time, median(d1, 2), '-');
            ylabel(['growth' tsd.wavelenghtData(w1).name]);
            xlabel('time [h]');
        end;
        
        %
        function h = plotSpecificGrowthW1VsW2(tsd, w1, w2, s, thresholdW1, useMedian)
            % default useMedian = 0
            if nargin == 5
                useMedian = 0;
            end
            d1 = tsd.getData(w1, s);
            % take only the median of w1
            if useMedian == 1
                d1 = median(d1, 2);
            end;
            %
            time = tsd.getTimes(w1);
            d1 = TimeSeriesData.calculateSpecificRate(time, d1);
            d2 = TimeSeriesData.filterData(tsd.getData(w2, s));
            % take only the median of w2
            if useMedian == 1
                d2 = median(d2, 2);
            end
            %
            % threshold for wvelength w1
            d1Filterd = TimeSeriesData.filterData(tsd.getData(w1, s));
            % take only the median of w1
            if useMedian == 1
                d1Filterd = median(d1Filterd, 2);
            end
            %
            d1(d1Filterd < thresholdW1) = [];
            d2(d1Filterd < thresholdW1) = [];
            d1Filterd(d1Filterd < thresholdW1) = [];
            h = plot(d1(:), d2(:)./d1Filterd(:), '.');
            ylabel(tsd.wavelenghtData(w2).name);
            % use second line to normalize gfp/od600
            %h = plot(d1(:), d2(:)./d1Filterd(:), '+');
            %ylabel([tsd.wavelenghtData(w2).name '/'...
            %    tsd.wavelenghtData(w1).name]);
            xlabel(['specific growth ' tsd.wavelenghtData(w1).name]);
        end;
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        function plotSpecificGrowthW1VsW2Many(tsd, w1, w2, samples, thresholdW1, useMedian)
            % default useMedian = 0
            if nargin == 5
                useMedian = 0;
            end
            
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            for i = 1:length(samples)
                hold on;
                h = tsd.plotSpecificGrowthW1VsW2(w1, w2, samples(i), thresholdW1, useMedian);
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
        end;
        
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        function plotManySamplesWithArea(tsd, wavelength, sampleNumbers)
            % construct a colormap
            labels = [];
            cmap = jet(length(sampleNumbers));
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotMedian(wavelength, sampleNumbers(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(sampleNumbers(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotRangesAsArea(wavelength, sampleNumbers(i));
                set(h, 'FaceColor', cmap(i, :));
            end
            set(gca, 'YScale', 'log', 'XLim', [0 24], 'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelenghtData(wavelength).name);
        end;
        
        % plot a matrix representing data from a plate
        function plotMatrix(tsd, wavelength, cycleNumber, column, line,...
                iOrders)
            if nargin > 5
                matrix2Plot =...
                    tsd.getPlateMatrix(wavelength, cycleNumber, iOrders);
            else
                matrix2Plot = tsd.getPlateMatrix(wavelength, cycleNumber);
            end
            imagesc(matrix2Plot);
            times = tsd.getTimes(wavelength);
            title(sprintf('%s at %0.1f h',...
                tsd.wavelenghtData(wavelength).name,...
                times(cycleNumber)));
            colorbar;
            set(gca, 'XTick', 1:size(matrix2Plot, 2),...
                'XTickLabel', line);
            set(gca, 'YTick', 1:size(matrix2Plot, 1),...
                'YTickLabel', column);
        end;
        
        function plotMatrixSmall(tsd, wavelength, cycleNumber, column,...
                line, iOrders)
            if nargin > 5
                matrix2Plot =...
                    tsd.getPlateMatrixSmall(wavelength, cycleNumber, iOrders);
            else
                matrix2Plot = tsd.getPlateMatrix(wavelength, cycleNumber);
            end
            imagesc(matrix2Plot);
            times = tsd.getTimes(wavelength);
            title(sprintf('%s at %0.1f h',...
                tsd.wavelenghtData(wavelength).name,...
                times(cycleNumber)));
            colorbar;
            set(gca, 'XTick', 1:size(matrix2Plot, 2),...
                'XTickLabel', line);
            set(gca, 'YTick', 1:size(matrix2Plot, 1),...
                'YTickLabel', column);
        end;
        
        % plot a surface plot representing data from a plate
        function plotSurface(tsd, wavelength, cycleNumber, column, line, iOrders)
            if nargin > 5
                matrix2Plot =...
                    tsd.getPlateMatrix(wavelength, cycleNumber, iOrders);
            else
                matrix2Plot = tsd.getPlateMatrix(wavelength, cycleNumber);
            end
            surf(matrix2Plot);
            times = tsd.getTimes(wavelength);
            title(sprintf('%s at %0.1f h',...
                tsd.wavelenghtData(wavelength).name,...
                times(cycleNumber)));
            colorbar;
            set(gca, 'XTick', 1:size(matrix2Plot, 2),...
                'XTickLabel', line);
            set(gca, 'YTick', 1:size(matrix2Plot, 1),...
                'YTickLabel', column);
        end;
        
        function plotSurfaceSmall(tsd, wavelength, cycleNumber, column,...
                line, iOrders)
            if nargin > 5
                matrix2Plot =...
                    tsd.getPlateMatrixSmall(wavelength, cycleNumber, iOrders);
            else
                matrix2Plot = tsd.getPlateMatrix(wavelength, cycleNumber);
            end
            surf(matrix2Plot);
            times = tsd.getTimes(wavelength);
            title(sprintf('%s at %0.1f h',...
                tsd.wavelenghtData(wavelength).name,...
                times(cycleNumber)));
            colorbar;
            set(gca, 'XTick', 1:size(matrix2Plot, 2),...
                'XTickLabel', line);
            set(gca, 'YTick', 1:size(matrix2Plot, 1),...
                'YTickLabel', column);
        end;
        
        % plot OD and GFP in the same plot foer every well in a plate
        function plotTimeSeriesMatrix(tsd)
            % od (first wavelength)
            odData = tsd.wavelenghtData(1).data;
            odData(odData <= 0) = NaN;
            odData = log10(odData);
            maxVal = max(odData(~isnan(odData)));
            %minVal = min(odData(~isnan(odData)));
            minVal = -2;
            odData = (odData - minVal) ./ (maxVal - minVal);
            % gfp (second wavelength)
            gfpData = tsd.wavelenghtData(2).data;
            gfpData(gfpData <= 0) = NaN;
            %gfpData = log10(gfpData);
            maxVal = max(gfpData(~isnan(gfpData)));
            %minVal = min(gfpData(~isnan(gfpData)));
            minVal = 1;
            gfpData = (gfpData - minVal) ./ (maxVal - minVal);
            % time
            times = tsd.getTimes(1);
            maxTime = max(times(~isnan(times)));
            minTime = min(times(~isnan(times)));
            k = 0;
            for j = 1:12,
                for i = 1:8,
                    k = k + 1;
                    subplot(8, 12, (i-1)*12 + j);
                    plot(times,...
                        odData(:, k),...
                        'k',...
                        'LineWidth', 1.5);
                    hold on;
                    plot(tsd.getTimes(2),...
                        gfpData(:, k),...
                        'g',...
                        'Color', [0 0.8 0],...
                        'LineWidth', 1.5);
                    hold off;
                    set(gca, 'YLim', [0 1],...
                        'XLim', [minTime maxTime],...
                        'YTickLabel', [],...
                        'Color', 'none',...
                        'LineWidth', 1.5);
                end;
            end;
        end;
        
        
        
        % search for the number of a sample with a given name
        function n = getSampleNumber(tsd, name)
            n = 0;
            for i = 1:length(tsd.samples)
                if strcmp(name, tsd.samples(i).name)
                    n = i;
                    return;
                end;
            end;
        end;
        
        % search for the column number that corresponds to a well
        function n = findColumnNumber(tsd, well)
            n = 0;
            for i = 1:length(tsd.map)
                if strcmp(well, tsd.map{i})
                    n = i;
                    return;
                end;
            end;
        end;
        
        % find the array of time delays between sample 'referenceSample'
        % and the samples in array sampleNumbers
        function tauArray = computeTimeDelays(tsd,...
                wavelength, referenceSample, sampleNumbers)
            time = tsd.getTimes(wavelength);
            tauArray = zeros(1, length(sampleNumbers));
            curve1 = tsd.getData(wavelength, referenceSample);
            for i = 1:length(sampleNumbers)
                s = sampleNumbers(i);
                curve2 = tsd.getData(wavelength, s);
                fToMin =...
                    @(tau)...
                    (TimeSeriesData.alignmentError(time,...
                    curve1, curve2, tau));
                tauArray(i) = fminbnd(fToMin, 0, 24,...
                    optimset('TolX', 1e-2));
            end
        end
        
        % find the array of time delays between sample 'referenceSample'
        % and the samples in array sampleNumbers
        function tauMatrix = computeTimeDelayMatrix(tsd,...
                wavelength, sampleNumbers)
            n = length(sampleNumbers);
            tauMatrix = zeros(n, n);
            for i = 1:n
                for j = (i+1):n
                    tauMatrix(i, j) = tsd.computeTimeDelays(...
                        wavelength, sampleNumbers(i), sampleNumbers(j),...
                        1e-2);
                end
            end
        end
        
        % calulate a total error of overlap for a given array of delays
        function errorVal = computeErrorMatrix(tsd,...
                wavelength, sampleNumbers, tauArray)
            tauArray = [0 tauArray];
            n = length(sampleNumbers);
            errorMatrix = zeros(n, n);
            time = tsd.getTimes(wavelength);
            for i = 1:n
                for j = (i+1):n
                    curve1 = tsd.getData(wavelength, sampleNumbers(i));
                    curve2 = tsd.getData(wavelength, sampleNumbers(j));
                    tau = tauArray(j) - tauArray(i);
                    errorMatrix(i, j) =...
                        TimeSeriesData.alignmentError(time, curve1, curve2, tau);
                end
            end
            errorVal = sum(errorMatrix(:));
        end
        
        
        % calculate an array of delays by minimizing the
        % total error of overlap between each pair of curves
        function tauArray = optimizeTauArray(tsd,...
                wavelength, sampleNumbers)
            tauArray = tsd.computeTimeDelays(wavelength,...
                sampleNumbers(1), sampleNumbers(2:end));
            fToMin =...
                @(tauArray)...
                (computeErrorMatrix(tsd,...
                wavelength, sampleNumbers, tauArray));
            tauArray = fminsearch(fToMin, tauArray);
            
        end
        
        
    end % methods
    
    
    methods(Static)
        % extract the text in cell 'i' in a line of the spreadsheet
        function t = textInCell(line, i)
            commas = strfind(line, ','); % find indexes of commas
            if i == 1
                if isempty(commas)
                    t = line;
                else
                    t = line(1:commas(i)-1);
                end
            else
                t = line(commas(i-1)+1:commas(i)-1);
            end
        end % end textInCell
        
        
        % extract the number in cell 'i' in a line of the spreadsheet
        function a = checkIfLineHasData(line)
            commas = strfind(line, ','); % find indexes of commas
            a = ~isempty(commas);
            % check if there is time data at entry 2
            try
                if isnan(TimeSeriesData.numberInCell(line, 2))
                    a = false;
                end
            catch
                a = false;
            end
        end % end textInCell
        
        
        % extract the number in cell 'i' in a line of the spreadsheet
        function n = numberInCell(line, i)
            n = str2double(TimeSeriesData.textInCell(line, i));
        end % end textInCell
        
        % extract all numbers from cell 'i' in a line of the spreadsheet
        function n = numberInLine(line, i)
            commas = strfind(line, ','); % find indexes of commas
            n = zeros(1, length(commas)-i);
            column = 1;
            for j = i:length(commas)
                str = line(commas(j-1)+1:commas(j)-1);
                if strcmp(str, 'OVER')
                    n(column) = Inf;
                else
                    n(column) = str2double(str);
                end;
                column = column + 1;
            end;
            % get the cells in last columns
            n(column) = str2double(line(commas(j)+1:end));
            % return only array of non NaN
            n = n(~isnan(n));
        end % end textInCell
        
        % filter data using to take out noise
        function dFiltered = filterData(d)
            windowSize = 5;
            dFiltered = filter(ones(1,windowSize)/windowSize,1,d, [], 1);
        end
        
        % calculate specific growth rate of a wavelength
        function r = calculateSpecificRate(time, data)
            % get the number of replicates
            nReps = size(data, 2);
            % log-transform data
            data = log(TimeSeriesData.filterData(data));
            % relicate time array
            timeDifferential = diff(time);
            timeDifferential = repmat(timeDifferential', 1, nReps);
            % filter data to reduce noice
            r = data;
            % calculate derivative
            r = diff(r, 1, 1)./timeDifferential;
            % replicate the first data point to make data same size
            % as original data
            r = [r(1, :); r];
        end;
        
        
        function err = alignmentError(time, curve1, curve2, tau)
            mCurve1 = median(curve1, 2);
            mCurve2 = median(curve2, 2);
            
            indexRef   = find(time<(time(end)-tau));
            indexOther = find(time>tau);
            ref        = mCurve1(indexRef);
            timeRef    = time(indexRef);
            other      = mCurve2(indexOther);
            timeOther  = time(indexOther)-tau;
            
            % take out first time point from timeRef to avoid NaN
            timeRef = timeRef(2:end);
            
            otherInterp = interp1(timeOther, other, timeRef);
            
            ref = ref(2:end);
            otherInterp = otherInterp(:);
            
            
            v = (ref - otherInterp).^2;
            err = sum(v);
        end;
        
        
        function str = num2strRound(n, i)
            
            % str = num2strRound(n, i) - works like num2str but rounds up if next
            % digit is >= 5.
            
            n = n.*10.^i;
            n = round(n);
            n = n./10.^i;
            str = num2str(n);
        end;
        
    end % end static methods
    
    
end % classdef