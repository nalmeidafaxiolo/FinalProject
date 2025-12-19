% ------------------------------------------------------------
%  Filter meteo data: keep daytime rows and convert humidity
% ------------------------------------------------------------

% Set input file path inside current folder (Preprocessing)
input_file = "mont_blanc_weather_2025.csv";

% Set output path to parent folder /data
output_file = "../mont_blanc_weather_2025_day.csv";

% Read the original CSV
data = readtable(input_file);

% Check required columns
requiredCols = ["is_day", "relative_humidity_2m"];
if ~all(ismember(requiredCols, data.Properties.VariableNames))
    error("Input file is missing required columns.");
end

% Filter daytime only
filtered = data(data.is_day == 1, :);

% Convert humidity from % to fraction
filtered.relative_humidity_2m = filtered.relative_humidity_2m / 100;

% Write filtered file into ../data
writetable(filtered, output_file);

% Display summary info
disp("Filtering complete.");
disp("Output saved to: " + output_file);
disp("Remaining rows: " + size(filtered,1));
