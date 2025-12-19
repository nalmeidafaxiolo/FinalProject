% ===============================================================
% Compile and run the C program producing "results_melt.csv"
% ===============================================================
% Compile the c script 
if ~isfile('snow_energy_balance')
    fprintf('>> Compiling the C program...\n');
    system('gcc snow_energy_balance.c -o snow_energy_balance -lm');
end

% Executes the c script
fprintf('>> Running the C program snow_energy_balance...\n');
system('./snow_energy_balance mont_blanc_weather_2025_day.csv');

melt_animation_slider_staticplots();




% ===============================================================
% MATLAB melt visualisation function
% ===============================================================

function melt_animation_slider_staticplots()

    % Read the CSV file obtained from the C file
    resultsFolder = '../results/';  % relative path to the results folder
    T = readtable(fullfile(resultsFolder, 'results_melt.csv'));


    % Extract variables
    lon  = T.lon; 
    lat  = T.lat; 
    alt  = T.altitude;
    melt_pure    = T.melt_pure_mmh;      
    melt_polluted = T.melt_polluted_mmh; 

    % convert the time values to time datatype
    T.time = datetime(T.time, 'InputFormat','yyyy-MM-dd''T''HH:mm');
    time = T.time;
    times  = unique(time);
    nTimes = numel(times);
    

  
    % Creates a fixed spatial grid to interpolate the data and avoid ...
    % gaps caused by irregular point locations
    nGrid = 200;
    LonGrid = linspace(min(lon), max(lon), nGrid);
    LatGrid = linspace(min(lat), max(lat), nGrid);
    [LonGrid, LatGrid] = meshgrid(LonGrid, LatGrid);

    % Compute time-averaged melt rates (static curves)
    melt_pure_mean    = zeros(nTimes,1);
    melt_polluted_mean = zeros(nTimes,1);
    melt_diff_mean     = zeros(nTimes,1);

    for k = 1:nTimes
        idx_k = (time == times(k));
        melt_pure_mean(k)    = mean(melt_pure(idx_k));
        melt_polluted_mean(k) = mean(melt_polluted(idx_k));
        melt_diff_mean(k)     = melt_polluted_mean(k) - melt_pure_mean(k);
    end


    % Create the main figure window used for all plots
    fig = figure('Name','Melt Visualization','Color',[1 1 1],...
                 'Position',[50 50 1500 950]);

    % Scaling factor for the Z-axis to make topography more visually...
    % pronounced
    zscale = 6;
    
   
    %% ========================= Left panel (pure snow) =================
    ax1 = axes(fig);

    % Position of the subplot (kept identical across panels for...
    % consistent layout)
    ax1.Position = [0.05 0.42 0.42 0.53];  
    
    % Create empty surface object (colors will be updated dynamically)
    hSurf1 = surf(ax1, LonGrid, LatGrid, zeros(size(LonGrid)),...
        zeros(size(LonGrid)));
    shading(ax1,'interp');
    colormap(ax1,'turbo');

    % Colorbar placed outside the axes to preserve identical panel dimensions
    c1 = colorbar(ax1);
    c1.Location = 'eastoutside';  

    % Color scaling for melt rate values
    clim(ax1,[0 6]);

    xlabel(ax1,'Longitude'); ylabel(ax1,'Latitude'); zlabel(ax1,'Altitude');
    title(ax1,'Pure snow melt (mm/h)');

    % View angle for 3D visualization
    view(ax1,40,35);
    grid(ax1,'on');

    % Enforce consistent aspect ratio between axes
    pbaspect(ax1, [1 1 0.35]);

    
    %% ========================= Right panel (polluted snow) ===============
    ax2 = axes(fig);

    % Position of the subplot (kept identical across panels for...
    % consistent layout)
    ax2.Position = [0.53 0.42 0.42 0.53]; 
    
    % Create empty surface object (colors will be updated dynamically)
    hSurf2 = surf(ax2, LonGrid, LatGrid, zeros(size(LonGrid)),...
        zeros(size(LonGrid)));
    shading(ax2,'interp');
    colormap(ax2,'turbo');

    % Colorbar placed outside the axes to preserve identical panel dimensions
    c2 = colorbar(ax2);
    c2.Location = 'eastoutside';

    % Color scaling for melt rate values
    clim(ax2,[0 6]);
    
    xlabel(ax2,'Longitude'); ylabel(ax2,'Latitude'); zlabel(ax2,'Altitude');
    title(ax2,'Polluted snow melt (mm/h)');

    % View angle for 3D visualization
    view(ax2,40,35);
    grid(ax2,'on');
    
    % Enforce consistent aspect ratio between axes
    pbaspect(ax2, [1 1 0.35]);


    %% ====================== Plotting time series of  ====================
    %% ============== mean melt rates for unpolluted & polluted ==========
    
    % bottom panel
    ax3 = subplot(3,1,3);

    % Positioning the figure to align with the 3D panels above
    ax3.Position = [0.10 0.10 0.82 0.22];

    % Keep existing plot when adding new lines
    hold(ax3, 'on');

    % Mean melt rate curves (pure, polluted, and their difference)
    plot(ax3, times, melt_pure_mean, 'b', 'LineWidth', 2);
    plot(ax3, times, melt_polluted_mean, 'r', 'LineWidth', 2);
    plot(ax3, times, melt_diff_mean, 'k', 'LineWidth', 2);

    ylabel(ax3,'Melt rate (mm/h)');
    xlabel(ax3,'Time');
    grid(ax3,'on');

    % Title and legend for the three time series
    title(ax3,'Mean melt rates (Pure, Polluted) + Difference');
    legend(ax3, 'Pure mean','Polluted mean','Polluted−Pure');

    % Format time axis in hours and minutes
    datetick(ax3,'x','HH:MM','keeplimits');

   
    %% ================= Time slider for the 3D models ==================

    % Create a label above the slider
    uicontrol('Style','text','String','Time index','Units','normalized',...
              'Position',[0.35 0.02 0.2 0.04],'BackgroundColor',[1 1 1]);

    % Create the slider to navigate through time steps
    hSlider = uicontrol('Style','slider','Min',1,'Max',nTimes,'Value',1,...
                        'Units','normalized','Position',[0.35 0.00 0.4 0.04],...
                        'SliderStep',[1/(nTimes-1) 1/(nTimes-1)],...
                        'Callback',@(src,evt) updateFrame(round(get(src,'Value'))));



   
    %% ========= Update frame function for 3D melt sufraces ===========

    function updateFrame(i)
        % Get the current time corresponding to the slider index
        t = times(i);
        idx = (time == t);

        % Extract data for the current time step
        lon_t  = lon(idx);
        lat_t  = lat(idx);
        alt_t  = alt(idx);
        mf_t   = melt_pure(idx);
        mp_t   = melt_polluted(idx);

        % Aggregate duplicate spatial points by averaging
        [uniquePts, ~, ic] = unique([lon_t, lat_t], 'rows');
        mf_u = accumarray(ic, mf_t, [], @mean);
        mp_u = accumarray(ic, mp_t, [], @mean);
        alt_u = accumarray(ic, alt_t, [], @mean);
    
        % Interpolate values onto the fixed grid for smooth 3D surfaces
        FG = griddata(uniquePts(:,1), uniquePts(:,2), mf_u,...
            LonGrid, LatGrid, 'natural');
        PG = griddata(uniquePts(:,1), uniquePts(:,2), mp_u,...
            LonGrid, LatGrid, 'natural');
        AG = griddata(uniquePts(:,1), uniquePts(:,2), alt_u,...
            LonGrid, LatGrid, 'natural');
        
        % Replace NaNs from interpolation with zeros
        FG(isnan(FG)) = 0;
        PG(isnan(PG)) = 0;
        AG(isnan(AG)) = 0;

        % Update 3D surfaces with altitude scaled and color mapped to ...
        % melt rates
        set(hSurf1, 'ZData', AG*zscale, 'CData', FG);
        set(hSurf2, 'ZData', AG*zscale, 'CData', PG);

        % Update titles with current time
        title(ax1, "Pure snow melt – " + string(t));
        title(ax2, "Polluted snow melt – " + string(t));
    end

    % Display initial frame at the first time step
    updateFrame(1);

     % Save the figure
    resultsFolder = '../results/';

    % Save in .fig format
    savefig(fig, fullfile(resultsFolder, 'melt_visualization.fig'));
end
