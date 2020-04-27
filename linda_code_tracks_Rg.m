% Linda - January 2020 
% CODE FOR RADIUS OF GYRATION.
% Cumulates cell means within condition. 
% 4 conditions applied to select for tracks (see below). 
% Here, my movies imported are all 300 ms/fr (originally or rendered).

% Use at your own discretion. Double-check content.

clear all, warning off, close all
set(0,'DefaultFigureWindowStyle','docked')
addpath('../functions')
addpath('/Applications/Fiji.app/scripts')

all_MeanRg = [];
all_percent_highRg = [];
all_nRg = [];
all_Rg = [];
%all_MeanRg_min = [];

% Call on all files with the desired condition
cells = dir('movie1_07022019_0.0_CONTROL2_30.08x_36.00y_4119.0a_4.0000c_6.0d.xml');
%cells = dir('*_phosphoTAU_*.xml');   % e.g. cells of lowTAU condition
%cells = [dir('*_FYN1_*.xml'); dir('*_FYN2_*.xml')];
% Other conditions are *_mediumTAU.xml, _highTAU, _lowTAUFYN, _mediumTAUFYN,
% _highTAUFYN, _FYN, _CONTROL
%OR
%cells = [dir('*_lowTAUFYN_*.xml'); dir('*_mediumTAUFYN_*.xml'); dir('*_highTAUFYN_*.xml')];


for d=1:length(cells)
    cell=cells(d, :);  % calling initially on 1st cell of this condition, and then so on.
    
file_path_tracks = cell.name;
clipZ = true;     % Remove Z coordinates, if you know you can.
scaleT = true;       % Use physical time for T.
tracks = importTrackMateTracks(file_path_tracks, clipZ, scaleT);
n_tracks = numel( tracks );

% To save new track coodinates with higher than the minimum number of spots

new_tracks = [];
for s = 1 : n_tracks
    if (numel(tracks{s}(:, 1)) > 5) && (tracks{s}(1, 1) == 0) && (tracks{s}(1, 2) ~= tracks{s}(2, 2)) && (tracks{s}(1, 3) ~= tracks{s}(2, 3))
        % 1st condition: Select only tracks with minimum number of points, now 5
        % 2nd condition: Select only trajectories of the cell that commence in the 1st frame of movie
        % 3rd and 4th condition: Getting rid of tracks that are not true trajectories, not moving at all in x and y
        new_tracks1 = tracks{s}(:, 1);
        i = find(new_tracks1 < 90);  % All time frame before frame 90s are included. I do this if I want to cut movie at certain frame.
        new_tracks1 = new_tracks1(i);
        
        new_tracks2 = tracks{s}(:, 2);
        new_tracks2 = new_tracks2(i);
        
        new_tracks3 = tracks{s}(:, 3);
        new_tracks3 = new_tracks3(i);
        
        new_tracks_final =[new_tracks1,new_tracks2,new_tracks3];
        new_tracks{s} = new_tracks_final;
    end
end
new_tracks = new_tracks';
tracks = new_tracks;
tracks(cellfun('isempty', tracks)) = [];   % getting rid of all cells without values
n_tracks = numel( tracks );        % number of tracks

% If I want to make a figure of the trackmate trajectories
figure(1)      
hold on
c = lines(n_tracks);
n_tracks = numel( tracks );
for s = 1 : n_tracks
x = tracks{s}(:, 2);
y = tracks{s}(:, 3);
plot(x, y, '.-', 'Color', c(s, :))
end
axis equal


% Radius of gyration Rg for each track
all_x = [];
all_y = [];
all_a = [];
all_b = [];
all_c = [];
mean_x_cum = zeros(n_tracks, 1);
mean_y_cum = zeros(n_tracks, 1);
Rg = zeros(n_tracks, 1);
tau_int = zeros(n_tracks, 1);
for s = 1 : n_tracks
    % x and y coordinates of every spot in a track
    all_x{s} = tracks{s}(:, 2);
    all_y{s} = tracks{s}(:, 3);
    mean_x = mean(all_x{s}(:, 1));
    mean_y = mean(all_y{s}(:, 1));
    % mean of x and y coordinates for each track
    mean_x_cum(s) = mean_x;
    mean_y_cum(s) = mean_y;
    % calculation of a = (x-X)^2 and b = (y-Y)^2 for each x and y
    for x_value = 1 : numel(all_x{s}(:, 1))
        a = ((all_x{s}(x_value, 1))-mean_x).^2;
        b = ((all_y{s}(x_value, 1))-mean_y).^2;
        c = a+b;
        all_a{s}(x_value) = a;
        all_b{s}(x_value) = b;
        all_c{s}(x_value) = c;   % c values then are all summed for a track
    end
    Rg(s) = sqrt((1/numel(all_x{s}(:, 1)))*sum(all_c{s}(1, :))); % Rg for a track
end
% You can make a histogram of hist(Rg) to see the distribution for Rg's of a cell

% If I want to have a threshold of Rg for the Rg mean calculation
Rg_min = [];
Rg_min_i = find(Rg > 0.5); % um
Rg_min = Rg(Rg_min_i);

% Looking at percent of higher_than_minimum_Rg cargo
percent_highRg = (numel(Rg_min)/numel(Rg))*100;


MeanRg = mean(Rg);     % The Rg arithmetic mean for one cell
%pd = fitdist(Rg, 'Exponential');     % Fit probability distribution - exponentially **DOUBLE-CHECK THIS
%MeanRg = pd.mu;                      % From above exponential fit, extract the mu = mean, here of Rg.
%MeanRg = mle(Rg, 'distribution', 'exp');  % Another way, usually this.
%MeanRg = median(Rg);    % Not really the mean, but the median here.

nRg = numel(Rg); % number of tracks per cell
%nRg_min = numel(Rg_min);

%Cumulate all Rg from all cells within condition, with its corresponding
%tau intensity value
tau_int(1:nRg,1) = str2double(regexp(cell.name, '\d+\.\d+', 'match', 'once'));
Rg_int = [Rg tau_int];
all_Rg = vertcat(all_Rg, Rg_int);

tau_int = str2double(regexp(cell.name, '\d+\.\d+', 'match', 'once'));
MeanRg_int = [MeanRg tau_int];
all_MeanRg = vertcat(all_MeanRg, MeanRg_int);   % The Rg means for each cell within condition and the cell's tau intensity value extracted from cell filename

% For Rg_min
%MeanRg_min_int = [MeanRg_min tau_int];
%all_MeanRg_min = vertcat(all_MeanRg_min, MeanRg_min_int);

percent_highRg_int = [percent_highRg tau_int];
all_percent_highRg = vertcat(all_percent_highRg, percent_highRg_int);

all_nRg = vertcat(all_nRg, nRg);            % The number of tracks for each cell within condition


end



