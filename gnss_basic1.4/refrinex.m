function [gps_time, prn, varargout] = refrinex(obs, obs_qual, clk_off);
% Function to reformat RINEX-II data.
%
% [gps_time, prn, data1, data2, data3, ..., datan] = refrinex(obs, obs_qual, clk_off);
%
% Function to reformat RINEX data from RD_RNX_IN (obs matrix) into separate matrices.
% This function takes the output from rd_rnx_o and reformats the observation matrix
% into a separate matrix for each observation type.  The observation type is defined
% in the type_str variable.
%
% Input:
%    obs      - observation data (nxm) (number_of_obs x 3+number_of_data_types)
%
%                    1         2     3     4        5        6
%                [GPS_week, GPS_sec, PRN, data(1), data(2), data(3), ...]
%                data fields depend on RINEX header file.  Use type_str
%                variable to identify each observation type.  Observations are
%                loaded in the data(n) variable in the same order as the typ_str.
%                Example obs matrix with 4 observation types (C1    L1    D1    P2)
%                   1      2     3        4            5           6         7
%                [1023  856792  28  20855400.126 109595864.844 -969.622 20855402.416]
%    obs_qual - observation quality matrix (nx2) [LLI Signal_strength]
%               LLI - loss of lock indicator (0-7) see RINEX format for meaning
%               Signal strngth (1-9, min-max), 0 - unknown/don't care
%    clk_off  - clock offset (nx1) seconds if available (optional)
%
% Output:
%    gps_time - GPS time tag for each time [GPS_week GPS_sec] (kx2). 
%               Each element corresponds to a row in thedata matrices.  
%               k is the number of obseration times in the data set.  
%               If the clock offset is provided in the input (3 input parameters), 
%               a 3-rd column is added to the time matrix that has the clock offset 
%               [GPS_week GPS_sec clk_off](kx3)
%    prn      - PRN numbers for each column in the satellite data matrices (jx1)
%    data1    - data matrix of observations corresponding to the first data type
%               in the obs matrix (column 4 of obs) (kxj) (e.g. L1, C1, D1, etc.)
%    data2    - data matrix of observations corresponding to the second data type
%               in the obs matrix (column 5 of obs) (kxj)
%    ....
%    datan    - data matrix of observations corresponding to the n-th data type
%               in the obs matrix (column 4+n-1 of obs) (kxj)
%
% Each data matrix has the following format
%           prn_1  prn_2 prn_3 ... prn_j
%     t0    L1_01  L1_02 L1_03 ... L1_0j
%     t1    L1_11  L1_12 L1_13 ... L1_1j
%     t2    L1_21  L1_22 L1_23 ... L1_2j
%     ...    ...    ...   ...  ...  ...
%     tk    L1_k1  L1_k2 L1_k3 ... L1_kj
%
% Note: 1) A variable number of output arguments are available based on the number of
%          observation types input in the obs matrix. 
%       2) GPS time is kept without rollovers (e.g. week 1025 is week 1 with a rollover).
%   
% See also READEPH, RD_RNX_O 

% Written by: Jimmy LaMance 10/17/99
% Copyright (c) 1999 by Constell, Inc.

%%%%% BEGIN VARIABLE CHECKING CODE %%%%%
% declare the global debug mode
global DEBUG_MODE

% Initialize the output variables
gps_time =[];
prn =[];
vargout = [];     % variable lenght output arguments, cell array i=1:num_output

% Check the number of input arguments and issues a message if invalid
msg = nargchk(2,3,nargin);
if ~isempty(msg)
  fprintf('%s  See help on REFRINEX for details.\n',msg);
  fprintf('Returning with empty outputs.\n\n');
  return
end

estruct.func_name = 'REFRINEX';

% Develop the error checking structure with required dimension, matching
% dimension flags, and input dimensions.
estruct.variable(1).name = 'obs';
estruct.variable(1).req_dim = [];
estruct.variable(1).var = obs;
  
estruct.variable(2).name = 'obs_qual';
estruct.variable(2).req_dim = [901 2];
estruct.variable(2).var = obs_qual;
  
if nargin==3,
  estruct.variable(3).name = 'clk_off';
  struct.variable(3).req_dim = [901 1];
  estruct.variable(3).var = clk_off;
end % if nargin==3
  
% Call the error checking function
stop_flag = err_chk(estruct);
  
if stop_flag == 1           
  fprintf('Invalid inputs to %s.  Returning with empty outputs.\n\n', ...
           estruct.func_name);
  return
end % if stop_flag == 1

%%%%% END VARIABLE CHECKING CODE %%%%%

%%%%% BEGIN ALGORITHM CODE %%%%%
% Get linear time from the obs matrix
t_lin = gpst2sec(obs(:,1:2));

% Find all of the unique time tags
t_lin_uniq = unique(t_lin);

% Sort all of the data in time order (it should already be this way, but 
% this insures no problems later)
[t_lin,I_sort] = sort(t_lin);
obs = obs(I_sort,:);
obs_qual = obs_qual(I_sort,:);
clk_off = clk_off(I_sort,:);

% Find the transition from one time to the next using intersect.  This
% gives the last time in the group that matches the others.
[t_int, I_all, I_uniq] = intersect(t_lin, t_lin_uniq);

% Because intersect gives the last time in the group, remove the last
% element of the output index because it lines up with the last
% observation.
I_all = I_all(1:end-1);

% Now use the cumsum function to fill in an indexing matrix that numbers the 
% observation times (e.g. [1 1 1 1 1 2 2 2 2... 2 3 3 ... 3 ... n]);
t_chng = zeros(size(t_lin));
t_chng(1) = 1;
t_chng(I_all+1) = ones(size(I_all));
I_uniq_times = find(t_chng == 1);
t_ndx = cumsum(t_chng);

% Now fill the output gps_time with the corresponding times (1:n)
gps_time = obs(I_uniq_times,1:2);

% Fill in the 3-rd column of the gps_time matrix if a clock offset is input
if nargin == 3
  gps_time(:,3) = clk_off(I_uniq_times);
end % if nargin == 3

% Find all of the unique prn numbers to size and order the output arrays
prn = unique(obs(:,3));

% Find the number of satellites, times, and observation types
num_sats = length(prn);
num_times = size(gps_time,1);
num_types = size(obs,2) - 3;

% Loop over the number of observation types and the number of satellites (don't
% loop over time or it will take forever).  This looping will be very small since
% it's loops over 10s of things, not thousands.
%keyboard
for i = 1:num_types
  % Size the output matrix (num_times x num_sats)
  clear data
  data = ones(num_times,num_sats) * NaN;
  
  % Loop over the satellites and fill in the data matrix
  for j = 1:num_sats
    I = find(obs(:,3) == prn(j));
    
    % Make sure we found some matching satellites
    if isempty(I)
      fprintf('Error in finding unique satellites in REFRINEX.\n');
      return
    end % if isempty(I)
    
    % Insert these observations
    data(t_ndx(I),j) = obs(I,i+3);
  end % for j = 1:num_sats
  
  % Copy the resulting data matrix into the vargour field
  varargout{i} = data;
  
end % for i = 1:num_types

%%%%% END ALGORITHM CODE %%%%%


