%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function produces three dimensional maps of the S1a, S1b, and S1a/S1b  means
%
%  Inputs:
%       - s1ab_x                  - S1a/S1b x position after golden selection
%       - s1ab_y                  - S1a/S1b y position after golden selection
%       - s1ab_z                  - S1a/S1b z position after golden selection
%       - det_edge                - Edge of the drift time histogram
%       - s1a_phe_both            - S1a pulse areas after golden selection
%       - s1b_phe_both            - S1b pulse areas after golden selection
%
%  Outputs:
%       - s1ab_xyz_xbins          - 1xN vector defining x binning in XYZ map
%       - s1ab_xyz_ybins          - 1xM vector defining y binning in XYZ map
%       - s1ab_xyz_zbins          - 1xM vector defining z binning in XYZ map
%       - s1a_xyz_mean            - MxNxS matrix of S1a means 
%       - s1a_xyz_mean_err        - MxNxS matrix of S1a mean errors
%       - s1b_xyz_mean            - MxNxS matrix of S1b means
%       - s1b_xyz_mean_err        - MxNxS matrix of S1b mean errors
%       - s1ab_xyz_mean           - MxNxS matrix of S1a/S1b means
%       - s1ab_xyz_mean_err       - MxNxS matrix of S1a/S1b mean errors
% 
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s1ab_xyz_xbins, s1ab_xyz_ybins, s1ab_xyz_zbins, s1a_xyz_mean, s1a_xyz_mean_err, s1b_xyz_mean, s1b_xyz_mean_err, s1ab_xyz_mean, s1ab_xyz_mean_err] = MeasureS1aS1b3DDep(s1ab_x, s1ab_y, s1ab_z, det_edge, s1a_phe_both, s1b_phe_both)

%Define bin sizes
s1ab_xyz_numbins     = floor((length(s1ab_z)/200)^(1/3)); %number of bins in one direction to get ~200 events per bin
s1ab_xyz_zstep       = det_edge/s1ab_xyz_numbins;
s1ab_xyz_xstep       = 50/s1ab_xyz_numbins;
s1ab_xyz_ystep       = 50/s1ab_xyz_numbins;
s1ab_xyz_xmax        = 25;
r_max                = 25;

%Define bins
s1ab_xyz_zbins       = 10+s1ab_xyz_zstep/2:s1ab_xyz_zstep:10+s1ab_xyz_zstep/2+s1ab_xyz_zstep*s1ab_xyz_numbins; %20 us
s1ab_xyz_xbins       = (-r_max+s1ab_xyz_xstep/2):s1ab_xyz_xstep:r_max;
s1ab_xyz_ybins       = (-r_max+s1ab_xyz_ystep/2):s1ab_xyz_ystep:r_max;

%Preallocate memory
s1a_xyz_mean         = zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_xyz_mean_err     = zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean         = zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean_err     = zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));  
Count_S1ab_3D        = zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            %Indicies for the map (x,y swapped due to matrix x-row, y-column definition)
            l = int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m = int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n = int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
            %Make a cut to select events in the bin
            q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    
            %Count the number of events per bin
            Count_S1ab_3D(l,m,n) = length(s1a_phe_both(q));
            
            %Calculate the map entry if there are 100 events in the bin, otherwise set to zero
            if (Count_S1ab_3D(l,m,n) >= 100) 
                s1a_z_hist_bins         = (0:2:400);
                s1a_z_hist_cut          = (q.' & inrange(s1a_phe_both,[0 400]));
                s1a_z_fit               = fit(s1a_z_hist_bins.',hist(s1a_phe_both(s1a_z_hist_cut),s1a_z_hist_bins).','gauss1');
                s1a_xyz_mean(l,m,n)     = s1a_z_fit.b1;
                s1a_xyz_mean_err(l,m,n) = s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both(s1a_z_hist_cut));
    
                s1b_z_hist_bin          = (0:2:300);
                s1b_z_hist_cut          = (q.' & inrange(s1b_phe_both,[0 300]) );
                s1b_z_fit               = fit(s1b_z_hist_bin.',hist(s1b_phe_both(s1b_z_hist_cut),s1b_z_hist_bin).','gauss1');
                s1b_xyz_mean(l,m,n)     = s1b_z_fit.b1;
                s1b_xyz_mean_err(l,m,n) = s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both(s1b_z_hist_cut));              
             else 
                s1a_xyz_mean(l,m,n)     = 0;
                s1a_xyz_mean_err(l,m,n) = 0;
                s1b_xyz_mean(l,m,n)     = 0;
                s1b_xyz_mean_err(l,m,n) = 0;
            end
                                      
        end
    end
k %#ok<NOPRT>
end

     s1ab_xyz_mean     = s1a_xyz_mean./s1b_xyz_mean;
     s1ab_xyz_mean_err = sqrt( (s1b_xyz_mean_err.*s1a_xyz_mean./(s1b_xyz_mean.^2)).^2 + (s1a_xyz_mean_err./s1b_xyz_mean).^2);
     
end
