%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function does 1D cubic interpolation within the bounds of the map,
%       and linear extrapolation based on the first/last 20% of the map 
%
%  Inputs:
%       - X             - X bins of the map
%       - Y             - Y bins of the map
%       - X0            - Value to interpolate/extrapolate at from the map
%                         Can be a vector specifying multiple locations                         
%
%  Outputs:
%       - first_interp  - The result of the interpolation/extrapolation at X0
%
%  Author:
%       - Richard Knoche
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ first_interp ] = RK_1DCubicInterp_LinearExtrap( X, Y, X0)

%Deal with multiple 
if length(X0)>1
    first_interp = interp1(X,Y,X0,'v5cubic');

    %find the working value of X0 closest to the NaN value of X0
    badX0        = X0(isnan(first_interp)); %#ok<NASGU>
    goodX0       = X0(~isnan(first_interp));
    goodY0       = first_interp(~isnan(first_interp));

    %sort to find the first and last 20% cubic values
    [sorted_goodX0, sorted_index] = sort(goodX0);
    sorted_goodY0                 = goodY0(sorted_index);
    first20_fit                   = polyfit(sorted_goodX0(1:floor(length(sorted_goodX0)*0.2)),sorted_goodY0(1:floor(length(sorted_goodX0)*0.2)),1);
    last20_fit                    = polyfit(sorted_goodX0(ceil(length(sorted_goodX0)*0.8):end),sorted_goodY0(ceil(length(sorted_goodX0)*0.8):end),1);

    %Doing linear extrapolation for events on the edges
    first_interp(isnan(first_interp) & X0<=sorted_goodX0(1))  = polyval(first20_fit,X0(isnan(first_interp) & X0<=sorted_goodX0(1)));
    first_interp(isnan(first_interp) & X0>=sorted_goodX0(end))= polyval(last20_fit,X0(isnan(first_interp) & X0>=sorted_goodX0(end)));

    %clean up for anything that is NaN in the middle of the interpolation
    %find the working value of X0 closest to the NaN value of X0
    badX0=X0(isnan(first_interp));
    if length(badX0)>= 1;
        goodX0 = X0(~isnan(first_interp));
        goodY0 = first_interp(~isnan(first_interp));

    for i=1:length(badX0);
       temp            = abs(badX0(i)-goodX0);
       [minx idx]      = min(temp); %#ok<ASGLU>
       FixedX0(i)      = goodX0(idx); %#ok<AGROW,NASGU>
       FixedX0Index(i) = idx; %#ok<AGROW>
    end

    first_interp(isnan(first_interp))=goodY0(FixedX0Index);
    end


else %Deal with the case of single X0 value
  
    first_interp = interp1(X,Y,X0,'v5cubic'); %normalize to right below the gate

    if isnan(first_interp)

        first_interp = interp1(X,Y,X,'v5cubic'); %normalize to right below the gate

        %find the working value of X0 closest to the NaN value of X0
        badX   = X(isnan(first_interp)); %#ok<NASGU>
        goodX  = X(~isnan(first_interp));
        goodY  = first_interp(~isnan(first_interp));

        %sort to find the first and last 20% cubic values
        [sorted_goodX, sorted_index] = sort(goodX);
        sorted_goodY                 = goodY(sorted_index);
        first20_fit                  = polyfit(sorted_goodX(1:floor(length(sorted_goodX)*0.2)),sorted_goodY(1:floor(length(sorted_goodX)*0.2)),1);
        last20_fit                   = polyfit(sorted_goodX(ceil(length(sorted_goodX)*0.8):end),sorted_goodY(ceil(length(sorted_goodX)*0.8):end),1);

        if X0 <= sorted_goodX(1)
            first_interp = polyval(first20_fit,X0);
        elseif X0 >= sorted_goodX(end)
            first_interp = polyval(last20_fit,X0);          
        else      
            first_interp = interp1(X,Y,X0,'nearest');
        end
    end
end

