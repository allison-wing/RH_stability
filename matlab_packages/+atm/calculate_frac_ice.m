function [fliq,fice,varargout] = calculate_frac_ice(T,varargin)
% Calculate the fraction of liquid and ice for a saturation adjustment
%
% [fliq,fice,[dfliqdT,dficedT]] = calculate_frac_ice(T,[type,ice,deltaT])
%

c = atm.load_constants(varargin{1:end});



if c.ice == 0 
    fliq = ones(size(T));
    fice = zeros(size(T));
    dfliqdT = 0;
    dficedT = 0;
else
    
    if c.deltaT == 0
        fice = (T<c.T0).*1.0;
    else
       fice = -( T-(c.T0) )./(c.deltaT);
       fice(fice<0) = 0;
       fice(fice>1) = 1;
    end
    fliq = 1-fice;
end

if nargout > 2
    
   dfliqdT = zeros(size(T));
   
   
   if c.ice ~= 0
      dfliqdT(T>c.T0-c.deltaT & T<=c.T0) = 1./c.deltaT;
   end
   
   varargout{1} = dfliqdT;
   if nargout >2; varargout{2} = -dfliqdT; end;
   
end



end
