function dists=CalcDists(x)

% Set radius of earth in m
r=6378.137*1000;

% Convert latitude and longitude from degrees to radians
lat=x.latitude/180*pi;
lon=x.longitude/180*pi;

lat_diff=abs(bsxfun(@minus,lat,lat'));
lon_diff=abs(bsxfun(@minus,lon,lon'));

dists=r*2*asin(sqrt((sin(lat_diff/2)).^2+cos(lat)*cos(lat').*(sin(lon_diff/2)).^2));
