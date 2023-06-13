% Постановка задачи 1
% Вычислите максимально точно геометрическое (евклидово) расстояние между шпилями Эйфелевой
% башни и Останкинской телебашни. В расчетах учтите превышение геоида над эллипсоидом WGS-84
clc

obj1_coords = struct();
obj2_coords = struct();

obj1_coords.name = 'Eifel Tower';
obj1_coords.lat = dms2rad([48 51 29]);
obj1_coords.lon = dms2rad([2 17 40]);
obj1_coords.altAboveSea = 38;
obj1_coords.height = 325;

obj2_coords.name = 'Ostankino Tower';
obj2_coords.lat = dms2rad([55 49 11]);
obj2_coords.lon = dms2rad([37 36 42]);
obj2_coords.altAboveSea = 160;
obj2_coords.height = 160;

% EGM2008
WGSh1 = 44.6404;
WGSh2 = 14.7992;

obj1_r = calcECEF(obj1_coords, WGSh1);
obj2_r = calcECEF(obj2_coords, WGSh2);

dist = norm(obj1_r - obj2_r);
disp(['Расстояние между шпилями ', num2str(dist), ' м'])
disp('----------------------------------------------')

mapWidth = 100; % см

xy1 = projectToMercatorMap(obj1_coords, mapWidth);
xy2 = projectToMercatorMap(obj2_coords, mapWidth);

disp(['Координаты точек на карте Меркатора шириной ', num2str(mapWidth), ' см'])
disp([obj1_coords.name, ':x = ', num2str(xy1(1)), ', y = ', num2str(xy1(2))])
disp([obj2_coords.name, ':x = ', num2str(xy2(1)), ', y = ', num2str(xy2(2))])

function posVec = calcECEF(coords, WGSh)
    posVec = [];

    meanEarthRadius = 6378137;
    heightAboveGeoid = coords.altAboveSea + coords.height + WGSh;
    meanEcc = 0.081819221456;

    meridRadius = meanEarthRadius / (sqrt(1 - meanEcc^2 * sin(coords.lat)^2));
    polarRadius = meanEarthRadius * (1 - meanEcc^2) / (sqrt(1 - meanEcc^2 * sin(coords.lat)^2));

    posVec(1) = (meridRadius + heightAboveGeoid) * cos(coords.lat) * cos(coords.lon);
    posVec(2) = (meridRadius + heightAboveGeoid) * cos(coords.lat) * sin(coords.lon);
    posVec(3) = (polarRadius + heightAboveGeoid) * sin(coords.lat);
end

function xy = projectToMercatorMap(coords, mapWidth)
    e = 0.081819221456;

    lat = coords.lat;
    lon = coords.lon;
    lon0 = 0;
    
    c = mapWidth / pi;

    x = c * (lon - lon0);
    y = c * (log(tan(lat/2 + pi/4)) - e/2 * log((1 + e * sin(lat))/(1 - e * sin(lat))));
    
    xy = [x,y];
end

function angleRad = dms2rad(dms)
    grads = dms(1);
    minutes = dms(2);
    secs = dms(2);
    
    angleRad = (grads + minutes/60 + secs/3600) * pi/180;
end