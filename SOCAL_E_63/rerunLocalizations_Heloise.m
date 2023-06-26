load('DetFile')
whale = loc3D_DOAintersect(DET, hydLoc, paramFile);

[whale, wfit] = weightedSplineFit(whale, smoothingParam)

sig_h = 3.722;
sig_H = .1852;
c = 1488;
sig_sml = sqrt( (sig_H/c)^2 + (1/c * sig_h)^2);
Asml = (2*pi*sig_sml^2)^(-12/2); 

for wn = 1:numel(whale)
    Iuse = find(~isnan(whale{wn}.wloc(:,1)));
    wloc = whale{wn}.wloc(Iuse, :);
    wlocSmooth = whale{wn}.wlocSmooth(Iuse, :);
    TDOA = whale{wn}.TDOA(Iuse, :);

    whale{wn}.CIx = nan(length(whale{wn}.TDet), 2);
    whale{wn}.CIy =  whale{wn}.CIx;
    whale{wn}.CIz =  whale{wn}.CIx;

    whale{wn}.CIxSmooth = nan(length(whale{wn}.TDet), 2);
    whale{wn}.CIySmooth = whale{wn}.CIxSmooth;
    whale{wn}.CIzSmooth = whale{wn}.CIxSmooth;

    for i = 1:length(Iuse)

        [CIx, CIy, CIz] = calcCI(TDOA, wloc, hydLoc, H1, H2, Asml, Alrg, LOC);

        whale{wn}.CIx(i, :) = CIx;
        whale{wn}.CIy(i, :) = CIy;
        whale{wn}.CIz(i, :) = CIz;


        [CIx, CIy, CIz] = calcCI(TDOA, wlocSmooth, hydLoc, H1, H2, Asml, Alrg, LOC);

        whale{wn}.CIxSmooth(i, :) = CIx;
        whale{wn}.CIySmooth(i, :) = CIy;
        whale{wn}.CIzSmooth(i, :) = CIz;
    end
end
