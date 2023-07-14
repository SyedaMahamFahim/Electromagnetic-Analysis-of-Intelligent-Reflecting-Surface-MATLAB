F0 = 5.8e9; 
lambda = physconst("LightSpeed")/F0;
k = 2*pi/lambda;

phi_inc = 0;
theta_inc = -45:1:45;
elevation_inc = 90 + theta_inc;
radius_inc = 100*lambda;
[x_inc,y_inc,z_inc] = sph2cart(phi_inc*pi/180,elevation_inc*pi/180,radius_inc);
inc_loc = [x_inc' y_inc' z_inc'];

phi_obs = 0;
theta_obs = 20;
elevation_obs = 90 - theta_obs;
radius_obs = 100*lambda;
[x_obs,y_obs,z_obs] = sph2cart(phi_obs*pi/180,elevation_obs*pi/180,radius_obs);
obs_loc = [x_obs, y_obs, z_obs];

rf = design(reflector,F0);
rf.Spacing = lambda/10;
rf.GroundPlaneLength = 0.5*lambda;
rf.GroundPlaneWidth = 0.5*lambda;
figure;
show(rf)

ifa = rectangularArray(Element==rf, Size==[10 10],...
    ColumnSpacing==rf.GroundPlaneLength,...
    RowSpacing==rf.GroundPlaneWidth);
figure;
show(ifa)
title("Finite IRS")

irs = infiniteArray(Element=rf);
irs.ScanAzimuth = phi_obs;
irs.ScanElevation = elevation_obs;
numSummationTerms(irs,20);
figure;
show(irs)
title("Infinite IRS")

% Compute the Array factor
AF=hArrayFactorCalc(irs,F0);

MagReflection=zeros(1,numel(theta_inc));
PhaseReflection=zeros(1,numel(theta_inc));

for mm = 1:numel(theta_inc)

    % Construct a planeWaveExcitation object with the IRS as element
    [pw,pol,dir] = hcalcPlaneWaveIncidence(theta_inc(mm),phi_inc,irs);

    % Compute incident field upon the IRS
    Ein = pol;
    Einc = dot(Ein,pol/norm(pol));

    % Compute the outward scattered field from the IRS
    [Eo,~] = EHfields(pw,F0,obs_loc');
    Eo = Eo*AF;
    Eobs = dot(Eo,pol/norm(pol));

    % Compute the reflection coefficient of the IRS
    MagReflection(1,mm) = abs(Eobs/Einc);
    PhaseReflection(1,mm) = (angle(Eobs/Einc))*180/pi;
end
figure;
plot(theta_inc,MagReflection,"b")
title("Reflection Coefficient Magnitude Against Incidence Angle")
xlabel("Incidence Angle (\theta_{inc} in degree)")
ylabel("Field Reflection Coefficient")

figure;
plot(theta_inc,PhaseReflection,"b")
title("Reflection Coefficient Phase Against Incidence Angle")
xlabel("Incidence angle (\theta_{inc} in degree)")
ylabel("Field Reflection Coefficient Phase")