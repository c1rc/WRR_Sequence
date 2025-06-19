%Code to plot stratigraphy from netCDF file
dtime =5.;
filename = ncfile;
age = ncread(filename,'at_layer:age'); %reading in age variable
age = (double(age))/1000.; %years -> ky 
model_timestep = age(1,2)-age(1,1); % difference between two timesteps
[width,time] = size(age); %getting the dimensions of age, i.e. x and layers

p_age = age(3,pl_tstep);
pstep = int64((p_age - age(3,1))/model_timestep+1);

dtime=5; % creating time boundaries (black lines)
deposit_thickness = ncread(filename,'at_layer:thickness');
t0 = ncread(filename,'at_layer:t0');
wd = ncread(filename,'at_layer:water_depth');
shore = ncread(filename,'at_grid:x_of_shore');
shelf = ncread(filename,'at_grid:x_of_shelf_edge');
%curve = ncread(filename,'at_cell:curvature');
sl = ncread(filename,'at_grid:sea_level__elevation');
perc_sand=ncread(filename,'at_layer:percent_sand');
porosity=ncread(filename,'at_layer:porosity');

psl = sl(pstep); % Present Sea Level
x = ncread(filename,'x_of_cell');
x = x/1000.; %meters -> km
dtime = dtime/model_timestep;  %draw every dtime layers
horzn = zeros(width,pstep+1);

% clf
for j = 1:width  %calculate depth to layers
    horzn(j,pstep+1) = -1*wd(j,pstep)+sl(pstep);
    for i= pstep:-1:1
            if isnan(deposit_thickness(j,i))
                dp = 0.;
            else 
                dp = deposit_thickness(j,i);
            end
            horzn(j,i) = horzn(j,i+1) - dp;
    end
end


% calculate outline of shoreface
dx = x(2)-x(1);
sf = 6 / dx; %shoreface width = 3/alpha

sh = shore/1000.; %distance to shore in km
xsh = sh/dx;      %number of grid points to shore
nsh = floor(xsh);    %interger number grid points before shore
for i = 1:pstep
    if nsh(i) >= width  
        nsh(i) = width - 2*sf;  
    end
end

nsh(1) = nsh(2);
top = nsh(pstep);
bot = nsh(2);

    %t1 = gradient(gradient (t0));
for i = 1:time-9  %zero curvatures landward of shore
    if nsh(i) > width-2  
        nsh(i) = width-2;  
    end
    a = nsh(i)+sf;
    if a > width  
        a = width;  
    end
    t0(1:a,i) = 0.; 
end
[c,ci] = max(t0.*wd); %shelf edge = thickest part of current deposity 
xw = zeros (width-nsh(pstep)+2,1); %plot light blue water
yw = zeros (width-nsh(pstep)+2,1);

xw(1:width-nsh(pstep)+1) = x(nsh(pstep):width);
yw(1:width-nsh(pstep)+1) = -wd(nsh(pstep):width,pstep)+sl(pstep);

xw(width-nsh(pstep)+2) = x(width);
yw(width-nsh(pstep)+2) = sl(pstep);

% fill(xw,yw,[.9 1 1],'LineStyle','none') %water = light blue
% hold on
xs = [x(top) x(width)];
ys = [sl(pstep) sl(pstep)];
% plot(xs,ys,'b')

% calculate outline of land
xg = zeros(pstep + top + bot-1,1); %creating array for xg
yg = zeros(pstep + top + bot-1,1); %creating array for yg

% go along bottom of land
xg(1:bot) = x (1:bot);
yg(1:bot) = horzn(1:bot,2);

% follow shoreline through layers
xg(bot+1:pstep+bot-1) = sh(2:pstep);pstep*2 -1 + sf*2;
for i = 1:pstep-1
    land = horzn(nsh(i+1),i+2);
    sea  = horzn(nsh(i+1)+1,i+2);
    frac = xsh(i+1)-nsh(i+1);
    yg(bot+i) = land*(1-frac) + sea*frac;
end

% back across top of land
xg(pstep+bot:pstep+bot+top-1) = x (top:-1:1);
yg(pstep+bot:pstep+bot+top-1) = horzn(top:-1:1,pstep+1);

%fill(xg,yg,[.5 .875 .5],'LineStyle','none') %land = green 

xsf = zeros(pstep*2 -1 + sf*2,1);
ysf = zeros(pstep*2 -1 + sf*2,1);

% follow shoreline up through layers
xsf(1:pstep-1) = sh(2:pstep);
for i = 1:pstep-1
    land = horzn(nsh(i+1),i+2);
    sea  = horzn(nsh(i+1)+1,i+2);
    frac = xsh(i+1)-nsh(i+1);
    ysf(i) = land*(1-frac) + sea*frac;
end

% across top shoreface
xsf(pstep:pstep+sf-1) = x(top+1:top+sf);
ysf(pstep:pstep+sf-1) = horzn(top+1:top+sf,pstep);

% down lower shoreface through layers
for i = 1:pstep-1
    xsf(pstep+sf+i-1) = sh(pstep-i+1)+sf*dx;
    land = horzn(nsh(pstep+1-i)+sf,pstep+2-i);
    sea  = horzn(nsh(pstep+1-i)+sf+1,pstep+2-i);
    frac = xsh(pstep+1-i)-nsh(pstep+1-i);
    ysf(pstep+sf+i-1) = land*(1-frac) + sea*frac;
end

%across bottom shoreface
xsf(2*pstep+sf-1:2*pstep+2*sf-1) = x(bot+sf:-1:bot);
ysf(2*pstep+sf-1:2*pstep+2*sf-1) = horzn(bot+sf:-1:bot,2);

%fill(xsf,ysf,[1 1 0],'LineStyle','none') %shoreface = yellow

% calculate outline of shelf
s_top = ci(pstep) - (nsh(pstep)+sf);
s_bot = ci(2) - (nsh(2)+sf);
%find amount of jumps in shelfedge and add length to arrays
dci = diff(ci(2:pstep));    
for i=1:pstep-2
    if abs(dci(i)) < 20 
        dci(i) = 0; 
    end
end
[a b]= size(nonzeros(dci));
plus = sum(abs(dci))-a;

xse = zeros(2*(pstep-1) + s_top + s_bot+plus-1,1);
yse = zeros(2*(pstep-1) + s_top + s_bot+plus-1,1);

%down lower shoreface
for i = 1:pstep-1
    xse(i) = sh(pstep-i+1)+sf*dx;
    land = horzn(nsh(pstep-i+1)+sf,pstep-i+2);
    sea  = horzn(nsh(pstep-i+1)+sf+1,pstep-i+2);
    frac = xsh(pstep-i+1)-nsh(pstep-i+1);
    yse(i) = land*(1-frac) + sea*frac;
end

%across lower shelf
xse(pstep:pstep+s_bot-1) = x(nsh(2)+sf+1:ci(2));
yse(pstep:pstep+s_bot-1) = horzn(nsh(2)+sf+1:ci(2),2);

xse(pstep+s_bot) = x(ci(2));
yse(pstep+s_bot) = horzn(ci(2),3);
extra = 0;
%up shelfedge
for i = 3:pstep-1
    if abs(dci(i-1)) > 1
        %p=pstep+s_bot+i-2
        sn = sign(dci(i-1));
        for j = 1:abs(dci(i-1))
            xse(pstep+s_bot+i-3+j+extra) = x(ci(i))+sn*dx*j;
            yse(pstep+s_bot+i-3+j+extra) = horzn(ci(i)+sn*j,i+1);
            %j
        end
    extra = extra + abs(dci(i-1))-1;
    end
    xse(pstep+s_bot+i-2+extra) = x(ci(i+1));
    yse(pstep+s_bot+i-2+extra) = horzn(ci(i+1),i+1);
    %time+s_bot+i-1+extra
end

%across top shelf
xse(2*pstep+s_bot+extra-2:2*(pstep-1)+s_top+s_bot+extra-1) = x(ci(pstep):-1:nsh(pstep)+sf+1);
yse(2*pstep+s_bot+extra-2:2*(pstep-1)+s_top+s_bot+extra-1)=horzn(ci(pstep):-1:nsh(pstep)+sf+1,pstep);

dtime = int64(dtime);
xb = zeros(width+2,1);
yb = zeros(width+2,1);

xb(1:width) = x(1:width);
yb(1:width) = horzn(1:width,2);
b = xlim;
d = ylim;
xb(width+1) = x(width);
xb(width+2) = b(1);
yb(width+1:width+2) = d(1);
% xb=[xb xb(1,1)]
% fill(xb,yb,[.9 .87 .84],'LineStyle','none') %below model = tan

erode = NaN(width,pstep+1); %plot erosion surfaces 
for i = 3:pstep+1
    for j = 1:width
        if ~isnan(deposit_thickness(j,i-1)) && isnan(deposit_thickness(j,i-2))
            erode(j,i) = 1;
        end
    end
end
dtime = int64(dtime);

