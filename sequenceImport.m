% Importing top and bottom surfaces for every from Sequence to CVFEM code.

seq_info=ncinfo(ncfile);
ncdisp(ncfile);

thickness=ncread(ncfile, 'at_layer:thickness')';
sed_dep_thickness=ncread(ncfile, 'at_node:sediment_deposit__thickness')';
bedrock_elev=ncread(ncfile, 'at_node:bedrock_surface__elevation')';
x_of_shore=ncread(ncfile, 'at_grid:x_of_shore')';
x_of_node=ncread(ncfile, 'x_of_node')';
x_of_cell=ncread(ncfile, 'x_of_cell')';
sediment_load=ncread(ncfile, 'at_grid:sediment_load')';
sea_level=ncread(ncfile, 'at_grid:sea_level__elevation')';
y_of_node=ncread(ncfile, 'y_of_node')';
t0=ncread(ncfile, 'at_layer:t0')';
perc_sand=ncread(ncfile, 'at_layer:percent_sand')';
porosity=ncread(ncfile, 'at_layer:porosity')';
wd = ncread(ncfile,'at_layer:water_depth')';

%Automation of input
sz_sl=size(sea_level);
sz_thickness=size(thickness);
sz_bed_elev=size(bedrock_elev);

tot_tstep=sz_sl(1,1);%total timesteps extracted from Sea level
seq_nodes=sz_bed_elev(1,2);%number of nodes from bed elevation 

csum_sdthickness=cumsum(sed_dep_thickness,1);
csum_thickness=cumsum(thickness,1);
sumthickness=sum(thickness,1);

bedrock_elev(:,1)=[];
bedrock_elev(:,sz_bed_elev(2)-1)=[];

sed_dep_thickness(:,1)=[];
sed_dep_thickness(:,sz_bed_elev(2)-1)=[];

for i=1:tot_tstep
%     figure(2)
    sl=sea_level(i,1);
    topsurface(i,:)=-1*wd(i,:)+sl;
    bedrock_position(i,:)=bedrock_elev(i,:);
    sed_surf(i,:)=bedrock_position(i,:)+thickness(1,:);
end



