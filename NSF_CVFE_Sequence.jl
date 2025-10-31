using CairoMakie
using Delaunay
using DelaunayTriangulation
using GeometryBasics
using Interpolations
using JLD2
using MAT
using Makie
using Makie.GeometryBasics
using NCDatasets
using Plots
using PolygonOps
using QHull
using SparseArrays
using StaticArrays
using Statistics
using VoronoiDelaunay

function generate_tri(x, y, ncol, cols)

    t1 = Int64[]
    first = 1
    for jj in 1:(cols-1)
        last = first + ncol[jj] + ncol[jj+1] - 1
        xcc = x[first:last]
        ycc = y[first:last]
        pts = hcat(xcc, ycc)
        mesh = delaunay(pts)
        tcc = mesh.simplices .+ (first - 1)
        append!(t1, tcc')

        first = first + ncol[jj]
    end

    return reshape(t1, 3, :)'
end


function stratigraphy_function(pl_tstep, model_timestep, age, sl, dtime, width, deposit_thickness, xc, shore)
    p_age = age[3, pl_tstep]
    pstep = Int64(round((p_age - age[3, 1]) / model_timestep + 1))
    psl = sl[pstep]


    x1 = xc ./ 1000.0
    dtime = dtime / model_timestep
    # Calculate depth to layers
    horzn = zeros(width, pstep + 1)

    for j in 1:width
        horzn[j, pstep+1] = -1 * wd[j, pstep] + sl[pstep]
        for i in pstep:-1:1
            dp = isnan(deposit_thickness[j, i]) ? 0.0 : deposit_thickness[j, i]
            horzn[j, i] = horzn[j, i+1] - dp
        end
    end

    # Calculate outline of shoreface
    dx = x1[2] - x1[1]
    sf = Int(6 / dx)  # Shoreface width = 3/alpha

    sh = shore ./ 1000  # Distance to shore in km
    xsh = sh / dx       # Number of grid points to shore
    nsh = floor.(Int, xsh)  # Integer number of grid points before shore

    for i in 1:pstep
        if nsh[i] >= width
            nsh[i] = width - 2 * sf
        end
    end

    nsh[1] = nsh[2]
    top = nsh[pstep]
    bot = nsh[2]


    # Zero curvatures landward of shore
    for i in 1:time-9
        if nsh[i] > width - 2
            nsh[i] = width - 2
        end
        a = nsh[i] + sf
        if a > width
            a = width
        end
        a = Int(a)
        t0[1:a, i] .= 0.0
    end

    # Shelf edge = thickest part of current deposit
    c, ci1 = findmax(t0 .* wd, dims=1)
    ci = [index[1] for index in ci1]


    # Calculate outline of water
    xw = zeros(width - nsh[pstep] + 2)
    yw = zeros(width - nsh[pstep] + 2)

    xw[1:width-nsh[pstep]+1] = x1[nsh[pstep]:width]
    yw[1:width-nsh[pstep]+1] = -wd[nsh[pstep]:width, pstep] .+ sl[pstep]

    xw[width-nsh[pstep]+2] = x1[width]
    yw[width-nsh[pstep]+2] = sl[pstep]

    # Calculate outline of land
    xg = zeros(pstep + top + bot - 1)
    yg = zeros(pstep + top + bot - 1)

    # Go along bottom of land
    xg[1:bot] = x1[1:bot]
    yg[1:bot] = horzn[1:bot, 2]

    # Follow shoreline through layers
    xg[bot+1:pstep+bot-1] = sh[2:pstep]
    for i in 1:pstep-1
        land = horzn[nsh[i+1], i+2]
        sea = horzn[nsh[i+1]+1, i+2]
        frac = xsh[i+1] - nsh[i+1]
        yg[bot+i] = land * (1 - frac) + sea * frac
    end

    # Back across top of land
    xg[pstep+bot:pstep+bot+top-1] = x1[top:-1:1]
    yg[pstep+bot:pstep+bot+top-1] = horzn[top:-1:1, pstep+1]


    # Calculate outline of shoreface
    xsf = zeros(pstep * 2 - 1 + sf * 2)
    ysf = zeros(pstep * 2 - 1 + sf * 2)

    # Follow shoreline up through layers
    for i in 1:pstep-1
        xsf[i] = sh[i+1]
        land = horzn[nsh[i+1], i+2]
        sea = horzn[nsh[i+1]+1, i+2]
        frac = xsh[i+1] - nsh[i+1]
        ysf[i] = land * (1 - frac) + sea * frac
    end

    # Across top shoreface
    xsf[pstep:pstep+sf-1] = x1[top+1:top+sf]
    ysf[pstep:pstep+sf-1] = horzn[top+1:top+sf, pstep]

    # Down lower shoreface through layers
    for i in 1:pstep-1
        idx = pstep + sf + i - 1
        xsf[idx] = sh[pstep-i+1] + sf * dx
        land = horzn[nsh[pstep+1-i]+sf, pstep+2-i]
        sea = horzn[nsh[pstep+1-i]+sf+1, pstep+2-i]
        frac = xsh[pstep+1-i] - nsh[pstep+1-i]
        ysf[idx] = land * (1 - frac) + sea * frac
    end

    # Across bottom shoreface
    xsf[2*pstep+sf-1:2*pstep+2*sf-1] = x1[bot+sf:-1:bot]
    ysf[2*pstep+sf-1:2*pstep+2*sf-1] = horzn[bot+sf:-1:bot, 2]


    # Calculate outline of shelf
    s_top = ci[pstep] - (nsh[pstep] + sf)
    s_bot = ci[2] - (nsh[2] + sf)

    # Find the amount of jumps in shelf edge and add length to arrays
    dci = diff(ci[2:pstep])
    for i in 1:pstep-2
        if abs(dci[i]) < 20
            dci[i] = 0
        end
    end

    # Calculate the size of non-zero elements in dci
    nonzero_dci = filter(x1 -> x1 != 0, dci)
    a = size(nonzero_dci, 1)

    # Calculate 'plus' based on the sum of absolute values of dci minus 'a'
    plus = sum(abs.(dci)) - a

    # Initialize xse and yse arrays
    xse = zeros(2 * (pstep - 1) + s_top + s_bot + plus - 1)
    yse = zeros(2 * (pstep - 1) + s_top + s_bot + plus - 1)

    # Down lower shoreface
    for i in 1:pstep-1
        xse[i] = sh[pstep-i+1] + sf * dx
        land = horzn[nsh[pstep-i+1]+sf, pstep-i+2]
        sea = horzn[nsh[pstep-i+1]+sf+1, pstep-i+2]
        frac = xsh[pstep-i+1] - nsh[pstep-i+1]
        yse[i] = land * (1 - frac) + sea * frac
    end

    # Across lower shelf
    xse[pstep:pstep+s_bot-1] = x1[nsh[2]+sf+1:ci[2]]
    yse[pstep:pstep+s_bot-1] = horzn[nsh[2]+sf+1:ci[2], 2]

    # Additional point at the end of the lower shelf
    xse[pstep+s_bot] = x1[ci[2]]
    yse[pstep+s_bot] = horzn[ci[2], 3]

    # Up shelf edge
    extra = 0
    for i in 3:pstep-1
        if abs(dci[i-1]) > 1
            sn = sign(dci[i-1])
            for j in 1:abs(dci[i-1])
                idx = pstep + s_bot + i - 3 + j + extra
                xse[idx] = x1[ci[i]] + sn * dx * j
                yse[idx] = horzn[ci[i]+sn*j, i+1]
            end
            extra += abs(dci[i-1]) - 1
        end
        xse[pstep+s_bot+i-2+extra] = x1[ci[i+1]]
        yse[pstep+s_bot+i-2+extra] = horzn[ci[i+1], i+1]
    end

    # Across top shelf
    range = 2*pstep+s_bot+extra-2:2*(pstep-1)+s_top+s_bot+extra-1
    xse[range] = x1[ci[pstep]:-1:nsh[pstep]+sf+1]
    yse[range] = horzn[ci[pstep]:-1:nsh[pstep]+sf+1, pstep]

    dtime = Int64(dtime)

    # Initialize arrays for plotting below the model
    xb = zeros(width + 2)
    yb = zeros(width + 2)

    xb[1:width] = x1[1:width]
    yb[1:width] = horzn[1:width, 2]

    # For demonstration, let's assume some values for b and d
    b = (0, maximum(x1))
    d = (minimum(yb), maximum(yb))

    xb[width+1] = x1[width]
    xb[width+2] = b[1]
    yb[width+1:width+2] .= d[1]

    # Initialize an array to plot erosion surfaces
    erode = fill(NaN, width, pstep + 1)

    for i in 3:pstep+1
        for j in 1:width
            if !isnan(deposit_thickness[j, i-1]) && isnan(deposit_thickness[j, i-2])
                erode[j, i] = 1
            end
        end
    end



    return xg, yg, xse, yse, xsf, ysf, xb, yb
end

# Create an output folder for the results
isdir("./Output") || mkdir("./Output")

ncfile = "./case54.nc"

# Read data from netCDF file
ds = NCDataset(ncfile);

thickness = (ds["at_layer:thickness"][:, :])';
sed_dep_thickness = (ds["at_node:sediment_deposit__thickness"][:, :])';
bedrock_elev = (ds["at_node:bedrock_surface__elevation"][:, :])';
x_of_shore = (ds["at_grid:x_of_shore"][:, :])';
x_of_node = (ds["x_of_node"][:])';
x_of_cell = (ds["x_of_cell"][:])';
sediment_load = (ds["at_grid:sediment_load"][:, :])';
sea_level = (ds["at_grid:sea_level__elevation"][:, :])';
y_of_node = (ds["y_of_node"][:])';
t0 = (ds["at_layer:t0"][:, :])';
perc_sand = (ds["at_layer:percent_sand"][:, :])';
porosity = (ds["at_layer:porosity"][:, :])';
wd = (ds["at_layer:water_depth"][:, :])';
shelf = ds["at_grid:x_of_shelf_edge"][:]
xc = ds["x_of_cell"][:]

# Convert age to ky (thousands of years)
age = (ds["at_layer:age"][:, :]) ./ 1000.0
model_timestep = round(age[1, 2] - age[1, 1])
width, time = size(age)

close(ds)

# Automation of input
sz_sl = size(sea_level);
sz_thickness = size(thickness);
sz_bed_elev = size(bedrock_elev);

tot_tstep = sz_sl[1];  # total timesteps extracted from Sea level
seq_nodes = sz_bed_elev[2];  # number of nodes from bed elevation 

bedrock_elev = bedrock_elev[:, 2:end-1];

topsurface = zeros(tot_tstep, length(x_of_cell));
bedrock_position = zeros(tot_tstep, length(x_of_cell));
sed_surf = zeros(tot_tstep, length(x_of_cell));


for i in 1:tot_tstep
    sl = sea_level[i]
    topsurface[i, :] = -1 .* wd[i, :] .+ sl
    bedrock_position[i, :] = bedrock_elev[i, :]
    sed_surf[i, :] = bedrock_position[i, :] .+ thickness[1, :]
end

# Define dtime for time boundaries
dtime = 50  # Creating time boundaries (black lines)
deposit_thickness = thickness'
t0 = t0'
wd = wd'
shore = x_of_shore'
sl = sea_level
perc_sand = perc_sand'
porosity = porosity'


Breath = (x_of_cell[2] - x_of_cell[1]) * (seq_nodes);
cols = seq_nodes - 2;
Delx = Breath / (cols - 1);
xcols = LinRange(0, Breath, cols);

AR = 100;
Dely = Delx / AR;

delt = 1000 * 365; #1000 Years

consea = 1;  #Salt concentration in sea 
rhorel = 1.025;#Relative density of saturated saline rhosat/rhowater
eps1 = 0.35; #Aquifer porosity
stov = 0.001; #0.001; #Aquifer storage


#Hydraulic conductivity m/day, othotropic
kxval = 0.0086 / 1000; #For Slope deposits NS
kyval = kxval / 10;

ksf = 8.6;
kg = ksf * 0.8;
kse = kxval * 2;

#diffusion and dispersion coefficients to be used in anisotropic dispersion
Dmol = 0.00003;#  molecular diffusion  m^2/day
aL = 50; # Longitudinal (m) 
aT = aL / 10; #Transverse 

## SEQUENCE IMPORT
totstep = tot_tstep;
sealevelin = sea_level[1, 1];#initial sealevel from sequence

new_columns = 200
xq = LinRange(0, Breath, new_columns)

eta = zeros(totstep, new_columns)
etabot = zeros(totstep, new_columns)
etabed = zeros(totstep, new_columns)

for i in 1:totstep
    interp_top = LinearInterpolation(xcols, topsurface[i, :], extrapolation_bc=Interpolations.Line())
    interp_bot = LinearInterpolation(xcols, bedrock_position[i, :], extrapolation_bc=Interpolations.Line())
    interp_bed = LinearInterpolation(xcols, sed_surf[i, :], extrapolation_bc=Interpolations.Line())

    eta[i, :] = interp_top.(xq)
    etabot[i, :] = interp_bot.(xq)
    etabed[i, :] = interp_bed.(xq)
end

cols = length(xq)
Delx = Breath / (cols - 1)
Dely = Delx / AR  # Assuming AR is already defined
xcols = xq;
Dely = 2;  # Overwrites the previous calculation of Dely


# Initial domain grid
ytop = eta[1, :]  # Initial elevation of top domain nodes
ybot = etabed[1, :] .- 4 * Dely  # Initial elevation of bottom domain nodes

# -2 * Dely creates a 2 Dely node basal thickness
if all(eta[1, :] .< eta[2, :])
    ytop = eta[1, :] .+ mean(eta[2, :] .- eta[1, :])
end

ncol = Int.(round.((ytop .- ybot) / Dely) .+ 1);  # Number of nodes in each column
N = Int(sum(round.((ytop .- ybot) / Dely) .+ 1));  # Number of nodes in domain

# Mutable struct
mutable struct DomainVariables
    phi::Vector{Float64}
    con::Vector{Float64}
    age_::Vector{Float64}
    # x::Array{Float64, 1}
    # y::Array{Float64, 1}
end



# Storage and initial values for head values
init_phi = maximum(ytop) * ones(N);
# phi_insert = zeros(1, cols)

# # Storage and initial values for concentration
# con_insert = zeros(1, cols)
# age_insert = zeros(1, cols)

# Set salt initial concentration of points vertically below sea level to value 1
init_con = Float64[]
init_age_ = Float64[]

domain_variables = DomainVariables(init_phi, init_con, init_age_);

for jj in 1:cols
    if ytop[jj] > sealevelin
        append!(domain_variables.con, zeros(ncol[jj]))
        append!(domain_variables.age_, zeros(ncol[jj]))
    else
        append!(domain_variables.con, ones(ncol[jj]))
        append!(domain_variables.age_, zeros(ncol[jj]))
    end
end


tolh = 1e-5;  # Convergence tolerance for head
tols = 1e-7; # For concentration
tola = 1e-4;  # For Age

printval = 1;  # Parameters to control plotting through time
printspace = 1;

# Printing timesteps
pt = [1, 100, 200, 300, 400, 500, 600, totstep];  # Assuming totstep is already defined

function run_groundwater_model(tot_tstep, sea_level, eta, etabed, age, ncol, sl, dtime, width, deposit_thickness, xc, shore, t0, wd, xq, domain_variables::DomainVariables)
    for tstep in 1:tot_tstep
        println("Running timestep $tstep")

        ## Adjustment of strat and set sea-level
        ytopold = copy(ytop)
        ybotold = copy(ybot)

        if tstep < totstep + 1
            global ytop = eta[tstep, :]
            global ybot = etabed[tstep, :] .- 4 * Dely
        end

        # Hardwire section to generate the spacing properly for the first timestep.
        # Not required for all of Sequence Simulations. Affect on final distribution is negligible
        if (tstep == 1) && all(eta[1, :] .< eta[2, :])
            # ytop = eta[1, :] .+ mean(eta[2, :] .- eta[1, :])
            # ytopold = ytop
            ytop = etabed[1, :]
        end

        pl_tstep = tstep + 1

        xg, yg, xse, yse, xsf, ysf, xb, yb = stratigraphy_function(pl_tstep, model_timestep, age, sl, dtime, width, deposit_thickness, xc, shore)

        # Number of nodes in each column
        ncolpre = ncol  # Store previous value
        global ncol = Int.(round.((ytop .- ybot) / Dely) .+ 1)  # Calculate new value

        Btop = Int[]  # List of nodes on top surface
        Btopsea = Int[]  # List of top domain nodes under sea
        x = Float64[]  # List of x locations of domain nodes
        y = Float64[] # List of y locations of domain nodes


        for jj in 1:cols
            append!(x, fill(Delx * (jj - 1), ncol[jj]))
            y_vals = LinRange(ybot[jj], ybot[jj] + (ncol[jj] - 2) * Dely, ncol[jj] - 1)
            append!(y, y_vals)
            push!(y, ytop[jj])

            # Find Btop, Btopsea: Nodes on top surface, top surface nodes under sea
            push!(Btop, length(y))

            if ytop[jj] < sea_level[tstep, 1]  # Assuming sea_level is indexed by tstep
                push!(Btopsea, length(y))
            end
        end

        Bbot = [1; Btop[1:end-1] .+ 1]  # List of bottom nodes

        N = length(x)  # Updated number of node points

        # Putting in the d_eta term by NS
        del_eta = (ytop .- ybot) .- (ytopold .- ybotold)  # Changes in elevation for each column at each timestep
        d_eta = Float64[]
        for i in 1:cols
            for j in 1:ncol[i]
                push!(d_eta, del_eta[i])
            end
        end

        phinew = domain_variables.phi
        connew = domain_variables.con
        agenew = domain_variables.age_

        phipre = domain_variables.phi
        conpre = domain_variables.con
        agepre = domain_variables.age_

        domain_variables.con = Float64[]
        domain_variables.phi = Float64[]
        domain_variables.age_ = Float64[]
        last = 0


        # Loop through columns
        for jj in 1:cols
            first = last + 1
            last = first + ncolpre[jj] - 1
            # phi, con and age if the nodes don't change
            if ncolpre[jj] == ncol[jj]
                append!(domain_variables.phi, phipre[first:last])
                append!(domain_variables.con, conpre[first:last])
                append!(domain_variables.age_, agepre[first:last])
            elseif ncolpre[jj] > ncol[jj]
                ndiff = ncolpre[jj] - ncol[jj]
                append!(domain_variables.phi, phipre[first:last-ndiff])
                append!(domain_variables.con, conpre[first:last-ndiff])
                append!(domain_variables.age_, agepre[first:last-ndiff])
            else
                ndiff = ncol[jj] - ncolpre[jj]
                append!(domain_variables.phi, phipre[first:last-1])
                append!(domain_variables.con, conpre[first:last-1])
                append!(domain_variables.age_, agepre[first:last-1])


                mnode = size(domain_variables.con, 1) # gives top node position in current column
                dy = y[mnode+ndiff+1] - y[mnode] # total change in y

                # interpolation ratio
                phirat = (phipre[last] - phipre[last-1]) / dy
                conrat = (conpre[last] - conpre[last-1]) / dy
                agerat = (agepre[last] - agepre[last-1]) / dy

                # Instead of For loop using append and fill 
                append!(domain_variables.phi, phipre[last-1] .+ phirat .* (Dely .* (1:ndiff)))
                append!(domain_variables.con, conpre[last-1] .+ conrat .* (Dely .* (1:ndiff)))
                append!(domain_variables.age_, agepre[last-1] .+ agerat .* (Dely .* (1:ndiff)))

                append!(domain_variables.phi, phinew[last])
                append!(domain_variables.con, connew[last])
                append!(domain_variables.age_, agenew[last])

            end
        end

        phinew = domain_variables.phi
        connew = domain_variables.con
        agenew = domain_variables.age_

        t = generate_tri(x, y, ncol, cols)

        Ntri = size(t, 1)  # The size of t is the number of triangle elements

        maxconnect = 0
        for ii in 1:N
            maxconnect = max(maxconnect, sum(t .== ii))
            # NOTE: sum(t .== ii) finds the number of elements connected to node ii
        end
        maxsup = 2 * maxconnect  # Double the nodes in support

        Volp = zeros(N)    # CV volume
        xmid = zeros(Ntri) # Element midpoint x-coordinates
        ymid = zeros(Ntri) # Element midpoint y-coordinates
        Nx = zeros(Ntri, 3)   # Derivatives of shape functions (x-component)
        Ny = zeros(Ntri, 3)   # Derivatives of shape functions (y-component)


        for itri in 1:Ntri
            k1 = t[itri, 1]  # Global number of 1st node in triangle itri
            k2 = t[itri, 2]  # 2nd node
            k3 = t[itri, 3]  # 3rd node

            # Element volume
            v = (x[k2] * y[k3] - x[k3] * y[k2] - x[k1] * y[k3] + x[k1] * y[k2] +
                 y[k1] * x[k3] - y[k1] * x[k2]) / 2

            # Contribution to control volume
            Volp[k1] += v / 3
            Volp[k2] += v / 3
            Volp[k3] += v / 3

            # Midpoint of element
            xmid[itri] = (x[k1] + x[k2] + x[k3]) / 3
            ymid[itri] = (y[k1] + y[k2] + y[k3]) / 3

            # Derivatives of shape functions
            Nx[itri, 1] = (y[k2] - y[k3]) / (2 * v)
            Nx[itri, 2] = (y[k3] - y[k1]) / (2 * v)
            Nx[itri, 3] = (y[k1] - y[k2]) / (2 * v)
            Ny[itri, 1] = -(x[k2] - x[k3]) / (2 * v)
            Ny[itri, 2] = -(x[k3] - x[k1]) / (2 * v)
            Ny[itri, 3] = -(x[k1] - x[k2]) / (2 * v)

        end

        points = vec(SVector.(xmid, ymid))
        npoints = vec(SVector.(x, y))

        xqq = []
        yqq = []
        append!(xqq, xq, xq[cols:-1:1], xq[1])
        append!(yqq, ybot .+ 2 * Dely, ybot[cols:-1:1], ybot[1] .+ 2 * Dely)
        push!(xg, xg[1])
        push!(yg, yg[1])
        push!(xse, xse[1])
        push!(yse, yse[1])
        push!(xsf, xsf[1])
        push!(ysf, ysf[1])
        # push!(xb,xb[1]); push!(yb,yb[1])
        if minimum(xg) > 0
            xg[xg .== minimum(xg)] .= 0
        end

        if maximum(yg) < maximum(y)
            yg[yg .== maximum(yg)] .= maximum(y)
        end



        polygon_g = SVector.(xg * 1000, yg)
        polygon_se = SVector.(xse * 1000, yse)
        polygon_sf = SVector.(xsf * 1000, ysf)
        polygon_bed = SVector.(xqq, yqq)

        in_mat_g = [inpolygon(p, polygon_g; in=true, on=false, out=false) for p in points]
        in_mat_se = [inpolygon(p, polygon_se; in=true, on=false, out=false) for p in points]
        in_mat_sf = [inpolygon(p, polygon_sf; in=true, on=false, out=false) for p in points]
        in_mat_bed = [inpolygon(p, polygon_bed; in=true, on=false, out=false) for p in points]
        ## Node In Mat     
        in_mat_ng = [inpolygon(p, polygon_g; in=true, on=false, out=false) for p in npoints]
        in_mat_nse = [inpolygon(p, polygon_se; in=true, on=false, out=false) for p in npoints]
        in_mat_nsf = [inpolygon(p, polygon_sf; in=true, on=false, out=false) for p in npoints]

        ## POROSITY Initialization by NS
        epsk = 0.5 * ones(size(x))

        eps_0s = 0.4
        bs = 1e-8
        eps_0c = 0.7
        bc = 1e-7

        k = 1 # current node
        l = 0 # last node of current column

        for i in 1:cols
            l = l + ncol[i]
            for j in 1:ncol[i]
                Z = y[l] - y[k]
                if in_mat_ng[k] == 1 || in_mat_nsf[k] == 1
                    if (Z * 2.3e4 - (domain_variables.phi[k] - y[k]) * 1e4) < 0
                        epsk[k] = eps_0s
                    else
                        epsk[k] = eps_0s * exp(-bs * (Z * 2.3e4 - (domain_variables.phi[k] - y[k]) * 1e4))
                    end
                else
                    if (Z * 2.3e4 - (domain_variables.phi[k] - y[k]) * 1e4) < 0
                        epsk[k] = eps_0c
                    else
                        epsk[k] = eps_0c * exp(-bc * (Z * 2.3e4 - (domain_variables.phi[k] - y[k]) * 1e4))
                    end
                end
                k = k + 1
            end
        end
        eps1 = Float64[]
        eps1 = epsk


        Big = 1e18
        BCh = zeros(N)  # Boundary coefficient values for head
        BCs = zeros(N)  # Boundary coefficient values for solute
        BBh = zeros(N)  # Fixed head values
        BBs = zeros(N)  # Fixed solute concentration values

        # Fixed head value on top boundary
        BCh[Btop] .= Big
        BBh[Btop] .= Big .* y[Btop]

        # Correction for sea water
        BBh[Btopsea] .= Big .* (y[Btopsea] .- (y[Btopsea] .- sea_level[tstep, 1]) .* rhorel)

        # Concentration fixed values: 0 above sea level, 1 below
        BCs[Btop] .= Big
        BBs[Btopsea] .= consea .* Big

        # Age boundary conditions
        BBa = zeros(N)
        BCa = zeros(N)
        BBa[Btopsea] .= Big
        BCa[Btop] .= Big
        BBa[Btopsea] .= Big

        BBvar = zeros(N)     # Variable density source
        kx = ones(Ntri) * kxval     # x-direction conductivity value
        ky = ones(Ntri) * kyval     # y-direction conductivity value

        AA = spzeros(N, N)
        BB = zeros(N, 1)


        for itri in 1:Ntri
            kx[itri] = kxval
            ky[itri] = kyval

            if xmid[itri] < shore[tstep] && ymid[itri] > sea_level[tstep]
                kx[itri] = kg
                ky[itri] = kg / 100
            end

            if in_mat_se[itri] == 1
                kx[itri] = kse
                ky[itri] = kse / 10
            elseif in_mat_sf[itri] == 1
                kx[itri] = ksf
                ky[itri] = ksf / 100
            elseif in_mat_g[itri] == 1
                kx[itri] = kg
                ky[itri] = kg / 100
            end


            # Node cycling order
            cyc = [1 2 3; 2 3 1; 3 1 2]

            for node in 1:3
                ii, jj, kk = cyc[node, :]


                k1 = t[itri, ii]  # Global node number
                k2 = t[itri, jj]
                k3 = t[itri, kk]

                Nx1, Nx2, Nx3 = Nx[itri, ii], Nx[itri, jj], Nx[itri, kk]
                Ny1, Ny2, Ny3 = Ny[itri, ii], Ny[itri, jj], Ny[itri, kk]

                # Face 1
                delx = (x[k1] + x[k2] + x[k3]) / 3 - (x[k1] + x[k2]) / 2
                dely = (y[k1] + y[k2] + y[k3]) / 3 - (y[k1] + y[k2]) / 2

                face1_k1 = kx[itri] * Nx1 * dely - ky[itri] * Ny1 * delx
                face1_k2 = kx[itri] * Nx2 * dely - ky[itri] * Ny2 * delx
                face1_k3 = kx[itri] * Nx3 * dely - ky[itri] * Ny3 * delx

                # Variable density source (face value)
                BBvar[k1] -= ky[itri] * ((rhorel - 1) / 12) * (5 * domain_variables.con[k1] + 5 * domain_variables.con[k2] + 2 * domain_variables.con[k3]) * delx

                # Face 2
                delx = -(x[k1] + x[k2] + x[k3]) / 3 + (x[k1] + x[k3]) / 2
                dely = -(y[k1] + y[k2] + y[k3]) / 3 + (y[k1] + y[k3]) / 2

                face2_k1 = kx[itri] * Nx1 * dely - ky[itri] * Ny1 * delx
                face2_k2 = kx[itri] * Nx2 * dely - ky[itri] * Ny2 * delx
                face2_k3 = kx[itri] * Nx3 * dely - ky[itri] * Ny3 * delx

                # Ax=b Solver
                AA[k1, k1] = AA[k1, k1] + delt * (face1_k1 + face2_k1)
                AA[k1, k2] = AA[k1, k2] + delt * (face1_k2 + face2_k2)
                AA[k1, k3] = AA[k1, k3] + delt * (face1_k3 + face2_k3)


                # Variable density source (face value)
                BBvar[k1] -= ky[itri] * ((rhorel - 1) / 12) * (5 * domain_variables.con[k1] + 2 * domain_variables.con[k2] + 5 * domain_variables.con[k3]) * delx
            end
        end

        # Storage set to 0 at first time step
        sto = stov
        if tstep == 1
            sto = 0
        end

        Badd = -BCh .- sto .* Volp
        AA = spdiagm(0 => Badd) + AA
        BB = -sto .* Volp .* domain_variables.phi - delt .* BBvar .- BBh - d_eta .* Volp * 1.3 * sto

        # Solve for domain_variables.phi
        phinew = AA \ BB

        conver = 1
        phipre = phinew


        # Calculate volume of flow stored per time step
        Vstore = sto .* (Volp .* (phinew - domain_variables.phi))

        domain_variables.phi = phinew


        AA = spzeros(N, N)
        BB = zeros(N, 1)

        qx = zeros(Ntri)
        qy = zeros(Ntri)

        qxval = []
        qyval = []

        Dxx = 0.0
        Dyy = 0.0
        Dxy = 0.0

        Netflux = zeros(N)

        eps2 = 0.35

        for itri in 1:Ntri
            cyc = [1 2 3; 2 3 1; 3 1 2]  # Node cycling order

            # Diffusion
            for node in 1:3
                ii, jj, kk = cyc[node, :]

                k1 = t[itri, ii]  # Global node number of element vertices
                k2 = t[itri, jj]
                k3 = t[itri, kk]

                Nx1, Nx2, Nx3 = Nx[itri, ii], Nx[itri, jj], Nx[itri, kk]
                Ny1, Ny2, Ny3 = Ny[itri, ii], Ny[itri, jj], Ny[itri, kk]

                # Dispersion for element
                if node == 1
                    # Contribution to discharge from fresh-water head
                    qxval = -kx[itri] * (Nx1 * domain_variables.phi[k1] + Nx2 * domain_variables.phi[k2] + Nx3 * domain_variables.phi[k3])
                    qyval = -ky[itri] * (Ny1 * domain_variables.phi[k1] + Ny2 * domain_variables.phi[k2] + Ny3 * domain_variables.phi[k3])

                    # Actual discharges at element midpoints
                    qxmid = qxval
                    qymid = qyval - ky[itri] * ((rhorel - 1) / 3) * (domain_variables.con[k1] + domain_variables.con[k2] + domain_variables.con[k3])

                    qx[itri] = qxmid
                    qy[itri] = qymid

                    qx2 = qxmid^2
                    qy2 = qymid^2
                    qabs = sqrt(qx2 + qy2)
                    eps2 = mean(eps1[[k1, k2, k3]])

                    # Dispersion
                    Dxx = aT * qabs + (aL - aT) * qx2 / qabs + Dmol * eps2
                    Dyy = aT * qabs + (aL - aT) * qy2 / qabs + Dmol * eps2
                    Dxy = (aL - aT) * qxmid * qymid / qabs
                end

                # Face 1

                qxface = qxval
                qyface = qyval - ky[itri] * ((rhorel - 1) / 12) * (5 * domain_variables.con[k1] + 5 * domain_variables.con[k2] + 2 * domain_variables.con[k3])
                delx = (x[k1] + x[k2] + x[k3]) / 3 - (x[k1] + x[k2]) / 2
                dely = (y[k1] + y[k2] + y[k3]) / 3 - (y[k1] + y[k2]) / 2

                # Contribution to flux across face due to dispersion
                face1_k1 = (Dxx * Nx1 + Dxy * Ny1) * dely - (Dyy * Ny1 + Dxy * Nx1) * delx
                face1_k2 = (Dxx * Nx2 + Dxy * Ny2) * dely - (Dyy * Ny2 + Dxy * Nx2) * delx
                face1_k3 = (Dxx * Nx3 + Dxy * Ny3) * dely - (Dyy * Ny3 + Dxy * Nx3) * delx

                # Contribution to flux across face due to flow (upwind scheme)
                qout = qxface * dely - qyface * delx  # Flow out of volume k
                if qout >= 0
                    face1_k1 -= qout
                else
                    face1_k2 -= qout
                end

                Netflux[k1] += qout

                qxface = qxval
                qyface = qyval - ky[itri] * ((rhorel - 1) / 12) * (5 * domain_variables.con[k1] + 2 * domain_variables.con[k2] + 5 * domain_variables.con[k3])

                delx = -(x[k1] + x[k2] + x[k3]) / 3 + (x[k1] + x[k3]) / 2
                dely = -(y[k1] + y[k2] + y[k3]) / 3 + (y[k1] + y[k3]) / 2

                # Contribution to flux across face due to dispersion
                face2_k1 = (Dxx * Nx1 + Dxy * Ny1) * dely - (Dyy * Ny1 + Dxy * Nx1) * delx
                face2_k2 = (Dxx * Nx2 + Dxy * Ny2) * dely - (Dyy * Ny2 + Dxy * Nx2) * delx
                face2_k3 = (Dxx * Nx3 + Dxy * Ny3) * dely - (Dyy * Ny3 + Dxy * Nx3) * delx

                # Contribution to flux across face due to flow (upwind scheme)
                qout = qxface * dely - qyface * delx  # Flow out of volume k
                if qout >= 0
                    face2_k1 -= qout
                else
                    face2_k3 -= qout
                end

                Netflux[k1] += qout

                # Ax=b Solver
                AA[k1, k1] = AA[k1, k1] + delt * (face1_k1 + face2_k1)
                AA[k1, k2] = AA[k1, k2] + delt * (face1_k2 + face2_k2)
                AA[k1, k3] = AA[k1, k3] + delt * (face1_k3 + face2_k3)

            end
        end


        BBs[Btopsea] .= Big .* consea

        # Alternative boundary condition for solute
        # nsub=size(Btopsea,2)
        # for ii=1:nsub
        #     bnode=Btopsea[ii]
        #     if domain_variables.phi[bnode-1,1] - domain_variables.phi[bnode,1] > 0
        #         BBs[bnode]=Big*domain_variables.con[bnode-1]
        #     end
        # end

        Badd = -BCs .- eps1 .* Volp .- Vstore
        AAcon = spdiagm(0 => Badd) + AA
        BBcon = -eps1 .* Volp .* domain_variables.con - BBs - d_eta .* Volp * 1.3 * sto .* domain_variables.con

        connew = AAcon \ BBcon #sparse solver
        # %connew=min(connew,1); %XXXX mass Limiter If needed
        domain_variables.con = connew

        for i = 1:N
            if domain_variables.con[i] > 1
                domain_variables.con[i] = 1
            end
        end


        BCa[Btop, 1] .= Big
        BBa[Btopsea, 1] .= 0

        Badd = -BCa .- eps1 .* Volp .- Vstore
        AAage = spdiagm(0 => Badd) + AA
        BBage = -eps1 .* Volp .* domain_variables.age_ - BBa - delt .* Volp .* eps1 - d_eta .* Volp .* 1.3 .* sto .* domain_variables.age_

        agenew = AAage \ BBage
        domain_variables.age_ = agenew

        # if any(tstep .== pt)
        #     filename_x = "./Output/case15_$(tstep).jld2"
        #     @save filename_x x y domain_variables.phi domain_variables.con domain_variables.age_ t ytop ybot xg yg xse yse xsf ysf xb yb tstep sea_level Breath cols tot_tstep xq kx ky
        # end

        xq = collect(xq)
        t = Matrix(t)

        if tstep > 0
            matfile = matopen("./Output/case54t1_$(tstep).mat", "w")

            write(matfile, "x", x)
            write(matfile, "y", y)
            write(matfile, "t", t)
            write(matfile, "con", domain_variables.con)
            write(matfile, "phi", domain_variables.phi)
            write(matfile, "age_", domain_variables.age_)
            write(matfile, "xq", xq)
            write(matfile, "ytop", ytop)
            write(matfile, "ybot", ybot)
            write(matfile, "kx", kx)
            write(matfile, "tstep", tstep)
            write(matfile, "sea_level", sea_level)
            write(matfile, "cols", cols)
            write(matfile, "tot_tstep", tot_tstep)
            write(matfile, "ky", ky)
            write(matfile, "qx", qx)
            write(matfile, "qy", qy)
            write(matfile, "eps1", eps1)

            close(matfile)
        end



    end
end

@time run_groundwater_model(
    500,
    sea_level,
    eta,
    etabed,
    age,
    ncol,
    sl,
    dtime,
    width,
    deposit_thickness,
    xc,
    shore,
    t0,
    wd,
    xq,
    domain_variables
)


