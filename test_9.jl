#import-------------------------
include("main.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-14
t_final = 10e-11
n_steps = round(Int, t_final / delta_t)

domain = Domain(4e-8, 4e-8, 8e-8)

animation_time = 15  # seconds of animation

n_He = 400
v_He = 789.45

n_Ar = 200
v_Ar = 249.88

n_bins = 32 #number of bins for the histograms

#auto-bins-defs fopr histograms------------------------------------------
z_bin_edges = range(-domain.z/2, domain.z/2, length = n_bins + 1)
z_bin_centers = [(z_bin_edges[i] + z_bin_edges[i+1]) / 2 for i in 1:n_bins]
z_bin_width   = domain.z / n_bins

#molecules initialisation helper ------------------------------------
function populate(molecule,n_molecules,initial_v)
    molecules = Vector{Molecule}(undef, n_molecules)
    for i in 1:n_molecules
        mol = copy(molecule)
        mol.position = (
            (rand() - 0.5) * domain.x,
            (rand() - 0.5) * domain.y,
            (rand() - 0.5) * domain.z
        )
        angle_vector_nonnormalised = (rand() - 0.5, rand() - 0.5, rand() - 0.5)
        angle_vector = angle_vector_nonnormalised ./ norm(angle_vector_nonnormalised)
        mol.speed = (
            angle_vector[1] * initial_v,
            angle_vector[2] * initial_v,
            angle_vector[3] * initial_v
        )
        molecules[i] = mol
    end
    return molecules
end

function get_color(m) 
    RGBf(
        (hash(m.formule_chimique * "redder") % 256 / 255),
        (hash(m.formule_chimique * "grent") % 256 / 255),
        (hash(m.formule_chimique * "bluey") % 256 / 255)
    )
end
#molecules initialisation ------------------------------------
Ar_molecules = populate(Ar,n_Ar,v_Ar)
He_molecules = populate(He,n_He,v_He)
molecules = vcat(Ar_molecules, He_molecules)

colors = [get_color(m) for m in molecules]

#simulation loop----------------------------------------------------
decouple = round(Int, n_steps / animation_time / 60)
states_history = [deepcopy(molecules)]
println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t, domain = domain)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
    end
    if i == round(n_steps / 2)
        print("50% done...\n")
    end
end
println("simulation done!")
#histogram helpers--------------------------------------
function get_bins(molecules)
    bins = [Molecule[] for _ in 1:n_bins]
    for mol in molecules
        z = mol.position[3]
        b = clamp(Int(floor((z + domain.z/2) / z_bin_width)) + 1, 1, n_bins)
        push!(bins[b], mol)
    end
    return bins
end    

function pressure_profile(molecules)
    bins = get_bins(molecules)

    pressures = zeros(Float64, n_bins)

    for i in 1:n_bins
        if !isempty(bins[i])
            pressures[i] = pression(bins[i], Domain(domain.x, domain.y, z_bin_width))
        end
    end

    return pressures
end

function cinetic_e_profile(molecules)
    bins = get_bins(molecules)

    cinetic_es = zeros(Float64, n_bins)

    for i in 1:n_bins
        if !isempty(bins[i])
            cinetic_es[i] = cineticEnergy(bins[i])
        end
    end

    return cinetic_es
end

function temperature_profile(molecules)
    bins = get_bins(molecules)

    temperatures = zeros(Float64, n_bins)

    for i in 1:n_bins
        if !isempty(bins[i])
            temperatures[i] = temperature(bins[i])
        end
    end

    return temperatures
end


#fig layout--------------------------------------
fig = Figure(size = (1800, 1200))

# Left: 3D scatter
ax3d = Axis3(fig[1:2, 1:2],
    title = "Molecule positions",
    limits = (
        -domain.x/2, domain.x/2,
        -domain.y/2, domain.y/2,
        -domain.z/2, domain.z/2
    )
)

# Right: z-position distribution
ax_hist = Axis(fig[2, 3],
    title = "Z-position distribution",
    xlabel = "Probability",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

#Down: histogrtams
ax_pres = Axis(fig[3, 1],
    title = "Z-pressure distribution",
    xlabel = "Pressure",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

ax_cin = Axis(fig[3, 2],
    title = "Z-cinetic energy distribution",
    xlabel = "cinetic energy",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

ax_temp = Axis(fig[3, 3],
    title = "Z-temperature distribution",
    xlabel = "temperature",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

#observables-----------------------------------------
obs_positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
obs_z_He = Observable([m.position[3] for m in He_molecules])
obs_z_Ar = Observable([m.position[3] for m in Ar_molecules])
obs_pres = Observable(pressure_profile(molecules))
obs_cin = Observable(cinetic_e_profile(molecules))
obs_temp = Observable(temperature_profile(molecules))

#plots-----------------
scatter!(ax3d, obs_positions, markersize = 8, color = colors)

hist!(ax_hist, obs_z_He,
    bins = z_bin_edges,
    direction = :x,
    normalization = :probability,
    color = (get_color(He), 0.5),
    label = "He"
)
hist!(ax_hist, obs_z_Ar,
    bins = z_bin_edges,
    direction = :x,
    normalization = :probability,
    color = (get_color(Ar), 0.5),
    label = "Ar"
)

barplot!(ax_pres,
    z_bin_centers,
    obs_pres,
    direction = :x,
    width     = z_bin_width,
)

barplot!(ax_cin,
    z_bin_centers,
    obs_cin,
    direction = :x,
    width     = z_bin_width,
)
barplot!(ax_temp,
    z_bin_centers,
    obs_temp,
    direction = :x,
    width     = z_bin_width,
)

#animation loop ----------------------------------------------------------------------------
running = Observable(true)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
        running[] = false
    end
end

screen = display(fig)
while running[]
    for state in states_history
        running[] || break
        obs_positions[] = [Point3f(m.position[1], m.position[2], m.position[3]) for m in state]
        obs_z_He[] = [m.position[3] for m in state if m.formule_chimique == "He"]
        obs_z_Ar[] = [m.position[3] for m in state if m.formule_chimique == "Ar"]
        obs_pres[] = pressure_profile(state)
        obs_cin[] = cinetic_e_profile(state)
        obs_temp[] = temperature_profile(state)
        sleep(1/60)  # 60 fps
    end
end
close(screen)