#import-------------------------
include("main.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-14
t_final = 1e-10
n_steps = round(Int, t_final / delta_t)

domain = Domain(4e-9, 4e-9, 1e-8)

animation_time = 30  # seconds of animation

n_He = 500
v_He = 1400

n_bins = 16 #number of bins for the histograms

#auto-bins-defs for histograms------------------------------------------
z_bin_edges = range(-domain.z/2, domain.z/2, length = n_bins + 1)
z_bin_centers = [(z_bin_edges[i] + z_bin_edges[i+1]) / 2 for i in 1:n_bins]
z_bin_width   = domain.z / n_bins

#molecules initialisation helper ------------------------------------
function populate(molecule, n_molecules, initial_v)
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
molecules = populate(He, n_He, v_He)
colors = [get_color(m) for m in molecules]

#simulation loop----------------------------------------------------
decouple = round(Int, n_steps / animation_time / 60)
states_history = [deepcopy(molecules)]
entropy_history = [entropy(molecules, domain,(10,10,10),200)]
println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t, domain = domain)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
        push!(entropy_history, entropy(molecules, domain,(10,10,10),200))
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
fig = Figure(size = (1800, 1800))

ax3d = Axis3(fig[1:2, 1:2],
    title = "Molecule positions",
    limits = (
        -domain.x/2, domain.x/2,
        -domain.y/2, domain.y/2,
        -domain.z/2, domain.z/2
    ),
    aspect = :data
)

ax_entropy = Axis(fig[1, 3],
    title = "Entropy over time",
    xlabel = "Frame",
    ylabel = "Entropy (bits)"
)

ax_temp = Axis(fig[2, 3],
    title = "Z-temperature distribution",
    xlabel = "temperature",
    ylabel = "z (m)",
    limits = (0, 1000, -domain.z/2, domain.z/2)
)

#observables-----------------------------------------
obs_positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
obs_entropy = Observable(entropy_history)
obs_temp = Observable(temperature_profile(molecules))

#plots-----------------
scatter!(ax3d, obs_positions, markersize = 8, color = colors)
lines!(ax_entropy, obs_entropy)
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
    for (frame, state) in enumerate(states_history)
        running[] || break
        obs_positions[] = [Point3f(m.position[1], m.position[2], m.position[3]) for m in state]
        obs_entropy[] = entropy_history[1:frame]
        obs_temp[] = temperature_profile(state)
        sleep(1/60)
    end
end
close(screen)