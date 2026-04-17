#import-------------------------
include("main.jl")
using GLMakie
using LinearAlgebra

#base variables-----------------------------
delta_t = 1e-14
t_final = 5e-11
n_steps = round(Int, t_final / delta_t)
domain = Domain(4e-8, 4e-8, 4e-8)

animation_time = 15  # seconds of animation
n_molecules = 500
initial_atom_speed = 1367

n_bins = 64 #number of bins for the histograms

#molecules initialisation------------------------------------
molecules = Vector{Molecule}(undef, n_molecules)
for i in 1:n_molecules
    mol = copy(He)
    mol.position = (
        (rand() - 0.5) * domain.x,
        (rand() - 0.5) * domain.y,
        (rand() - 0.5) * domain.z
    )
    angle_vector_nonnormalised = (rand() - 0.5, rand() - 0.5, rand() - 0.5)
    angle_vector = angle_vector_nonnormalised ./ norm(angle_vector_nonnormalised)
    mol.speed = (
        angle_vector[1] * initial_atom_speed,
        angle_vector[2] * initial_atom_speed,
        angle_vector[3] * initial_atom_speed
    )
    molecules[i] = mol
end

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

#pressure vs z helper--------------------------------------------
z_bin_edges = range(-domain.z/2, domain.z/2, length = n_bins + 1)
z_bin_centers = [(z_bin_edges[i] + z_bin_edges[i+1]) / 2 for i in 1:n_bins]
z_bin_width   = domain.z / n_bins


function pressure_profile(molecules)
    bins = [Molecule[] for _ in 1:n_bins]
    for mol in molecules
        z = mol.position[3]
        b = clamp(Int(floor((z + domain.z/2) / z_bin_width)) + 1, 1, n_bins)
        push!(bins[b], mol)
    end

    pressures = zeros(Float64, n_bins)

    for i in 1:n_bins
        if !isempty(bins[i])
            pressures[i] = pression(bins[i], Domain(domain.x, domain.y, z_bin_width))
        end
    end

    return pressures
end

#figure layout------------------------------------------------------
fig = Figure(size = (1800, 900))

# Left: 3D scatter
ax3d = Axis3(fig[1, 1],
    title = "Molecule positions",
    limits = (
        -domain.x/2, domain.x/2,
        -domain.y/2, domain.y/2,
        -domain.z/2, domain.z/2
    )
)

# Middle: z-position distribution
ax_hist = Axis(fig[1, 2],
    title = "Z-position distribution",
    xlabel = "Probability",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

# Right: pressure vs z (manual barplot)
ax_pres = Axis(fig[1, 3],
    title  = "Local pressure vs z",
    xlabel = "Pressure (Pa)",
    ylabel = "z (m)",
    limits = (0, nothing, -domain.z/2, domain.z/2)
)

colsize!(fig.layout, 1, Auto())
colsize!(fig.layout, 2, Fixed(350))
colsize!(fig.layout, 3, Fixed(350))

#observables--------------------------------------------------------
obs_positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
obs_z         = Observable([m.position[3] for m in molecules])
obs_pres      = Observable(pressure_profile(molecules))

scatter!(ax3d, obs_positions, markersize = 8, color = :dodgerblue)

hist!(ax_hist, obs_z,
    bins          = z_bin_edges,
    direction     = :x,
    normalization = :probability,
)

barplot!(ax_pres,
    z_bin_centers,
    obs_pres,
    direction = :x,
    width     = z_bin_width,
)

#animation loop-----------------------------------------------------
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
        obs_z[]         = [m.position[3] for m in state]
        obs_pres[]      = pressure_profile(state)
        reset_limits!(ax_hist, yauto = false)
        reset_limits!(ax_pres, yauto = false)
        sleep(1/60)
    end
end
close(screen)
