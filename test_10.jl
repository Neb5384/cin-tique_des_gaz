#import-------------------------
include("main.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-14
t_final = 5e-11
n_steps = round(Int, t_final / delta_t)

domain = Domain(1e-8, 1e-8, 1e-8)

animation_time = 8  # seconds of animation

n_He = 400
v_He = 1400

#molecules initialisation helper ------------------------------------
function populate(molecule, n_molecules, initial_v)
    molecules = Vector{Molecule}(undef, n_molecules)
    for i in 1:n_molecules
        mol = copy(molecule)
        mol.position = (
            - (rand()) * domain.x/2,
            - (rand()) * domain.y/2,
            - (rand()) * domain.z/2
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

println("expanding walls")
domain = Domain(2e-8, 1e-8, 1e-8)

@time for i in n_steps: n_steps + 5000
    step(molecules, delta_t = delta_t, domain = domain)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
        push!(entropy_history, entropy(molecules, domain,(20,10,10),200))
    end
    if i == round(n_steps / 2)
        print("50% done...\n")
    end
end

println("simulation done!")

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

#observables-----------------------------------------
obs_positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
obs_entropy = Observable(entropy_history)

#plots-----------------
scatter!(ax3d, obs_positions, markersize = 8, color = colors)
lines!(ax_entropy, obs_entropy)

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
        sleep(1/60)
    end
end
close(screen)