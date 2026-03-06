#import----------------
include("main.jl")
using GLMakie



#base variables------------------
delta_t = 1e-13
t_final = 100e-9
n_steps = round(Int, t_final / delta_t)

domain = Domain(2e-6,2e-6,2e-6)

animation_time = 30 #secondes d'animation

n_molecules = 100

individual_color_va = 0.1



#molecules initialization---------------------------------------
base_molecules = [He, Ne, N2, O2]


thermal_speeds = Dict("He" => 1200, "Ne" => 560, "N2" => 470, "O2" => 440)

molecules = Vector{Molecule}(undef, n_molecules)
for i in 1:n_molecules
    base = rand(base_molecules)
    mol = copy(base)

    mol.position = (
        (rand() -0.5) * domain.x,
        (rand() -0.5) * domain.y,
        (rand() -0.5) * domain.z
    )

    v_max = thermal_speeds[mol.formule_chimique]
    mol.speed = (
        (rand() * 2 - 1) * v_max,
        (rand() * 2 - 1) * v_max,
        (rand() * 2 - 1) * v_max
    )

    molecules[i] = mol
end

function get_color(name, seed, variation)
    va = variation
    nva = 1- va
    RGBf(
        nva * (hash(name * "redder") % 256 / 255) + va * (hash(seed + 17) % 256 / 255),
        nva * (hash(name * "grent") % 256 / 255) + va * (hash(seed + 42) % 256 / 255),
        nva * (hash(name * "bluey") % 256 / 255) + va * (hash(seed + 5) % 256 / 255)
    )
end

molecule_colors = [get_color(m.formule_chimique, i,individual_color_va) for (i, m) in enumerate(molecules)]


#animation parameters-------------
fig = Figure()
ax = Axis3(fig[1,1],limits=(-1e-6,1e-6,-1e-6,1e-6,-1e-6,1e-6))

positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
colors = Observable(molecule_colors)
scatterplot = scatter!(ax, positions,color=colors,markersize=20)

decouple = round(Int, n_steps/animation_time/60)


#simulation loop which copies states that will be shown imn animation later--------------
states_history = [deepcopy(molecules)]
println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t,domain = domain)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
    end
end
println("simulation done !")

#testing ifcinetic enrgy have been conserved-----------------------
first_cinetic_energy = cineticEnergy(states_history[1])
last_cinetic_energy = cineticEnergy(states_history[end])

println("difference energie cinetique : ", (last_cinetic_energy -first_cinetic_energy) / first_cinetic_energy, " %" )

#animation loop ----------------------------------------------
running = Observable(true)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
        running[] = false
        close(screen)
    end
end

screen = display(fig)
while running[]
    for state in states_history
        running[] || break
        positions[] = [Point3f(m.position[1], m.position[2],m.position[3]) for m in state]
        sleep(1/60)  # 60 fps
    end
end
close(screen)


