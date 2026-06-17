#import-------------------------
include("main2D.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-15
t_final = 5e-11
n_steps = round(Int, t_final / delta_t)

domain = Domain(1e-8, 1e-8) 

animation_time = 5  # seconds of animation

n_Ne = 100
t_ref = 10.0
min_dist_factor = 0.9   # particules séparées d'au moins 0.9*sigma à l'initialisation

#molecules initialisation helpers ------------------------------------

# chaque composante de vitesse suit une gaussienne d'écart-type sqrt(k_B T / m)
function maxwellBoltzmannSpeed(masse::Float64, temp::Float64)
    c_boltzmann = 1.380649e-23
    sigma_v = sqrt(c_boltzmann * temp / masse)
    return (randn() * sigma_v, randn() * sigma_v)
end

# place n_molecules copies de `molecule`, vitesses tirées à la température `temp`,
# positions aléatoires dans `domain` sans chevauchement (distance mini min_dist_factor*sigma)
function populate(molecule::Molecule, n_molecules::Int, temp::Float64, domain::Domain, min_dist_factor::Float64)
    molecules = Vector{Molecule}(undef, n_molecules)
    for i in 1:n_molecules
        mol = copy(molecule)
        mol.speed = maxwellBoltzmannSpeed(mol.masse, temp)

        placed = false
        attempts = 0
        while !placed
            attempts += 1
            candidate = ((rand() - 0.5) * domain.x, (rand() - 0.5) * domain.y)
            ok = true
            for j in 1:(i-1)
                if norm(candidate .- molecules[j].position) < min_dist_factor * mol.sigma
                    ok = false
                    break
                end
            end
            if ok
                mol.position = candidate
                placed = true
            elseif attempts > 100_000
                error("Impossible de placer la molécule $i sans chevauchement après $attempts tentatives : le domaine est probablement trop petit pour $n_molecules atomes.")
            end
        end
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
molecules = populate(Ne, n_Ne, t_ref, domain, min_dist_factor)
colors = [get_color(m) for m in molecules]

#simulation loop----------------------------------------------------
decouple = max(1, round(Int, n_steps / animation_time / 60))
states_history = [deepcopy(molecules)]
temperature_history = [temperature(molecules)]
println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t, domain = domain, t_ref = t_ref)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
        push!(temperature_history, temperature(molecules))
    end
    if i == round(n_steps / 2)
        print("50% done...\n")
    end
end

println("simulation done!")

#fig layout--------------------------------------
fig = Figure(size = (1200, 600))

ax2d = Axis(fig[1, 1],
    title = "Positions des atomes (Ne)",
    limits = (-domain.x/2, domain.x/2, -domain.y/2, domain.y/2),
    aspect = DataAspect()
)

ax_temp = Axis(fig[1, 2],
    title = "Température au cours du temps",
    xlabel = "Frame",
    ylabel = "Température (K)"
)

#observables-----------------------------------------
obs_positions   = Observable([Point2f(m.position[1], m.position[2]) for m in molecules])
obs_temperature = Observable(temperature_history)

#plots-----------------
scatter!(ax2d, obs_positions, markersize = 10, color = colors)
lines!(ax_temp, obs_temperature, label = "T(t)")
hlines!(ax_temp, [t_ref], color = :red, linestyle = :dash, label = "T_ref")
axislegend(ax_temp)

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
        obs_positions[]   = [Point2f(m.position[1], m.position[2]) for m in state]
        obs_temperature[] = temperature_history[1:frame]
        sleep(1/60)
    end
end
close(screen)