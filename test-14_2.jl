#import-------------------------
include("main2D.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-15
t_final = 5e-11
n_steps = round(Int, t_final / delta_t)

domain = Domain(1e-8, 1e-8) 

animation_time = 10  # seconds of animation

n_Ne = 100
min_dist_factor = 0.9   # particules séparées d'au moins 0.9*sigma à l'initialisation

t_ref_1 = 10.0
t_ref_2 = 40.0

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

#molecules initialisation (une population par température) ------------------------------------
molecules_1 = populate(Ne, n_Ne, t_ref_1, domain, min_dist_factor)
molecules_2 = populate(Ne, n_Ne, t_ref_2, domain, min_dist_factor)
colors_1 = [get_color(m) for m in molecules_1]
colors_2 = [get_color(m) for m in molecules_2]

#simulation loop (les deux systèmes avancent en parallèle, indépendamment l'un de l'autre)----------------------------------------------------
decouple = max(1, round(Int, n_steps / animation_time / 60))
states_history_1 = [deepcopy(molecules_1)]
states_history_2 = [deepcopy(molecules_2)]

println("simulating.....")
@time for i in 1:n_steps
    step(molecules_1, delta_t = delta_t, domain = domain, t_ref = t_ref_1)
    step(molecules_2, delta_t = delta_t, domain = domain, t_ref = t_ref_2)
    if i % decouple == 0
        push!(states_history_1, deepcopy(molecules_1))
        push!(states_history_2, deepcopy(molecules_2))
    end
    if i == round(n_steps / 2)
        print("50% done...\n")
    end
end

println("simulation done!")

#fig layout--------------------------------------
fig = Figure(size = (1200, 600))

ax2d_1 = Axis(fig[1, 1],
    title = "T_ref = $(t_ref_1) K",
    limits = (-domain.x/2, domain.x/2, -domain.y/2, domain.y/2),
    aspect = DataAspect()
)

ax2d_2 = Axis(fig[1, 2],
    title = "T_ref = $(t_ref_2) K",
    limits = (-domain.x/2, domain.x/2, -domain.y/2, domain.y/2),
    aspect = DataAspect()
)

#observables-----------------------------------------
obs_positions_1 = Observable([Point2f(m.position[1], m.position[2]) for m in molecules_1])
obs_positions_2 = Observable([Point2f(m.position[1], m.position[2]) for m in molecules_2])

# diamètre du marqueur = sigma (composition fixe, donc calculé une seule fois)
sizes_1 = [m.sigma for m in molecules_1]
sizes_2 = [m.sigma for m in molecules_2]

#plots-----------------
scatter!(ax2d_1, obs_positions_1, markersize = sizes_1, markerspace = :data, color = colors_1)
scatter!(ax2d_2, obs_positions_2, markersize = sizes_2, markerspace = :data, color = colors_2)

#animation loop ----------------------------------------------------------------------------
running = Observable(true)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
        running[] = false
    end
end

screen = display(fig)
while running[]
    for frame in 1:length(states_history_1)
        running[] || break
        obs_positions_1[] = [Point2f(m.position[1], m.position[2]) for m in states_history_1[frame]]
        obs_positions_2[] = [Point2f(m.position[1], m.position[2]) for m in states_history_2[frame]]
        sleep(1/60)
    end
end
close(screen)