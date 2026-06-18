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

T1 = 40.0   # premiere temperature
N1 = 15000 # nbre de pas pour premiere temperature

T2 = 10.0 # 2 eme temperature
N2 = 15000 # nbre de pas de changement

n_bins = 25   # résolution de la grille pour la distribution de présence

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

#temperature change helpers ------------------------------------

# retourne la température cible au pas `step` selon le profil :
# constante T1 pendant N1 pas, decroissance lineaire pendant N2 pas, puis constante T2
function movingTemperature(step::Int)
    if step <= N1
        return T1
    elseif step <= N1 + N2
        frac = (step - N1) / N2
        return T1 + (T2 - T1) * frac
    else
        return T2
    end
end

# retourne (titre, couleur de fond) pour le pas `step`, selon la même découpe en 3 phases
function phaseInfo(step::Int)
    if step <= N1
        return ("Phase 1 : palier à T1 = $(T1) K", RGBf(1.0, 0.93, 0.85))     # léger orange (chaud)
    elseif step <= N1 + N2
        return ("Phase 2 : refroidissement $(T1) K → $(T2) K", RGBf(1.0, 1.0, 0.85))  # léger jaune (transition)
    else
        return ("Phase 3 : palier à T2 = $(T2) K", RGBf(0.85, 0.93, 1.0))     # léger bleu (froid)
    end
end

#distribution de présence spatiale helpers ------------------------------------

# convertit une coordonnée continue dans [-extent/2, extent/2] en indice de bin (1..n_bins)
function positionToBin(pos::Float64, extent::Float64, n_bins::Int)
    frac = (pos + extent / 2) / extent
    idx = floor(Int, frac * n_bins) + 1
    return clamp(idx, 1, n_bins)
end

# incrémente l'histogramme 2D avec les positions courantes des molécules
function accumulateHistogram!(hist::Matrix{Float64}, molecules::Vector{Molecule}, domain::Domain, n_bins::Int)
    for m in molecules
        ix = positionToBin(m.position[1], domain.x, n_bins)
        iy = positionToBin(m.position[2], domain.y, n_bins)
        hist[ix, iy] += 1.0
    end
end

# normalise un histogramme de comptages en probabilité sans dimension (somme sur tous les bins = 1)
function toProbability(hist::Matrix{Float64})
    total = sum(hist)
    return total > 0 ? hist ./ total : hist
end

#molecules initialisation ------------------------------------
molecules = populate(Ne, n_Ne, T1, domain, min_dist_factor)
colors = [get_color(m) for m in molecules]

#histogrammes de présence, un par phase ------------------------------------
hist_phase1 = zeros(n_bins, n_bins)
hist_phase2 = zeros(n_bins, n_bins)
hist_phase3 = zeros(n_bins, n_bins)

#simulation loop ----------------------------------------------------
decouple = max(1, round(Int, n_steps / animation_time / 60))
states_history = [deepcopy(molecules)]
steps_history = [1]   # pas de simulation correspondant à chaque état stocké, pour retrouver la phase à l'animation

println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t, domain = domain, t_ref = movingTemperature(i))

    if i <= N1
        accumulateHistogram!(hist_phase1, molecules, domain, n_bins)
    elseif i <= N1 + N2
        accumulateHistogram!(hist_phase2, molecules, domain, n_bins)
    else
        accumulateHistogram!(hist_phase3, molecules, domain, n_bins)
    end

    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
        push!(steps_history, i)
    end
    if i == round(n_steps / 2)
        print("50% done...\n")
    end
end

println("simulation done!")

#normalisation des distributions de présence en probabilités sans dimension------------
prob_phase1 = toProbability(hist_phase1)
prob_phase2 = toProbability(hist_phase2)
prob_phase3 = toProbability(hist_phase3)

bin_centers_x = range(-domain.x/2 + domain.x/(2*n_bins), domain.x/2 - domain.x/(2*n_bins), length = n_bins)
bin_centers_y = range(-domain.y/2 + domain.y/(2*n_bins), domain.y/2 - domain.y/(2*n_bins), length = n_bins)

vmax = max(maximum(prob_phase1), maximum(prob_phase2), maximum(prob_phase3))

#fig layout : une seule figure/fenêtre --------------------------------------
# - l'animation occupe un bloc de 2 lignes x 3 colonnes en haut
# - chaque heatmap occupe une case 1x1 sur la ligne du bas (+ une colonne pour la colorbar commune)
fig = Figure(size = (1000, 1000))

init_title, init_bg = phaseInfo(1)

ax2d = Axis(fig[1:2, 1:3],
    title = init_title,
    backgroundcolor = init_bg,
    limits = (-domain.x/2, domain.x/2, -domain.y/2, domain.y/2),
    aspect = DataAspect()
)

ax_d1 = Axis(fig[3, 1], title = "Phase 1 : T = $(T1) K", aspect = DataAspect())
ax_d2 = Axis(fig[3, 2], title = "Phase 2 : $(T1) K → $(T2) K", aspect = DataAspect())
ax_d3 = Axis(fig[3, 3], title = "Phase 3 : T = $(T2) K", aspect = DataAspect())

heatmap!(ax_d1, bin_centers_x, bin_centers_y, prob_phase1, colorrange = (0, vmax), colormap = :viridis)
heatmap!(ax_d2, bin_centers_x, bin_centers_y, prob_phase2, colorrange = (0, vmax), colormap = :viridis)
hm3 = heatmap!(ax_d3, bin_centers_x, bin_centers_y, prob_phase3, colorrange = (0, vmax), colormap = :viridis)

Colorbar(fig[3, 4], hm3, label = "probabilité de présence")

#observables-----------------------------------------
obs_positions = Observable([Point2f(m.position[1], m.position[2]) for m in molecules])

# diamètre du marqueur = sigma (composition fixe, donc calculé une seule fois)
sizes = [m.sigma for m in molecules]

#plots-----------------
scatter!(ax2d, obs_positions, markersize = sizes, markerspace = :data, color = colors)

#animation loop ----------------------------------------------------------------------------
running = Observable(true)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
        running[] = false
    end
end

screen = display(fig)
while running[]
    for frame in 1:length(states_history)
        running[] || break
        obs_positions[] = [Point2f(m.position[1], m.position[2]) for m in states_history[frame]]

        # marquage de la phase courante : titre + léger fond de couleur
        title, bgcolor = phaseInfo(steps_history[frame])
        ax2d.title = title
        ax2d.backgroundcolor = bgcolor

        sleep(1/60)
    end
end
close(screen)