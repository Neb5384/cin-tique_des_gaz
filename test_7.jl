#import-------------------------
include("main.jl")
using GLMakie

#base variables-----------------------------
delta_t = 1e-14
t_final = 20e-11
n_steps = round(Int, t_final / delta_t)

domain = Domain(10e-9,10e-9,10e-9)

animation_time = 5 #secondes d'animation

n_molecules = 1000. 

initial_atom_speed = 1400

#molecules initialisation------------------------------------
molecules = Vector{Molecule}(undef, n_molecules)
for i in 1:n_molecules

    mol = copy(He)

    mol.position = (
        (rand() -0.5) * domain.x,
        (rand() -0.5) * domain.y,
        (rand() -0.5) * domain.z
    )

    #uniformally distributed angles with same speed
    angle_vector_nonormalised = (rand()-0.5,rand()-0.5,rand()-0.5)
    angle_vector = angle_vector_nonormalised./ norm(angle_vector_nonormalised)
    mol.speed = (
        angle_vector[1] * initial_atom_speed,
        angle_vector[2] * initial_atom_speed,
        angle_vector[3] * initial_atom_speed
    )

    molecules[i] = mol
end



#animation parameters--------------------------------------
fig = Figure(size=(1200, 1200))
ax = Axis3(fig[1,1],limits=(-domain.x/2,domain.x/2,-domain.y/2,domain.y/2,-domain.z/2,domain.z/2))

positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
scatterplot = scatter!(ax, positions,markersize=15)

decouple = round(Int, n_steps/animation_time/60)

#analytics-tracking--------------------------------------------------
speeds = zeros(Float64, n_steps)
speeds[1] = speedMagnitude(molecules)

temperatures = zeros(Float64, n_steps)
temperatures[1] = temperature(molecules)

pressions = zeros(Float64, n_steps)
pressions[1] = pression(molecules,domain)


#simulation loop which copies states that will be shown imn animation later-----------------
states_history = [deepcopy(molecules)]
println("simulating.....")
@time for i in 1:n_steps
    step(molecules, delta_t = delta_t,domain = domain)

    speeds[i] = speedMagnitude(molecules)
    temperatures[i] = temperature(molecules)
    pressions[i] = pression(molecules,domain)

    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
    end
    if i == round(n_steps / 2)
        print("50% done... \n")
    end    
end
println("simulation done !")



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
        positions[] = [Point3f(m.position[1], m.position[2],m.position[3]) for m in state]
        sleep(1/60)  # 60 fps
    end
end
close(screen)


#analytics plots and results------------------------------------------------------------------
fig_analytics = Figure(size=(1200, 1200))
times = (1:n_steps) .* delta_t

#mean speed plot----
ax2 = Axis(fig_analytics[1,1], xlabel = "Time (s)", ylabel = "Speed (m/s)", title = "Mean Speed over Time")
lines!(ax2, times, speeds)

#distribution of maghitudesd of speeds at last timestep---
precision_steps = 40

magnitudes = [norm(m.speed) for m in molecules]
speed_precision = (maximum(magnitudes)) /precision_steps

magnitude_distribution = zeros(Int,precision_steps)

for s in magnitudes
    magnitude_idx = Int(round(s/speed_precision))
    magnitude_distribution[magnitude_idx] += 1
end    

ax3 = Axis(fig_analytics[1,2], xlabel = "speed (m/s)", ylabel = "number", title = "speed magnitude distribution")
speeds3 = (1:precision_steps)*speed_precision
barplot!(ax3,speeds3,magnitude_distribution)

#temperature plot 
ax4 = Axis(fig_analytics[2,1], xlabel = "Time (s)", ylabel = "temperature (K)", title = "Temperature over Time")
lines!(ax4, times, temperatures)

#pression plot
ax5 = Axis(fig_analytics[2,2], xlabel = "Time (s)", ylabel = "Pression (Pa)", title = "Pression over Time")
lines!(ax5, times, pressions)


screen2 = display(fig_analytics)
