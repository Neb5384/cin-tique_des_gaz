include("main.jl")
using GLMakie

molecule1 = copy(O2)
molecule1.position = (1e-7,0,0)
molecule1.speed = (-300,0,0)

molecule2 = copy(O2)
molecule2.position = (-1e-7,0,0)
molecule2.speed = (300,0,0)

molecule3 = copy(O2)
molecule3.position = (1e-7,5e-7,0)
molecule3.speed = (-100,0,0)

molecule4 = copy(O2)
molecule4.position = (-1e-7,5e-7,0)
molecule4.speed = (700,0,0)

molecule5 = copy(He)
molecule5.position = (1e-7,-5e-7,0)
molecule5.speed = (-300,0,0)

molecule6 = copy(O2)
molecule6.position = (-1e-7,-5e-7,0)
molecule6.speed = (300,0,0)

molecule7 = copy(O2)
molecule7.position = (1e-7,10e-7,0)
molecule7.speed = (-300,0,500)

molecule8 = copy(O2)
molecule8.position = (-1e-7,10e-7,0)
molecule8.speed = (300,0,500)



molecules::Vector = [molecule1,molecule2,molecule3,molecule4,molecule5,molecule6,molecule7,molecule8]

delta_t = 1e-12
t_final = 1e-9
n_steps = round(Int, t_final / delta_t)


fig = Figure()
ax = Axis3(fig[1,1],limits=(-1e-6,1e-6,-1e-6,1e-6,-1e-6,1e-6))

positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
scatterplot = scatter!(ax, positions, markersize=20)

animation_time = 2 #secondes d'animation
decouple = round(Int, n_steps/animation_time/60)


#simulation loop which copies states that will be shown imn animation later
states_history = [deepcopy(molecules)]
for i in 1:n_steps
    step(molecules, delta_t = delta_t)
    if i % decouple == 0
        push!(states_history, deepcopy(molecules))
    end
end

#testing if momentum and cinetic enrgy have been conserved
first_momentum = momentum(states_history[1])
last_momentum = momentum(states_history[end])

println("difference quantite de mouvment : ", (last_momentum -first_momentum) / first_momentum, " %" )

first_cinetic_energy = cineticEnergy(states_history[1])
last_cinetic_energy = cineticEnergy(states_history[end])

println("difference energie cinetique : ", (last_cinetic_energy -first_cinetic_energy) / first_cinetic_energy, " %" )

#animation loop 
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


