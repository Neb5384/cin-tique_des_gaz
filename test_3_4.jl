include("main.jl")
using GLMakie

#instance

O2.speed = (500.0, 300.0, 100.0)

delta_t = 1e-12
t_final = 1e-9
n_steps = round(Int, t_final / delta_t)

positions = Vector{NTuple{3, Float64}}(undef, n_steps + 1)
positions[1] = O2.position

for i in 2:n_steps + 1
    computeNextPosition(O2,acceleration =(0.0,0.0,500.0), delta_t = delta_t)
    positions[i] = O2.position
end

# Animation GIF
fig = Figure(resolution = (800, 600))
ax  = Axis3(fig[1, 1],
    xlabel = "x (m)",
    ylabel = "y (m)",
    zlabel = "z (m)",
    title  = "Trajectoire de O₂"
)

# Observable pour le tracé progressif
idx = Observable(1)
head = @lift(Point3f.(
    [positions[$idx][1]],
    [positions[$idx][2]],
    [positions[$idx][3]]
))
scatter!(ax, head, color = :red, markersize = 12)

max_dist = maximum(maximum(p) for p in positions)
limits!(ax, 0, max_dist, 0, max_dist, 0, max_dist)

screen = display(fig)

animation_time = 2
decouple = round(Int, n_steps/animation_time/60)

running = Observable(true)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
        running[] = false
        close(screen)
    end
end

while running[]
    for i in 1:decouple:n_steps+1
        running[] || break
        idx[] = i
        sleep(1/60)  # 60 fps
    end
end