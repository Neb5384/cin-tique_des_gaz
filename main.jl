#AAAAAARGHH AAAARGH AAAAAAGH

using LinearAlgebra

#structure de molecule ---------------------------------------------------------------------------------
Base.@kwdef mutable struct Molecule
    formule_chimique::String
    masse::Float64
    rayon::Float64
    position::NTuple{3,Float64}
    speed::NTuple{3,Float64}
end

function copy(molecule::Molecule)
    Molecule(
        formule_chimique = molecule.formule_chimique,
        masse    = molecule.masse,
        rayon    = molecule.rayon,
        position = molecule.position,
        speed    = molecule.speed
    )
end

#variables (here for retrofitting...)--------------------------
#g = -9.81e13
g = 0.0

edge_temp = true
(t_plus, t_minus) = (700, 300)

#simulation computations---------------------------------------------------------------------------------------


function computeNextPosition(molecule::Molecule; acceleration::NTuple{3,Float64} = (0.0,0.0,g), delta_t::Float64)
    molecule.speed = molecule.speed .+ acceleration .* delta_t
    molecule.position = molecule.position .+ molecule.speed .* delta_t
end

function computeCollisionSpeed(molecule1::Molecule, molecule2::Molecule)
    r1,v1,m1 = molecule1.position, molecule1.speed, molecule1.masse
    r2,v2,m2 = molecule2.position, molecule2.speed, molecule2.masse

    molecule1.speed = v1 .- 2*m2/(m1+m2) *  dot((v1.-v2),(r1.-r2)) / dot((r1.-r2),(r1.-r2)) .* (r1.-r2)
    molecule2.speed = v2 .+ 2*m1/(m1+m2) *  dot((v1.-v2),(r1.-r2)) / dot((r1.-r2),(r1.-r2)) .* (r1.-r2)
end

function isCollision(molecule1::Molecule, molecule2::Molecule)
    vector_dist = molecule1.position .- molecule2.position
    dist = sqrt(vector_dist[1]^2 + vector_dist[2]^2 + vector_dist[3]^2)
    dist <= molecule1.rayon + molecule2.rayon
end

#Bounding box-----------------------------------------------------------------------------------------

struct Domain
    x:: Float64
    y:: Float64
    z:: Float64
end

function domainVolume(domain::Domain)
    domain.x * domain.y * domain.z
end

function reflect1D(position::Float64,domain::Float64)
    if position > 0
        domain - position
    else
        -domain - position 
    end    
end            

function specularReflection(molecule::Molecule, domain::Domain)
    position = molecule.position
    speed = molecule.speed
    if abs(position[1]) > domain.x/2
        speed = (-speed[1],speed[2],speed[3])
        position = (reflect1D(position[1],domain.x),position[2],position[3])
    end
    if abs(position[2]) > domain.y/2
        speed = (speed[1],-speed[2],speed[3])
        position = (position[1],reflect1D(position[2],domain.y),position[3])
    end
    if abs(position[3]) > domain.z/2
        speed = (speed[1],speed[2],-speed[3])
        position = (position[1],position[2],reflect1D(position[3],domain.z))    
    end
    molecule.speed = speed
    molecule.position = position
end

function temperatureReflection(molecule::Molecule, (t_plus, t_minus)::Tuple{Int,Int})
    c_boltzmann = 1.380649e-23

    position = molecule.position
    speed = molecule.speed
    
    if position[3] > domain.z/2
        sigma = sqrt(c_boltzmann * t_plus / molecule.masse)
        speed = (randn() * sigma, randn() * sigma, - sigma * sqrt(-2*log(rand())))
        position = (position[1],position[2],reflect1D(position[3],domain.z))    
    end
    if position[3] < - domain.z/2
        sigma = sqrt(c_boltzmann * t_minus / molecule.masse)
        speed = (randn() * sigma, randn() * sigma, + sigma * sqrt(-2*log(rand())))
        position = (position[1],position[2],reflect1D(position[3],domain.z))    
    end

    molecule.speed = speed
    molecule.position = position
    return
end    

#step function, does all that is needed to be done in a step--------------------------------------------------------------------

function step(molecules::Vector{Molecule};delta_t::Float64, domain::Domain)
    for (i, molecule) in enumerate(molecules)
        for other_molecule in molecules[i+1:end]
            if isCollision(molecule, other_molecule) computeCollisionSpeed(molecule, other_molecule);break end
        end
    end
    for molecule in molecules
        computeNextPosition(molecule, delta_t = delta_t)
        if edge_temp
            temperatureReflection(molecule, (t_plus, t_minus))
        end  
        specularReflection(molecule,domain)  
    end
end

#energy, movement quantity, speed magnitude, temperature, pression computations-------------------------------------------------------------------------

function momentum(molecule::Molecule)
    molecule.masse .* molecule.speed 
end
function momentum(molecules::Vector{Molecule})
    mom = (0.0,0.0,0.0)
    for molecule in molecules
        mom = mom.+ momentum(molecule)
    end
    mom
end

function cineticEnergy(molecule::Molecule)
    0.5 * molecule.masse * dot(molecule.speed,molecule.speed) 
end
function cineticEnergy(molecules::Vector{Molecule})
    cin_e = 0.0
    for molecule in molecules
        cin_e += cineticEnergy(molecule)
    end
    cin_e
end

function speedMagnitude(molecule::Molecule)
    norm(molecule.speed)
end
function speedMagnitude(molecules::Vector{Molecule})
    sum(speedMagnitude(m) for m in molecules) / length(molecules)
end  

function temperature(molecules::Vector{Molecule})
    c_boltzmann = 1.380649e-23
    sum(m.masse * dot(m.speed,m.speed) for m in molecules)/length(molecules)  / (3 * c_boltzmann)
end

function pression(molecules::Vector{Molecule},domain::Domain)
    sum(m.masse * dot(m.speed,m.speed) for m in molecules) / (3 * domainVolume(domain))
end    


#entropy computations---------------------------------------------------------------------------------------

function slice_volume(molecules::Vector{Molecule}, domain::Domain, num_slices::Int, dimension::Int)
    sliced_molecules = [Molecule[] for _ in 1:num_slices]
    domain_size = getfield(domain, dimension)

    for mol in molecules
        position = mol.position[dimension]
        bin = floor(Int, (position / domain_size + 0.5) * num_slices) + 1
        bin = clamp(bin, 1, num_slices)
        push!(sliced_molecules[bin], mol)
    end
    return sliced_molecules
end

function slice_sq_speed(molecules::Vector{Molecule}, num_slices::Int, max_sq_speed::Int)
    sliced_molecules = [Molecule[] for _ in 1:num_slices]

    for mol in molecules
        bin = floor(Int, norm(mol.speed)^2 / max_sq_speed * num_slices) + 1
        bin = clamp(bin, 1, num_slices)
        push!(sliced_molecules[bin], mol)
    end    
    return sliced_molecules
end    

function entropy_on_slices(sliced_molecules::Vector{Vector{Molecule}}, tot_mol::Int)
    entropy = 0
    for slice in sliced_molecules
        proba = length(slice) / tot_mol
        if proba != 0
            entropy -= proba * log2(proba)
        end    
    end    
    return entropy
end


function position_entropy(molecules::Vector{Molecule}, domain::Domain, num_slices::NTuple{3, Int})
    tot_mol = length(molecules)
    x_sliced = slice_volume(molecules,domain,num_slices[1],1)
    y_sliced = slice_volume(molecules,domain,num_slices[2],2)
    z_sliced = slice_volume(molecules,domain,num_slices[3],3)

    entropy = entropy_on_slices(x_sliced, tot_mol) + entropy_on_slices(y_sliced, tot_mol) + entropy_on_slices(z_sliced, tot_mol)
    return entropy
end    

function speed_entropy(molecules::Vector{Molecule}, num_slices::Int, max_sq_speed::Int)
    tot_mol = length(molecules)
    v_sliced = slice_sq_speed(molecules,num_slices,max_sq_speed)

    entropy = entropy_on_slices(v_sliced, tot_mol)
    return entropy
end    

function entropy(molecules::Vector{Molecule}, domain::Domain, pos_slices::NTuple{3, Int}, sq_speed_slices::Int)
    return position_entropy(molecules,domain,pos_slices) + speed_entropy(molecules,sq_speed_slices,200000)
end    


#some Base molecules--------------------------------------------
#rayon de Van der Waals correspond à la moitié de la distance minimale entre les noyaux de deux atomes non liés, lorsqu'ils s'approchent au plus près sans former de liaison chimique

He = Molecule(
    formule_chimique = "He",
    masse    = 6.646e-27,
    rayon    = 1.10e-10,            #hard sphere model radius 
    position = (0.0, 0.0, 0.0),
    speed    = (0.0, 0.0, 0.0)
)
Ne = Molecule(
    formule_chimique = "Ne",
    masse    = 3.351e-26,
    rayon    = 1.54e-10,
    position = (0.0, 0.0, 0.0),
    speed    = (0.0, 0.0, 0.0)
)
N2 = Molecule(
    formule_chimique = "N2",
    masse    = 4.652e-26,
    rayon    = 1.55e-10,
    position = (0.0, 0.0, 0.0),
    speed    = (0.0, 0.0, 0.0)
)
O2 = Molecule(
    formule_chimique = "O2",
    masse    = 5.314e-26,
    rayon    = 1.52e-10,
    position = (0.0, 0.0, 0.0),
    speed    = (0.0, 0.0, 0.0)
)
Ar = Molecule(
    formule_chimique = "Ar",
    masse    = 6.634e-26,
    rayon    = 1.88e-10,
    position = (0.0, 0.0, 0.0),
    speed    = (0.0, 0.0, 0.0)
)


