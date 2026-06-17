#AAAAAARGHH AAAARGH AAAAAAGH

using LinearAlgebra

#structure de molecule ---------------------------------------------------------------------------------
Base.@kwdef mutable struct Molecule
    formule_chimique::String
    masse::Float64
    sigma::Float64      # Lennard-Jones sigma
    epsilon::Float64     # Lennard-Jones well depth (J)
    position::NTuple{2,Float64}
    speed::NTuple{2,Float64}
end
 
function copy(molecule::Molecule)
    Molecule(
        formule_chimique = molecule.formule_chimique,
        masse    = molecule.masse,
        sigma    = molecule.sigma,
        epsilon  = molecule.epsilon,
        position = molecule.position,
        speed    = molecule.speed
    )
end

#Bounding box-----------------------------------------------------------------------------------------

struct Domain
    x:: Float64
    y:: Float64
end

function domainArea(domain::Domain)
    domain.x * domain.y
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
        speed = (-speed[1],speed[2])
        position = (reflect1D(position[1],domain.x),position[2])
    end
    if abs(position[2]) > domain.y/2
        speed = (speed[1],-speed[2])
        position = (position[1],reflect1D(position[2],domain.y))
    end
    molecule.speed = speed
    molecule.position = position
end

#simulation computations---------------------------------------------------------------------------------------

function computeNextPosition(molecule::Molecule; acceleration::NTuple{2,Float64} = (0.0,0.0), delta_t::Float64)
    molecule.speed = molecule.speed .+ acceleration .* delta_t
    molecule.position = molecule.position .+ molecule.speed .* delta_t
end
  
function isInRange(molecule1::Molecule, molecule2::Molecule)
    vector_dist = molecule1.position .- molecule2.position
    dist = sqrt(vector_dist[1]^2 + vector_dist[2]^2)
    dist <= (molecule1.sigma + molecule2.sigma) / 2 * 3
end

function lennardJonesForce(molecule1::Molecule, molecule2::Molecule)
    sigma   = (molecule1.sigma + molecule2.sigma) / 2
    epsilon = sqrt(molecule1.epsilon * molecule2.epsilon)
 
    r_ij = molecule2.position .- molecule1.position
    r = norm(r_ij)
 
 
    return -24*epsilon * (2*sigma^12/r^12 - sigma^6/r^6) / r^2 .* r_ij
end

function velocityRescale(molecules::Vector{Molecule}, t_ref)
    factor = sqrt(t_ref / temperature(molecules))
    for molecule in molecules
        molecule.speed = molecule.speed .* factor
    end    
end


#step function, does all that is needed to be done in a step--------------------------------------------------------------------

function step(molecules::Vector{Molecule}; delta_t::Float64, domain::Domain, t_ref::Float64)
    n = length(molecules)
    forces = [(0.0, 0.0) for _ in 1:n]

    # calculate Lennard-Jones forces
    for (i, molecule) in enumerate(molecules)
        for j in (i+1):n
            other_molecule = molecules[j]
            if isInRange(molecule, other_molecule)
                f = lennardJonesForce(molecule, other_molecule)
                forces[i] = forces[i] .+ f
                forces[j] = forces[j] .- f
            end
        end
    end

    # semi-implicit euler and specular reflexion
    for (i, molecule) in enumerate(molecules)
        computeNextPosition(molecule, acceleration = forces[i] ./ molecule.masse, delta_t = delta_t)
        specularReflection(molecule, domain)
    end

    # velocity rescaling
    velocityRescale(molecules, t_ref)
end


#energy, movement quantity, speed magnitude, temperature, pression computations-------------------------------------------------------------------------

function momentum(molecule::Molecule)
    molecule.masse .* molecule.speed 
end
function momentum(molecules::Vector{Molecule})
    mom = (0.0,0.0)
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
    sum(m.masse * dot(m.speed,m.speed) for m in molecules)/length(molecules)  / (2 * c_boltzmann)
end

function pression(molecules::Vector{Molecule},domain::Domain)
    sum(m.masse * dot(m.speed,m.speed) for m in molecules) / (2 * domainArea(domain))
end    


#some Base molecules--------------------------------------------

Ne = Molecule(
    formule_chimique = "Ne",
    masse    = 3.351e-26,
    sigma    = 2.74e-10,
    epsilon  = 4.91511044e-22,                 
    position = (0.0, 0.0),
    speed    = (0.0, 0.0)
)
