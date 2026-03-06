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

#simulation computations---------------------------------------------------------------------------------------

function computeNextPosition(molecule::Molecule; acceleration::NTuple{3,Float64} = (0.0,0.0,0.0), delta_t::Float64)
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


#step function, does all that is needed to be done in a step--------------------------------------------------------------------

function step(molecules::Vector{Molecule};delta_t::Float64, domain::Domain)
    for (i, molecule) in enumerate(molecules)
        for other_molecule in molecules[i+1:end]
            if isCollision(molecule, other_molecule) computeCollisionSpeed(molecule, other_molecule);break end
        end
    end
    for molecule in molecules
        computeNextPosition(molecule, delta_t = delta_t)
        specularReflection(molecule,domain)
    end
end

#energy and movement quantity computations-------------------------------------------------------------------------

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
    cin_e = 0
    for molecule in molecules
        cin_e += cineticEnergy(molecule)
    end
    cin_e
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



#some Base molecules--------------------------------------------
#rayon de Van der Waals correspond à la moitié de la distance minimale entre les noyaux de deux atomes non liés, lorsqu'ils s'approchent au plus près sans former de liaison chimique

He = Molecule(
    formule_chimique = "He",
    masse    = 6.647e-27,
    rayon    = 1.40e-10,
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

