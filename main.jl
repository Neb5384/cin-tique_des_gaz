#AAAAAARGHH AAAARGH AAAAAAGH

#rayon de Van der Waals correspond à la moitié de la distance minimale entre les noyaux de deux atomes non liés, lorsqu'ils s'approchent au plus près sans former de liaison chimique
using LinearAlgebra

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

function computeNextPosition(molecule::Molecule; acceleration::NTuple{3,Float64} = (0.0,0.0,0.0), delta_t::Float64)
    molecule.speed = molecule.speed .+ acceleration .* delta_t
    molecule.position = molecule.position .+ molecule.speed .* delta_t
end

function computeCollisionSpeed(molecule1::Molecule, molecule2::Molecule)
    r1,v1,m1 = molecule1.position, molecule1.speed, molecule1.masse
    r2,v2,m2 = molecule2.position, molecule2.speed, molecule2.masse

    molecule1.speed = v1 .- 2*m1/(m1+m2) *  dot((v1.-v2),(r1.-r2)) / dot((r1.-r2),(r1.-r2)) .* (r1.-r2)
    molecule2.speed = v2 .+ 2*m2/(m1+m2) *  dot((v1.-v2),(r1.-r2)) / dot((r1.-r2),(r1.-r2)) .* (r1.-r2)
end

function isCollision(molecule1::Molecule, molecule2::Molecule)
    vector_dist = molecule1.position .- molecule2.position
    dist = sqrt(vector_dist[1]^2 + vector_dist[2]^2 + vector_dist[3]^2)
    dist <= molecule1.rayon + molecule2.rayon
end

function step(molecules::Vector{Molecule};delta_t::Float64)
    for (i, molecule) in enumerate(molecules)
        for other_molecule in molecules[i+1:end]
            if isCollision(molecule, other_molecule) computeCollisionSpeed(molecule, other_molecule);break end
        end
    end
    for molecule in molecules
        computeNextPosition(molecule, delta_t = delta_t)
    end
end

function momentum(molecule::Molecule)
    molecule.masse * norm(molecule.speed) 
end
function momentum(molecules::Vector{Molecule})
    mom = 0
    for molecule in molecules
        mom += momentum(molecule)
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

