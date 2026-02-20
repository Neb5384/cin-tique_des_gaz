#AAAAAARGHH AAAARGH AAAAAAGH

#rayon de Van der Waals correspond à la moitié de la distance minimale entre les noyaux de deux atomes non liés, lorsqu'ils s'approchent au plus près sans former de liaison chimique

Base.@kwdef mutable struct Molecule
    formule_chimique::String
    masse::Float64
    rayon::Float64
    position::NTuple{3,Float64}
    speed::NTuple{3,Float64}
end

function computeNextPosition(molecule::Molecule; acceleration::NTuple{3,Float64} = (0.0,0.0,0.0), delta_t::Float64)
    molecule.speed = molecule.speed .+ acceleration .* delta_t
    molecule.position = molecule.position .+ molecule.speed .* delta_t
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