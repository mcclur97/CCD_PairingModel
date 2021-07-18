##Packages
using OMEinsum

## Defining Functions
function init_pairing_v(g, pnum, hnum)
    """
    returns potential matrices of the pairing model in three relevant channels
    
    param g: strength of the pairing interaction
    param pnum: number of particle states
    param hnum: number of hole states
    
    return v_pppp, v_pphh, v_hhhh: np.array(pnum,pnum,pnum,pnum), 
                                   np.array(pnum,pnum,hnum,hnum), 
                                   np.array(hnum,hnum,hnum,hnum), 
                                   The interaction as a 4-indexed tensor in three channels.
    """
    v_pppp = zeros((pnum,pnum,pnum,pnum))
    v_pphh = zeros((pnum,pnum,hnum,hnum))
    v_hhhh = zeros((hnum,hnum,hnum,hnum))

    gval = -0.5 * g
    for a in (1:2:pnum)
        for b in (1:2:pnum)
            v_pppp[a, a+1, b, b+1] = gval
            v_pppp[a+1, a, b, b+1] = -gval
            v_pppp[a, a+1, b+1, b] = -gval
            v_pppp[a+1, a, b+1, b] = gval
        end
    end

    for a in (1:2:pnum)
        for i in (1:2:hnum)
            v_pphh[a, a+1, i, i+1] = gval
            v_pphh[a+1, a, i, i+1] = -gval
            v_pphh[a, a+1, i+1, i] = -gval
            v_pphh[a+1, a, i+1, i] = gval
        end
    end
            
    for j in (1:2:hnum)
        for i in (1:2:hnum)
            v_hhhh[j, j+1, i, i+1] = gval
            v_hhhh[j+1, j, i, i+1] = -gval
            v_hhhh[j, j+1, i+1, i] = -gval
            v_hhhh[j+1, j, i+1, i] = gval
        end
    end

    return v_pppp, v_pphh, v_hhhh
end

function init_pairing_fock(delta, g, pnum, hnum)
    """
    initializes the Fock matrix of the pairing model
    
    param delta: Single-particle spacing, as in Eq. (8.41)
    param g: pairing strength, as in eq. (8.42)
    param pnum: number of particle states
    param hnum: number of hole states
    
    return f_pp, f_hh: The Fock matrix in two channels as numpy arrays np.array(pnum,pnum), np.array(hnum,hnum). 
    """
    deltaval = 0.5*delta
    gval = -0.5*g
    f_pp = zeros((pnum,pnum))
    f_hh = zeros((hnum,hnum))

    for i in (1:2:hnum)
        f_hh[i, i] = (deltaval*(i-1))+gval
        f_hh[i+1, i+1] = (deltaval*(i-1))+gval
    end
        
    for a in (1:2:pnum)
        f_pp[a, a] = deltaval*(hnum+(a-1))
        f_pp[a+1, a+1] = deltaval*(hnum+(a-1))
    end
    
    return f_pp, f_hh
end

function init_t2(v_pphh,f_pp,f_hh)
    """
    Initializes t2 amlitudes as in MBPT2, see first equation on page 345
    
    param v_pphh: pairing tensor in pphh channel
    param f_pp:   Fock matrix in pp channel
    param f_hh:   Fock matrix in hh channel
    
    return t2: numpy array in pphh format, 4-indices tensor
    """
    pnum = size(f_pp)[1]
    hnum = size(f_hh)[1]
    t2_new = zeros((pnum,pnum,hnum,hnum))
    for i in (1:1:hnum)
        for j in (1:1:hnum)
            for a in (1:1:pnum)
                for b in (1:1:pnum)
                    denom = (f_hh[i, i] + f_hh[j, j] - f_pp[a, a] - f_pp[b, b])
                    t2_new[a, b, i, j] = v_pphh[a, b, i, j] / denom
                end
            end
        end
    end
    return t2_new
end

# CCD equations. Note that the "->abij" assignment is redundant, because indices are ordered alphabetically.
# Nevertheless, we retain it for transparency.
function ccd_iter(v_pppp, v_pphh, v_hhhh, f_pp, f_hh, t2)
    """
    Performs one iteration of the CCD equations (8.34), using also intermediates for the nonliniar terms
    
    param v_pppp: pppp-channel pairing tensor, numpy array
    param v_pphh: pphh-channel pairing tensor, numpy array
    param v_hhhh: hhhh-channel pairing tensor, numpy array
    param f_pp: Fock matrix in pp channel
    param f_hh: Fock matrix in hh channel
    param t2: Initial t2 amplitude, tensor in form of pphh channel
    
    return t2_new: new t2 amplitude, tensor in form of pphh channel
    """
    pnum = size(f_pp)[1]
    hnum = size(f_hh)[1]
    Hbar_pphh = (  v_pphh 
                 + ein"bc,acij->abij"(f_pp,t2) 
                 - ein"ac,bcij->abij"(f_pp,t2) 
                 - ein"abik,kj->abij"(t2,f_hh)
                 + ein"abjk,ki->abij"(t2,f_hh)
                 + 0.5*ein"abcd,cdij->abij"(v_pppp,t2) 
                 + 0.5*ein"abkl,klij->abij"(t2,v_hhhh)
                )

    # hh intermediate, see (8.47)
    chi_hh = 0.5* ein"cdkl,cdjl->kj"(v_pphh,t2)

    Hbar_pphh = Hbar_pphh - (  ein"abik,kj->abij"(t2,chi_hh) 
                             - ein"abik,kj->abji"(t2,chi_hh) )

    # pp intermediate, see (8.46)
    chi_pp = -0.5* ein"cdkl,bdkl->cb"(v_pphh,t2)

    Hbar_pphh = Hbar_pphh + (  ein"acij,cb->abij"(t2,chi_pp) 
                             - ein"acij,cb->baij"(t2,chi_pp) )

    # hhhh intermediate, see (8.48)
    chi_hhhh = 0.5 * ein"cdkl,cdij->klij"(v_pphh,t2)

    Hbar_pphh = Hbar_pphh + 0.5 * ein"abkl,klij->abij"(t2,chi_hhhh)

    # phph intermediate, see (8.49)
    chi_phph= + 0.5 * ein"cdkl,dblj->bkcj"(v_pphh,t2)


    Hbar_pphh = Hbar_pphh + (  ein"bkcj,acik->abij"(chi_phph,t2)
                             - ein"bkcj,acik->baij"(chi_phph,t2)
                             - ein"bkcj,acik->abji"(chi_phph,t2)
                             + ein"bkcj,acik->baji"(chi_phph,t2) )
                 
    t2_new = zeros((pnum,pnum,hnum,hnum))
    for i in (1:1:hnum)
        for j in (1:1:hnum)
            for a in (1:1:pnum)
                for b in (1:1:pnum)
                    denom = (f_hh[i, i] + f_hh[j, j] - f_pp[a, a] - f_pp[b, b])
                    t2_new[a,b,i,j] = (t2[a, b, i, j] + Hbar_pphh[a, b, i, j] / denom )
                end
            end
        end
    end
    return t2_new
end

function ccd_energy(v_pphh,t2)
    """
    Computes CCD energy. Call as 
    energy = ccd_energy(v_pphh,t2)
    
    param v_pphh: pphh-channel pairing tensor, numpy array
    param t2: t2 amplitude, tensor in form of pphh channel
    
    return energy: CCD correlation energy
    """
    
    erg = 0.25*sum(v_pphh .* t2)
    return erg 
end
## Main Program 

pnum = 4 # number of particle states
hnum = 4 # number of hole states
delta = 1.0
g = 0.5

print("parameters: ")
print("delta =", delta, ", g =", g)
print("
")

v_pppp, v_pphh, v_hhhh = init_pairing_v(g,pnum,hnum)
f_pp, f_hh = init_pairing_fock(delta,g,pnum,hnum)

# Initialize T2 amplitudes from MBPT2
t2 = init_t2(v_pphh,f_pp,f_hh)
erg = ccd_energy(v_pphh,t2)



# Exact MBPT2 for comparison, see last equation on page 365 
exact_mbpt2 = -0.25*g^2*(1.0/(2.0+g) + 2.0/(4.0+g) + 1.0/(6.0+g))

## Actual Calculation
niter=60
time = @elapsed for iter in (1:1:niter)
    t2_new = ccd_iter(v_pppp,v_pphh,v_hhhh,f_pp,f_hh,t2)
    erg = ccd_energy(v_pphh,t2_new)
    t2 = 0.5 * (t2_new + t2)
end
#Prints out final answer
print("TIME:",time)
print("
")
print("ENERGY:",erg)
print("
")
print("ITERATIONS:",niter)
