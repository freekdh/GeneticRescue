BDtestpars <- function(AB = 0, Ab = 1, aB = 10, ab = 0) {
    list(
    AB0 = AB,
    Ab0 = Ab,
    aB0 = aB,
    ab0 = ab,
    bA = 1.1,
    ba = 1.1,
    dA = 1.05,
    da = 1.15,
    r = 0.5
    )
}

IBStestpars <- function() {
    list(
    b=1.1,
    dA=1,
    da=0.8,
    nloci=100,
    ninit0=20,
    ninit1=2,
    ngen=50,
    nrep=1000,
    rec=0.5,
    k=20)
}