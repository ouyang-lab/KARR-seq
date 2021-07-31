#/usr/bin/env python3

import scipy.sparse as sps
import numpy as np
import sys

# Scripts were retrieved from the follow git repository
# https://github.com/ay-lab/HiCKRy

##FUNCTION DESCRIPTION
# knighRuizAlg is an implementation of the matrix balancing algorithm
#  developed by Knight and Ruiz. The goal is to take a matrix A and
#  find a vector x such that, diag(x)*A*diag(x) returns a doubly
#  stochastic matrix

##PARAMETERS
#A is a given numpy array
#tol is error tolerance
#f1 boolean indicating if the intermediate convergance statistics
# should also be outputted
def knightRuizAlg(A, tol=1e-6, f1 = False):
    n = A.shape[0]
    e = np.ones((n,1), dtype = np.float64)
    res = []


    Delta = 3
    delta = 0.1
    x0 = np.copy(e)
    g = 0.9

    etamax = eta = 0.1
    stop_tol = tol*0.5
    x = np.copy(x0)

    rt = tol**2.0
    v = x * (A.dot(x))
    rk = 1.0 - v
    #rho_km1 = np.dot(rk.T, rk)[0, 0]
    rho_km1 = ((rk.transpose()).dot(rk))[0,0]
    rho_km2 = rho_km1
    rout = rold = rho_km1

    MVP = 0 #we'll count matrix vector products
    i = 0 #outer iteration count

    if f1:
        print(("it        in. it      res\n"), end=' ')

    while rout > rt: #outer iteration
        i += 1

        if i > 30:
            break

        k = 0
        y = np.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)

        while rho_km1 > innertol: #inner iteration by CG
            k += 1
            if k == 1:
                Z = rk / v
                p = np.copy(Z)
                #rho_km1 = np.dot(rk.T, Z)
                rho_km1 = (rk.transpose()).dot(Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            if k > 10:
                break


            #update search direction efficiently
            w = x * A.dot(x * p) + v * p
            #alpha = rho_km1 / np.dot(p.T, w)[0,0]
            alpha = rho_km1 / (((p.transpose()).dot(w))[0,0])
            ap = alpha * p
            #test distance to boundary of cone
            ynew = y + ap

            if np.amin(ynew) <= delta:

                if delta == 0:
                    break

                ind = np.where(ap < 0.0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = np.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            #rho_km1 = np.dot(rk.T, Z)[0,0]
            rho_km1 = ((rk.transpose()).dot(Z))[0,0]
        x *= y
        v = x * (A.dot(x))
        rk = 1.0 - v
        #rho_km1 = np.dot(rk.T, rk)[0,0]
        rho_km1 = ((rk.transpose()).dot(rk))[0,0]
        rout = rho_km1
        MVP += k + 1

        #update inner iteration stopping criterion
        rat = rout/rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if f1:
            print(("%03i %06i %03.3f %e %e \n") % \
                (i, k, res_norm, rt, rout), end=' ')
            res.append(res_norm)
    if f1:
        print(("Matrix - vector products = %06i\n") % \
            (MVP), end=' ')

    #X = np.diag(x[:,0])
    #x = X.dot(A.dot(X))
    return [x,i,k]

def removeZeroDiagonalCSR(mtx, i=0, toRemovePre=None):
    iteration = 0
    toRemove = []
    ctr = 0

    if toRemovePre is not None:
        for items in toRemovePre:
            toRemove.append(items)

    if i == 0:
        diagonal = mtx.diagonal()
        #print diagonal
        for values in diagonal:
            if values == 0:
                toRemove.append(ctr)
            ctr += 1

    else:
        rowSums = mtx.sum(axis=0)
        rowSums = list(np.array(rowSums).reshape(-1,))
        rowSums = list(enumerate(rowSums))
        for value in rowSums:
            if int(value[1]) == 0:
                toRemove.append(value[0])
                rowSums.remove(value)
                rowSums.sort(key=lambda tup: tup[1])
                size = len(rowSums)
                perc = i/100.0
                rem = int(perc * size)
                while ctr < rem:
                    toRemove.append(rowSums[ctr][0])
                    ctr += 1
    list(set(toRemove))
    toRemove.sort()
    #print toRemove
    mtx = dropcols_coo(mtx, toRemove)
    for num in toRemove:
        if iteration != 0:
            num -= iteration
        removeRowCSR(mtx,num)
        iteration +=1
    return [mtx, toRemove]

def dropcols_coo(M, idx_to_drop):
    idx_to_drop = np.unique(idx_to_drop)
    C = M.tocoo()
    keep = ~np.in1d(C.col, idx_to_drop)
    C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    C.col -= idx_to_drop.searchsorted(C.col) # decrement column indices
    C._shape = (C.shape[0], C.shape[1] - len(idx_to_drop))
    return C.tocsr()

def removeRowCSR(mat, i):
    if not isinstance(mat, sps.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

def addZeroes(mtx, toAdd):
    indices = np.array(toAdd)
    i = indices - np.arange(len(indices))
    mtx = np.insert(np.insert(mtx,i,0,axis=1),i,0,axis=0)
    return mtx

def main():
    #test = sps.rand(5,5,density=0.5,format='csr')
    test = np.random.randint(0,20+1, size=(5,5))
    test = (test + test.T)/2

    test[2,:] = 0.0
    test[:,2] = 0.0
    test[3,:] = 0.0
    test[:,3] = 0.0

    test = sps.csr_matrix(test)
    print(test.toarray())

    percentOfSparseToRemove = 0 #20
    mtxAndRemoved = removeZeroDiagonalCSR(test, percentOfSparseToRemove)
    initialSize = test.shape[0]
    test, removed = mtxAndRemoved[0], mtxAndRemoved[1]
    newSize = test.shape[0]

    print(removed)
    print(initialSize, newSize)

    CC = ((test.sum())/(test.shape[0]*2))
    CCother = test.sum()/test.size
    result = knightRuizAlg(test)
    col = result[0]
    x = sps.diags(col.flatten(), 0, format='csr')
    mtx = x.dot(test.dot(x))

    print(mtx.toarray())
    CCmtx = CC * mtx
    CCothermtx = CCother * mtx
    print(CCmtx.toarray())


    CCmtx = addZeroes(CCmtx.toarray(), removed)
    print(CCmtx)

    assert False

    print()
    print(CCothermtx.toarray())
    print((test.toarray().sum(axis=0)))
    print((CCmtx.toarray().sum(axis=0)))
    print((CCothermtx.toarray().sum(axis=0)))
    #print(result.sum(axis=0))

if __name__=="__main__":
    main()

