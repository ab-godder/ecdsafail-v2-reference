"""Merged raw shift/add and a single quotient/remainder network."""
from .arith import add, add_const, gt_const_into, gt_into, cswap, _selected_chain, _selected_carry

def decode(C,k):
    out=C.alloc(5)
    (tmp,)=C.alloc(1)
    for v,t in enumerate(out):
        for i,b in enumerate(k):
            if not (v>>i&1):C.x(b)
        C.ccx(k[0],k[1],tmp);C.ccx(tmp,k[2],t);C.ccx(k[0],k[1],tmp)
        for i,b in enumerate(k):
            if not (v>>i&1):C.x(b)
    C.free([tmp])
    return out

def undec(C,k,out):
    (tmp,)=C.alloc(1)
    for v,t in enumerate(out):
        for i,b in enumerate(k):
            if not (v>>i&1):C.x(b)
        C.ccx(k[0],k[1],tmp);C.ccx(tmp,k[2],t);C.ccx(k[0],k[1],tmp)
        for i,b in enumerate(k):
            if not (v>>i&1):C.x(b)
    C.free([tmp]);C.free(out)

def reduce_fast(C,raw,q):
    n=q.bit_length();f=(1<<n)-q;w=(5*f).bit_length()
    low=raw[:n];k=raw[n:];L,H=low[:w],low[w:]
    prefix=[None]+H[:1]
    for bit in H[1:]:
        (p,)=C.alloc(1);C.ccx(prefix[-1],bit,p);prefix.append(p)
    sels=decode(C,k);masks=[(i+1)*f for i in range(5)]
    chain=_selected_chain(C,L,sels,masks);(e,)=C.alloc(1)
    if H:C.ccx(prefix[-1],chain[-1],e)
    else:C.cx(chain[-1],e)
    for i in reversed(range(w)):
        _selected_carry(C,L[i],sels,masks,i,chain[i],chain[i+1]);C.free([chain[i+1]])
    undec(C,k,sels)
    (tmp,)=C.alloc(1)
    C.ccx(k[0],k[1],tmp);C.ccx(e,tmp,k[2]);C.ccx(k[0],k[1],tmp)
    C.free([tmp]);C.ccx(e,k[0],k[1]);C.cx(e,k[0])
    sels=decode(C,k);masks=[i*f for i in range(5)]
    chain=_selected_chain(C,L,sels,masks)
    if H:C.ccx(prefix[-1],chain[-1],e)
    else:C.cx(chain[-1],e)
    C.free([e])
    for i in reversed(range(1,len(H))):
        C.ccx(prefix[i],H[i],prefix[i+1]);C.ccx(chain[-1],prefix[i],H[i]);C.free([prefix[i+1]])
    if H:C.cx(chain[-1],H[0])
    for i in reversed(range(w)):
        _selected_carry(C,L[i],sels,masks,i,chain[i],chain[i+1]);C.free([chain[i+1]])
        for sel,mask in zip(sels,masks):
            if mask>>i&1:C.cx(sel,L[i])
        if chain[i] is not None:C.cx(chain[i],L[i])
    undec(C,k,sels)
    return k,[]

def reduce_general(C,raw,q):
    # Binary restoring division on the same raw sum, never modular doubling.
    k=C.alloc(3)
    for i in reversed(range(3)):
        C.x(k[i]);gt_const_into(C,q<<i,raw,k[i])
        add_const(C,raw,-(q<<i),ctrl=k[i])
    return k,raw[q.bit_length():]

def clear_quotient(C,k,g,j,a,y,q):
    # High quotient bit is g*[R>Y] only when the low j bits are zero.
    v,u,z=C.alloc(3)
    C.ccx(j,k[1],v)
    C.x(k[0]);C.ccx(g,k[0],u);C.x(k[0])
    C.x(v);C.ccx(u,v,z);C.x(v)
    C.x(j);cswap(C,j,[k[1]],[k[2]]);C.x(j)
    gt_into(C,a,y,k[2],ctrl=z)
    C.x(j);cswap(C,j,[k[1]],[k[2]]);C.x(j)
    C.x(v);C.ccx(u,v,z);C.x(v)
    C.x(k[0]);C.ccx(g,k[0],u);C.x(k[0])
    C.ccx(j,k[1],v);C.free([v,u,z])
    # k mod4 = (g*R-Y)*q^{-1} mod4, for j=2.
    b0,b1,borrow=C.alloc(3)
    C.ccx(g,a[0],b0);C.ccx(g,a[1],b1)
    C.x(b0);C.ccx(b0,y[0],borrow);C.x(b0)
    C.ccx(j,b1,k[1]);C.ccx(j,y[1],k[1]);C.ccx(j,borrow,k[1])
    if q&2:C.ccx(j,k[0],k[1])
    C.x(b0);C.ccx(b0,y[0],borrow);C.x(b0)
    C.ccx(g,a[1],b1);C.ccx(g,a[0],b0);C.free([b0,b1,borrow])
    C.ccx(g,a[0],k[0]);C.cx(y[0],k[0])

def tick(C,mod,g,R,S,m):
    assert m==2 and len(g)==3
    n=mod.n;l0,ext,top,zero=C.alloc(4)
    shifted=S+[ext]
    for i in reversed(range(n)):cswap(C,g[2],[shifted[i]],[shifted[i+1]])
    add(C,R[1:]+[l0,zero],shifted,top,ctrl=g[0]);C.free([zero])
    C.ccx(g[0],R[0],l0)
    raw=[l0]+shifted+[top]
    if 5*((1<<n)-mod.q)<1<<n:k,extra=reduce_fast(C,raw,mod.q)
    else:k,extra=reduce_general(C,raw,mod.q)
    low=raw[:n]
    clear_quotient(C,k,g[0],g[2],R,low,mod.q)
    # Put the low output bit on the old high S wire; all other surplus wires clean.
    C.cx(l0,S[-1]);C.cx(S[-1],l0);C.cx(l0,S[-1])
    if extra:C.free(k)
    C.free([l0,ext,top]);S[:]=S[-1:]+S[:-1]
    cswap(C,g[1],R,S)
