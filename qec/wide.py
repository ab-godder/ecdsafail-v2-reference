"""Carry trees that spend clean scratch wires for depth.

Every routine here is a permutation built from X, CX and CCX that restores each wire it does not
name as an output, and leaves every scratch wire it allocates clean.
"""


def target(C, *inputs):
    """A clean wire for a gate that reads `inputs`: the circuit's alloc_after when it has one, which
    picks a wire that is free by the time the inputs are, else its next wire."""
    pick = getattr(C, "alloc_after", None)
    return pick(*inputs) if pick else C.alloc(1)[0]


class Tape:
    """Records gates so a compute can be replayed in reverse as its uncompute."""

    def __init__(self):
        self.ops = []

    def x(self, a):
        self.ops.append(("x", a))

    def cx(self, a, b):
        self.ops.append(("cx", a, b))

    def ccx(self, a, b, c):
        self.ops.append(("ccx", a, b, c))

    def play(self, C, reverse=False):
        for name, *wires in reversed(self.ops) if reverse else self.ops:
            getattr(C, name)(*wires)


def and_tree(T, bits, alloc):
    """AND of bits on a fresh wire, by a balanced tree of len(bits) - 1 CCX. One bit is returned as
    it is."""
    level = list(bits)
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            t = alloc()
            T.ccx(level[i], level[i + 1], t)
            nxt.append(t)
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    return level[0]


def _spans(n, block):
    """[lo, hi) blocks covering n bits: `block` bits each, or, when it is a tuple, its sizes from
    the bottom with the last one repeated."""
    sizes = block if isinstance(block, tuple) else (block,)
    out, lo, k = [], 0, 0
    while lo < n:
        hi = min(lo + sizes[min(k, len(sizes) - 1)], n)
        out.append((lo, hi))
        lo, k = hi, k + 1
    return out


def _carry_out(T, A, B, alloc, block):
    """The carry out of A + B, no carry in, onto A's top wire. Blocks of `block` bits ripple in
    place as CDKM majorities, each from its own clean carry-in, which leaves A[i] ^ B[i] on B[i].
    A binary tree then joins the blocks: G_hi ^= P_hi G_lo, with P the AND of a block's B wires,
    made only for a node that is some join's upper half."""
    n = len(A)
    nodes = []
    for lo, hi in _spans(n, block):
        cin = alloc()
        T.cx(A[lo], B[lo])
        T.cx(A[lo], cin)
        T.ccx(cin, B[lo], A[lo])
        for i in range(lo + 1, hi):
            T.cx(A[i], B[i])
            T.cx(A[i], A[i - 1])
            T.ccx(A[i - 1], B[i], A[i])
        nodes.append({"g": A[hi - 1], "bits": B[lo:hi], "kids": None})
    level = nodes
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            nxt.append({"g": level[i + 1]["g"], "kids": (level[i], level[i + 1])})
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    root = level[0]

    def mark(node, need):
        node["need"] = need
        if node["kids"]:
            mark(node["kids"][0], need)
            mark(node["kids"][1], True)

    mark(root, False)

    def build(node):
        if node["kids"] is None:
            node["p"] = and_tree(T, node["bits"], alloc) if node["need"] else None
            return
        lo, hi = node["kids"]
        build(lo)
        build(hi)
        # The product first: it is the next join's P, and both gates hold hi's P.
        node["p"] = None
        if node["need"]:
            node["p"] = alloc()
            T.ccx(hi["p"], lo["p"], node["p"])
        T.ccx(hi["p"], lo["g"], hi["g"])

    build(root)
    return root["g"]


def gt_tree(C, A, B, t, ctrl, block=4):
    """t ^= ctrl * [A > B], with A, B and ctrl restored: the carry out of A + ~B."""
    assert len(A) == len(B) >= 1 and ctrl not in (*A, *B, t)
    taken = []

    def alloc():
        (w,) = C.alloc(1)
        taken.append(w)
        return w

    T = Tape()
    for x in B:
        T.x(x)
    top = _carry_out(T, list(A), list(B), alloc, block)
    T.play(C)
    C.ccx(ctrl, top, t)
    T.play(C, reverse=True)
    C.free(taken[::-1])


class PrefixTree:
    """AND(bits) on .top, and write(t): bits += t mod 2^len(bits), on clean scratch.

    Blocks of `block` bits each hold the chain q_k = AND(bits[lo .. lo+k]), and a balanced tree
    over the chain tops holds each subtree's AND. write sends the carry t down the tree
    (c_hi = c * AND(lo)), adds within each block (bit_k ^= cin * q_(k-1)), and clears the scratch
    in two halves. The old chains and tree clear against each block's new bits XOR its cin: a block
    whose carry in is 1 held all ones up to its first zero and now holds zeros there, so that XOR
    has the old chain, and a one-bit block's old bit. The carries clear against the same chains and tree rebuilt on the
    complemented new bits, because c * AND(old lo) = c * AND(~new lo): with c = 1 the lower half
    was all ones exactly when it is all zeros now.
    """

    def __init__(self, C, bits, block=4):
        self.C, self.bits = C, list(bits)
        self.taken, self.blocks, self.top = [], [], None
        h = len(self.bits)
        if h == 0:
            return
        self.chains, self.tree = Tape(), Tape()
        for lo in range(0, h, block):
            seg = self.bits[lo:lo + block]
            q = [seg[0]]
            for x in seg[1:]:
                t = self._alloc()
                self.chains.ccx(q[-1], x, t)
                q.append(t)
            self.blocks.append({"seg": seg, "q": q, "a": q[-1], "kids": None})
        level = self.blocks
        while len(level) > 1:
            nxt = []
            for i in range(0, len(level) - 1, 2):
                a = self._alloc()
                self.tree.ccx(level[i]["a"], level[i + 1]["a"], a)
                nxt.append({"a": a, "kids": (level[i], level[i + 1])})
            if len(level) % 2:
                nxt.append(level[-1])
            level = nxt
        self.root = level[0]
        self.top = self.root["a"]
        self.chains.play(C)
        self.tree.play(C)

    def _alloc(self):
        (w,) = self.C.alloc(1)
        self.taken.append(w)
        return w

    def release(self):
        pass

    def restore(self):
        pass

    def write(self, t, mid=None):
        """mid(z), when given, runs while z holds AND(~new bits), before t's last use: t * z is
        the carry out of the bits, t * AND(old bits)."""
        C = self.C
        if not self.bits:
            return
        carries = Tape()
        made = []

        def down(node, c):
            if node["kids"] is None:
                node["cin"] = c
                return
            lo, hi = node["kids"]
            ch = target(C, c, lo["a"])
            made.append(ch)
            carries.ccx(c, lo["a"], ch)
            C.ccx(c, lo["a"], ch)
            down(lo, c)
            down(hi, ch)

        down(self.root, t)
        for blk in self.blocks:
            seg, q, cin = blk["seg"], blk["q"], blk["cin"]
            for k in range(len(seg) - 1, 0, -1):
                C.ccx(cin, q[k - 1], seg[k])
            C.cx(cin, seg[0])
        for blk in self.blocks:
            for x in blk["seg"]:
                C.cx(blk["cin"], x)
        self.tree.play(C, reverse=True)
        self.chains.play(C, reverse=True)
        for blk in self.blocks:
            for x in blk["seg"]:
                C.cx(blk["cin"], x)
        for x in self.bits:
            C.x(x)
        self.chains.play(C)
        self.tree.play(C)
        if mid is not None:
            mid(self.top)
        carries.play(C, reverse=True)
        self.tree.play(C, reverse=True)
        self.chains.play(C, reverse=True)
        for x in self.bits:
            C.x(x)
        C.free(made[::-1])
        C.free(self.taken[::-1])


class LowAdd:
    """The carries of L + K, where K is the XOR of the masks whose selectors are set, and then
    either undo() or write(): L[1:] += K's carries and bits (bit 0 is left alone, as the tick's
    chain leaves it), the quiet run through quiet_write.

    L's bits a .. b-1 are the quiet run: no mask has a bit there, so the carry crosses it only when
    it is all ones, which mtop holds. K's other bits go on clean wires, and each stretch of them
    ripples in place as CDKM majorities in blocks of `block` bits, every block from its own clean
    carry-in. A block's products pi_k = p_lo ... p_(lo+k-1), p = L ^ K, turn its carry out into the
    true one: G ^= P * C. .top holds the carry out of L's top bit.
    """

    def __init__(self, C, L, selectors, masks, run, mtop, block):
        self.C, self.L = C, list(L)
        self.selectors, self.masks = list(selectors), list(masks)
        a, b = run
        w = len(self.L)
        self.a, self.b, self.mtop = a, b, mtop
        self.fan, self.up = Tape(), Tape()
        self.taken = []
        active = [i for i in range(w) if not a <= i < b]
        self.K = {i: self._alloc() for i in active}
        for i in active:
            for s in self._on(i):
                self.fan.cx(s, self.K[i])
        self.blocks = []
        spans = [(lo, min(lo + block, a)) for lo in range(0, a, block)]
        spans += [(lo, min(lo + block, w)) for lo in range(b, w, block)]
        L_, K = self.L, self.K
        for lo, hi in spans:
            z = self._alloc()
            if self._on(lo):                            # else K[lo] is 0, and so is its carry out
                self.up.cx(K[lo], L_[lo])
                self.up.cx(K[lo], z)
                self.up.ccx(z, L_[lo], K[lo])
            for i in range(lo + 1, hi):
                self.up.cx(K[i], L_[i])
                self.up.cx(K[i], K[i - 1])
                self.up.ccx(K[i - 1], L_[i], K[i])
            self.blocks.append({"lo": lo, "hi": hi, "z": z, "pi": None, "cin": None})
        self.pis, self.join = Tape(), Tape()
        self.skip = None
        carry = None
        for blk in self.blocks:
            lo, hi = blk["lo"], blk["hi"]
            if lo == b and b > a:
                carry = self._skip(carry)
            if carry is not None:
                pi = [L_[lo]]
                for i in range(lo + 1, hi):
                    t = self._alloc()
                    self.pis.ccx(pi[-1], L_[i], t)
                    pi.append(t)
                blk["pi"], blk["cin"] = pi, carry
                self.join.ccx(carry, pi[-1], K[hi - 1])
            carry = K[hi - 1]
        if b == w and b > a:
            carry = self._skip(carry)
        self.top = carry
        self.fan.play(C)
        self.up.play(C)
        self.pis.play(C)
        self.join.play(C)

    def _alloc(self):
        (w,) = self.C.alloc(1)
        self.taken.append(w)
        return w

    def _skip(self, carry):
        self.skip = self._alloc()
        self.join.ccx(carry, self.mtop, self.skip)
        return self.skip

    def _on(self, i):
        return [s for s, mask in zip(self.selectors, self.masks) if mask >> i & 1]

    def _free(self):
        self.C.free(self.taken[::-1])

    def undo(self):
        C = self.C
        self.join.play(C, reverse=True)
        self.pis.play(C, reverse=True)
        self.up.play(C, reverse=True)
        self.fan.play(C, reverse=True)
        self._free()

    def write(self, quiet_write):
        """quiet_write(c, mid) adds the carry wire c into the quiet run, and calls mid(z) once,
        before c's last use, with a wire z for which c * z = c * [the run was all ones]: the old
        mtop before the run changes, or AND(~run) after. The quiet run's write does not wait for
        the top stretch, whose first fall waits for whoever still reads .top."""
        C, L, K = self.C, self.L, self.K
        # Make each later block's state the single chain's: true inner carries, its products
        # cleared, and its carry-in wire traded for the wire below it.
        for blk in self.blocks:
            cin, pi, lo, hi = blk["cin"], blk["pi"], blk["lo"], blk["hi"]
            if cin is None:
                continue
            for k in range(1, hi - lo):
                C.ccx(cin, pi[k - 1], K[lo + k - 1])
        self.pis.play(C, reverse=True)
        for blk in self.blocks:
            if blk["cin"] is None:
                continue
            if self._on(blk["lo"]):
                C.cx(blk["z"], blk["cin"])
            for s in self._on(blk["lo"]):
                C.cx(s, blk["z"])
        first = self.blocks[0]

        def below(i):
            for blk in self.blocks:
                if blk["lo"] == i:
                    return blk["z"] if blk["cin"] is None else blk["cin"]
            return K[i - 1]

        def fall(i):
            c = below(i)
            C.ccx(c, L[i], K[i])
            C.cx(K[i], c)
            if i:
                C.cx(c, L[i])
            else:
                C.cx(K[i], L[i])

        def top(z):
            for i in reversed(range(self.b, len(L))):
                fall(i)
            C.ccx(K[self.a - 1], z, self.skip)

        if self.skip is not None:
            quiet_write(K[self.a - 1], top)
        else:
            for i in reversed(range(self.b, len(L))):
                fall(i)
        for i in reversed(range(self.a)):
            fall(i)
        self.fan.play(C, reverse=True)
        self._free()
