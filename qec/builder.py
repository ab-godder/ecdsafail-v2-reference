"""Reversible circuit builder over the v2 arithmetic IR's permutation subset.

An op is four int32 values (opcode, a, b, c). Unused operand slots repeat the
last used wire so wire remapping can treat every slot uniformly.

In hierarchical mode a body is its primitive op array plus a list of Calls,
each run just before primitive op `pos`. `wires` holds the body's own wire for
each of sorted(block.wires), so applying a block never remaps the calls inside
it.
"""
from array import array
from typing import NamedTuple

X, CX, CCX, SWAP = 0, 1, 2, 3
CCX_TICKS = 8
OPNAMES = {X: "X", CX: "CX", CCX: "CCX", SWAP: "SWAP"}


class Call(NamedTuple):
    pos: int
    block: "Block"
    wires: list
    inverse: bool


class Block:
    """A recorded op range that can be replayed forward or reversed onto any
    register layout of the same shape. Wires outside the layout are zero at the
    block's start and end and are mapped onto fresh ancillas."""

    def __init__(self, ops, flat_in, flat_out, calls=()):
        self.ops = ops
        self.flat_in = flat_in
        self.flat_out = flat_out
        self.calls = calls
        self.wires = set(ops[1::4]) | set(ops[2::4]) | set(ops[3::4])
        for call in calls:
            self.wires.update(call.wires)


def flatten(regs):
    return [w for r in regs for w in r]


def unflatten_into(regs, flat):
    i = 0
    for r in regs:
        r[:] = flat[i:i + len(r)]
        i += len(r)


class Circuit:
    def __init__(self, hierarchical=False):
        self.prims = array("i")
        self.calls = [] if hierarchical else None
        self.num_wires = 0
        self.free_stack = []
        self.ready = {}

    # allocation

    def alloc(self, k=1):
        out = []
        for _ in range(k):
            if self.free_stack:
                out.append(self.free_stack.pop())
            else:
                out.append(self.num_wires)
                self.num_wires += 1
        return out

    def free(self, wires):
        self.free_stack.extend(reversed(wires))

    # gates

    def x(self, a):
        self.prims.extend((X, a, a, a))
        self.ready[a] = self.ready.get(a, 0) + 1

    def _touch(self, wires, ticks):
        ready = self.ready
        t = max(ready.get(w, 0) for w in wires) + ticks
        for w in wires:
            ready[w] = t

    def cx(self, a, b):
        self.prims.extend((CX, a, b, b))
        self._touch((a, b), 1)

    def ccx(self, a, b, c):
        if a == b:
            self.cx(a, c)
            return
        # Lowering lets a chain through the first control overlap more, so the
        # control that comes free last goes first. The order changes no logic.
        ready = self.ready
        if ready.get(b, 0) > ready.get(a, 0):
            a, b = b, a
        self.prims.extend((CCX, a, b, c))
        self._touch((a, b, c), CCX_TICKS)

    def swap(self, a, b):
        self.prims.extend((SWAP, a, b, b))
        self._touch((a, b), 3)

    @property
    def ops(self):
        if self.calls is not None:
            raise TypeError("a hierarchical Circuit holds only its top-level primitives; expand it with qec.flatten")
        return self.prims

    def __len__(self):
        return len(self.ops) // 4

    # blocks

    def record(self, fn, regs):
        """Run fn(*regs) without emitting it. Returns a Block and leaves the
        circuit, the allocator, and regs exactly as they were."""
        flat_in = flatten(regs)
        start = len(self.prims)
        first_call = len(self.calls) if self.calls is not None else 0
        saved_stack, saved_hw, saved_ready = list(self.free_stack), self.num_wires, self.ready
        self.ready = {}
        fn(*regs)
        self.ready = saved_ready
        flat_out = flatten(regs)
        ops = self.prims[start:]
        del self.prims[start:]
        calls = ()
        if self.calls is not None:
            calls = [call._replace(pos=call.pos - start // 4) for call in self.calls[first_call:]]
            del self.calls[first_call:]
        self.free_stack, self.num_wires = saved_stack, saved_hw
        unflatten_into(regs, flat_in)
        return Block(ops, flat_in, flat_out, calls)

    def apply(self, block, regs, inverse=False):
        src_live, dst_live = (block.flat_out, block.flat_in) if inverse else (block.flat_in, block.flat_out)
        cur = flatten(regs)
        assert len(cur) == len(src_live)
        if not block.wires:
            unflatten_into(regs, [dict(zip(src_live, cur))[w] for w in dst_live])
            return
        size = max(max(block.wires), max(src_live, default=0), max(dst_live, default=0)) + 1
        table = [-1] * size
        for w, v in zip(src_live, cur):
            table[w] = v
        scratch = sorted(w for w in block.wires if table[w] < 0)
        assert set(dst_live) <= set(src_live) | block.wires
        fresh = self.alloc(len(scratch))
        for w, v in zip(scratch, fresh):
            table[w] = v
        if self.calls is None:
            self._replay(block, table, inverse)
        else:
            wires = [table[w] for w in sorted(block.wires)]
            if len(set(wires)) != len(wires):
                raise ValueError(f"call wires are not distinct: {wires}")
            self.calls.append(Call(len(self.prims) // 4, block, wires, inverse))
        new_live = [table[w] for w in dst_live]
        live = set(new_live)
        self.free([w for w in list(cur) + fresh if w not in live])
        unflatten_into(regs, new_live)

    def _replay(self, block, table, inverse):
        src = block.ops
        out = array("i", bytes(len(src) * src.itemsize))
        if inverse:
            rev = src[::-1]
            out[0::4] = rev[3::4]
            slots = (rev[2::4], rev[1::4], rev[0::4])
        else:
            out[0::4] = src[0::4]
            slots = (src[1::4], src[2::4], src[3::4])
        get = table.__getitem__
        for k, slot in enumerate(slots, 1):
            out[k::4] = array("i", map(get, slot))
        self.prims.extend(out)

    def record_shape(self, fn, sizes):
        """Record fn on fresh placeholder registers of the given sizes."""
        saved_stack, saved_hw = list(self.free_stack), self.num_wires
        regs = [self.alloc(k) for k in sizes]
        block = self.record(fn, regs)
        self.free_stack, self.num_wires = saved_stack, saved_hw
        return block

    def run(self, fn, regs):
        block = self.record(fn, regs)
        self.apply(block, regs)
        return block

