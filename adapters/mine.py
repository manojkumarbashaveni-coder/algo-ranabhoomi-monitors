"""An anytime local-search solver for the monitor placement problem.

The important bit here is that all decisions use ``dist_table``: candidate
locations are deliberately not confused with demand locations.
"""

import time

from adapter import Solver, dist_table


class MySolver(Solver):
    def solve(self, instance, submit_candidate):
        table = dist_table(instance)
        weights = instance.weights
        n, m, k = instance.size, len(instance.sites), instance.k
        started = time.perf_counter()
        # Leave a little room for the harness and for returning the final plan.
        deadline = started + 4.80

        def cost_and_values(sites):
            values = [min(row[s] for s in sites) for row in table]
            return sum(w * d for w, d in zip(weights, values)), values

        # A valid result right away matters to the early scoring checkpoints.
        # This small deterministic shuffle does not depend on Python's salted hash.
        seed = 2166136261
        for ch in instance.digest.encode("ascii"):
            seed = ((seed ^ ch) * 16777619) & 0xFFFFFFFF
        ids = list(range(m))
        for i in range(m - 1, 0, -1):
            seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF
            j = seed % (i + 1)
            ids[i], ids[j] = ids[j], ids[i]
        submit_candidate({"sites": ids[:k]})

        # Greedy construction.  Keeping the per-point best charged distance makes
        # this O(k*m*n), and it is normally complete well before the first check.
        chosen = []
        chosen_mask = [False] * m
        values = [instance.penalty] * n
        for _ in range(k):
            best_s, best_cost = -1, None
            for s in range(m):
                if chosen_mask[s]:
                    continue
                c = 0
                for i, row in enumerate(table):
                    d = row[s]
                    if d > values[i]:
                        d = values[i]
                    c += weights[i] * d
                if best_cost is None or c < best_cost:
                    best_s, best_cost = s, c
            chosen.append(best_s)
            chosen_mask[best_s] = True
            for i, row in enumerate(table):
                if row[best_s] < values[i]:
                    values[i] = row[best_s]

        current_cost = sum(w * d for w, d in zip(weights, values))
        global_sites, global_cost = list(chosen), current_cost
        submit_candidate({"sites": list(chosen)})

        def best_swap(sites, mask):
            """Return the best single exchange, or None.

            For each removal we precompute its exact residual assignment, then
            evaluate every closed candidate against it.  This is substantially
            stronger than a single pass because the caller repeats it to a local
            optimum.
            """
            incumbent = None
            incumbent_cost = None
            for pos in range(k):
                if time.perf_counter() >= deadline:
                    return None
                rest = sites[:pos] + sites[pos + 1:]
                base = [min(row[s] for s in rest) for row in table]
                for add in range(m):
                    if mask[add]:
                        continue
                    c = 0
                    for i, row in enumerate(table):
                        d = row[add]
                        if d > base[i]:
                            d = base[i]
                        c += weights[i] * d
                    if incumbent_cost is None or c < incumbent_cost:
                        incumbent_cost = c
                        incumbent = (pos, add)
            return incumbent_cost, incumbent

        # A compact deterministic RNG for perturbations/restarts.
        def next_rand():
            nonlocal seed
            seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF
            return seed

        def best_pair_repair(sites, mask):
            """Try a small large-neighborhood move: remove two, greedily add two.

            The two removed facilities are deliberately barred from the repair.
            Otherwise greedy often simply restores the original local optimum and
            this neighborhood adds no information beyond a 1-swap.  We enumerate
            every pair when k is small and sample pairs for the large-k family.
            """
            pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
            for i in range(len(pairs) - 1, 0, -1):
                j = next_rand() % (i + 1)
                pairs[i], pairs[j] = pairs[j], pairs[i]
            if k > 10:
                pairs = pairs[:48]
            best_cost, best_sites = None, None
            for left, right in pairs:
                if time.perf_counter() >= deadline:
                    break
                old_left, old_right = sites[left], sites[right]
                kept = [s for p, s in enumerate(sites) if p != left and p != right]
                base = [min(row[s] for s in kept) for row in table]

                first, first_cost = -1, None
                for add in range(m):
                    if mask[add] or add == old_left or add == old_right:
                        continue
                    c = 0
                    for i, row in enumerate(table):
                        d = row[add]
                        if d > base[i]:
                            d = base[i]
                        c += weights[i] * d
                    if first_cost is None or c < first_cost:
                        first, first_cost = add, c
                after_first = [min(base[i], row[first]) for i, row in enumerate(table)]

                second, second_cost = -1, None
                for add in range(m):
                    if mask[add] or add == old_left or add == old_right or add == first:
                        continue
                    c = 0
                    for i, row in enumerate(table):
                        d = row[add]
                        if d > after_first[i]:
                            d = after_first[i]
                        c += weights[i] * d
                    if second_cost is None or c < second_cost:
                        second, second_cost = add, c
                if best_cost is None or second_cost < best_cost:
                    candidate = list(sites)
                    candidate[left], candidate[right] = first, second
                    best_cost, best_sites = second_cost, candidate
            return best_cost, best_sites

        # Search a basin to a 1-swap local optimum, then remove a few monitors
        # and greedily repair.  The latter often discovers useful coordinated
        # exchanges on clustered instances where a pure 1-swap search is stuck.
        sites, mask = list(chosen), chosen_mask[:]
        while time.perf_counter() < deadline:
            while time.perf_counter() < deadline:
                move = best_swap(sites, mask)
                if move is None:
                    break
                new_cost, pair = move
                old_cost, _ = cost_and_values(sites)
                if pair is None or new_cost >= old_cost:
                    break
                pos, add = pair
                mask[sites[pos]] = False
                sites[pos] = add
                mask[add] = True
                if new_cost < global_cost:
                    global_cost, global_sites = new_cost, list(sites)
                    submit_candidate({"sites": list(global_sites)})
            if time.perf_counter() >= deadline:
                break

            pair_cost, pair_sites = best_pair_repair(sites, mask)
            if pair_sites is not None and pair_cost < global_cost:
                global_cost, global_sites = pair_cost, pair_sites
                sites = list(pair_sites)
                mask = [False] * m
                for s in sites:
                    mask[s] = True
                submit_candidate({"sites": list(global_sites)})
                # Newly opened sites can expose further improving 1-swaps.
                continue
            if time.perf_counter() >= deadline:
                break

            # Start the next basin from the best known plan, replacing 2--4
            # randomly selected monitors with a greedy repair.  It is fine for
            # this trial plan to be worse; only global improvements are submitted.
            sites = list(global_sites)
            mask = [False] * m
            for s in sites:
                mask[s] = True
            destroy = min(k - 1, 2 + next_rand() % 3)
            positions = list(range(k))
            for i in range(destroy):
                j = i + next_rand() % (k - i)
                positions[i], positions[j] = positions[j], positions[i]
            removed = set(positions[:destroy])
            kept = [s for i, s in enumerate(sites) if i not in removed]
            for i in removed:
                mask[sites[i]] = False
            base = [min((row[s] for s in kept), default=instance.penalty) for row in table]
            sites = kept
            for _ in range(destroy):
                if time.perf_counter() >= deadline:
                    break
                best_s, best_c = -1, None
                for s in range(m):
                    if mask[s]:
                        continue
                    c = sum(w * min(base[i], row[s])
                            for i, (w, row) in enumerate(zip(weights, table)))
                    if best_c is None or c < best_c:
                        best_s, best_c = s, c
                sites.append(best_s)
                mask[best_s] = True
                for i, row in enumerate(table):
                    if row[best_s] < base[i]:
                        base[i] = row[best_s]
            if len(sites) != k:
                break

        return {"sites": global_sites}
