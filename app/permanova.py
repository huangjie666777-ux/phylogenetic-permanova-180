"""Batch-constrained PERMANOVA over stored KR distance matrices.

Only group labels are permuted; the distance matrix and every group's size
stay fixed.  When batches are supplied, labels are swapped *within* each
batch so that the number of slots per group inside every batch is preserved
(no cross-batch shuffling).

Sums of squares follow McArdle & Anderson (2001) on squared distances:
    SST  = sum_{i<j} d_ij^2 / n
    SSW  = sum_g (1/n_g) * sum_{i<j, both in g} d_ij^2
    SSB  = SST - SSW
    F    = (SSB/(a-1)) / (SSW/(n-a))
    R2   = SSB / SST
"""

from __future__ import annotations

import random

from .community import CommunityError


def validate_permanova_request(payload) -> list[dict]:
    if not isinstance(payload, dict) or "samples" not in payload:
        raise CommunityError("request body must contain 'samples'")
    raw_samples = payload["samples"]
    if not isinstance(raw_samples, list):
        raise CommunityError("'samples' must be a list")
    if not 4 <= len(raw_samples) <= 10:
        raise CommunityError("expected 4..10 samples, got "
                             + str(len(raw_samples)))

    seen: set[str] = set()
    have_batch: list[bool] = []
    for idx, sample in enumerate(raw_samples):
        where = "sample[" + str(idx) + "]"
        if not isinstance(sample, dict):
            raise CommunityError(where + " must be an object")
        sid = sample.get("sample_id")
        if not isinstance(sid, str) or not sid:
            raise CommunityError(where + ".sample_id must be a non-empty string")
        if sid in seen:
            raise CommunityError("duplicate sample_id '" + sid + "'")
        seen.add(sid)
        job_id = sample.get("job_id")
        if not isinstance(job_id, str) or not job_id:
            raise CommunityError(where + " ('" + sid + "').job_id must be a "
                                 "non-empty string")
        counts = sample.get("counts")
        if not isinstance(counts, dict):
            raise CommunityError(where + " ('" + sid + "').counts must be an "
                                 "object mapping query id to read count")
        for qid, count in counts.items():
            if not isinstance(qid, str) or not qid:
                raise CommunityError(where + " ('" + sid + "'): query id keys "
                                     "must be non-empty strings")
            if isinstance(count, bool) or not isinstance(count, int) or count < 0:
                raise CommunityError(where + " ('" + sid + "') query '" + qid +
                                     "': count must be a non-negative integer, "
                                     "got " + repr(count))
        group = sample.get("group")
        if not isinstance(group, str) or not group:
            raise CommunityError(where + " ('" + sid + "').group must be a "
                                 "non-empty string")
        batch = sample.get("batch")
        if batch is not None and (not isinstance(batch, str) or not batch):
            raise CommunityError(where + " ('" + sid + "').batch must be a "
                                 "non-empty string when provided")
        have_batch.append(batch is not None)

    if any(have_batch) and not all(have_batch):
        missing = [raw_samples[i]["sample_id"] for i, h in enumerate(have_batch)
                   if not h]
        raise CommunityError("batch is optional, but once any sample supplies "
                             "it every sample must; missing on: "
                             + ", ".join(missing))

    groups = [s["group"] for s in raw_samples]
    sizes: dict[str, int] = {}
    for g in groups:
        sizes[g] = sizes.get(g, 0) + 1
    if len(sizes) < 2:
        raise CommunityError("at least two groups required, got one group ('"
                             + groups[0] + "')")
    small = sorted(g for g, k in sizes.items() if k < 2)
    if small:
        raise CommunityError("every group needs at least 2 samples; group(s) "
                             + ", ".join(repr(g) for g in small) + " have "
                             "fewer than 2")

    for key, lo, hi in (("permutations", 1, 9999), ("seed", None, None)):
        if key in payload and payload[key] is not None:
            val = payload[key]
            if isinstance(val, bool) or not isinstance(val, int):
                raise CommunityError("'" + key + "' must be an integer, got "
                                     + repr(val))
            if lo is not None and not lo <= val <= hi:
                raise CommunityError("'" + key + "' must be in " + str(lo)
                                     + ".." + str(hi) + ", got " + str(val))
    return raw_samples


def sums_of_squares(d2, labels: list[str], group_names: list[str]):
    """Return (SST, SSW, SSB) from squared distances d2[i][j]."""
    n = len(labels)
    members: dict[str, list[int]] = {g: [] for g in group_names}
    for i, lab in enumerate(labels):
        members[lab].append(i)
    sst = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            sst += d2[i][j]
    sst /= n
    ssw = 0.0
    for g in group_names:
        idx = members[g]
        pair_sum = 0.0
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                pair_sum += d2[idx[a]][idx[b]]
        ssw += pair_sum / len(idx)
    return sst, ssw, sst - ssw


def pseudo_f(sst: float, ssw: float, n: int, a: int) -> tuple[float, float]:
    ssb = sst - ssw
    df_between = a - 1
    df_within = n - a
    f_stat = (ssb / df_between) / (ssw / df_within)
    return f_stat, ssb / sst


def _batch_plans(positions: list[int], group_order: list[str],
                 counts: dict[str, int]):
    """All distinct within-batch label vectors with fixed per-group counts.

    Each plan is a tuple aligned to the batch positions sorted ascending.
    """
    positions = sorted(positions)
    index = {p: k for k, p in enumerate(positions)}
    vec: list[str | None] = [None] * len(positions)

    def fill(k: int):
        if k == len(group_order):
            yield tuple(vec)
            return
        g = group_order[k]
        need = counts[g]
        free = [p for p in positions if vec[index[p]] is None]
        chosen: list[tuple[int, ...]] = []

        def combos(start: int, picked: list[int]) -> None:
            if len(picked) == need:
                chosen.append(tuple(picked))
                return
            for i in range(start, len(free)):
                combos(i + 1, picked + [free[i]])

        combos(0, [])
        for combo in chosen:
            for p in combo:
                vec[index[p]] = g
            yield from fill(k + 1)
            for p in combo:
                vec[index[p]] = None

    yield from fill(0)


def legal_labelings(labels: list[str], batches: list[str] | None):
    """Enumerate distinct label vectors via per-batch swaps only.

    The original labeling is included.  Raises CommunityError when no legal
    swap can change the grouping.
    """
    n = len(labels)
    batch_ids = ["__all__"] * n if batches is None else batches
    batch_names = sorted(set(batch_ids))
    group_names = sorted(set(labels))

    positions_by_batch: dict[str, list[int]] = {}
    for i, b in enumerate(batch_ids):
        positions_by_batch.setdefault(b, []).append(i)

    plans_by_batch: dict[str, list[tuple]] = {}
    for b in batch_names:
        pos = positions_by_batch[b]
        counts = {g: 0 for g in group_names}
        for p in pos:
            counts[labels[p]] += 1
        plans_by_batch[b] = list(_batch_plans(pos, group_names, counts))

    total = 1
    for plans in plans_by_batch.values():
        total *= len(plans)
    if total <= 1:
        raise CommunityError(
            "no legal relabeling can change the grouping: within every batch "
            "the group slots are fixed (each batch holds one group on every "
            "slot pattern); cross-batch shuffling is not a valid test")

    combined: list[tuple] = []
    current: list[str | None] = [None] * n
    batch_positions = {b: sorted(positions_by_batch[b]) for b in batch_names}

    def assemble(k: int) -> None:
        if k == len(batch_names):
            combined.append(tuple(current))
            return
        b = batch_names[k]
        pos = batch_positions[b]
        for plan in plans_by_batch[b]:
            for p, lab in zip(pos, plan):
                current[p] = lab
            assemble(k + 1)

    assemble(0)
    return combined


def run_permanova(matrix: list[list[float]], samples: list[dict],
                  permutations: int, seed: int) -> dict:
    n = len(samples)
    labels = [s["group"] for s in samples]
    batches_in = [s.get("batch") for s in samples]
    batches = batches_in if any(b is not None for b in batches_in) else None
    group_names = sorted(set(labels))
    a = len(group_names)

    d2 = [[matrix[i][j] ** 2 for j in range(n)] for i in range(n)]

    sst, ssw, _ = sums_of_squares(d2, labels, group_names)
    if sst <= 0.0:
        raise CommunityError("total sum of squares is 0: all pairwise "
                             "squared distances are 0, so the design cannot "
                             "be tested")
    if ssw <= 0.0:
        raise CommunityError("within-group sum of squares is 0: samples "
                             "within every group coincide exactly, so "
                             "pseudo-F is not finite; redesign the comparison")

    f_obs, r2 = pseudo_f(sst, ssw, n, a)
    labelings = legal_labelings(labels, batches)

    if len(labelings) <= permutations:
        method = "exact enumeration of all distinct legal relabelings"
        hits = 0
        for vec in labelings:
            t, w, _ = sums_of_squares(d2, list(vec), group_names)
            f_stat, _ = pseudo_f(t, w, n, a)
            if f_stat + 1e-12 >= f_obs:
                hits += 1
        p_value = hits / len(labelings)
        draws = len(labelings)
        used_seed = None
    else:
        method = ("random permutations (group labels swapped within batches, "
                  "with replacement)")
        rng = random.Random(seed)
        hits = 0
        for _ in range(permutations):
            vec = list(labelings[rng.randrange(len(labelings))])
            t, w, _ = sums_of_squares(d2, vec, group_names)
            f_stat, _ = pseudo_f(t, w, n, a)
            if f_stat + 1e-12 >= f_obs:
                hits += 1
        p_value = (hits + 1) / (permutations + 1)
        draws = permutations
        used_seed = seed

    return {
        "sample_order": [s["sample_id"] for s in samples],
        "groups": labels,
        "batches": batches,
        "sum_of_squares": {"total": sst, "within": ssw,
                           "between": sst - ssw},
        "degrees_of_freedom": {"between": a - 1, "within": n - a},
        "pseudo_F": f_obs,
        "R2": r2,
        "p_value": p_value,
        "method": method,
        "permutations": draws,
        "seed": used_seed,
    }

