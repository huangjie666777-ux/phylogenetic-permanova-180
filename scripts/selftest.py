"""Self-test: run a placement in-process and sanity-check the results."""

import json
import pathlib

from fastapi.testclient import TestClient

from app.tree import parse_annotated_newick, quote_name
from app.community import pair_distance
from app.main import app

ROOT = pathlib.Path(__file__).resolve().parent.parent
EX = ROOT / "examples"


def main() -> None:
    client = TestClient(app)
    payload = {
        "reference_fasta": (EX / "reference.fasta").read_text(),
        "query_fasta": (EX / "queries.fasta").read_text(),
        "newick": (EX / "tree.nwk").read_text(),
    }
    resp = client.post("/place", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "q3" in data["unplaceable"], data
    for qid in ("q1", "q2"):
        places = data["placements"][qid]
        assert len(places) == data["edge_count"] == 7
        weights = [p["like_weight_ratio"] for p in places]
        assert abs(sum(weights) - 1.0) < 1e-9, weights
        lls = [p["likelihood"] for p in places]
        assert lls == sorted(lls, reverse=True)
        for p in places:
            assert 0.0 <= p["pendant_length"] <= 2.0
            assert p["distal_length"] >= 0.0
    # q1 matches refA exactly: its best edge should touch refA's pendant edge
    best_q1 = data["placements"]["q1"][0]
    print("q1 best edge:", best_q1)

    dl = client.get("/place/" + data["job_id"] + "/jplace")
    assert dl.status_code == 200
    doc = json.loads(dl.content)
    assert doc["version"] == 3
    assert doc["fields"] == ["edge_num", "likelihood", "like_weight_ratio",
                             "distal_length", "pendant_length"]
    assert "{" in doc["tree"] and "}" in doc["tree"]
    assert len(doc["placements"]) == 2

    # rejection path: illegal character must be located and reject the batch
    bad = dict(payload)
    bad["query_fasta"] = ">qX\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGX\n"
    resp = client.post("/place", json=bad)
    assert resp.status_code == 422
    problems = resp.json()["detail"]["problems"]
    assert any("qX" in p and "column 40" in p for p in problems), problems
    print("rejection OK:", problems[0])

    # RNA uracil must not be silently accepted as DNA
    bad_u = dict(payload)
    bad_u["query_fasta"] = ">qU\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGU\n"
    resp = client.post("/place", json=bad_u)
    assert resp.status_code == 422
    problems_u = resp.json()["detail"]["problems"]
    assert any("qU" in p and "column 40" in p and "'U'" in p
               for p in problems_u), problems_u
    print("U rejection OK:", problems_u[0])

    # ---- community KR distances on the stored job ----
    tree_text = doc["tree"]
    geo = parse_annotated_newick(tree_text)
    assert sorted(geo["lengths"]) == list(range(data["edge_count"]))
    assert all(v > 0 for v in geo["lengths"].values())

    # Newick labels containing commas/parens/quotes must survive quoting
    weird = ("(('we,ird':0.1{0},'x''y':0.12{1}):0.08{2},"
             "('r c':0.11{3},d:0.09{4}):0.07{5},e:0.2{6});")
    g2 = parse_annotated_newick(weird)
    assert g2["leaf"][0] == "we,ird", g2["leaf"]
    assert g2["leaf"][1] == "x'y", g2["leaf"]
    assert g2["lengths"][6] == 0.2
    assert quote_name("a'b,c") == "'a''b,c'"

    job_id = data["job_id"]
    comp = {
        "samples": [
            {"sample_id": "site1", "job_id": job_id,
             "counts": {"q1": 8, "q3": 3}},
            {"sample_id": "site2", "job_id": job_id,
             "counts": {"q1": 2, "q2": 6, "q3": 1}},
            {"sample_id": "site3", "job_id": job_id,
             "counts": {"q1": 5, "q2": 3}},
        ]
    }
    resp = client.post("/community/distances", json=comp)
    assert resp.status_code == 200, resp.text
    cdata = resp.json()
    assert cdata["sample_ids"] == ["site1", "site2", "site3"]
    mat = cdata["matrix"]
    assert len(mat) == 3 and all(len(row) == 3 for row in mat)
    for i in range(3):
        assert mat[i][i] == 0.0
        for j in range(3):
            assert mat[i][j] == mat[j][i]
            assert mat[i][j] >= 0.0
    assert mat[0][1] > 0.0 and mat[0][2] > 0.0
    for rep in cdata["samples"]:
        assert rep["valid_total"] > 0
    assert cdata["samples"][0]["excluded"] == [
        {"id": "q3", "count": 3,
         "reason": data["unplaceable"]["q3"]}]
    for pair in cdata["pairs"]:
        i, j = pair["i"], pair["j"]
        edges = [c["edge_num"] for c in pair["edge_contributions"]]
        assert edges == list(range(data["edge_count"]))
        contribs = [c["contribution"] for c in pair["edge_contributions"]]
        assert all(c >= 0.0 for c in contribs)
        assert abs(sum(contribs) - pair["distance"]) < 1e-9
        assert abs(pair["distance"] - mat[i][j]) < 1e-12
    print("community matrix:", [[round(v, 6) for v in row] for row in mat])

    # identical compositions (up to scaling) -> distance 0
    same = {"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 3, "q2": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 6, "q2": 2}}]}
    resp = client.post("/community/distances", json=same)
    assert resp.status_code == 200, resp.text
    assert resp.json()["matrix"][0][1] == 0.0

    def reject(body, fragment):
        r = client.post("/community/distances", json=body)
        assert r.status_code == 422, r.text
        assert fragment in r.json()["detail"]["problem"], r.text
        print("community rejected:", r.json()["detail"]["problem"])

    reject({"samples": [{"sample_id": "a", "job_id": job_id,
                         "counts": {"q1": 1}}]}, "2..10 samples")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}}]},
        "duplicate sample_id")
    reject({"samples": [
        {"sample_id": "a", "job_id": "nope", "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "unknown job_id")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"qX": 1}}]},
        "unknown query id")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": -1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "non-negative integer")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1.5}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "non-negative integer")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q3": 2}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "valid total must be positive")

    # a zero count for an *unknown* query id must be located and rejected,
    # not silently skipped
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1, "qZ": 0}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "unknown query id 'qZ'")

    # deep-tree regression: distal mass must be accumulated through *immediate*
    # child edges only, otherwise mass is double-counted at every ancestor
    deep = ("((((L:0.1{0},x:0.1{1}):0.1{2},y:0.12{3}):0.13{4},"
            "z:0.14{5}):0.15{6},w:0.16{7});")
    dgeo = parse_annotated_newick(deep)
    at_tip = {0: [(0.1, 1.0)]}
    x_tip = {1: [(0.1, 1.0)]}
    empty = {}
    d_lx = pair_distance(dgeo, at_tip, x_tip)["distance"]
    d_ll = pair_distance(dgeo, at_tip, at_tip)["distance"]
    d_lw = pair_distance(dgeo, at_tip,
                         {7: [(0.16, 1.0)]})["distance"]
    assert abs(d_ll) < 1e-12, d_ll
    assert abs(d_lx - 0.2) < 1e-9, (d_lx, "L<->x must be 0.1+0.1")
    assert abs(d_lw - (0.1 + 0.1 + 0.13 + 0.15 + 0.16)) < 1e-9, d_lw
    print("deep-tree KR OK:", d_lx, d_lw)

    # ---- batch-constrained PERMANOVA over the KR matrix ----
    from app.permanova import run_permanova, sums_of_squares, pseudo_f

    # within-group compositions differ slightly so SSW > 0 (SSW == 0 is a
    # rejection, tested separately below)
    comps = [{"q1": 8, "q2": 2}, {"q1": 7, "q2": 3},
             {"q1": 2, "q2": 8}, {"q1": 3, "q2": 7}]
    design4 = {"samples": [
        {"sample_id": "a1", "job_id": job_id, "counts": dict(comps[0]),
         "group": "A"},
        {"sample_id": "a2", "job_id": job_id, "counts": dict(comps[1]),
         "group": "A"},
        {"sample_id": "b1", "job_id": job_id, "counts": dict(comps[2]),
         "group": "B"},
        {"sample_id": "b2", "job_id": job_id, "counts": dict(comps[3]),
         "group": "B"}],
        "permutations": 99, "seed": 123}

    def permanova(body):
        return client.post("/community/permanova", json=body)

    resp = permanova(design4)
    assert resp.status_code == 200, resp.text
    pdata = resp.json()
    assert pdata["sample_order"] == ["a1", "a2", "b1", "b2"]
    assert pdata["groups"] == ["A", "A", "B", "B"]
    assert pdata["batches"] is None
    assert pdata["degrees_of_freedom"] == {"between": 1, "within": 2}
    soq = pdata["sum_of_squares"]
    assert abs(soq["between"] - (soq["total"] - soq["within"])) < 1e-12
    n = 4
    f_manual = ((soq["between"] / 1) / (soq["within"] / 2))
    assert abs(pdata["pseudo_F"] - f_manual) < 1e-9
    assert abs(pdata["R2"] - soq["between"] / soq["total"]) < 1e-12
    # 4 labels with fixed sizes 2/2 give C(4,2)=6 distinct assignments;
    # requested 99 >= 6 -> exact enumeration including the observed labeling
    assert pdata["permutations"] == 6, pdata["permutations"]
    assert "exact" in pdata["method"]
    assert abs((pdata["p_value"] * 6) - round(pdata["p_value"] * 6)) < 1e-9
    assert pdata["p_value"] >= 1.0 / 6.0 - 1e-12
    assert pdata["seed"] is None
    print("PERMANOVA exact:", {k: pdata[k] for k in
          ("pseudo_F", "R2", "p_value", "permutations")})

    # same request/seed reproduces; random mode needs fewer draws than the
    # 6? no — 6 is tiny, force random path by enlarging the design to n=8
    # (C(8,4)=70 > requested 20)
    design8 = {"samples": [
        {"sample_id": "a" + str(i + 1), "job_id": job_id,
         "counts": dict(comps[i % 2]), "group": "A"} for i in range(4)] + [
        {"sample_id": "b" + str(i + 1), "job_id": job_id,
         "counts": dict(comps[2 + i % 2]), "group": "B"} for i in range(4)],
        "permutations": 20, "seed": 123}
    r1 = permanova(design8).json()
    r2 = permanova(design8).json()
    assert r1 == r2, "same request + seed must reproduce"
    assert r1["permutations"] == 20 and "random" in r1["method"]
    assert r1["seed"] == 123
    assert abs(r1["p_value"] - round(r1["p_value"] * 21) / 21) < 1e-9
    design8["seed"] = 124
    r3 = permanova(design8).json()
    assert r3["seed"] == 124
    print("PERMANOVA random reproducible p =", r1["p_value"])

    # batches: balanced swaps within batches (one A + one B each) are legal
    batched = {"samples": [
        {"sample_id": "a1", "job_id": job_id, "counts": dict(comps[0]),
         "group": "A", "batch": "s1"},
        {"sample_id": "b1", "job_id": job_id, "counts": dict(comps[2]),
         "group": "B", "batch": "s1"},
        {"sample_id": "a2", "job_id": job_id, "counts": dict(comps[1]),
         "group": "A", "batch": "s2"},
        {"sample_id": "b2", "job_id": job_id, "counts": dict(comps[3]),
         "group": "B", "batch": "s2"}],
        "permutations": 999, "seed": 7}
    rb = permanova(batched)
    assert rb.status_code == 200, rb.text
    bdata = rb.json()
    # 2 plans per batch independently -> 4 distinct legal labelings total
    assert bdata["permutations"] == 4, bdata["permutations"]
    assert bdata["batches"] == ["s1", "s1", "s2", "s2"]
    print("batch PERMANOVA p =", bdata["p_value"])

    def reject_pm(body, fragment):
        r = permanova(body)
        assert r.status_code == 422, r.text
        assert fragment in r.json()["detail"]["problem"], r.text
        print("permanova rejected:", r.json()["detail"]["problem"])

    # batch completely confounded with group: no within-batch swap exists,
    # and cross-batch shuffling must not be treated as a valid test
    confounded = json.loads(json.dumps(batched))
    confounded["samples"][2]["batch"] = "s1"
    confounded["samples"][1]["batch"] = "s2"
    reject_pm(confounded, "no legal relabeling")

    reject_pm({"samples": design4["samples"][:3]}, "4..10 samples")
    one_group = json.loads(json.dumps(design4))
    for s in one_group["samples"]:
        s["group"] = "X"
    reject_pm(one_group, "at least two groups")
    singleton = json.loads(json.dumps(design4))
    singleton["samples"][3]["group"] = "C"
    reject_pm(singleton, "at least 2 samples")
    no_group = json.loads(json.dumps(design4))
    del no_group["samples"][0]["group"]
    reject_pm(no_group, ".group must be a non-empty string")
    partial_batch = json.loads(json.dumps(batched))
    del partial_batch["samples"][0]["batch"]
    reject_pm(partial_batch, "every sample must")
    bad_perm = json.loads(json.dumps(design4))
    bad_perm["permutations"] = 10000
    reject_pm(bad_perm, "permutations")
    bad_perm["permutations"] = 0
    reject_pm(bad_perm, "permutations")
    bad_seed = json.loads(json.dumps(design4))
    bad_seed["seed"] = 1.5
    reject_pm(bad_seed, "'seed' must be an integer")

    # SST == 0: all samples identical -> reject, never a degenerate test
    identical = json.loads(json.dumps(design4))
    for s in identical["samples"]:
        s["counts"] = {"q1": 1}
    reject_pm(identical, "total sum of squares is 0")

    # SSW == 0: within each group samples coincide exactly -> reject
    ssw0 = json.loads(json.dumps(design4))
    ssw0["samples"][1]["counts"] = dict(ssw0["samples"][0]["counts"])
    ssw0["samples"][3]["counts"] = dict(ssw0["samples"][2]["counts"])
    reject_pm(ssw0, "within-group sum of squares is 0")

    # statistics cross-check: hand-computed sums of squares from the matrix
    resp4 = permanova(design4).json()
    dm = client.post("/community/distances", json={
        "samples": [{k: v for k, v in s.items() if k in ("sample_id",
                     "job_id", "counts")} for s in design4["samples"]]}).json()
    mat = dm["matrix"]
    labs = ["A", "A", "B", "B"]
    d2 = [[mat[i][j] ** 2 for j in range(4)] for i in range(4)]
    sst_h = sum(d2[i][j] for i in range(4) for j in range(i + 1, 4)) / 4
    ssw_h = (d2[0][1] / 2) + (d2[2][3] / 2)
    assert abs(resp4["sum_of_squares"]["total"] - sst_h) < 1e-12
    assert abs(resp4["sum_of_squares"]["within"] - ssw_h) < 1e-12

    # a second job on a different reference tree must be rejected pairwise
    other = dict(payload)
    other["newick"] = ("((refA:0.1,refB:0.12):0.08,"
                       "(refC:0.11,refD:0.09):0.07,refE:0.31);")
    r2 = client.post("/place", json=other)
    assert r2.status_code == 200, r2.text
    job2 = r2.json()["job_id"]
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job2, "counts": {"q1": 1}}]},
        "different annotated reference tree")
    print("SELFTEST PASSED")


if __name__ == "__main__":
    main()
