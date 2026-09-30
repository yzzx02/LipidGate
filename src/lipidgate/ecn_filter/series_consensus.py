"""Observed-RT consensus and ordered homolog-series support, without moving peaks."""

from itertools import combinations

import numpy as np

from .workflow import _Curve


def _predict(model, x):
    return np.polynomial.polynomial.polyval(x, model[0].coefficients)


def _consensus(group, tolerance):
    # Repeated spectra at one retention cluster must not win by repetition.
    pool = group.copy()
    pool["_bin"] = (pool._rt_minutes / 0.05).round()
    pool = pool.sort_values(["_score", "_rank"], ascending=[False, True], kind="stable")
    pool = (
        pool.drop_duplicates(["total_C", "_bin"]).groupby("total_C", sort=True).head(20)
    )
    x = pool.total_C.to_numpy(float)
    y = pool._rt_minutes.to_numpy(float)
    scores = pool._score.fillna(0).to_numpy(float)
    carbons = np.unique(x)
    if len(carbons) < 3:
        return None
    slots = [np.flatnonzero(x == c) for c in carbons]
    rng = np.random.default_rng(20260914)
    hypotheses = []
    for degree, trials in ((1, 350), (2, 150)):
        if degree == 2 and len(carbons) < 6:
            continue
        for _ in range(trials):
            chosen = rng.choice(len(carbons), degree + 1, replace=False)
            selected = [rng.choice(slots[i]) for i in chosen]
            coef = np.polynomial.polynomial.polyfit(x[selected], y[selected], degree)
            grid = np.linspace(x.min(), x.max(), 64)
            if (
                np.min(
                    np.polynomial.polynomial.polyval(
                        grid, np.polynomial.polynomial.polyder(coef)
                    )
                )
                <= 0
            ):
                continue
            residual = abs(y - np.polynomial.polynomial.polyval(x, coef))
            nearest = np.array(
                [idx[np.lexsort((-scores[idx], residual[idx]))[0]] for idx in slots]
            )
            good = nearest[residual[nearest] <= tolerance]
            minimum = 3 if degree == 1 else 6
            if len(good) < minimum:
                continue
            quality = (
                len(good) - (degree - 1),
                -np.median(residual[good]),
                np.mean(scores[good]),
            )
            hypotheses.append((quality, degree, good))
    if not hypotheses:
        return None
    _, degree, selected = max(hypotheses, key=lambda h: h[0])
    # Refine only using actual, distinct-carbon observations inside consensus.
    for _ in range(3):
        coef = np.polynomial.polynomial.polyfit(x[selected], y[selected], degree)
        residual = abs(y - np.polynomial.polynomial.polyval(x, coef))
        nearest = np.array(
            [idx[np.lexsort((-scores[idx], residual[idx]))[0]] for idx in slots]
        )
        updated = nearest[residual[nearest] <= tolerance]
        if len(updated) < (3 if degree == 1 else 6):
            break
        selected = updated
    points = pool.iloc[selected].drop(columns="_bin")
    coef = np.polynomial.polynomial.polyfit(
        points.total_C.to_numpy(float), points._rt_minutes.to_numpy(float), degree
    )
    bounds = (float(points.total_C.min()), float(points.total_C.max()))
    derivative = np.polynomial.polynomial.polyval(
        np.linspace(*bounds, 64), np.polynomial.polynomial.polyder(coef)
    )
    if derivative.min() <= 0:
        return None
    curve = _Curve(
        "quadratic" if degree == 2 else "linear", tuple(coef), *bounds, len(points)
    )
    return curve, points, "independent_consensus"


def _crosses(a, b):
    left = max(a[0].carbon_min, b[0].carbon_min)
    right = min(a[0].carbon_max, b[0].carbon_max)
    if left > right:
        return False
    ca = np.pad(a[0].coefficients, (0, 3 - len(a[0].coefficients)))
    cb = np.pad(b[0].coefficients, (0, 3 - len(b[0].coefficients)))
    diff = ca - cb
    x = [left, right]
    if abs(diff[2]) > 1e-14:
        vertex = -diff[1] / (2 * diff[2])
        if left < vertex < right:
            x.append(vertex)
    x = np.array(x)
    return bool(np.any(_predict(a, x) < _predict(b, x) - 1e-8))


def _observed_pair(group, db, accepted, threshold):
    """Fit two observed carbon points; neighboring models only validate them."""
    if group.total_C.nunique() != 2:
        return None
    bounds = (float(group.total_C.min()), float(group.total_C.max()))
    covering = {
        k: m
        for k, m in accepted.items()
        if m[0].carbon_min <= bounds[0] and m[0].carbon_max >= bounds[1]
    }
    brackets = [
        (hi - lo, lo, hi)
        for lo, hi in combinations(sorted(covering), 2)
        if lo < db < hi
    ]
    if brackets:
        _, low, high = min(brackets)
    else:
        nearest = sorted(covering, key=lambda k: (abs(k - db), k))[:2]
        if len(nearest) < 2 or abs(nearest[0] - db) != 1:
            return None
        low, high = sorted(nearest)
    weight = (db - low) / (high - low)
    expected = (1 - weight) * _predict(
        covering[low], np.array(bounds)
    ) + weight * _predict(covering[high], np.array(bounds))
    pool = group.copy()
    pool["_bin"] = (pool._rt_minutes / 0.05).round()
    pool = pool.sort_values(
        ["_score", "_rank"], ascending=[False, True], kind="stable"
    ).drop_duplicates(["total_C", "_bin"])
    choices = [
        pool.loc[pool.total_C.eq(c) & (pool._rt_minutes - t).abs().le(threshold)].head(
            20
        )
        for c, t in zip(bounds, expected)
    ]
    best = None
    for i in choices[0].index:
        for j in choices[1].index:
            points = pool.loc[[i, j]].drop(columns="_bin")
            y = points._rt_minutes.to_numpy(float)
            slope = (y[1] - y[0]) / (bounds[1] - bounds[0])
            if slope <= 0:
                continue
            curve = _Curve("linear", (y[0] - slope * bounds[0], slope), *bounds, 2)
            model = (curve, points, f"observed_two_point_DB{low}_DB{high}")
            if any(
                _crosses(other, model) if k < db else _crosses(model, other)
                for k, other in accepted.items()
            ):
                continue
            quality = (
                float(points._score.fillna(0).mean()),
                -float(abs(y - expected).mean()),
            )
            if best is None or quality > best[0]:
                best = (quality, model)
    return None if best is None else best[1]


def fit_ordered_family(candidates, threshold=0.5):
    """Return validated DB models, using only supplied high-confidence rows.

    Two actual carbon points define their own line, validated by two DB
    neighbors. Such observed lines may validate adjacent two-point series.
    Borrowed single-point models cannot recursively supply new references.
    """
    groups = {int(db): g for db, g in candidates.groupby("total_DB")}
    independent = {
        db: m
        for db, g in groups.items()
        if (m := _consensus(g, threshold / 2)) is not None
    }
    # Resolve crossings by support, not by cosmetically shifting coefficients.
    accepted = {}
    priority = sorted(independent, key=lambda db: (-len(independent[db][1]), db))
    for db in priority:
        model = independent[db]
        if len(model[1]) == 3:
            neighbors = [
                other
                for other_db, other in independent.items()
                if other_db != db
                and abs(other_db - db) <= 2
                and max(model[0].carbon_min, other[0].carbon_min)
                < min(model[0].carbon_max, other[0].carbon_max)
                and not (
                    _crosses(other, model) if other_db < db else _crosses(model, other)
                )
            ]
            if len(neighbors) < 2:
                continue
        if any(
            _crosses(other, model) if other_db < db else _crosses(model, other)
            for other_db, other in accepted.items()
        ):
            continue
        accepted[db] = model
    # Extend in either DB direction until no additional observed pair qualifies.
    # Each accepted extension contributes two real high-confidence anchors,
    # rather than recursively manufacturing coefficients from borrowed curves.
    while True:
        added = False
        for db, group in groups.items():
            if db not in accepted:
                model = _observed_pair(group, db, accepted, threshold)
                if model is not None:
                    accepted[db] = model
                    added = True
        if not added:
            break
    references = dict(accepted)
    for db, group in groups.items():
        if db in accepted or len(references) < 2:
            continue
        brackets = []
        for low, high in combinations(sorted(references), 2):
            if low < db < high:
                brackets.append((high - low, low, high))
        if not brackets:
            # At a family edge, permit at most one DB step of borrowing.
            nearest = sorted(references, key=lambda k: abs(k - db))[:2]
            if abs(nearest[0] - db) > 1:
                continue
            low, high = sorted(nearest)
        else:
            _, low, high = min(brackets)
        a, b = references[low], references[high]
        left = max(a[0].carbon_min, b[0].carbon_min)
        right = min(a[0].carbon_max, b[0].carbon_max)
        pool = group.loc[group.total_C.between(left, right)].copy()
        if pool.empty:
            continue
        weight = (db - low) / (high - low)
        ca = np.pad(a[0].coefficients, (0, 3 - len(a[0].coefficients)))
        cb = np.pad(b[0].coefficients, (0, 3 - len(b[0].coefficients)))
        coef = (1 - weight) * ca + weight * cb
        pool["_residual"] = abs(
            pool._rt_minutes
            - np.polynomial.polynomial.polyval(pool.total_C.to_numpy(float), coef)
        )
        points = (
            pool.loc[pool._residual <= threshold]
            .sort_values(
                ["_residual", "_score"], ascending=[True, False], kind="stable"
            )
            .drop_duplicates("total_C")
            .drop(columns="_residual")
        )
        if points.empty:
            continue
        bounds = (float(points.total_C.min()), float(points.total_C.max()))
        support = f"neighbor_supported_DB{low}_DB{high}"
        if len(points) == 2:
            points = points.sort_values("total_C")
            y = points._rt_minutes.to_numpy(float)
            slope = (y[1] - y[0]) / (bounds[1] - bounds[0])
            coef = np.array([y[0] - slope * bounds[0], slope, 0.0])
            support = f"observed_two_point_DB{low}_DB{high}"
        derivative = np.polynomial.polynomial.polyval(
            np.linspace(*bounds, 64), np.polynomial.polynomial.polyder(coef)
        )
        if derivative.min() <= 0:
            continue
        linear = abs(coef[2]) < 1e-12
        curve = _Curve(
            "linear" if linear else "quadratic",
            tuple(coef[:2] if linear else coef),
            *bounds,
            len(points),
        )
        model = (curve, points, support)
        if any(
            _crosses(other, model) if other_db < db else _crosses(model, other)
            for other_db, other in accepted.items()
        ):
            continue
        accepted[db] = model
    return accepted
