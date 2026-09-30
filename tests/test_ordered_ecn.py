import numpy as np
import pandas as pd

from lipidgate.ecn_filter.series_consensus import fit_ordered_family


def rows(db, carbon, rt):
    return [
        {"total_DB": db, "total_C": c, "_rt_minutes": t, "_score": 90.0, "_rank": 1}
        for c, t in zip(carbon, rt)
    ]


def test_alternate_observed_rt_replaces_high_outlier():
    data = rows(0, [16, 18, 20, 22, 24], [4, 5, 6, 7, 8]) + rows(
        1, [16, 18, 20, 22, 24], [3, 4, 5, 12, 7]
    )
    data += rows(1, [22], [6])
    result = fit_ordered_family(pd.DataFrame(data))
    curve, points, _ = result[1]
    assert points.loc[points.total_C.eq(22), "_rt_minutes"].iloc[0] == 6
    assert np.polynomial.polynomial.polyval(22, curve.coefficients) < 7


def test_sparse_neighbors_support_one_two_and_three_points():
    data = rows(0, [30, 32, 34, 36, 38], [8, 9, 10, 11, 12]) + rows(
        4, [30, 32, 34, 36, 38], [4, 5, 6, 7, 8]
    )
    data += (
        rows(1, [34], [9])
        + rows(2, [32, 36], [7, 9])
        + rows(3, [30, 34, 38], [5, 7, 9])
    )
    result = fit_ordered_family(pd.DataFrame(data))
    assert set(result) == {0, 1, 2, 3, 4}
    assert result[1][2].startswith("neighbor_supported")
    assert result[1][0].carbon_min == result[1][0].carbon_max == 34
    assert len(result[2][1]) == 2
    assert result[3][0].model_type == "linear"


def test_inconsistent_sparse_point_is_not_forced_onto_curve():
    data = (
        rows(0, [30, 32, 34, 36], [8, 9, 10, 11])
        + rows(2, [30, 32, 34, 36], [6, 7, 8, 9])
        + rows(1, [34], [20])
    )
    assert 1 not in fit_ordered_family(pd.DataFrame(data))


def test_overlapping_accepted_curves_do_not_cross():
    data = rows(0, [30, 32, 34, 36, 38], [8, 9, 10, 11, 12]) + rows(
        1, [30, 32, 34], [6, 9, 12]
    )
    result = fit_ordered_family(pd.DataFrame(data))
    if 1 in result:
        x = np.linspace(result[1][0].carbon_min, result[1][0].carbon_max, 100)
        assert np.all(
            np.polynomial.polynomial.polyval(x, result[0][0].coefficients)
            >= np.polynomial.polynomial.polyval(x, result[1][0].coefficients) - 1e-8
        )


def test_isolated_three_point_series_without_db_context_is_not_confirmed():
    data = rows(0, [30, 32, 34, 36], [8, 9, 10, 11]) + rows(
        1, [30, 32, 34, 36], [7, 8, 9, 10]
    )
    data += rows(20, [37, 41, 49], [3, 3.5, 4.7])
    assert 20 not in fit_ordered_family(pd.DataFrame(data))


def test_low_confidence_rows_cannot_change_ordered_models():
    from lipidgate.ecn_filter.confidence_filter import fit_high_confidence_ecn

    data = []
    for db in (0, 2):
        for c in (30, 32, 34, 36, 38):
            data.append(
                {
                    "matched_name": f"PC({c}:{db})",
                    "compound_class": "PC",
                    "rt_minutes": c / 2 - db,
                    "final_score": 90,
                    "result_rank": 1,
                    "置信度": "高",
                }
            )
    high = pd.DataFrame(data)
    _, before, _ = fit_high_confidence_ecn(high, ordered_series=True)
    extra = high.iloc[[2]].copy()
    extra["置信度"] = "低"
    extra["final_score"] = 100
    extra["rt_minutes"] = 50
    _, after, anchors = fit_high_confidence_ecn(
        pd.concat([high, extra], ignore_index=True), ordered_series=True
    )
    pd.testing.assert_frame_equal(before, after)
    assert anchors["置信度"].eq("高").all()


def test_observed_two_point_series_can_validate_next_db():
    data = rows(8, [54, 56, 58, 60], [14.5, 15, 15.5, 16])
    data += rows(9, [54, 56, 58, 60], [14.1, 14.6, 15.1, 15.6])
    data += rows(10, [58, 60], [14.8, 15.4])
    data += rows(11, [58, 60], [14.63, 15.22])
    result = fit_ordered_family(pd.DataFrame(data))
    for db, expected in [(10, [14.8, 15.4]), (11, [14.63, 15.22])]:
        assert result[db][2].startswith("observed_two_point")
        np.testing.assert_allclose(
            np.polynomial.polynomial.polyval([58, 60], result[db][0].coefficients),
            expected,
        )


def test_two_point_crossing_and_remote_pair_are_not_confirmed():
    data = rows(8, [54, 56, 58, 60], [14.5, 15, 15.5, 16])
    data += rows(9, [54, 56, 58, 60], [14.1, 14.6, 15.1, 15.6])
    data += rows(10, [58, 60], [14.9, 15.7])
    data += rows(20, [58, 60], [3, 4])
    result = fit_ordered_family(pd.DataFrame(data))
    assert not result.get(10, (None, None, ""))[2].startswith("observed_two_point")
    assert 20 not in result
