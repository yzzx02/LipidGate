from __future__ import annotations

import pytest


def test_pyopenms_peak_height_uses_apex_not_integrated_area():
    oms = pytest.importorskip("pyopenms")
    from lipidbench.runners.run_pyopenms import filter_feature_map_by_height

    feature_map = oms.FeatureMap()
    for mz, area, height in ((500.0, 1e9, 50.0), (600.0, 1e3, 500.0)):
        feature = oms.Feature()
        feature.setMZ(mz)
        feature.setIntensity(area)
        feature.setMetaValue("max_height", height)
        feature_map.push_back(feature)
    filtered = filter_feature_map_by_height(feature_map, 100.0)
    assert filtered.size() == 1
    assert filtered[0].getMZ() == 600.0
    assert feature_map.size() == 2


def test_xcms_minimum_samples_and_peak_height_reach_r_runner(monkeypatch):
    from lipidbench.runners.run_xcms import extract_xcms_params, run_xcms

    parameters = extract_xcms_params({"parameters": {"xcms": {"peak_picking": {
        "peakwidth": [6.0, 48.0], "noise": 900.0,
        "min_peak_height": 2500.0, "minFraction": 0.4, "minSamples": 3,
    }}}})
    assert parameters["minwidth"] == 6.0 and parameters["maxwidth"] == 48.0
    assert parameters["noise"] == 900.0
    assert parameters["min_maxo"] == 2500.0
    assert parameters["frac"] == 0.4 and parameters["min_samples"] == 3

    commands = []
    monkeypatch.setattr("subprocess.run", lambda command, check: commands.append(command))
    run_xcms(input_dir="input", output_file="output.csv", **parameters)
    command = commands[0]
    assert command[command.index("--min_samples") + 1] == "3"
    assert command[command.index("--min_maxo") + 1] == "2500.0"
