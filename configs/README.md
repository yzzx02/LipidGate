# Project settings

Create a project and save parameters in the GUI; `lipidgate.project.json` is
the portable settings format. Input files are referenced by path and each run
snapshots its parameters. Re-run saved settings with:

```powershell
python -m lipidgate project-run --project "path/to/project"
```

The previous benchmark YAML / peak-truth model settings are retired and are
not accepted as a replacement for a saved desktop project.
