# AI assistance with LPERFECT metrics

LPERFECT can emit **GPT-friendly JSON metrics** so you can hand results to an AI
assistant for configuration review, usage guidance, or tuning suggestions. Two
complementary reports are available:

- **`metrics.assistant` / `--ai-metrics`**: compact hydrology + compute summary
  with wall-clock time, key parameters, and budget checks.
- **`metrics.parallelization` / `--parallel-metrics`**: detailed performance
  sampling (wall-clock per step, migration ratios, GPU/CPU timing).

Both reports are logged on rank 0; if you provide an output path they are also
written to disk as JSON.

### Choose detailed vs compact JSON

Both GPT-friendly reports can be emitted as **detailed** (pretty-printed with
indentation) or **compact** (minified, no extra whitespace) JSON. Use the CLI
flags `--ai-metrics-format {detailed|compact}` and `--parallel-metrics-format
{detailed|compact}` or set `metrics.assistant.format` / `metrics.parallelization.format`
in the config file to match your workflow (pretty for debugging, compact for
token-efficient copy/paste into an AI assistant).

## Enabling the AI-friendly hydrology + compute report

The assistant-oriented report focuses on **what you ran** (grid size, travel
times, device, MPI layout) and **how it behaved** (rain/runoff/outflow volumes,
mass balance, hydrology flags). Enable it via CLI:

```bash
python main.py --config config.json --ai-metrics --ai-metrics-output runs/ai_metrics.json
```

…or via JSON:

```json
{
  "metrics": {
    "assistant": {
      "enabled": true,
      "output": "runs/ai_metrics.json"
    }
  }
}
```

### What the assistant metrics contain

- `domain`: domain label, grid shape, active/total cells, and cell-area summary
  (scalar or min/mean/max for grids).
- `time`: configured start time, `dt_s`, total simulation window, executed steps,
  and wall-clock duration (MPI max when distributed).
- `compute`: device, MPI ranks/world size, shared-memory workers, travel-time
  summary (fixed values or active-cell medians/ranges), particle volume.
- `hydrology`: initial abstraction ratio, outflow sink flag, rain/runoff/outflow/
  storage volumes, runoff/rain ratio, mass-error stats, hydrology status flag.
- `status` and `notes`: quick booleans plus reminders that the structure is
  optimized for AI ingestion.

Example (truncated):

```json
{
  "domain": { "name": "campania", "shape_rc": [200, 300], "active_cells": 52000 },
  "time": { "dt_s": 5.0, "simulation_window_s": 7200.0, "elapsed_wall_clock_s": 182.4 },
  "compute": { "device": "cpu", "mpi_ranks": 4, "shared_memory_workers": 16 },
  "hydrology": { "rain_m3": 2.3e7, "runoff_m3": 1.1e7, "runoff_to_rain_ratio": 0.48 },
  "status": { "numerical_ok": true, "hydrological_ok": true }
}
```

Use this JSON as context for an AI assistant to ask questions like:

- “Is the runoff/rain ratio plausible for these Curve Numbers and rainfall totals?”
- “How should I adjust travel-time parameters to reduce wall-clock time while
  preserving hydrological consistency?”

## Using the performance sampling report

For step-level performance diagnostics (e.g., migration ratios, hop throughput,
GPU vs. CPU timing), enable:

```bash
python main.py --config config.json --parallel-metrics --parallel-metrics-output runs/perf_metrics.json
```

Details of the performance schema live in `docs/parallelization_schema.md`. You
can combine both reports in the same run when you need **both** hydrological
assessment and performance tuning input for an AI helper.
