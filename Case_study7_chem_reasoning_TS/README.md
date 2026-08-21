# Case Study 7: Competing TDSO pathways

This directory contains the Case Study 7 benchmark input and artifacts reported
for GPT-5.4, GPT-5.2, and Claude Opus 4.5. Outputs are preserved as generated,
including incomplete or chemically incorrect runs.

## Input

- `input.json` is the exact benchmark question. Its dependency file is embedded
  as base64 so the input is self-contained.
- `inputs/TS.png` is the same dependency extracted for convenient inspection.
- `full_ts_prompt.txt` is the exact full prompt transmitted to GPT-5.4 for the
  transition-state benchmark.
- `endpoint_only_prompt.txt` is the exact prompt used for the GPT-5.2 and Claude
  Opus 4.5 endpoint-only baselines.

## Runs

| Directory | Model | Effort | Task | Repetitions |
| --- | --- | --- | --- | --- |
| `runs/gpt-5.4_max` | `gpt-5.4` | max | Full endpoint construction, direct TS guess, Sella optTS, frequency analysis, and relative energies | 5 |
| `runs/gpt-5.2_medium_endpoint_only` | `gpt-5.2` | medium | Endpoint structures only | 5 |
| `runs/claude-opus-4-5_high_endpoint_only` | `claude-opus-4-5-20251101` | high | Endpoint structures only | 5 |

Each `repeat_N/ts-study` directory is the corresponding agent-generated output
tree. Temporary visualization sidecars and `.xtboptok` marker files were omitted;
the molecular structures, reports, scripts, trajectories, and calculation outputs
were retained.

## Session IDs

| Model | Repeat | Session ID |
| --- | ---: | --- |
| GPT-5.4 | 1 | `Estructural_Case7__gpt-5.4__max__repeat_1__21af4fe2` |
| GPT-5.4 | 2 | `Estructural_Case7__gpt-5.4__max__repeat_2__75d6a9a0` |
| GPT-5.4 | 3 | `Estructural_Case7__gpt-5.4__max__repeat_3__cf14d77d` |
| GPT-5.4 | 4 | `Estructural_Case7__gpt-5.4__max__repeat_4__216019aa` |
| GPT-5.4 | 5 | `Estructural_Case7__gpt-5.4__max__repeat_5__327da67c` |
| GPT-5.2 | 1 | `Estructural_Case7__gpt-5.2__medium__repeat_1__933acf8f` |
| GPT-5.2 | 2 | `Estructural_Case7__gpt-5.2__medium__repeat_2__75438750` |
| GPT-5.2 | 3 | `Estructural_Case7__gpt-5.2__medium__repeat_3__9e5f9e23` |
| GPT-5.2 | 4 | `Estructural_Case7__gpt-5.2__medium__repeat_4__e0882c1b` |
| GPT-5.2 | 5 | `Estructural_Case7__gpt-5.2__medium__repeat_5__c4a8b142` |
| Claude Opus 4.5 | 1 | `Estructural_Case7__claude-opus-4-5-20251101__high__repeat_1__44ddc27f` |
| Claude Opus 4.5 | 2 | `Estructural_Case7__claude-opus-4-5-20251101__high__repeat_2__f1af81a1` |
| Claude Opus 4.5 | 3 | `Estructural_Case7__claude-opus-4-5-20251101__high__repeat_3__069562cd` |
| Claude Opus 4.5 | 4 | `Estructural_Case7__claude-opus-4-5-20251101__high__repeat_4__c7884835` |
| Claude Opus 4.5 | 5 | `Estructural_Case7__claude-opus-4-5-20251101__high__repeat_5__5cd99d5b` |

## Metadata

- `metadata/gpt-5.4_max.json` contains per-run responses, generated-file lists,
  token usage, timing, and process telemetry from the native benchmark logger.
- `metadata/gpt-5.4_scientific_evaluation.json` contains artifact hashes,
  workflow auditing, ground-truth comparisons, frequency signals, and verdicts.
- The GPT-5.2 and Claude Opus 4.5 JSON files contain their native benchmark
  telemetry, including responses, failures, generated files, token usage,
  timing, and process data.
- `SHA256SUMS` provides checksums for the published input, metadata, and run
  artifacts.

The cost fields are logger estimates. A zero value does not imply a free run:
the baseline manifests mark cost coverage as incomplete where the provider/model
price was unavailable to the logger.
