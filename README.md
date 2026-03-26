# reaction-smarts — Phase 34

Validate 6 synthetic reaction SMARTS patterns against a compound library.
Reports which reactions are applicable to each compound and which families
are most synthetically tractable.

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv
```

## Outputs

| File | Description |
|---|---|
| `output/reaction_applicability.csv` | Compound × reaction — 1 if applicable |
| `output/reaction_heatmap.png` | Heatmap with family color bar |

## Reactions

| Reaction | Description |
|---|---|
| amide_coupling | Carboxylic acid → amide bond |
| ester_formation | Carboxylic acid → ester |
| reductive_amination | Aldehyde/ketone → amine |
| sulfonamide_form | Sulfonyl chloride → sulfonamide |
| buchwald_amine | Aryl halide → aryl amine (Buchwald-Hartwig) |
| suzuki_coupling | Aryl halide → biaryl (Suzuki) |
