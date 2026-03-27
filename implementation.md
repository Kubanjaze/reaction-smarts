# Phase 34 — Reaction SMARTS Validator

**Version:** 1.1 | **Tier:** Micro | **Date:** 2026-03-26

## Goal
Validate a library of synthetic reaction SMARTS against the compound library.
For each reaction, test whether it can transform any compound into a valid product.
Report which reactions are applicable, how many compounds they can process, and
which compounds are transformed by each reaction.

CLI: `python main.py --input data/compounds.csv`

Outputs: reaction_applicability.csv, reaction_heatmap.png

## Reaction Library (6 reactions)
| Reaction | SMARTS | Description |
|---|---|---|
| amide_coupling | `[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[N:2]` | Acid + amine → amide |
| ester_formation | `[C:1](=O)[OH].[OH:2][C,c]>>[C:1](=O)[O:2]` | Acid + alcohol → ester |
| reductive_amination | `[C:1]=O.[NH2:2]>>[C:1][N:2]` | Aldehyde/ketone + amine → amine |
| sulfonamide_formation | `[S:1](=O)(=O)[Cl].[NH2:2]>>[S:1](=O)(=O)[N:2]` | Sulfonyl chloride + amine |
| buchwald | `[c:1][Br,Cl].[NH2:2]>>[c:1][N:2]` | Aryl halide + amine → aryl amine (Pd) |
| suzuki | `[c:1][Br,Cl].[c:2][B]>>[c:1][c:2]` | Aryl halide + boronic acid (Pd) |

## Logic
- `AllChem.ReactionFromSmarts(rxn_smarts)` to parse each reaction
- For each compound, test whether it contains reactant substructures (reactant matching)
- Use `rxn.RunReactants((mol,))` for single-reactant reactions, or flag as multi-reactant
- Count: n_applicable = number of compounds where reaction produces valid products
- Valid product = non-None, sanitizable mol

## Outputs
- reaction_applicability.csv: compound × reaction — 1 if applicable, 0 if not
- reaction_heatmap.png: heatmap compound × reaction with family color bar

## Key Concepts
- RDKit reaction SMARTS (`AllChem.ReactionFromSmarts`, `RunReactants`)
- Reactant substructure matching for synthetic feasibility
- Coverage of 6 common medicinal chemistry reactions (amide, ester, reductive amination, sulfonamide, Buchwald, Suzuki)

## Verification Checklist
- [x] 45/45 compounds parsed without SMILES failures
- [x] reaction_applicability.csv contains compound x reaction binary matrix
- [x] reaction_heatmap.png saved to output/
- [x] Only Buchwald and Suzuki applicable (8/45 compounds with Cl/Br handles)
- [x] No false positives: 4 non-applicable reactions correctly return 0/45

## Risks
- Reaction SMARTS only check for reactant substructure presence, not synthetic yield or selectivity
- Multi-reactant reactions tested in single-reactant mode may miss valid transformations requiring two substrates

## Actual Results (v1.1)

| Reaction | Applicable | Notes |
|---|---|---|
| amide_coupling | 0/45 | No free acids in library |
| ester_formation | 0/45 | No free acids in library |
| reductive_amination | 0/45 | No free aldehydes/ketones |
| sulfonamide_form | 0/45 | No sulfonyl chlorides |
| buchwald_amine | 8/45 | Cl/Br-substituted compounds |
| suzuki_coupling | 8/45 | Same Cl/Br-substituted compounds |

**Key insight:** This CETP library consists of finished compounds — no reactive acid/aldehyde handles. Only aryl halide substituents (Cl, Br) allow further functionalization via Pd-catalyzed cross-coupling.
**Most versatile:** benz_002_Cl, benz_003_Br, benz_012_dichloro (2 applicable reactions each)
**45/45 valid**
