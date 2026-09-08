# TS study from `TS.png`

## Structural inference used
- The left scaffold in `TS.png` was interpreted as an **NH-TDSO six-membered S,N,N heterocycle** consistent with a **1,2,3-thiadiazine 1-oxide-like ring**.
- The carbon **alpha to sulfur** was set to **Me**.
- The two other peripheral ring substituents were set to **H**.
- The **N-alkyl group** was set to **Me**.
- The mechanistic transformation inferred from the figure is **SO extrusion / ring contraction to an N-methyl pyrazole**.

## Pathways built
### Pathway A
- **Substrate:** N-methyl on the sulfur-adjacent ring N.
- **Product:** **1-methyl-5-methylpyrazole + SO**.

### Pathway B
- **Substrate:** N-methyl on the other ring N tautomer/regioisomer.
- **Product:** **2-methyl-5-methylpyrazole + SO**.

## What was done
- Built and optimized both substrates at **GFN2-xTB**.
- Built both products as optimized pyrazole isomers plus separated SO.
- Generated TS guesses by direct geometry editing of each substrate.
- Ran **Sella optTS at xTB level** for both pathways.
- Verified the Sella-relaxed stationary points by numerical frequency analysis.

## Outcome of TS verification
- For **both pathways**, the Sella-relaxed stationary points gave **0 imaginary frequencies** at GFN2-xTB.
- Additional frequency checks on constrained TS-like geometries also gave **0 imaginary frequencies**.
- Therefore, **no validated first-order saddle point was located for either pathway** on this closed-shell xTB surface.
- The saved `pathway_A_TS_optimized.xyz` and `pathway_B_TS_optimized.xyz` are therefore the **closest Sella-relaxed stationary points obtained**, not confirmed TSs.

## Relative-energy summary
Reference = **Pathway A substrate = 0.00 kcal/mol**.

| Pathway | State | Rel. energy vs A substrate (kcal/mol) | Rel. energy vs own substrate (kcal/mol) | Validation |
|---|---:|---:|---:|---|
| A | substrate | 0.00 | 0.00 | minimum |
| A | Sella-relaxed TS candidate | 23.09 | 23.09 | **not validated** (0 imag) |
| A | product (pyrazole + SO) | 16.72 | 16.72 | fragment energy sum |
| B | substrate | 10.18 | 0.00 | minimum |
| B | Sella-relaxed TS candidate | 20.44 | 10.26 | **not validated** (0 imag) |
| B | product (pyrazole + SO) | 17.49 | 7.32 | fragment energy sum |

## Mechanistic takeaways
- **Pathway A substrate** is lower in energy than **Pathway B substrate** by **10.18 kcal/mol**.
- On this xTB surface, **Product A** is lower than **Product B** by **0.78 kcal/mol**.
- The closed-shell GFN2-xTB/Sella treatment did **not** locate a verified first-order saddle for the SO-extrusion/ring-contraction step, so the saved TS-optimized structures should be treated as **TS candidates only**.

## Files
- Structures: `ts-study/structures/`
- Logs, trajectories, and frequency outputs: `ts-study/calculations/pathway_A/` and `ts-study/calculations/pathway_B/`
- Energy table: `ts-study/energies.csv`
