# TS study: tautomer-dependent N-methylation from the TS.png substitution pattern

## Inference from `TS.png`
The figure was used only to define the heteroaromatic substitution pattern after simplification:
- the ring carbon corresponding to the carbon **alpha to sulfur** in the precursor was set to **Me**,
- all other peripheral ring substituents were set to **H**,
- the N-alkyl group was set to **Me**.

Because the prompt explicitly invokes **pre-alkylation N–H tautomerism**, the two neutral starting structures were taken as the two N–H tautomers of methylpyrazole:

- **Pathway A substrate:** 3-methyl-1H-pyrazole + diazomethane
- **Pathway B substrate:** 5-methyl-1H-pyrazole + diazomethane

To make the alkylation step a closed-shell neutral PES suitable for direct xTB/Sella TS optimization, **diazomethane** was used as the methyl donor:

- **Pathway A product:** 1,3-dimethylpyrazole + N2
- **Pathway B product:** 1,5-dimethylpyrazole + N2

## Modeled transformation
Concerted/near-concerted **N-methyl transfer from diazomethane to the pyrazole N–H tautomer**, with concurrent H transfer to the methylene carbon and N2 extrusion.

## Calculations performed
- Fragment geometries optimized at **GFN2-xTB**.
- Pre-reactive substrate complexes built explicitly and optimized with xTB.
- TS guesses built by direct editing of the substrate complexes.
- TS optimization performed with **Sella + xTB (GFN2-xTB)**.
- Frequency checks performed with **xTB Hessians** on the optimized TS candidates.

## Energetic summary
Reference values are in `energies.csv`.

### Pathway A
- Substrate complex energy: **-26.3669564243 Eh**
- TS energy: **-26.3668267771 Eh**
- Barrier from substrate complex: **+0.081 kcal/mol**
- Product energy (optimized separated products): **-26.3876057021 Eh**
- Reaction energy from substrate complex: **-12.958 kcal/mol**

### Pathway B
- Substrate complex energy: **-26.3684158034 Eh**
- TS energy: **-26.3683118925 Eh**
- Barrier from substrate complex: **+0.065 kcal/mol**
- Product energy (optimized separated products): **-26.3888423209 Eh**
- Reaction energy from substrate complex: **-12.818 kcal/mol**

### Cross-path comparison
- Pathway B substrate complex is lower than pathway A by **0.916 kcal/mol**.
- Pathway B also gives the slightly lower TS.
- Both modeled methylations are strongly exergonic at this xTB level.

## Frequency verification
### Pathway A TS candidate
Projected vibrational frequencies include:
- **-78.6 cm-1** (dominant imaginary mode)
- additional soft negatives: **-50.7, -14.6, -5.1 cm-1**

Interpretation: the dominant unstable mode is consistent with methyl transfer, but the xTB/Sella stationary point retains extra low-curvature directions and should be treated as a **soft, partially converged TS candidate rather than a rigorously clean first-order saddle**.

### Pathway B TS candidate
Projected vibrational frequencies include:
- **-50.8 cm-1** (dominant imaginary mode)
- extra very soft negatives: **-21.5, -11.7 cm-1**

Interpretation: this is the **cleaner** of the two xTB TS candidates, though the extra very low negatives indicate residual flat intermolecular degrees of freedom typical of loose reagent/product complexes.

## Files saved
Top-level requested structures are in `ts-study/structures/`.
Raw xTB/Sella outputs, fragment optimizations, Hessians, and logs are in `ts-study/calculations/`.
