# TS study: Me-substituted NH-TDSO to N-methyl pyrazoles

## Inference from `TS.png`
The left scaffold is an **NH 1,2,3-thiadiazine 1-oxide (TDSO)** and the depicted transformation is **sulfur monoxide (SO) extrusion / ring contraction** to an **N-alkyl pyrazole**.  
Using the requested substitutions:
- ring carbon **alpha to sulfur** = **Me**
- all other peripheral ring substituents = **H**
- N-alkyl group = **Me**

The pre-alkylation NH tautomerism of the TDSO ring implies two neutral N-methylated TDSO regioisomers after methylation, so two pathways were modeled.

## Modeled pathways
### Pathway A
- **Substrate:** N-methyl TDSO regioisomer corresponding to alkylation at the NH-adjacent ring N  
  SMILES used: `CC1=CC=NN(C)S1=O`
- **Product:** **1,5-dimethylpyrazole + SO**  
  pyrazole SMILES: `CC1=CC=NN1C`

### Pathway B
- **Substrate:** N-methyl TDSO regioisomer corresponding to alkylation at the other ring N (modeled in its zwitterionic resonance form)  
  SMILES used: `CC1=CC=[N+](C)[N-]S1=O`
- **Product:** **1,3-dimethylpyrazole + SO**  
  pyrazole SMILES: `CC1=NN(C)C=C1`

## Computational protocol
- Initial 3D structures: Open Babel generation from chemically inferred SMILES
- Minima: **GFN2-xTB** geometry optimization
- TS guesses: made by **directly editing the optimized substrates** toward C–N bond formation and S–N / S–C cleavage
- TS refinement: **Sella optTS at GFN2-xTB level**
- TS verification: xTB analytical Hessian / frequency analysis
- Energies reported in `energies.csv`: **electronic energies on the same closed-shell GFN2-xTB surface**

## TS verification
- **Pathway A TS:** one imaginary frequency, **-1270.4 cm-1**
- **Pathway B TS:** one imaginary frequency, **-381.5 cm-1**

Both optimized TS structures are first-order saddle points on the modeled singlet surface.

## Relative energies
See `energies.csv` for the full table. Key values (kcal/mol):

| Pathway | Substrate | TS vs own substrate | Product vs own substrate |
|---|---:|---:|---:|
| A | 0.0 | 31.1 | 16.6 |
| B | 0.0 | 10.2 | 7.3 |

Using the **lowest substrate overall** as zero:
- substrate A = 0.0
- substrate B = +10.1
- TS A = +31.1
- TS B = +20.4
- product A = +16.6
- product B = +17.4

## Mechanistic conclusion
- The reaction is a **cheletropic SO extrusion / ring-contraction** from N-methyl TDSO to N-methyl pyrazole.
- **Pathway B is kinetically favored** on the modeled surface because its barrier from its own substrate is much lower (**10.2 kcal/mol**) than Pathway A (**31.1 kcal/mol**).
- **Pathway A is thermodynamically favored at the substrate level** and gives the slightly lower-energy product (**product A lower than product B by ~0.8 kcal/mol**).
- Therefore, within this GFN2-xTB singlet-surface model, the two pathways are **not isoenergetic**: A is lower before extrusion, whereas B extrudes SO more easily once the corresponding N-methyl TDSO is formed.

## Files produced
- `structures/pathway_A_substrate.xyz`
- `structures/pathway_A_TS_Initial_Guess.xyz`
- `structures/pathway_A_TS_optimized.xyz`
- `structures/pathway_A_product.xyz`
- `structures/pathway_B_substrate.xyz`
- `structures/pathway_B_TS_Initial_Guess.xyz`
- `structures/pathway_B_TS_optimized.xyz`
- `structures/pathway_B_product.xyz`
- pathway-specific calculation outputs under `calculations/pathway_A/` and `calculations/pathway_B/`
- `energies.csv`
