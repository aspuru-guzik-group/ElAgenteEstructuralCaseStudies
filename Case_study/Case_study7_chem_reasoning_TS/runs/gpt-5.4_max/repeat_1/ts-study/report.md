# TS study: NH-TDSO to N-methyl pyrazoles

## Structural inference from `TS.png`
- Core substrate class: **1,2,6-thiadiazine 1-oxide (TDSO)**.
- Requested substitution pattern applied here:
  - ring carbon **alpha to sulfur = Me**
  - all other peripheral ring substituents = H
  - **N-alkyl = Me**
- Because the exchangeable NH can tautomerize between the two ring nitrogens **before alkylation**, two neutral N-methylated starting structures are possible:
  - **Pathway A substrate**: methyl on the ring N adjacent to S=O; this maps to **1,5-dimethylpyrazole** after SO extrusion.
  - **Pathway B substrate**: methyl on the distal ring N (best represented as a neutral zwitterionic valence isomer); this maps to **1,3-dimethylpyrazole** after SO extrusion.

## Transformation modeled
A **concerted cheletropic SO extrusion / ring contraction** from the N-methyl TDSO to the corresponding N-methyl pyrazole was modeled on the **GFN2-xTB singlet surface**. Product files contain the pyrazole plus a far-separated SO fragment so atom counts remain conserved.

## Transition-state verification
- **Pathway A TS**: one significant imaginary frequency, **-433.00 cm^-1**.
- **Pathway B TS**: one significant imaginary frequency, **-389.26 cm^-1**.
- In both optimized TS structures, the forming N–C bond and the breaking S–N / S–C contacts are all in the expected partially formed/broken range.

## Relative energies (kcal/mol)
Reference = lowest substrate energy (Pathway A substrate = 0.0 kcal/mol).

| Pathway | Species | Model | Rel. E | Local ΔE / ΔE‡ |
|---|---|---:|---:|---:|
| A | substrate | singlet | 0.00 | 0.00 |
| A | TS | singlet | 22.95 | barrier = 22.95 |
| A | product + SO | singlet SO (same-surface) | 16.70 | reaction = 16.70 |
| B | substrate | singlet | 10.14 | 0.00 |
| B | TS | singlet | 20.54 | barrier = 10.39 |
| B | product + SO | singlet SO (same-surface) | 17.47 | reaction = 7.33 |

## Mechanistic interpretation
- **Pathway A** leads to **1,5-dimethylpyrazole + SO**.
- **Pathway B** leads to **1,3-dimethylpyrazole + SO**.
- Pathway B starts from a **higher-energy substrate** (**10.14 kcal/mol** above A), but its TS is **lower in absolute energy** than Pathway A by **2.42 kcal/mol**.
- Consequently, Pathway B has the **lower local barrier** (**10.39 vs 22.95 kcal/mol**) and also the slightly **lower overall TS energy** relative to the global minimum substrate (**20.54 vs 22.95 kcal/mol**).

## Important model caveat
Isolated SO has a triplet ground state experimentally, while this mechanistic comparison was carried out on a **closed-shell singlet surface** to keep the TS search and both pathways internally consistent. `energies.csv` therefore includes both:
1. the **same-surface singlet-SO product energies** used for pathway comparison, and
2. a **triplet-SO reference** for the separated products.
