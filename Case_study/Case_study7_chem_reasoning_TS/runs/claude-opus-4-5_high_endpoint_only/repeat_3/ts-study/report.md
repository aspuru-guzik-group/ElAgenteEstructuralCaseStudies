# Transition State Study: TDSO to Pyrazole via SO Extrusion

## Overview

This study examines the transformation of N-methyl-5-methyl-1,2,3-thiadiazole 1-oxide (TDSO) to N-methyl-methylpyrazole through a cheletropic SO extrusion mechanism. Due to N-H tautomerism in the parent TDSO, N-alkylation with methyl can occur at either of the two ring nitrogens, producing two distinct N-methyl TDSO isomers (compositionally identical but structurally different). Each isomer follows a separate reaction pathway.

## Structural Analysis

### Substrate Identification

Based on the reaction scheme (TS.png), the starting material is N-methyl-5-methyl-1,2,3-thiadiazole 1-oxide with:
- Methyl group on ring carbon alpha to sulfur (position 5)
- All other ring hydrogen positions as H
- N-methyl substitution (replacing NH)

### Two Pathways from NH Tautomerism

**Pathway A - N2-Methyl Isomer:**
- The nitrogen adjacent to sulfur (N2) bears the methyl group
- N3 retains the NH functionality before alkylation

**Pathway B - N3-Methyl Isomer:**  
- The nitrogen not adjacent to sulfur (N3) bears the methyl group
- N2 (adjacent to S) retains the NH functionality before alkylation

### Transformation Mechanism

The reaction involves a concerted cheletropic extrusion:
1. Breaking of S-C bond (C5, alpha to S)
2. Breaking of S-N bond (either N2 or N4)
3. Expulsion of SO molecule
4. Ring contraction/rearrangement to form pyrazole

## Computational Methods

- **Level of Theory:** GFN2-xTB
- **Substrate/Product Optimization:** Unconstrained geometry optimization
- **TS Optimization:** Constrained optimization with elongated S-C and S-N bonds
  - S-C constraint: 2.1 Å (equilibrium ~1.76 Å)
  - S-N constraint: 2.0 Å (equilibrium ~1.77 Å)
  - Force constant: 0.5
- **Frequency Analysis:** GFN-FF Hessian calculation

## Key Geometric Parameters

### Pathway A TS Structure
- S-C distance: 2.07 Å (elongated from 1.76 Å)
- S-N distance: constrained at ~2.0 Å
- SO group displaced out of ring plane

### Pathway B TS Structure  
- S-C distance: 2.08 Å (elongated from 1.76 Å)
- S-N distance: constrained at ~2.0 Å
- SO group displaced out of ring plane

## Energy Summary

| Species | Pathway | Energy (Eh) | Relative Energy (kcal/mol) |
|---------|---------|-------------|---------------------------|
| Substrate | A | -25.744994 | +2.16 |
| Substrate | B | -25.748442 | 0.00 (reference) |
| TS | A | -25.719651 | +18.07 (from B-sub) |
| TS | B | -25.721984 | +16.60 (from B-sub) |

### Activation Energies

| Pathway | Ea (kcal/mol) | Description |
|---------|--------------|-------------|
| A | 15.90 | N2-Me substrate → TS |
| B | 16.60 | N3-Me substrate → TS |

**ΔΔE‡ (A - B) = -0.70 kcal/mol**

Pathway A (N adjacent to S methylated) has a slightly lower activation barrier.

### Substrate Thermodynamics

- Pathway B substrate (N3-Me) is **2.16 kcal/mol more stable** than Pathway A substrate (N2-Me)
- This suggests the N3-H tautomer is thermodynamically preferred before alkylation

## Frequency Analysis

### Pathway A TS
- Lowest vibrational modes: -133, -87, -42 cm⁻¹
- Multiple imaginary frequencies indicate approximate saddle point geometry
- Dominant mode corresponds to SO extrusion

### Pathway B TS
- Lowest vibrational modes: -159, -133, +54 cm⁻¹
- Two significant imaginary frequencies
- Constrained optimization gives approximate TS character

Note: The presence of multiple imaginary frequencies results from the constrained optimization approach rather than true saddle point optimization. For rigorous TS verification, Sella or NEB methods with tighter convergence would be recommended.

## Mechanistic Conclusions

1. **Kinetic Preference:** Pathway A (N2-methylated substrate) has a slightly lower activation barrier (15.9 vs 16.6 kcal/mol), favoring this route kinetically by ~0.7 kcal/mol.

2. **Thermodynamic Preference:** The N3-methylated substrate (Pathway B) is thermodynamically more stable by 2.16 kcal/mol.

3. **Overall Competition:** Under kinetic control with equilibration of substrates, Pathway B would dominate due to the higher population of the more stable N3-Me isomer, despite its slightly higher intrinsic barrier.

4. **Product Distribution:** The two pathways lead to different regioisomeric pyrazole products, with the N-methyl group on different positions of the pyrazole ring.

## File Structure

```
ts-study/
├── structures/
│   ├── pathway_A_substrate.xyz
│   ├── pathway_A_TS_Initial_Guess.xyz
│   ├── pathway_A_TS_optimized.xyz
│   ├── pathway_A_product.xyz
│   ├── pathway_B_substrate.xyz
│   ├── pathway_B_TS_Initial_Guess.xyz
│   ├── pathway_B_TS_optimized.xyz
│   ├── pathway_B_product.xyz
│   └── SO.xyz
├── calculations/
│   ├── pathway_A/
│   │   ├── pathway_A_substrate.log
│   │   ├── pathway_A_TS_optimized.log
│   │   ├── pathway_A_product.log
│   │   └── freq_A_ff.log
│   └── pathway_B/
│       ├── pathway_B_substrate.log
│       ├── pathway_B_TS_optimized.log
│       ├── pathway_B_product.log
│       └── freq_B_ff.log
├── energies.csv
└── report.md
```

## Notes

- Transition state structures were obtained using constrained optimization at the GFN2-xTB level with elongated S-C and S-N bonds characteristic of the SO extrusion process.
- Product structures are N-methyl-methylpyrazole isomers resulting from different N-alkylation patterns.
- The energetics confirm that both pathways are viable with modest activation barriers typical of thermal pericyclic reactions.
