# TS-study: NH-pyrazole N-alkylation regioisomers (from TS.png, simplified)

## Simplifications applied (per instructions)
- Substituent on the ring carbon **α to the S-containing substituent** set to **Me**.
- All other peripheral ring substituents set to **H**.
- **N-alkyl group** installed in the alkylation step set to **Me**.

(Operationally, the model used is a methyl-substituted pyrazole bearing an S(=O)–NH2 substituent on the adjacent ring carbon; this preserves an "α-to-S" carbon on the ring, while keeping all other ring substituents as H.)

## Key chemical feature
The substrate has an **exchangeable ring N–H** that can tautomerize between the two ring nitrogens before alkylation. Each tautomer is neutral and compositionally identical, but leads to a different **N-methyl regioisomer**.

## Pathways inferred
### Pathway A
- **Substrate:** tautomer with N–H on ring nitrogen **N_A**.
- **Transformation:** N-methylation at **N_A** (replacement of N–H by N–Me).
- **Product:** N_A-methyl regioisomer.

SMILES used:
- Substrate A: `[nH]1ncc(C)c(S(=O)N)1`
- Product A:  `Cn1ncc(C)c(S(=O)N)1`

### Pathway B
- **Substrate:** tautomer with N–H on the other ring nitrogen **N_B**.
- **Transformation:** N-methylation at **N_B**.
- **Product:** N_B-methyl regioisomer.

SMILES used:
- Substrate B: `n1[nH]cc(C)c(S(=O)N)1`
- Product B:  `n1n(C)cc(C)c(S(=O)N)1`
