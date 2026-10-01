# Conversation summary: BindCraft2 properties, pH handling, and new pH-sensitive loss

Summary of a Claude Code session in this repository. Sections 1–2 come from reading the source. Section 3 describes code that was added. No campaign was run.

## 1. What `--humanize` and `--protease-stable` change

Each flag expands to `<name>=true` (`bindcraft/cli.py`) and loads a preset from `settings/property/`. The preset turns on extra loss terms for the gradient design stage and adds final acceptance filters.

### `--humanize` (`settings/property/humanize.json`)
- Settings: `weights_humanization=1.0`, `humanization_species=human`, `humanization_mhc2_weight=1.0`, `humanization_coupling_weight=0.5`, `max_mhc_anchor_score_final=1.95`.
- Loss `humanization` (`bindcraft/loss.py`): slides a 9-residue window over the binder's amino-acid probabilities. It scores each window against MHC anchor-preference panels in `bindcraft/developability.py`.
  - Class I: 13 HLA-A/B alleles. Mouse H-2 alleles are available through `humanization_species`.
  - Class II: 7 HLA-DRB1 alleles.
  - Each window combines anchor match, anchor-pair coupling (`humanization_coupling_weight`) and an optional groove-burial hydrophobicity term (`humanization_hydro_weight`, default 0).
  - Soft maximum over windows and alleles, so the worst epitope dominates. Total = class I + `mhc_class_ii_weight` × class II.
- Only redesignable residues are scored (`redesignable_residue_mask`).
- Filter `MHC_Anchor_Score` (`bindcraft/filters.py`): same score on the argmax sequence, capped at 1.95.
- It does not pull the sequence toward human framework or germline sequences. It is an MHC-presentation proxy only. Some docs describe it as favouring "human-like sequence features", which is broader than the code. The docs themselves mark it "in development".

### `--protease-stable` (`settings/property/protease_stable.json`)
Three loss terms and three filters:

| Loss | Weight | Filter | Ceiling |
|---|---|---|---|
| `protease_sites` | 1.0 | `Protease_Site_Score` | 0.5 |
| `exposed_loops` | 0.5 | `Exposed_Loop_Fraction` | 0.7 |
| `exposed_termini` | 0.5 | `Terminus_Exposure` | 0.6 |

- `protease_sites`: per-position probability of being a protease P1 residue, weighted over trypsin (K/R), chymotrypsin (F/Y/W, L/M at 0.5), elastase (A/V/G/S etc., weight 0.4) and pepsin (F/L/W/Y, weight 0.5). A following proline blocks cleavage for the first three. Mean over redesignable residues.
- `exposed_loops`: penalises residues that are both solvent-exposed (soft neighbour count) and loop-like (low pLDDT, combined with a geometry test). β-strands are not counted as loops (`exposed_loops_distinguish_sheets`).
- `exposed_termini`: penalises exposure of the first and last 3 residues.
- The filters recompute the sequence score from the argmax sequence, and the structural scores from solvent-accessible area and secondary structure.

### Pipeline effect
- Both properties act only on the design-stage loss and the accept/reject filters.
- No protease or MHC logic appears in `MPNN_stage.py`, `proteinmpnn.py` or `bindcraft/mpnn/`. ProteinMPNN redesign is not biased by either flag, but its outputs must still pass the property filters.
- Both are proxies. Neither measures immunogenicity or serum stability. Extra objectives trade off against interface quality, so stack properties sparingly.

## 2. pH and protonation in BindCraft2 (before our changes)

- There was no setting for pH or protonation state, and no protonation or titration code.
- The only pH-related pieces were reporting and sequence priors:
  - `Binder_pI` and `Binder_Net_Charge` (`bindcraft/filters.py`) use Henderson–Hasselbalch pKa values (free His pKa 6.0). The charge is reported at a hard-coded `REPORTED_PH = 7.4`. They are report and filter metrics only, not losses.
  - `mpnn_variant` (`neutral`, `negative`, `positive`) picks the ProteinMPNN weights. Default `negative`.
  - `aa_bias` can favour or exclude residues, for example H.
  - `Interface_<AA>_Count` filter metrics (4 Å cutoff) count residues of one amino acid at the interface, for example `Interface_H_Count`.
- Changing `REPORTED_PH` to 6.4 would only change the single reported value. It would not give a 6.4-versus-7.5 comparison, so it is of limited use for a pH switch.
- No existing loss covered pH-dependent binding.

## 3. New `--ph-sensitive` property (implemented)

A new property was added following the same pattern as `--humanize` and `--protease-stable`. Four files were changed:

### Files changed

| File | Change |
|---|---|
| `bindcraft/loss.py` | Added `AMINO_ACID_INDEX` import; added `interface_his_acid_pairing` loss |
| `bindcraft/filters.py` | Added `Interface_His_Acid_Pairs` filter metric |
| `bindcraft/settings.py` | Added `min_interface_his_acid_pairs_final` to `CAMPAIGN_SETTING_NAMES`; added `CampaignFeature('pH sensitivity', ...)` |
| `settings/property/ph_sensitive.json` | New property preset |

### Loss: `interface_his_acid_pairing` (`loss.py`)

Registered as `@loss('interface_his_acid_pairing', target_weighting='binds_target')`. Enabled by `weights_interface_his_acid_pairing`.

- **Binder side (differentiable):** reads the soft probability of histidine at each redesignable position via `amino_acid_probabilities()`.
- **Target side (fixed):** reads D and E from the target's sequence. Since the target is not redesigned, `softmax(target.sequence)` collapses to one-hot, giving a 0/1 acid mask.
- **Proximity (differentiable):** slices the AF2 distogram for binder × target pairs and computes the fraction of probability mass within `cutoff` (default 8 Å).
- **Score per binder position:** `P(H at i) × Σⱼ [acid(j) × P(close to j)]`.
- Returns the **negative mean** over redesignable residues — minimising the loss rewards interface histidines near target acids.

### Filter: `Interface_His_Acid_Pairs` (`filters.py`)

Counts binder His residues within 4 Å (Cβ, or Cα for Gly) of a target D or E in the hard predicted structure. This is the acceptance gate.

### Settings integration (`settings.py`)

- `CampaignFeature('pH sensitivity', ...)` detects `weights_interface_his_acid_pairing` and wires a `ModalityCheck` with `higher=True` on `Interface_His_Acid_Pairs` keyed to `min_interface_his_acid_pairs_final`.
- `min_interface_his_acid_pairs_final` was added to `CAMPAIGN_SETTING_NAMES` so it passes the setting-validation check at campaign startup.

### Preset (`settings/property/ph_sensitive.json`)

```json
{
  "description": "Favour interface histidines paired with target acidic residues (D/E)...",
  "weights_interface_his_acid_pairing": 1.0,
  "min_interface_his_acid_pairs_final": 1
}
```

Use as `--ph-sensitive` on the CLI or `"ph_sensitive": true` in a campaign JSON.

### Bugs caught and fixed during review

| Bug | Would have caused | Fix |
|---|---|---|
| `if not acid_prob.any()` — Python conditional on a traced JAX array inside `jax.jit` | `ConcretizationTypeError` during gradient design stage | Removed the early exit; the math already returns 0 when `acid_prob` is all zeros |
| `min_interface_his_acid_pairs_final` missing from `CAMPAIGN_SETTING_NAMES` | "unrecognized campaign settings" error at startup | Added to the frozenset |

### Caveats

- This enriches the binder for the geometric ingredients of a pH switch (His near D/E). It does **not** predict pH-dependent affinity.
- The loss does not model protonation states, pKa shifts, or electrostatic complementarity.
- Candidates should be validated with external tools (PROPKA, APBS, constant-pH MD) or in the lab at both pH values.
- Like all properties, the extra objective trades off against interface quality.
- ProteinMPNN redesign is not biased by this loss, but its outputs must pass the `Interface_His_Acid_Pairs` filter.

## 4. Testing

- All changed files pass Python AST parsing.
- The preset appears in `settings/property/` and would be listed by `bindcraft design --list-properties`.
- No conda environment with both `jax` and `biotite` was available on this machine, so import-level and runtime tests could not be run. The implementation has not been tested in a real campaign.
