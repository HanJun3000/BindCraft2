# AtomWorks annotation-driven conditioning for BindCraft2 — Tier 1 plan

**Status:** proposed, not implemented.
**Drafted:** 2026-10-29, against `main` at `b6fb3ae`.
**Audience:** whoever implements this, with no prior context on the investigation.

---

## 1. Goal

Let a campaign take its binding-site conditioning from **annotations carried inside the target
structure file** instead of from residue numbers typed into the campaign JSON — the way
RFD4-Proteina reads an annotated mmCIF.

Today:

```json
{"targets": [{"name": "hPDL1", "target_path": "hPDL1.pdb", "chains": "A", "hotspots": "54,56,66,115"}]}
```

After this change, either of:

```bash
bindcraft design design.json --structure-annotations     # auto-detect, opt-in
```
```json
{"targets": [{"name": "hPDL1", "target_path": "hPDL1_annotated.cif"}],
 "structure_annotations": true}
```
```json
{"targets": [{"name": "hPDL1", "target_path": "x.cif", "hotspots": "@condition_hotspot_1_atom"}]}
```

The explicit `@column` form needs no flag; only auto-detection is gated.

Scope is **targets only** (`HOTSPOT`, `COLDSPOT`). Binder scaffolds are Tier 2, see §8.

---

## 2. Why this is a small change

BC2 already has a per-residue annotation layer. `ResidueFlags` (`bindcraft/protein.py:42`) is an
`IntFlag`, and every `Protein` carries a per-residue `uint8` `flags` array:

```
DESIGN  TEMPLATE  SEQUENCE  CONTACT  HOTSPOT  COLDSPOT  CYCLIC  PADDING
```

**Nothing downstream of ingest ever sees a residue number.** Every consumer reads the flags:

| Consumer | Location |
| --- | --- |
| hotspot pLDDT loss | `loss.py:196` |
| `binder_coldspot` loss | `loss.py:468` |
| `coldspot_repel` loss | `loss.py:505` |
| interface-contact loss | `loss.py:551` |
| `Hotspot_Contact_Fraction` | `filters.py:209` |
| `Coldspot_Contact_Fraction` | `filters.py:219` |
| off-epitope contact | `filters.py:417` |
| forced-targeting epitope | `epitope_targeting.py:53` |

Residue numbers exist in exactly one place — a string compiler that runs at ingest:

```
"hotspots": "54,56,66,115"
  -> target_binding_site_flags()      protein_preparation.py:21   "A54+HOTSPOT,A56+HOTSPOT,..."
  -> annotated_residue_flags()        protein.py:712
  -> scaffold_edit_spans() / _EDIT_RE protein.py:675 / :657       -> flags array
```

So this task is **a second compiler into the same array**. Losses, filters and metrics are
untouched. The flag model is already close to RFD4's `C_HOT` / `C_SEQ` / `C_CRD`; it is simply
populated from text today rather than from the file.

---

## 3. Findings (all verified, see §3.4 to reproduce)

### 3.1 RFD4's annotated-CIF column convention

Taken from `RFD4-Proteina/data/benchmarks/ppi/ppi_all_atom/pdl1.cif.zst`. The `atom_site` category
carries the 21 standard mmCIF columns plus:

```
_atom_site.annotation_expsegmin_1_residue     "-1" x924, "50" x1
_atom_site.annotation_expsegmax_1_residue
_atom_site.condition_sequence_1_residue       three-letter codes (identity reveal)
_atom_site.condition_index_1_residue          "True" x925
_atom_site.annotation_src-component_1_atom
_atom_site.condition_hotspot_1_atom           "True" x6, "False" x919
_atom_site.condition_hbond_1_atom
_atom_site.condition_index-group_1_token
_atom_site.atom_array_index
```

Four consequences for the design:

1. **The hotspot column is `condition_hotspot_1_atom`, not `is_hotspot`.** The pattern is
   `condition_<name>_<n>_<level>` / `annotation_<name>_<n>_<level>`, with level in
   `{atom, residue, token}` and `<n>` an instance counter. Auto-detection must **prefix-match**,
   not exact-match.
2. **Values are the Python `bool` string form** `"True"` / `"False"` — what AtomWorks' writer
   produces. A tolerant parser is still wanted (`1`/`0`, `.`/`?` for CIF-absent).
3. **Hotspots are atom-level and sparse** — 6 atoms over 3 residues (B40, B99, B107). The
   atom-to-residue reduction (a residue is flagged when *any* of its atoms is) is therefore
   mandatory, not a convenience. This matches RFD4's training definition: a sparse interface patch.
4. **RFD4 has no coldspot condition.** `is_coldspot` / `condition_coldspot_` is a BC2 convention we
   are free to define.

### 3.2 BC2 already ingests an RFD4 annotated CIF unmodified

BC2's real ingest chain (`read_structure_atoms` -> `polymer_atom_records` ->
`atom37_arrays_from_records`) was run verbatim against that file:

```
chain labels  : ['A', 'B']
polymer chains: ['B']
atoms read    : 925
atoms kept    : 924
  chain B: 115 residues, resid 1-115, seq starts NAFTVTVPKDLY
```

Chain A is RFD4's **collapsed expandable-segment placeholder**: one `HETATM` whose `auth_comp_id`
is literally the string `"50-120"` (the 50-120 residue binder to be designed). BC2's existing
hetero/no-backbone filter at `protein.py:630` drops it silently, *before* the unsupported-residue
check at `:631` can see it. `structure_chain_names` then reports only `['B']`, so the user does not
even need to supply `"chains"`.

**No change to `polymer_atom_records` is required.** This was initially expected to be a blocker and
is not.

### 3.3 Custom annotations survive BC2's exact reader call

An annotated CIF written by AtomWorks 3.0.4 (`to_cif_file(..., extra_fields=["is_hotspot"])`) read
back through the same call BC2 makes at `protein.py:339`:

```python
pdbx.get_structure(block, model=1, extra_fields=["b_factor", "is_hotspot"])
# annotations: [..., 'b_factor', 'is_hotspot']
# hotspot res_id round-tripped: [54, 56, 66, 115]
```

Biotite 1.6 reads **arbitrary** `atom_site` columns as string annotations — see the final loop of
`biotite.structure.io.pdbx.convert._fill_annotations`:

```python
# Handle all remaining custom fields
for field in extra_fields:
    array.set_annotation(field, atom_site[field].as_array(str))
```

**Therefore BC2 needs no AtomWorks dependency to read annotated CIFs.** An annotated CIF is just an
mmCIF with extra columns, and BC2 already uses biotite's pdbx reader. AtomWorks is only needed to
*author* such files.

Caveat: `atom_site[field]` raises `KeyError` on an absent column, so membership must be checked
first. Minor: AtomWorks' writer emits no `B_iso_or_equiv`, so b-factors come back `NaN`; only
`score.py:76` consumes b-factors, and only for scored designs, so this is harmless here.

### 3.4 Reproducing the findings

```bash
# 1. RFD4 column names and values
zstd -dc RFD4-Proteina/data/benchmarks/ppi/ppi_all_atom/pdl1.cif.zst > /tmp/pdl1.cif
grep '^_atom_site\.' /tmp/pdl1.cif

# 2. BC2 ingest of that file - AST-extract the jax-free ingest functions from protein.py
#    (read_structure_atoms, polymer_atom_records, atom37_arrays_from_records, ...) and run them.
#    They depend only on biotite + numpy, so any env with biotite will do.

# 3. Annotation round trip - needs atomworks
#    load_any() -> aa[0] -> set_annotation("is_hotspot", ...) -> to_cif_file(extra_fields=[...])
#    -> pdbx.CIFFile.read -> pdbx.get_structure(..., extra_fields=[...])
```

---

## 4. Settings surface

Both forms live in the **existing** per-target `hotspots` / `coldspots` fields. No new target key.

| `hotspots` value | Behaviour |
| --- | --- |
| *(empty)* | auto-detect convention columns, **only when `structure_annotations` is on**; nothing found = nothing happens, silently |
| `"54,56"` | residue numbers win. If the file also carries a convention column, print that it was ignored |
| `"@condition_hotspot_1_atom"` | that exact column. Error if absent |
| `"@col,A200"` | union of both |

`@` is safe as a sigil: `_EDIT_RE` (`protein.py:657`) is `(?P<chain>[A-Za-z]+)(?P<start>\d+)...`, so
`@` can never begin a valid residue span, and `cleaned_setting` (`preflight.py:59`) leaves it intact
(`hotspots` is in `COMMA_SEPARATED_SETTINGS`, which only strips whitespace within spans).

Mixing works for free: both forms compile into the same flag array and are OR-ed.

### Auto-detection is a property preset (decided)

Auto-detection is **off by default** and turned on by a flag. BC2 derives its boolean CLI flags
automatically from the files in `settings/property/` — `design_property_flags()` at `cli.py:60` is
literally `{'--' + name.replace('_', '-'): name for name in shipped_preset_names('property')}`. So
adding one file:

```json
// settings/property/structure_annotations.json
{
  "description": "Read hotspots and coldspots from annotations carried in the target structure",
  "structure_annotations": true
}
```

gives `--structure-annotations` **and** the JSON entry `"structure_annotations": true` with no new
CLI code. `settings/property/forced_targeting.json` is the direct precedent: a property preset that
changes how hotspots are used.

Turning it off is omitting the flag. There is no `--no-*` machinery in BC2 and none is added.

An explicit `@col` is a declaration, not auto-detection, so it works regardless of the flag.

**Precedence rule, stated once:** an explicit declaration always wins over the file, and the
override is announced rather than silent. There is no path by which a file silently overrides what
the campaign JSON says.

---

## 5. Implementation

Roughly 140 lines across 6 files. No new dependency. No change to `Protein`'s shape, the losses,
the filters, or `TARGET_SETTING_NAMES`.

### 5.1 `bindcraft/protein.py` — the substance (~90 lines)

New module constants:

```python
ANNOTATION_COLUMN_PREFIX = '@'

#RFD4 writes condition_<name>_<n>_<level>; a plain AtomWorks annotation is usually is_<name>
CONVENTION_ANNOTATION_COLUMNS: dict[ResidueFlags, tuple[str, ...]] = {
    ResidueFlags.HOTSPOT:  ('condition_hotspot_', 'is_hotspot'),
    ResidueFlags.COLDSPOT: ('condition_coldspot_', 'is_coldspot'),
}

#every atom_site column a coordinate reader already names; the rest are annotations
STANDARD_ATOM_SITE_COLUMNS = frozenset({
    'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
    'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code',
    'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
    'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'})

ANNOTATION_TRUE  = frozenset({'true', '1', 'yes', 'y', 't'})
ANNOTATION_FALSE = frozenset({'false', '0', 'no', 'n', 'f', '', '.', '?'})
```

Verify `STANDARD_ATOM_SITE_COLUMNS` against the RFD4 fixture: excluding it must leave exactly the
8 `annotation_*`/`condition_*` columns plus `atom_array_index`.

1. **`read_structure_atoms` (`:332`)** — read every `atom_site` column not in
   `STANDARD_ATOM_SITE_COLUMNS` as a string array and attach them to each atom record under a new
   `'annotations'` key (`{}` for PDB input, which cannot carry columns).

   Read **all** custom columns unconditionally rather than adding an `extra_fields` parameter: this
   keeps the `functools.lru_cache` key as `(source,)`, so asking for a different column later does
   not re-read the file. Cost on the RFD4 fixture is 9 string columns over 925 atoms.

2. **New `annotation_residue_flags(chain_records, residue_index, column_flags) -> list[int]`** —
   sibling of `annotated_residue_flags` (`:712`). Reduces atom-level values to residue level by OR
   (any flagged atom flags the residue), keyed on `res_id` exactly as `atom37_arrays_from_records`
   (`:645`) does. Raises on a value in neither `ANNOTATION_TRUE` nor `ANNOTATION_FALSE`, naming the
   column, the value and the offending atom.

3. **New `structure_annotation_columns(source) -> tuple[str, ...]`** — the custom columns a file
   carries. Used for auto-detection and for error messages.

4. **New `resolved_annotation_columns(source, selection, flag) -> tuple[str, dict[str, ResidueFlags]]`**
   — turns one `hotspots`/`coldspots` value plus the file into `(residue_spans_text, {column: flag})`,
   applying the §4 precedence table. Prefix-matches the convention tuples; exact-matches a `@name`.

5. **`Protein.from_structure` (`:215`)** — new keyword
   `annotations: dict[str, ResidueFlags] | None = None`. OR its per-residue result into the existing
   `annotated_residue_flags(...)` output. The edit-string path is otherwise untouched, so the
   current behaviour cannot regress.

Keep 2-4 generic over `ResidueFlags` so Tier 2 is a table entry rather than a rewrite.

### 5.2 `bindcraft/protein_preparation.py` (~25 lines)

- **`target_binding_site_flags` (`:21`)** — skip spans beginning with `@`. Without this,
  `chain_qualified_residue_span` prefixes a chain name (`A@is_hotspot`), `_EDIT_RE` fails to match
  it, and `scaffold_edit_spans` raises "carry ... outside any span". Filtering here keeps that raise
  as a genuine consistency check that nothing leaked through.
- **New `target_annotation_columns(target, settings)`** — wraps `resolved_annotation_columns` for
  both flags and returns the merged `{column: flag}` mapping.
- **`prepare_targets` (`:76`)** — pass `annotations=` through to `from_structure`. Raise if a FASTA
  target names `@col` (a sequence carries no columns). Extend the existing coldspot `print` to name
  the source and the resolved residue count for both flags:

  ```
  target=hPDL1 hotspots=@condition_hotspot_1_atom residues=3 (B40,B99,B107) coldspots=none
  ```

  Note this print already fires once per trajectory for coldspots; match that rather than changing it.

### 5.3 `bindcraft/settings.py` + `settings/property/` (3 lines, 1 new file)

- add `structure_annotations` to `CAMPAIGN_SETTING_NAMES` (`:37`) so it is a recognised setting
- add `settings/property/structure_annotations.json` (contents in §4); this creates the
  `--structure-annotations` flag automatically
- read it as `settings.get('structure_annotations', False)` — **no** entry in
  `settings/core/default.json`, so the default stays off without a line to maintain

`TARGET_SETTING_NAMES` (`:39`) needs **no** change — `hotspots` and `coldspots` are already there.
`PreparedTargetState` keeps carrying the declaration string unchanged, so `campaign_facts`
(`campaign_output.py:510`) and `campaign_log.py:162` record what was *asked for*
(`hPDL1:@condition_hotspot_1_atom`) with no change. The resolved residues are recorded separately
by §5.4.

### 5.4 `bindcraft/campaign_output.py` (~15 lines) — closing the round trip

**This is the non-obvious bug in a naive implementation.** `rank.py:214` re-scores finished designs:

```python
score_design(structure, hotspots=target.get('hotspots', ''), coldspots=target.get('coldspots', ''))
```

Hand it `@condition_hotspot_1_atom` and it looks for that column on a **BC2 output** structure,
which has none — so `Hotspot_Contact_Fraction` silently becomes `None`.

BC2 already has the pattern for this. `designed_span_stamp` (`:539`) writes binder `DESIGN`/`CONTACT`
into the CIF `bindcraft` category as residue spans, and `score.py:29` (`redesigned_binder_flags`)
reads them back. So:

- generalize **`flagged_residue_spans` (`:529`)** to take `start` / `stop` / `residue_offset`, so a
  merged multi-chain receptor can be split per `WrittenChain` (`protein.py:written_chains`) — a
  merged receptor is one `Protein` written as several chains, and the spans must follow that split;
- extend **`designed_span_stamp` (`:539`)** to also stamp target `hotspot_residues` /
  `coldspot_residues` from the non-binder chains.

Single caller, `campaign.py:216`, so the change is local.

This also gives reproducibility for free: the campaign records the *declaration*, each design
records the *resolved residues*.

### 5.5 `bindcraft/score.py` (~5 lines)

Mirror `redesigned_binder_flags` (`:29`): in `design_states` (`:33`), when `hotspots`/`coldspots` are
not given, read `hotspot_residues` / `coldspot_residues` from the structure stamp via
`read_structure_metadata`. `bindcraft score` then works on an annotation-driven design with no
arguments, and `rank.py` is fixed without being touched.

---

## 6. Validation and error messages

Most of this is free: `design_input_problems` (`preflight.py:155`) already calls `prepare_targets`
inside `try/except (ValueError, KeyError, OSError)`, so every raise below surfaces as a preflight
error before the campaign starts.

Raise on:

- `@col` against a `.pdb` — PDB format cannot carry extra columns
- `@col` against a FASTA target
- `@col` absent from the file — **list the custom columns it does carry**
- a value in neither `ANNOTATION_TRUE` nor `ANNOTATION_FALSE`
- **a resolved column matching zero residues in the selected chains** — this is the check that
  catches a real mistake, e.g. selecting chain A of the RFD4 fixture (the placeholder) instead of B

Do **not** raise when auto-detection finds no column: that is the normal case for an ordinary PDB
and must stay silent.

Print (not raise) when an explicit numeric declaration causes a present convention column to be
ignored.

---

## 7. Docs

- `README.md` — the "Binding site" bullet under *Design against your own target*
- `docs/source/reference.md#targeting-options`
- `examples/pdl1_annotated.json` plus its annotated CIF

Generate the example fixture from the shipped `settings/target/structures/hPDL1.pdb` rather than
copying from the RFD4 checkout, so the repository stays self-contained.

---

## 8. Out of scope, and where Tier 2 goes

Not in this change:

- **Binder-scaffold annotations** (Tier 2). The scaffold presets carry edit strings like
  `"A26!,A27-33,A34!,A35*,...,A52-57(5-7),..."` (`settings/modality/VHH.json`). Annotating
  `scaffolds/VHH.cif` directly would retire them, and §5.1's generic design makes the flag mapping
  trivial — but the `(5-7)` **length sampling cannot be expressed per residue**. RFD4 solves exactly
  this with `annotation_expsegmin_1_residue` / `annotation_expsegmax_1_residue`; copy that mechanism
  when tackling Tier 2.
- **Any AtomWorks dependency.** Reading needs only biotite (§3.3). AtomWorks is for authoring.
- **`condition_sequence_*` -> `ResidueFlags.SEQUENCE`**, and the other RFD4 conditions.
- **Emitting annotation columns on BC2 output.** The stamp (§5.4) closes the round trip without them.
  Writing real columns would make BC2 output valid RFD4/AtomWorks input, which is a nice follow-up.
- **Replacing BC2's structure reader with AtomWorks `parse()`** (assemblies, CCD chemistry,
  missing-atom modelling, `sel()` queries). Separable; do it for those features, not for annotations.

---

## 9. Decisions taken

1. **Auto-detect is opt-in, off by default,** enabled by the `--structure-annotations` flag or the
   equivalent JSON entry, realised as a property preset (§4). This avoids a stray annotated CIF
   silently conditioning a campaign, and costs no CLI code.
2. **Setting name** — `structure_annotations`.

Nothing is left open. Implement as written.

---

## 10. Verification checklist

The repository has **no test suite** (CI only builds docs), so verify by hand:

- [ ] Existing numeric `hotspots` campaigns are byte-identical in behaviour — run
      `examples/pdl1.json` preflight before and after and compare resolved settings
- [ ] `examples/il2_receptor.json` still resolves chain-prefixed hotspots (`A67,...,B229-231`)
      across a merged two-chain receptor
- [ ] RFD4 fixture, chain B auto-detected: resolves to 3 hotspot residues (B40, B99, B107)
- [ ] RFD4 fixture, chain A selected: fails preflight with the zero-residue message
- [ ] AtomWorks-written `is_hotspot` fixture auto-detects
- [ ] `@condition_hotspot_1_atom` explicit form resolves identically to auto-detection
- [ ] `@nonexistent_column` fails preflight and lists the file's real columns
- [ ] `@col` on a `.pdb` and on a FASTA target both fail preflight
- [ ] numeric hotspots beside a present convention column: numbers win, note printed
- [ ] auto-detection does nothing without the flag; `--structure-annotations` enables it
- [ ] `@col` resolves with the flag absent (explicit declaration is never gated)
- [ ] `bindcraft design --help` lists `--structure-annotations` among the property flags
- [ ] a completed design's CIF stamp carries `hotspot_residues`
- [ ] `bindcraft score <design.cif>` with no `--hotspots` reports `Hotspot_Contact_Fraction`
- [ ] `rank.py` re-scoring produces a non-`None` `Hotspot_Contact_Fraction` for an
      annotation-driven campaign
