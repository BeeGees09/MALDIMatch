# MALDImatch

An R package for matching MALDI m/z peaks to LC-MS peptide identifications. Provides tools to compute peptide ion masses and isotopic distributions, build METASPACE-compatible databases, match MALDI peaks within a user-defined ppm tolerance, and visualise the overlap with Venn diagrams.

## Installation

```r
# Install dependencies first
install.packages(c("OrgMassSpecR", "data.table", "VennDiagram", "devtools"))


```r
devtools::install_github("YOUR_USERNAME/maldimatch")
```

---

## Quick Start

```r
library(maldimatch)

# --- Load data ---
lcms <- read.csv(
  "~/combined_modified_peptide.csv"
)

maldi <- read.table(
  "~/T-ReX_MALDI.csv",
  header = TRUE,
  sep    = ";",
  skip   = 8,
  fill   = TRUE
)
```

---

## Functions

### 1. `ComputeIon()`

Computes monoisotopic adduct masses, molecular formula, and isotopic distribution (M+0 to M+4) for each peptide in the LC-MS data frame.

```r
lcms_ann <- ComputeIon(lcms)
```

**Custom adducts** — pass any named numeric vector of mass offsets (Da):

```r
my_adducts <- c(
  "[M+H]+"   = 1.00728,
  "[M+Na]+"  = 22.98922,
  "[M+K]+"   = 38.96316,
  "[M+NH4]+" = 18.034
)
lcms_ann <- ComputeIon(lcms, adducts = my_adducts)
```

**Custom modifications** — extend the built-in PTM list:

```r
my_mods <- c(
  maldimatch:::.mod_adjust,
  list("Y_79.9663" = c(C=0, H=0, N=0, O=3, S=0))  # phosphotyrosine
)
lcms_ann <- ComputeIon(lcms, mod_adjust = my_mods)
```

**Save a simplified CSV:**

```r
lcms_ann <- ComputeIon(
  lcms,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/peptide_MH_formula_isotopes.csv"
)
```

Columns added to the returned data frame:

| Column | Description |
|---|---|
| `MH_plus` | Monoisotopic [M+H]+ m/z |
| `_M_Na__plus`, etc. | Additional adduct m/z columns (named from `adducts`) |
| `Formula` | Molecular formula string |
| `Max_Peak_mz` | m/z of the most-abundant isotope peak |
| `Max_Peak_Label` | Label of the most-abundant isotope peak (e.g. `M+2`) |
| `M0_mz` - `M4_mz` | m/z for isotope peaks M+0 to M+4 |
| `M0_%` - `M4_%` | Relative intensity (%) for isotope peaks M+0 to M+4 |

---

### 2. `DbLCMS()`

Builds a METASPACE-compatible tab-delimited database from the annotated LC-MS data frame.

```r
db <- DbLCMS(
  lcms_ann,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/Kidney_ECM_LCMS_DB.tsv"
)
```

Requires columns `Gene`, `Modified.Sequence`, and `Formula` (all present after `ComputeIon()`).  
Output columns: `id`, `name` (`Gene|ModifiedSequence`), `formula`.

---

### 3. `MALDIMatches()`

Matches a vector of MALDI-observed m/z values to LC-MS peptides using both the monoisotopic (M+0) and the most-abundant isotope peak, within a ppm tolerance.

```r
matches <- MALDIMatches(
  maldi_mz    = maldi$m.z,
  lcms        = lcms_ann,
  tolerance   = 10,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/Matched_peptides_KPMP_ECM.csv"
)
```

| Argument | Description |
|---|---|
| `maldi_mz` | Numeric vector of MALDI observed m/z values |
| `lcms` | Annotated data frame from `ComputeIon()` |
| `tolerance` | ppm threshold (default `10`) |
| `output_path` | Optional path to save results as CSV |

Returns a `data.table` with: `mz_MALDI`, `closest_value_LCMS`, `ppm_error`, `iso_peak_matched`, plus all original LC-MS annotation columns.

---

### 4. `PlotVenn()`

Draws a pairwise Venn diagram showing the overlap between MALDI m/z values and LC-MS peptide ions.

```r
PlotVenn(
  maldi_mz    = maldi$m.z,
  lcms        = lcms_ann,
  matches     = matches,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/mz_overlap_venn.png"
)
```

| Argument | Description |
|---|---|
| `maldi_mz` | Numeric vector of MALDI m/z values (sets total MALDI count) |
| `lcms` | Annotated data frame from `ComputeIon()` (sets total LC-MS count) |
| `matches` | `data.table` returned by `MALDIMatches()` (sets intersection count) |
| `output_path` | Optional path to save PNG at 300 dpi |
| `colors` | Fill colours for the two circles (default: `c("lightblue", "lightgreen")`) |
| `...` | Additional arguments passed to `VennDiagram::draw.pairwise.venn()` |

---

## Full Workflow

```r
library(maldimatch)

# Load raw data
lcms <- read.csv(
  "C:/Users/gorm567/OneDrive - PNNL/Documents/PNNL/Data/KPMP/ECM/output/combined_modified_peptide.csv"
)
maldi <- read.table(
  "F:/Data/KPMP/ECM/ECM_TISAC/T-ReX_v1.csv",
  header = TRUE, sep = ";", skip = 8, fill = TRUE
)

# 1. Compute ion masses
lcms_ann <- ComputeIon(
  lcms,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/peptide_MH_formula_isotopes.csv"
)

# 2. Build METASPACE database
db <- DbLCMS(
  lcms_ann,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/Kidney_ECM_LCMS_DB.tsv"
)

# 3. Match MALDI to LC-MS (10 ppm)
matches <- MALDIMatches(
  maldi_mz    = maldi$m.z,
  lcms        = lcms_ann,
  tolerance   = 10,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/Matched_peptides_KPMP_ECM.csv"
)

# 4. Venn diagram
PlotVenn(
  maldi_mz    = maldi$m.z,
  lcms        = lcms_ann,
  matches     = matches,
  output_path = "F:/Data/KPMP/ECM/ECM_TISAC/mz_overlap_venn.png"
)
```

---

## Built-in Modifications

| Key | Modification | Elements added |
|---|---|---|
| `C_57.0215` | Carbamidomethyl (C) | +C2H3N1O1 |
| `M_15.9949` | Oxidation (M) | +O1 |
| `P_15.9949` | Oxidation (P) | +O1 |
| `n_42.0106` | Acetylation (N-term) | +C2H2O1 |

Pass a custom `mod_adjust` list to `ComputeIon()` to extend or override these.

---

## Dependencies

- [`OrgMassSpecR`](https://cran.r-project.org/package=OrgMassSpecR) — monoisotopic mass and isotopic distribution calculations
- [`data.table`](https://cran.r-project.org/package=data.table) — fast merge of match results
- [`VennDiagram`](https://cran.r-project.org/package=VennDiagram) — Venn diagram rendering
- [`grid`](https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid-package.html) — graphics display

---

## Author

Brittney Gorman — brittney.gorman@pnnl.gov

