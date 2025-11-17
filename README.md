# autapomorphicLoci
Script to get autapomorphic loci / SNPs from multiple FASTA alignments.

# Autapomorphy Scanner for Population-Assigned FASTA Loci

**Detect autapomorphic SNPs per population** from a directory of aligned DNA FASTA files (one locus per file), with strict per-site coverage requirements and an optional **permutation (null) mode**.

---

## ✨ What this tool does

A site is called **autapomorphic** for population *P* when:

1. it is **polymorphic** among `A/T/C/G`;
2. **exclusivity**: the candidate base occurs in **exactly one** population;
3. **fixation**: within that population, all valid bases at the site are **identical** to the candidate base;
4. **coverage**: every population has at least `min_samples_per_pop` valid bases at that site.

It supports:

* **Regex mode** — extract population labels from sequence IDs.
* **Lookup mode** — map sequence IDs to populations via a two-column file.

**Both modes** apply a **locus-level filter**: a FASTA is kept only if *every* population present in that file has at least `min_samples_per_pop` sequences overall.

---

## 🧠 Why it matters

Autapomorphic SNPs are useful **diagnostic molecular phenotypes** in taxonomic inference. Sparse or uneven sampling can inflate “unique” states. The **permutation mode** randomizes population labels by site to estimate how many autapomorphies you would expect **by chance** given the same data structure. A small (or near-zero) null count strengthens the claim that observed diagnostics are **not** sampling artifacts.

---

## 📦 Requirements

* **Perl 5**
* **BioPerl** (`Bio::SeqIO`)
* **List::Util** (core; for `shuffle`)

> **Important:** All sequences within a FASTA must be **aligned** and have **equal length**. Only `A/T/C/G` are considered valid (gaps/ambiguities are ignored).

---

## 🚀 Usage

```bash
perl getAutopomorphic.pl <fastas_directory> <pop_regex | lookup_file> <min_samples_per_pop> <shuffle_reps>
```

**Arguments**

* `<fastas_directory>`
  Directory with aligned DNA FASTA files. Files matched by `^.+\.fa[sta]*$`.

* `<pop_regex | lookup_file>`

  * **Regex mode:** if this is **not** a readable file, it is treated as a Perl regex applied to each sequence ID.
    The population label is built from **`$1$2`**, so your regex should define **at least one** capture group (the second is optional).
    Example: `'(POP\d+)'` extracts `POP1` from `POP1_indivA`.
  * **Lookup mode:** if this **is** a readable file, it must contain two whitespace-separated columns:
    `sequence_id  population_label`
    Every sequence ID must appear in this file.

* `<min_samples_per_pop>`
  Integer ≥ 1. Per-site **coverage threshold** (valid `A/T/C/G`) required **in every population**.

* `<shuffle_reps>`

  * `0` → deterministic report.
  * `R > 0` → **permutation mode** with `R` replicates.

---

## 🧪 Inclusion filters

* **Locus-level (both modes):** keep a FASTA only if **every population present in that file** has **≥ `min_samples_per_pop` sequences** overall.
* **Site-level (both modes):** test a site only if **every population** has **≥ `min_samples_per_pop` valid bases (`A/T/C/G`)** at that position.

---

## 📤 Output

### Deterministic report (`shuffle_reps = 0`)

Tab-delimited table to **STDOUT**:

```
Locus    count_<pop1>    pos_<pop1>    base_<pop1>    count_<pop2>    pos_<pop2>    base_<pop2> ...
```

* One row per locus with ≥ 1 autapomorphy.
* `count_<pop>` — number of autapomorphic sites for that population in that locus.
* `pos_<pop>` — **1-based** site positions, **pipe-separated** (e.g., `5|17|20`).
* `base_<pop>` — concatenated bases in the same order (e.g., `ACT`).

> Progress and diagnostics are printed to **STDERR**.

### Permutation (null) mode (`shuffle_reps > 0`)

* Prints **one line** with `R` space-separated integers: the number of autapomorphic **events** per replicate.
* **Note:** the statistic is **events**, not “loci with ≥1 event”. Multi-allelic sites can contribute >1 event.

---

## 📘 Examples

**Regex mode (deterministic)**

```bash
perl getAutopomorphic.pl loci/ '(POP\d+)' 3 0 > autapomorphies.tsv
```

**Lookup mode (permutation, R=1000)**

```bash
perl getAutopomorphic.pl loci/ samples2pop.txt 2 1000 > null_counts.txt
```

---

## 📈 Interpreting the permutation output

Summarize the `R` values (mean/median/IQR/95% CI) and compute an **empirical p-value** comparing the observed statistic (from the deterministic run) with the null distribution:

[
p = \frac{1 + #{\text{replicates with } T_{\text{null}} \ge T_{\text{obs}}}}{1 + R}
]

* **Observed ≫ null** (small *p*): diagnostics unlikely by chance → supports taxonomic diagnosis.
* **Observed ~ null**: pattern may arise from sampling/labeling; interpret with caution.

> If you prefer a **per-locus (0/1)** null statistic instead of “events”, adapt the permutation loop to mark presence/absence per locus per replicate.

---

## ⚠️ Caveats & tips

* Keep only **DNA** FASTAs in the directory; the filename pattern is permissive.
* Ambiguities/gaps are ignored and reduce per-site coverage.
* LD between sites is not modeled; permutations are done **per site**.
* In regex mode, the script dies if the regex can’t extract a population label; in lookup mode, it dies if a sequence ID is missing from the mapping.

---

## 📚 Citation

If this tool contributes to your work, please cite this repository and also:

da Cunha Nunes, J.C., Polo, É., Capurucho, J.M.G. et al. **Cryptic diversity in the Amazonian Red-shouldered Tanager (Tachyphonus phoenicius) Swainson 1838 with a neotype designation and revalidation of Tachyphonus saucius Strickland 1844**. J Ornithol (2025). https://doi.org/10.1007/s10336-025-02344-7

---

## 🛠️ License

Specify your preferred license (e.g., MIT, GPL-3.0) here.

---

*Questions or feature requests?* Open an issue or PR.
