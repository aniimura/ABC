/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.TorsionIntegralGood
import ABC3.Meta.Claim

/-!
# 第 1337 ブロック —— **一意化の族は選択で取れる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——6 個の仮説が 1 本の選択に潰れる

`htFalt_veluQuotientFull_le`（第 704）と `lemma_3_5_velu`（在庫）は

    `P`・`Cv`・`hΔ`・`hPC`・`hell1`・`hell2`

の 6 個を仮説として受けていた。★これらは**すべて無条件に取れる**:

| 仮説 | 出どころ |
|---|---|
| `hell1` | `Δ` は単射準同型で移る |
| `P`・`Cv`・`hPC` | `exists_periodPair_of_isElliptic`（在庫、無条件） |
| `hΔ` | `latticeDisc_ne_zero`（在庫、無条件） |
| `hell2` | mathlib の `(C • W).IsElliptic` インスタンス |

☆したがって `lemma_3_5_velu` の入力は `hmin`・`hint`・`hss`・`hdeg` の 4 個だけになる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★**体の埋め込みで楕円性は移る**——★**無条件**（第 1337）。 -/
theorem isElliptic_map_of_ringHom (E : WeierstrassCurve L) [E.IsElliptic] (σ : L →+* ℂ) :
    (E.map σ).IsElliptic :=
  ⟨isUnit_iff_ne_zero.2 (by
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff σ σ.injective).2 E.isUnit_Δ.ne_zero)⟩

/-- ★★★★★★★★★★★★★★★★★★★★
**一意化の族は無条件に取れる**——★**無条件**（第 1337）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆各 `σ : L →+* ℂ` について `E ⊗_σ ℂ` は楕円だから、
`exists_periodPair_of_isElliptic`（在庫）が `P σ`・`Cv σ` を与える。
★選択公理で族にまとめるだけである。 -/
theorem exists_uniformization_family (E : WeierstrassCurve L) [E.IsElliptic] :
    ∃ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) ∧
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) ∧
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) ∧
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) := by
  have hell1 : ∀ σ : L →+* ℂ, (E.map σ).IsElliptic :=
    fun σ => isElliptic_map_of_ringHom E σ
  choose P Cv hPC using fun σ : L →+* ℂ =>
    exists_periodPair_of_isElliptic (E.map σ) (hell1 σ)
  refine ⟨P, Cv, fun σ => latticeDisc_ne_zero (P σ), hPC, hell1, fun σ => ?_⟩
  haveI := hell1 σ
  infer_instance

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`ht^Falt(E′) ≤ ht^Falt(E) + 2·log l`——一意化の族を外した形**（第 1337）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これで第 704 の入力は `hmin`（`E` は至る所 Néron 極小）と
`hint`（`E′` は整）の**2 個だけ**になった。 -/
theorem htFalt_veluQuotientFull_le' (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (l : ℕ) (hl : 0 < l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hmin : ∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0)
    (hint : ∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_uniformization_family E
  exact htFalt_veluQuotientFull_le E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 一意化の族を外した形**（第 1337）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★受ける入力は `hmin`・`hint`・`hss`・`hdeg` の **4 個だけ**である。 -/
theorem lemma_3_5_velu' (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), 0 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E') →
      degInfOf L E' = (l : ℝ) * degInfOf L E →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu eps heps
  refine ⟨C, ?_⟩
  intro L _ _ E E' _ _ l hl Q hQ hE' hmin hint hss hdeg
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_uniformization_family E
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin hint hss hdeg

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`ht^Falt(E′) ≤ ht^Falt(E) + 2·log l`——★無条件**（第 1338）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これで同種写像の高さ評価から仮説がすべて消えた**:

| かつての仮説 | 消し方 |
|---|---|
| `P`・`Cv`・`hΔ`・`hPC`・`hell1`・`hell2` | `exists_uniformization_family`（第 1337） |
| `hmin`・`hint` | `hfin`（第 1049 の欠損版）に弱めてある |
| `hfin` | `hfin_of_veluQuotientFull`（第 1083・1149、無条件） |

☆受けるのは `l` が素数・`Q` の位数が `l`・`E′` が Vélu の商、という原文の設定だけである。 -/
theorem htFalt_veluQuotientFull_le_uncond [inst : DecidableEq L] (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] {l : ℕ} (hl : l.Prime)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    htFaltOf L E' ≤ htFaltOf L E + 2 * Real.log l := by
  -- ★`Decidable` のインスタンスは一意なので古典的なものに揃える（配管）
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_uniformization_family E
  exact htFalt_veluQuotientFull_le_of_defect E E' l hl.pos Q hQ hE' P Cv hΔ hPC hell1 hell2
    (hfin_of_veluQuotientFull E E' hl Q hQ hE')

/-! ## ★出典の紐付け(`.src`) -/

def htFalt_veluQuotientFull_le_uncond.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ht^Falt(E′) ≤ ht^Falt(E) + 2 log l。★無条件)",
    sectionId := "genell-lemma-3-5" }

def htFalt_veluQuotientFull_le_uncond.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_veluQuotientFull_le_of_defect(第 1049、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_of_defect") 1,
    .citation "[ABC3]" "hfin_of_veluQuotientFull(第 1083・1149、無条件)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hfin_of_veluQuotientFull") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1338）**——同種写像の高さ評価から" ++
       "**仮説がすべて消えた**。" ++
       "☆これで witness の `faltingsHeight_quotLCyclic` は配管だけになる。") 3 ]

def isElliptic_map_of_ringHom.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(体の埋め込みで楕円性は移る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_uniformization_family.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一意化の族は無条件に取れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def htFalt_veluQuotientFull_le'.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ht^Falt の同種評価——一意化の族を外した形)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu'.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一意化の族を外した形——入力は hmin・hint・hss・hdeg の 4 個)",
    sectionId := "genell-lemma-3-5" }

def exists_uniformization_family.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_periodPair_of_isElliptic(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_periodPair_of_isElliptic") 1,
    .citation "[ABC3]" "latticeDisc_ne_zero(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeDisc_ne_zero") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1337）**——第 704 と `lemma_3_5_velu` が受けていた " ++
       "6 個の仮説（`P`・`Cv`・`hΔ`・`hPC`・`hell1`・`hell2`）が" ++
       "**1 本の選択に潰れた**。☆残る入力は `hmin`・`hint`・`hss`・`hdeg` である。") 2 ]

end ABC3.Found.GaloisRep
