/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# `ht^Falt` は `j` だけで決まる（同じ体の上で）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`EllModuliData` の witness で `EllClass := ℂ`（`j` 不変量）を取るには

> **`cls E = cls E′` なら `ht^Falt(E) = ht^Falt(E′)`**

が要る（`ResearchPaper/ellmoduli-witness-status.json` の `designChoice`）。
★本ファイルはその**同じ体の上での形**を無条件で証明する。

## ★★★★★捻りの議論は要らなかった

☆着手前は「同じ `j` の 2 曲線は捻りだから、`√d` を添加した体へ上げて同型にする」
という筋を見込んでいた。★★実際には**両側とも `j` だけで書ける**ので、その必要はない:

| 側 | 量 | `j` だけで決まる理由 |
|---|---|---|
| アルキメデス | `archNorm E σ = curveArchInv (E×σℂ)` | `curveArchInv_congr_j`（在庫、`CurveArchInv.lean:154`） |
| 有限素点 | `minDeltaExp p E = max(0, −v_p(j))` | `minDeltaExp_eq_maxJ_of_semistable`（`§9-1165`、第 738） |

★有限素点の側だけが半安定性を要る（捻りで還元型が変わっても `j` は変わらないから、
半安定でない曲線では `v_p(Δ_min)` は `j` だけでは決まらない）。

## ★本ファイルで取れるもの

| 定理 | 内容 |
|---|---|
| `archNorm_congr_j` / `archSum_congr_j` | ★★★アルキメデス側は `j` だけで決まる（★無条件） |
| `jExp_congr_j` | ★`v_p(j)` は当然 `j` だけで決まる |
| `minDeltaExp_congr_j_of_semistable` | ★★★★半安定なら `v_p(Δ_min)` も |
| `degInfOf_congr_j_of_semistable` | ★★★★★`deg∞` は `j` だけで決まる |
| `htFaltOf_congr_j_of_semistable` | ★★★★★★★★**`ht^Falt` は `j` だけで決まる** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll IsDedekindDomain NumberField WeierstrassCurve

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★アルキメデス側 -/

/-- ★★★★★**`archNorm` は `j` だけで決まる**——★無条件。

★`archNorm E σ = curveArchInv (E×σℂ)` であり、`ℂ` 上の曲線の
アルキメデス不変量は同型不変量（`curveArchInv_congr_j`）だからである。 -/
theorem archNorm_congr_j (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hj : E.j = E'.j) (σ : L →+* ℂ) : archNorm E σ = archNorm E' σ :=
  curveArchInv_congr_j inferInstance inferInstance
    (by rw [WeierstrassCurve.map_j, WeierstrassCurve.map_j, hj])

/-- ★★★★★★**`archSum` は `j` だけで決まる**——★無条件。 -/
theorem archSum_congr_j (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hj : E.j = E'.j) : archSum L E = archSum L E' := by
  rw [archSum, archSum]
  exact Finset.sum_congr rfl (fun σ _ => by rw [archNorm_congr_j E E' hj σ])

/-! ## ★★★★★★有限素点側 -/

/-- ★★`v_p(j)` は当然 `j` だけで決まる。 -/
theorem jExp_congr_j (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic] (hj : E.j = E'.j) : jExp p E = jExp p E' := by
  rw [jExp, jExp]
  by_cases h : E.j = 0
  · have h' : E'.j = 0 := by rw [← hj]; exact h
    rw [dif_pos h, dif_pos h']
  · have h' : E'.j ≠ 0 := by rw [← hj]; exact h
    rw [dif_neg h, dif_neg h']
    congr 1
    exact Units.ext hj

/-- ★★★★★★**半安定なら `v_p(Δ_min)` は `j` だけで決まる**。 -/
theorem minDeltaExp_congr_j_of_semistable (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : SemistableAt p E) (hss' : SemistableAt p E') (hj : E.j = E'.j) :
    minDeltaExp p E = minDeltaExp p E' := by
  rw [minDeltaExp_eq_maxJ_of_semistable p E hss, minDeltaExp_eq_maxJ_of_semistable p E' hss',
    jExp_congr_j p E E' hj]

/-- ★★★★★★★**半安定なら `deg∞` は `j` だけで決まる**。 -/
theorem degInfOf_congr_j_of_semistable (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') (hj : E.j = E'.j) :
    degInfOf L E = degInfOf L E' := by
  rw [degInfOf, degInfOf]
  congr 1
  exact finsum_congr (fun p => by
    rw [minDeltaExp_congr_j_of_semistable p E E' (hss p) (hss' p) hj])

/-! ## ★★★★★★★★★★★★★★★★到達点 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `ht^Falt` は `j` だけで決まる**（同じ体の上で）——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これが `EllModuliData` の `faltingsHeight : EllClass → ℝ` を
`EllClass := ℂ`（`j` 不変量）の上で定義できる根拠である
（`ResearchPaper/ellmoduli-witness-status.json` の `designChoice`）。

☆残るのは**体が違う場合**——`emb₁(j₁) = emb₂(j₂)` の 2 つの `SSCurve` を
共通の体へ上げる段である（`htFaltOf_baseChange_of_semistable`、`§9-1165`、第 738 が受け皿）。 -/
theorem htFaltOf_congr_j_of_semistable (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') (hj : E.j = E'.j) :
    htFaltOf L E = htFaltOf L E' := by
  rw [htFaltOf, htFaltOf, degInfOf_congr_j_of_semistable E E' hss hss' hj,
    archSum_congr_j E E' hj]

/-! ## ★出典の紐付け(`.src`) -/

def archSum_congr_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(archSum は j だけで決まる。★無条件)",
    sectionId := "genell-prop-3-4" }

def htFaltOf_congr_j_of_semistable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(半安定な曲線の ht^Falt は j だけで決まる——同じ体の上で。★無条件)",
    sectionId := "genell-prop-3-4" }

def htFaltOf_congr_j_of_semistable.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは体が違う場合——emb₁(j₁) = emb₂(j₂) の 2 つの SSCurve を" ++
       "共通の数体へ上げる段である(第 738 の htFaltOf_baseChange_of_semistable が受け皿)") 5 ]

end ABC3.Found.GaloisRep
