/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Uniformization
import ABC3.Meta.Claim

/-!
# 第 1330 ブロック —— **格子曲線の Vélu の商は楕円**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——**Vélu の定理が在庫から出た**

第 1328-1329 で「`VeluQuotOK` の根は Vélu の定理 1 本」と測ったが、
★**解析側の在庫にそれがあった**:

| 在庫 | 内容 |
|---|---|
| `exists_isogeny_periodPair`（証明済み） | 位数 `l` の点から大きい格子 `Λ′` を作る |
| `exists_velu_rep`（証明済み） | 代表系 `T`（`0 ∈ T`、`|T| = l`） |
| `g₂_isogeny`・`g₃_isogeny`（証明済み） | `g₂(Λ′)`・`g₃(Λ′)` が Vélu の和で書ける |
| `latticeCurve_eq_veluQuotientFull`（証明済み） | `latticeCurve Λ′ = veluQuotientFull (latticeCurve Λ) S` |
| `isElliptic_latticeCurve′`（証明済み） | 格子曲線はつねに楕円 |

☆これらを繋ぐと **`veluQuotientFull (latticeCurve P) S` は楕円**である。

★★★これで `VeluQuotOK` の楕円性は、`ℂ` へ埋め込めば出る道が見えた
——数体上の曲線は `ℂ` に埋め込め（`veluQuotientFull_map` で `Δ` が移る）、
`ℂ` 上の楕円曲線は格子曲線に変数変換できる（`exists_periodPair_of_isElliptic`、在庫）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**格子曲線の Vélu の商は大きい格子の曲線**——★**無条件**（第 1330）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆在庫の 4 本（`exists_isogeny_periodPair`・`exists_velu_rep`・`g₂_isogeny`・`g₃_isogeny`）を
`latticeCurve_eq_veluQuotientFull` に渡すだけである。 -/
theorem exists_latticeCurve_eq_veluQuotientFull (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (P' : PeriodPair) (T : Finset ℂ), (0 : ℂ) ∈ T ∧ T.card = l ∧
      latticeCurve P' = veluQuotientFull (latticeCurve P)
        ((T.erase 0).image (fun w => (latticePointX P w, latticePointY P w))) := by
  classical
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  exact ⟨P', T, h0T, hcard,
    latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**格子曲線の Vélu の商は楕円**——★**無条件**（第 1330）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが「Vélu の定理」の、本プロジェクトで必要な形である。 -/
theorem exists_isElliptic_veluQuotientFull_lattice (P : PeriodPair)
    (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ T : Finset ℂ, (0 : ℂ) ∈ T ∧ T.card = l ∧
      (veluQuotientFull (latticeCurve P)
        ((T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)))).IsElliptic := by
  obtain ⟨P', T, h0T, hcard, heq⟩ :=
    exists_latticeCurve_eq_veluQuotientFull P hΔ hl hQ
  refine ⟨T, h0T, hcard, ?_⟩
  rw [← heq]
  exact isElliptic_latticeCurve' P'

/-! ## ★出典の紐付け(`.src`) -/

def exists_latticeCurve_eq_veluQuotientFull.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子曲線の Vélu の商は大きい格子の曲線。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isElliptic_veluQuotientFull_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子曲線の Vélu の商は楕円——Vélu の定理の必要な形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isElliptic_veluQuotientFull_lattice.needs : List ProofObligation :=
  [ .citation "[ABC3]" "latticeCurve_eq_veluQuotientFull(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeCurve_eq_veluQuotientFull") 1,
    .citation "[ABC3]" "exists_isogeny_periodPair・exists_velu_rep・g₂_isogeny・g₃_isogeny(在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_isogeny_periodPair") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1330）**——第 1328-1329 で「根は Vélu の定理 1 本」と測ったが、" ++
       "**解析側の在庫にそれがあった**。☆これで `VeluQuotOK` の楕円性は" ++
       "`ℂ` へ埋め込めば出る道が見えた。") 3 ]

end ABC3.Found.GenEll
