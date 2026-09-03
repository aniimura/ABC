/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVelu

/-!
# 第 962 ブロック —— **★★★★★★★★★★★★★★★★Vélu の商の `Δ` は `l¹²` 倍**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——(D3) の (d5)

`tateParam_quot_velu_of_torsion`（第 948）は
`(veluCurve (tateCurveAt q hq) v w).map` の**楼円性**を受けている。

★それは `c₄`・`c₆` の帳簿から**そのまま出る**:

    `c₄(Vélu の商) = l⁴·c₄(E_{q^l})`、`c₆(Vélu の商) = l⁶·c₆(E_{q^l})`

なので、`1728Δ = c₄³ - c₆²` より

    `Δ(Vélu の商) = l¹²·Δ(E_{q^l})`

☆`l` が単元で `E_{q^l}` が楼円なら、Vélu の商も楼円である。

★★**同種写像の理論は一切要らない**——`c₄`・`c₆` の等式だけで済む。

| 定理 | 内容 |
|---|---|
| `Delta_velu_tate_eq` | ★★★★★★★★**`Δ = l¹²·Δ`** |
| `isElliptic_veluCurve_tate_map` | ★★★★★★★★★★★★★★★★**Vélu の商は楼円** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
  [IsAdicComplete I R]

/-- ★★★★★★★★**Vélu の商の `Δ` は `l¹²` 倍**。

☆`1728Δ = c₄³ - c₆²` と `c₄ ↦ l⁴c₄`・`c₆ ↦ l⁶c₆` だけである。
★`R` は `CharZero` の整域なので `1728` は消せる。 -/
theorem Delta_velu_tate_eq (q : R) (hq : q ∈ I) (l : ℕ) (hql : q ^ l ∈ I) (v w : R)
    (h4 : (tateCurveAt q hq).c₄ + 240 * v = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄)
    (h6 : (tateCurveAt q hq).c₆ + 504 * v + 6048 * w
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆) :
    (veluCurve (tateCurveAt q hq) v w).Δ = (l : R) ^ 12 * (tateCurveAt (q ^ l) hql).Δ := by
  have hc4 : (veluCurve (tateCurveAt q hq) v w).c₄
      = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄ := by
    rw [veluCurve_c₄]; exact h4
  have hc6 : (veluCurve (tateCurveAt q hq) v w).c₆
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆ := by
    rw [veluCurve_c₆, tateCurveAt_b₂, mul_one]; exact h6
  have e1 := WeierstrassCurve.c_relation (veluCurve (tateCurveAt q hq) v w)
  have e2 := WeierstrassCurve.c_relation (tateCurveAt (q ^ l) hql)
  have h1728ne : (1728 : R) ≠ 0 := by norm_num
  refine mul_left_cancel₀ h1728ne ?_
  calc (1728 : R) * (veluCurve (tateCurveAt q hq) v w).Δ
      = (veluCurve (tateCurveAt q hq) v w).c₄ ^ 3
        - (veluCurve (tateCurveAt q hq) v w).c₆ ^ 2 := e1
    _ = ((l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄) ^ 3
        - ((l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆) ^ 2 := by rw [hc4, hc6]
    _ = (l : R) ^ 12 * ((tateCurveAt (q ^ l) hql).c₄ ^ 3
        - (tateCurveAt (q ^ l) hql).c₆ ^ 2) := by ring
    _ = (l : R) ^ 12 * ((1728 : R) * (tateCurveAt (q ^ l) hql).Δ) := by rw [← e2]
    _ = (1728 : R) * ((l : R) ^ 12 * (tateCurveAt (q ^ l) hql).Δ) := by ring

/-- ★★★★★★★★★★★★★★★★**Vélu の商は楼円である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 962）**——これが (D3) の (d5) である。
`Δ(Vélu の商) = l¹²·Δ(E_{q^l})` で、`l` が単元で `E_{q^l}` が楼円なら、
右辺は単元である。 -/
theorem isElliptic_veluCurve_tate_map {K : Type} [Field K] [Algebra R K]
    (q : R) (hq : q ∈ I) (l : ℕ) (hql : q ^ l ∈ I) (v w : R)
    (hlu : IsUnit ((l : R)))
    (h4 : (tateCurveAt q hq).c₄ + 240 * v = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄)
    (h6 : (tateCurveAt q hq).c₆ + 504 * v + 6048 * w
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆)
    (hell : ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic) :
    ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic := by
  refine ⟨?_⟩
  rw [WeierstrassCurve.map_Δ, Delta_velu_tate_eq q hq l hql v w h4 h6, map_mul, map_pow]
  refine IsUnit.mul ((hlu.map (algebraMap R K)).pow 12) ?_
  have hΔ : ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).Δ
      = algebraMap R K ((tateCurveAt (q ^ l) hql).Δ) := WeierstrassCurve.map_Δ _ _
  rw [← hΔ]
  exact hell.isUnit

/-! ## ★出典の紐付け(`.src`) -/

def Delta_velu_tate_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の Δ は l¹² 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluCurve_tate_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は楼円である。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
