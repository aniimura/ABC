/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VeluImage
import ABC3.Found.GenEll.PointVariableChange

/-!
# 第 949 ブロック —— **★★★★★★★★★★★★★★★★Vélu の点集合の帳簿**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——**(D2) の 2 枚**

`j_veluQuotientFull_variableChange`（第 915）は 2 つのものを欲しがる:

1. `hS`——**点集合が反転で安定**であること
2. 変数変換先の点集合が `S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2))` の形であること

★本ブロックはその 2 枚を `⟨Q⟩∖{O}` について作る。

☆(1) は `-(k·Q) = (l-k)·Q` だけである——`l·Q = 0` なので。
☆(2) は `image_pointCoords_rhPoint_nsmul`（体準同型の場合）の**変数変換版**である。

| 定理 | 内容 |
|---|---|
| `pointCoords_neg` | ★`O` でなければ `-P` の座標は `(x, negY x y)` |
| `nsmul_eq_neg_nsmul_of_addOrderOf` | ★`(l-k)·Q = -(k·Q)` |
| `pointCoords_image_negY_stable` | ★★★★★★★★★★★★**`hS`——点集合は反転で安定** |
| `image_pointCoords_vcPoint_nsmul` | ★★★★★★★★★★★★★★★★**変数変換先の点集合** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine

section PointSet

variable {F : Type} [Field F]

/-- ★**`O` でなければ `-P` の座標は `(x, negY x y)`**。 -/
theorem pointCoords_neg {W : WeierstrassCurve F} {P : W.toAffine.Point} (hP : P ≠ 0) :
    pointCoords (-P)
      = ((pointCoords P).1,
         W.toAffine.negY (pointCoords P).1 (pointCoords P).2) := by
  rcases P with _ | ⟨x, y, h⟩
  · exact absurd rfl hP
  · rfl

/-- ★**`l·Q = 0` なら `(l-k)·Q = -(k·Q)`**（`k ≤ l`）。 -/
theorem nsmul_eq_neg_nsmul_of_addOrderOf {G : Type} [AddGroup G] {Q : G} {l k : ℕ}
    (hl : l • Q = 0) (hkl : k ≤ l) : (l - k) • Q = -(k • Q) := by
  refine eq_neg_of_add_eq_zero_left ?_
  rw [← add_nsmul, Nat.sub_add_cancel hkl, hl]

open scoped Classical in
/-- ★★★★★★★★★★★★**`⟨Q⟩∖{O}` の座標集合は反転で安定である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが `j_veluQuotientFull_variableChange`（第 915）の `hS` である。
☆`-(k·Q) = (l-k)·Q` で、`0 < k < l` なら `0 < l-k < l` だから。 -/
theorem pointCoords_image_negY_stable {W : WeierstrassCurve F}
    {Q : W.toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l) :
    ∀ z ∈ ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)),
      ((z.1, W.toAffine.negY z.1 z.2) : F × F)
        ∈ ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) := by
  intro z hz
  obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hz
  rw [Finset.mem_erase, Finset.mem_range] at hk
  have hkne : k • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hk.1 hk.2
  have hk0 : 0 < k := Nat.pos_of_ne_zero hk.1
  have hl0 : l • Q = 0 := by
    rw [← hQ]
    exact addOrderOf_nsmul_eq_zero Q
  refine Finset.mem_image.2 ⟨l - k, ?_, ?_⟩
  · rw [Finset.mem_erase, Finset.mem_range]
    omega
  · rw [nsmul_eq_neg_nsmul_of_addOrderOf hl0 hk.2.le, pointCoords_neg hkne]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**変数変換先の `⟨Q⟩∖{O}` の座標集合**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`image_pointCoords_rhPoint_nsmul`（体準同型版）の**変数変換版**である。
☆これが `j_veluQuotientFull_variableChange`（第 915）の右辺の形を与える。 -/
theorem image_pointCoords_vcPoint_nsmul (C : VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic]
    {Q : W.toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l) :
    ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C W Q))
      = (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
          (fun z : F × F => (vcX C z.1, vcY C z.1 z.2)) := by
  rw [Finset.image_image]
  refine Finset.image_congr ?_
  intro k hk
  rw [Finset.mem_coe, Finset.mem_erase, Finset.mem_range] at hk
  have hkne : k • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hk.1 hk.2
  simp only [Function.comp_apply]
  rw [← vcPoint_nsmul, pointCoords_vcPoint_of_ne C W hkne]

end PointSet

/-! ## ★★★★★★★★★★★★★★★★★★★★(D2)——商を Tate モデルの側へ移す -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**`⟨Q⟩` による Vélu の商の `j` は
変数変換先で計算しても同じ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 950）**——これが `Lemma 3.5` の **(D2)** である。
`E ⊗ Lv` の側で与えられた商を、Tate モデル `C • (E ⊗ Lv)` の側の商に移す。

☆道は 3 段:

1. 点集合は反転で安定（第 949）なので `j_veluQuotientFull_variableChange`（第 915）が使える
2. 変数変換先の点集合は `⟨vcPoint C W Q⟩∖{O}` の座標集合（第 949）
3. `j_congr_curve` で繋ぐ -/
theorem j_veluQuotientFull_nsmul_variableChange {F : Type} [Field F]
    (C : VariableChange F) (W E' : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic]
    {Q : W.toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l) (h2 : (2 : F) ≠ 0)
    (hE' : E' = veluQuotientFull W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    [E'.IsElliptic]
    [(veluQuotientFull (C • W)
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C W Q)))).IsElliptic] :
    E'.j = (veluQuotientFull (C • W)
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C W Q)))).j := by
  have himg := image_pointCoords_vcPoint_nsmul C W hQ
  have hcurve : veluQuotientFull (C • W)
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C W Q)))
      = veluQuotientFull (C • W)
        ((((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
          (fun z : F × F => (vcX C z.1, vcY C z.1 z.2))) := by
    rw [himg]
  haveI : (veluQuotientFull (C • W)
      ((((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
        (fun z : F × F => (vcX C z.1, vcY C z.1 z.2)))).IsElliptic := by
    rw [← hcurve]; infer_instance
  rw [j_veluQuotientFull_variableChange C W E' _ (pointCoords_image_negY_stable hQ) h2 hE']
  exact (j_congr_curve hcurve).symm

def j_veluQuotientFull_nsmul_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(⟨Q⟩ による Vélu の商の j は変数変換先で計算しても同じ。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★第 967 —— 大域の商を Tate モデル側の商の `j` に繋ぐ -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**`L` の上で与えた Vélu の商の `j` を、
底変換して変数変換した先の商の `j` に繋ぐ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 967）**——これが `isMuAtBadPrimes_of_veluQuotient` の
`hW′j` を作る段である。☆道は 2 段:

1. `veluQuotientFull_baseChange`——`E′ ⊗ K` は `E ⊗ K` の `⟨Q ⊗ K⟩` による商
2. `j_veluQuotientFull_nsmul_variableChange`（第 950）——それを `C • (E ⊗ K)` 側へ

★点の位数は体拡大で保たれる（`addOrderOf_rhPoint`、在庫）ので、
第 950 の `hQ` はそのまま満たされる。 -/
theorem j_map_velu_vcPoint {L K : Type} [Field L] [Field K] (φ : L →+* K)
    (C : WeierstrassCurve.VariableChange K)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [(E.map φ).IsElliptic]
    [(C • (E.map φ)).IsElliptic] [(E'.map φ).IsElliptic]
    {l : ℕ} {Q : E.toAffine.Point} (hQ : addOrderOf Q = l) (h2 : (2 : K) ≠ 0)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    [(veluQuotientFull (C • (E.map φ))
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C (E.map φ) (rhPoint φ E Q))))).IsElliptic] :
    (E'.map φ).j = (veluQuotientFull (C • (E.map φ))
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C (E.map φ) (rhPoint φ E Q))))).j := by
  have hQ' : addOrderOf (rhPoint φ E Q) = l := by rw [addOrderOf_rhPoint, hQ]
  exact j_veluQuotientFull_nsmul_variableChange C (E.map φ) (E'.map φ) hQ' h2
    (veluQuotientFull_baseChange φ E E' hQ hE')

def j_map_velu_vcPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(大域の Vélu の商の j を底変換・変数変換先の商の j に繋ぐ。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(O でなければ -P の座標は (x, negY x y)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def nsmul_eq_neg_nsmul_of_addOrderOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5((l-k)·Q = -(k·Q)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_image_negY_stable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(⟨Q⟩∖{O} の座標集合は反転で安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def image_pointCoords_vcPoint_nsmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換先の ⟨Q⟩∖{O} の座標集合。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
