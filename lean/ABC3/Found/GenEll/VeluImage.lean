/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Velu

/-!
# 第 885 ブロック —— **★★★★★★★★★★番号で書いた部分群の Vélu の商**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

`veluQuotientFull W S` の `S` は**座標の対の Finset**であるが、
`c4_velu_tate`・`c6_velu_tate`（第 853・867）の側は
**`(range l).erase 0` にわたる和**で書かれている。

★本ファイルはその両者を繋ぐ:

    `veluQuotientFull (W ⊗ K) (s.image (i ↦ (X i, Y i))) = (veluCurve W v w) ⊗ K`

だだし `v = ∑_{i∈s} v₂(X i, Y i)`、`2w = ∑_{i∈s} (u + 2·v₂·X)`。

☆鍵は `Finset.sum_image`（単射なら像の上の和は元の上の和）と
`veluV2_map`・`veluU_map`（環準同型と可換）だけである。

| 定理 | 内容 |
|---|---|
| `veluVFull_image` | ★★`v` の側 |
| `veluWFull_image` | ★★`2w` の側（`/2` を避けるため 2 倍で述べる） |
| `veluQuotientFull_image_eq` | ★★★★★★★★★★**商の一致** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve Finset
open scoped Classical

variable {R : Type} [CommRing R] {K : Type} [Field K] [Algebra R K]

/-- ★★**`v` の側**——番号にわたる和に直す。 -/
theorem veluVFull_image (W : WeierstrassCurve R) (s : Finset ℕ) (X Y : ℕ → R)
    (hinj : ∀ i ∈ s, ∀ j ∈ s,
      ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)
        = ((algebraMap R K (X j), algebraMap R K (Y j)) : K × K) → i = j) :
    veluVFull (W.map (algebraMap R K))
        (s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)))
      = algebraMap R K (∑ i ∈ s, veluV2 W (X i) (Y i)) := by
  rw [veluVFull, Finset.sum_image hinj, map_sum]
  exact Finset.sum_congr rfl (fun i _ => (veluV2_map (algebraMap R K) W (X i) (Y i)).symm)

/-- ★★**`2w` の側**——`veluWFull` の `/2` を避けるため 2 倍で述べる。 -/
theorem veluWFull_image (W : WeierstrassCurve R) (s : Finset ℕ) (X Y : ℕ → R)
    (hinj : ∀ i ∈ s, ∀ j ∈ s,
      ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)
        = ((algebraMap R K (X j), algebraMap R K (Y j)) : K × K) → i = j)
    (h2 : (2 : K) ≠ 0) :
    2 * veluWFull (W.map (algebraMap R K))
        (s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)))
      = algebraMap R K
          (∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) := by
  rw [veluWFull, Finset.sum_image hinj, map_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [map_add, map_mul, map_mul, map_ofNat, ← veluU_map, ← veluV2_map]
  have hu : (2 : K) * (algebraMap R K (veluU W (X i) (Y i)) / 2)
      = algebraMap R K (veluU W (X i) (Y i)) := by
    field_simp
  linear_combination hu

/-- ★★★★★★★★★★**座標の対の像で取った Vélu の商は
`veluCurve W v w` の底変換である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが大域の `veluQuotientFull` と局所の `veluCurve` を繋ぐ一本である。 -/
theorem veluQuotientFull_image_eq (W : WeierstrassCurve R) (s : Finset ℕ) (X Y : ℕ → R)
    (hinj : ∀ i ∈ s, ∀ j ∈ s,
      ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)
        = ((algebraMap R K (X j), algebraMap R K (Y j)) : K × K) → i = j)
    (v w : R) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ s, veluV2 W (X i) (Y i))
    (hw : 2 * w = ∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) :
    veluQuotientFull (W.map (algebraMap R K))
        (s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)))
      = (veluCurve W v w).map (algebraMap R K) := by
  have hV : veluVFull (W.map (algebraMap R K))
      (s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)))
      = algebraMap R K v := by
    rw [veluVFull_image W s X Y hinj, hv]
  have hW : veluWFull (W.map (algebraMap R K))
      (s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)))
      = algebraMap R K w := by
    have h := veluWFull_image W s X Y hinj h2
    rw [← hw, map_mul, map_ofNat] at h
    exact mul_left_cancel₀ h2 h
  rw [veluQuotientFull, veluCurve_map, hV, hW]

/-! ## ★出典の紐付け(`.src`) -/

def veluVFull_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の対の像にわたる v の和を番号の和に直す。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_image_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の対の像で取った Vélu の商 = veluCurve W v w の底変換。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
