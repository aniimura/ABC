/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Found.GenEll.JScale

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

/-! ## ★★★★★★★★★★大域の Vélu の商は底変換で局所の Vélu の商になる -/

/-- ★★★★★★★★★★**`E′ = E/⟨Q⟩` を底変換すると `E′ ⊗ A = (E ⊗ A)/⟨Q ⊗ A⟩`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが「大域の `E′ = E/H` を各悪い素点で完備化に落とす」段である。
☆`veluQuotientFull_map`（底変換と可換）と
`image_pointCoords_rhPoint_nsmul`（点の像の座標）を並べるだけである。 -/
theorem veluQuotientFull_baseChange {F A : Type} [Field F] [Field A] (f : F →+* A)
    (E E' : WeierstrassCurve F) [E.IsElliptic] [(E.map f).IsElliptic]
    {l : ℕ} {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    E'.map f = veluQuotientFull (E.map f)
      (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • rhPoint f E Q))) := by
  rw [hE', veluQuotientFull_map, image_pointCoords_rhPoint_nsmul f E hQ]

/-! ## ★★★★★★★★反転で安定な点集合の上では奇関数の和は消える -/

/-- ★★**反転で安定な集合の上での奇関数の和は `0`**。

☆`∑ g = ∑ g∘σ = −∑ g` なので `2∑g = 0`、体で `2 ≠ 0` なら `∑g = 0`。 -/
theorem sum_eq_zero_of_neg_involution {F : Type} [Field F] (S : Finset (F × F))
    (σ : F × F → F × F) (hσS : ∀ z ∈ S, σ z ∈ S) (hσσ : ∀ z ∈ S, σ (σ z) = z)
    (g : F × F → F) (hg : ∀ z ∈ S, g (σ z) = - g z) (h2 : (2 : F) ≠ 0) :
    ∑ z ∈ S, g z = 0 := by
  have hswap : ∑ z ∈ S, g z = ∑ z ∈ S, g (σ z) := by
    refine Finset.sum_nbij' (i := fun z => σ z) (j := fun z => σ z) ?_ ?_ ?_ ?_ ?_
    · exact fun a ha => hσS a ha
    · exact fun a ha => hσS a ha
    · exact fun a ha => hσσ a ha
    · exact fun a ha => hσσ a ha
    · intro a ha
      exact congrArg g (hσσ a ha).symm
  have hneg : ∑ z ∈ S, g (σ z) = - ∑ z ∈ S, g z := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl hg
  have h0 : (2 : F) * ∑ z ∈ S, g z = 0 := by
    have hself := hswap.trans hneg
    linear_combination hself
  exact (mul_eq_zero.1 h0).resolve_left h2

/-- ★★★★★★**`veluQuotientFull_variableChange` の仮説 `hB`**——
反転で安定な点集合なら自動的に成り立つ。

☆`2y + a₁x + a₃` は `y ↦ negY x y` で符号が反転する。 -/
theorem sum_negY_eq_zero {F : Type} [Field F] (W : WeierstrassCurve F) (S : Finset (F × F))
    (hS : ∀ z ∈ S, ((z.1, W.toAffine.negY z.1 z.2) : F × F) ∈ S) (h2 : (2 : F) ≠ 0) :
    ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) = 0 := by
  refine sum_eq_zero_of_neg_involution S (fun z => (z.1, W.toAffine.negY z.1 z.2)) hS ?_ _ ?_ h2
  · intro z _
    have hnn : W.toAffine.negY z.1 (W.toAffine.negY z.1 z.2) = z.2 := by
      simp only [WeierstrassCurve.Affine.negY]; ring
    show ((z.1, W.toAffine.negY z.1 (W.toAffine.negY z.1 z.2)) : F × F) = z
    rw [hnn]
  · intro z _
    show 2 * (W.toAffine.negY z.1 z.2) + W.a₁ * z.1 + W.a₃ = -(2 * z.2 + W.a₁ * z.1 + W.a₃)
    rw [WeierstrassCurve.Affine.negY]
    show 2 * (-z.2 - W.toAffine.a₁ * z.1 - W.toAffine.a₃) + W.a₁ * z.1 + W.a₃
      = -(2 * z.2 + W.a₁ * z.1 + W.a₃)
    ring

/-- ★★★★★★**`veluQuotientFull_variableChange` の仮説 `hBx`**——同じ理由。

☆`x` は `σ` で不変なので、奇関数に偶関数を掛けても奇のままである。 -/
theorem sum_negY_mul_x_eq_zero {F : Type} [Field F] (W : WeierstrassCurve F)
    (S : Finset (F × F))
    (hS : ∀ z ∈ S, ((z.1, W.toAffine.negY z.1 z.2) : F × F) ∈ S) (h2 : (2 : F) ≠ 0) :
    ∑ Q ∈ S, (2 * Q.2 + W.a₁ * Q.1 + W.a₃) * Q.1 = 0 := by
  refine sum_eq_zero_of_neg_involution S (fun z => (z.1, W.toAffine.negY z.1 z.2)) hS ?_ _ ?_ h2
  · intro z _
    have hnn : W.toAffine.negY z.1 (W.toAffine.negY z.1 z.2) = z.2 := by
      simp only [WeierstrassCurve.Affine.negY]; ring
    show ((z.1, W.toAffine.negY z.1 (W.toAffine.negY z.1 z.2)) : F × F) = z
    rw [hnn]
  · intro z _
    show (2 * (W.toAffine.negY z.1 z.2) + W.a₁ * z.1 + W.a₃) * z.1
      = -((2 * z.2 + W.a₁ * z.1 + W.a₃) * z.1)
    rw [WeierstrassCurve.Affine.negY]
    show (2 * (-z.2 - W.toAffine.a₁ * z.1 - W.toAffine.a₃) + W.a₁ * z.1 + W.a₃) * z.1
      = -((2 * z.2 + W.a₁ * z.1 + W.a₃) * z.1)
    ring

/-! ## ★★★★★★★★★★Vélu の商を変数変換先のモデルへ移す -/

/-- ★★★★★★★★★★**Vélu の商の `j` は変数変換先で計算しても同じ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが「`E ⊗ Lv` の Vélu の商を Tate モデルへ移す」段である。
☆`veluQuotientFull_variableChange`（証明済み）の仮説 `hB`・`hBx` は、
点集合が反転で安定なら `sum_negY_eq_zero`・`sum_negY_mul_x_eq_zero`（第 912）が与える。 -/
theorem j_veluQuotientFull_variableChange {F : Type} [Field F]
    (C : WeierstrassCurve.VariableChange F) (W E' : WeierstrassCurve F)
    (S : Finset (F × F))
    (hS : ∀ z ∈ S, ((z.1, W.toAffine.negY z.1 z.2) : F × F) ∈ S) (h2 : (2 : F) ≠ 0)
    (hE' : E' = veluQuotientFull W S) [E'.IsElliptic]
    [(veluQuotientFull (C • W)
      (S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)))).IsElliptic] :
    E'.j = (veluQuotientFull (C • W)
      (S.image (fun Q => (vcX C Q.1, vcY C Q.1 Q.2)))).j := by
  have hVC := veluQuotientFull_variableChange C W S
    (sum_negY_eq_zero W S hS h2) (sum_negY_mul_x_eq_zero W S hS h2)
  refine ABC3.Found.GenEll.j_eq_of_smul_eq C E' _ ?_
  rw [hE']
  exact hVC.symm

/-! ## ★★★★★★★★★★`veluCurve` の形での移行 -/

/-- ★★★★★★★★★★**`veluCurve` の `j` は変数変換先で計算しても同じ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 915 は `veluQuotientFull`（点集合の形）だったが、
捧りを使うと商は `veluCurve` の形で出てくる（第 925）。
☆こちらが**非分裂の場合に使える形**である。

☆`v` は重み 4、`w` は重み 6（ただし `r` の分だけ `w − r·v` にずれる）。 -/
theorem j_veluCurve_variableChange {F : Type} [Field F]
    (C : WeierstrassCurve.VariableChange F) (W E' : WeierstrassCurve F) (v w : F)
    (hE' : E' = veluCurve W v w) [E'.IsElliptic]
    [(veluCurve (C • W) (((C.u⁻¹ : Fˣ) : F) ^ 4 * v)
      (((C.u⁻¹ : Fˣ) : F) ^ 6 * (w - C.r * v))).IsElliptic] :
    E'.j = (veluCurve (C • W) (((C.u⁻¹ : Fˣ) : F) ^ 4 * v)
      (((C.u⁻¹ : Fˣ) : F) ^ 6 * (w - C.r * v))).j := by
  refine ABC3.Found.GenEll.j_eq_of_smul_eq C E' _ ?_
  rw [hE']
  exact (veluCurve_variableChange C W v w).symm

/-! ## ★出典の紐付け(`.src`) -/

def j_veluCurve_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluCurve の j は変数変換先で計算しても同じ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def j_veluQuotientFull_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の j は変数変換先で計算しても同じ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_negY_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(反転で安定な点集合なら ∑(2y + a₁x + a₃) = 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(大域の E/⟨Q⟩ は底変換で局所の (E ⊗ A)/⟨Q ⊗ A⟩ になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluVFull_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の対の像にわたる v の和を番号の和に直す。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_image_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の対の像で取った Vélu の商 = veluCurve W v w の底変換。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
