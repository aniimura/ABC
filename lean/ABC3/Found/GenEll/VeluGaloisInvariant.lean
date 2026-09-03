/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Velu
import ABC3.Found.GaloisRep.SemiLinear
import ABC3.Meta.Claim

/-!
# 第 1153 ブロック —— **Vélu の和は `Gal` 不変**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`LCyclicReading` の節点 1 の核

第 1152 の枠（`Skeleton/GenEll/LCyclicReading.lean`）の節点 1 は

> `Gal`-安定な位数 `l` の部分群 `H ⊆ E(L̄)` に対し、Vélu の和
> `v = Σ veluV2`・`w` は `Gal` 不変なので `L` に降りる

であった。★本ファイルは**その「`Gal` 不変」の側を無条件に取る**。

## ★★★★★★★★機構

係数を固定する体自己同型 `Φ`（`FixesCoeffs`）に対し、`veluGx`・`veluGy` は
**係数と座標の多項式**なので

    `Φ(veluGx W x y) = veluGx W (Φx) (Φy)`,   `Φ(veluGy W x y) = veluGy W (Φx) (Φy)`

である。☆したがって座標の集合 `S` が `Φ` で保たれるなら、和

    `veluVFull W S = Σ_{Q ∈ S} veluV2(Q)`,   `veluWFull W S = Σ_{Q ∈ S} (veluU(Q)/2 + veluV2(Q)·x_Q)`

は `Φ` で不変である——添字を `Φ` で並べ替えるだけだからである。

★★**したがって `veluQuotientFull W S` の係数はすべて `Φ` で固定される**。

## ☆残る 1 歩（本ファイルでは取らない）

`L̄^{Gal(L̄/L)} = L` を使って「すべての `σ` で不変 ⟹ `L` の元」と言う段。
★これは mathlib の無限 Galois 理論の navigation であり、別の節点である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve

variable {M : Type} [Field M] (W : WeierstrassCurve M) {Φ : M ≃+* M}

/-! ## ★★★★★★点ごとの量の半線型同変性 -/

/-- ★★★★★**`veluGx` は半線型同変**——★**無条件**。 -/
theorem veluGx_semi (h : FixesCoeffs W.toAffine Φ) (x y : M) :
    Φ (veluGx W x y) = veluGx W (Φ x) (Φ y) := by
  simp only [veluGx, map_sub, map_add, map_mul, map_pow, map_ofNat, h.a₁, h.a₂, h.a₄]

/-- ★★★★★**`veluGy` は半線型同変**——★**無条件**。 -/
theorem veluGy_semi (h : FixesCoeffs W.toAffine Φ) (x y : M) :
    Φ (veluGy W x y) = veluGy W (Φ x) (Φ y) := by
  simp only [veluGy, map_sub, map_neg, map_mul, map_ofNat, h.a₁, h.a₃]

/-- ★★★★★**`veluV2` は半線型同変**——★**無条件**。 -/
theorem veluV2_semi (h : FixesCoeffs W.toAffine Φ) (x y : M) :
    Φ (veluV2 W x y) = veluV2 W (Φ x) (Φ y) :=
  veluGx_semi W h x y

/-- ★★★★★**`veluU` は半線型同変**——★**無条件**。 -/
theorem veluU_semi (h : FixesCoeffs W.toAffine Φ) (x y : M) :
    Φ (veluU W x y) = veluU W (Φ x) (Φ y) := by
  rw [veluU, veluU, map_pow, veluGy_semi W h]

/-! ## ★★★★★★★★★★★★座標の集合が保たれるなら和は不変 -/

/-- ☆座標の対に `Φ` を当てる写像。 -/
def semiPair (Φ : M ≃+* M) (z : M × M) : M × M := (Φ z.1, Φ z.2)

theorem semiPair_injective (Φ : M ≃+* M) : Function.Injective (semiPair Φ) := by
  rintro ⟨a, b⟩ ⟨c, d⟩ hcd
  simp only [semiPair, Prod.mk.injEq] at hcd
  rw [Prod.mk.injEq]
  exact ⟨Φ.injective hcd.1, Φ.injective hcd.2⟩

/-- ☆`Φ` で保たれる有限集合は `Φ` の像でちょうど自分自身になる（単射だから）。 -/
theorem image_semiPair_eq (Φ : M ≃+* M) (S : Finset (M × M))
    (hS : ∀ z ∈ S, semiPair Φ z ∈ S) [DecidableEq (M × M)] :
    S.image (semiPair Φ) = S :=
  Finset.eq_of_subset_of_card_le
    (fun z hz => by
      obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hz
      exact hS y hy)
    (le_of_eq (Finset.card_image_of_injective S (semiPair_injective Φ)).symm)

/-- ★★★★★★★★★★★★**`S` が `Φ` で保たれるなら `veluVFull` は `Φ` で不変**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆和の添字を `Φ` で並べ替えるだけである。 -/
theorem veluVFull_semi (h : FixesCoeffs W.toAffine Φ) (S : Finset (M × M))
    (hS : ∀ z ∈ S, semiPair Φ z ∈ S) :
    Φ (veluVFull W S) = veluVFull W S := by
  classical
  have hinj : ∀ x ∈ S, ∀ y ∈ S, semiPair Φ x = semiPair Φ y → x = y :=
    fun _ _ _ _ hxy => semiPair_injective Φ hxy
  have hsum : (∑ z ∈ S.image (semiPair Φ), veluV2 W z.1 z.2)
      = ∑ z ∈ S, veluV2 W (semiPair Φ z).1 (semiPair Φ z).2 := Finset.sum_image hinj
  rw [image_semiPair_eq Φ S hS] at hsum
  rw [veluVFull, map_sum]
  calc ∑ z ∈ S, Φ (veluV2 W z.1 z.2)
      = ∑ z ∈ S, veluV2 W (semiPair Φ z).1 (semiPair Φ z).2 :=
        Finset.sum_congr rfl (fun z _ => veluV2_semi W h z.1 z.2)
    _ = ∑ z ∈ S, veluV2 W z.1 z.2 := hsum.symm

/-- ★★★★★★★★★★★★**`S` が `Φ` で保たれるなら `veluWFull` は `Φ` で不変**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`/2` は `Φ 2 = 2` なので害が無い。 -/
theorem veluWFull_semi (h : FixesCoeffs W.toAffine Φ) (S : Finset (M × M))
    (hS : ∀ z ∈ S, semiPair Φ z ∈ S) :
    Φ (veluWFull W S) = veluWFull W S := by
  classical
  have hinj : ∀ x ∈ S, ∀ y ∈ S, semiPair Φ x = semiPair Φ y → x = y :=
    fun _ _ _ _ hxy => semiPair_injective Φ hxy
  have hsum : (∑ z ∈ S.image (semiPair Φ), (veluU W z.1 z.2 / 2 + veluV2 W z.1 z.2 * z.1))
      = ∑ z ∈ S, (veluU W (semiPair Φ z).1 (semiPair Φ z).2 / 2
          + veluV2 W (semiPair Φ z).1 (semiPair Φ z).2 * (semiPair Φ z).1) :=
    Finset.sum_image hinj
  rw [image_semiPair_eq Φ S hS] at hsum
  rw [veluWFull, map_sum]
  calc ∑ z ∈ S, Φ (veluU W z.1 z.2 / 2 + veluV2 W z.1 z.2 * z.1)
      = ∑ z ∈ S, (veluU W (semiPair Φ z).1 (semiPair Φ z).2 / 2
            + veluV2 W (semiPair Φ z).1 (semiPair Φ z).2 * (semiPair Φ z).1) :=
        Finset.sum_congr rfl (fun z _ => by
          rw [map_add, map_div₀, map_mul, veluU_semi W h, veluV2_semi W h, map_ofNat]
          rfl)
    _ = ∑ z ∈ S, (veluU W z.1 z.2 / 2 + veluV2 W z.1 z.2 * z.1) := hsum.symm

/-! ## ★★★★★★★★★★★★★★★★商曲線の係数はすべて固定される -/

/-- ★★★★★★★★★★★★★★★★
**`S` が `Φ` で保たれるなら Vélu の商の係数は `Φ` で固定される**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の節点 1 の核である
——`H` が `Gal`-安定なら `E/H` の係数は `Gal` で固定され、
`L̄^{Gal} = L` により `L` に降りる。☆最後の 1 歩は別の節点である。 -/
theorem fixesCoeffs_veluQuotientFull (h : FixesCoeffs W.toAffine Φ) (S : Finset (M × M))
    (hS : ∀ z ∈ S, semiPair Φ z ∈ S) :
    FixesCoeffs (veluQuotientFull W S).toAffine Φ := by
  refine ⟨h.a₁, h.a₂, h.a₃, ?_, ?_⟩
  · show Φ ((veluCurve W (veluVFull W S) (veluWFull W S)).a₄) = _
    rw [veluCurve_a₄, map_sub, map_mul, map_ofNat, veluVFull_semi W h S hS, h.a₄]
    rfl
  · show Φ ((veluCurve W (veluVFull W S) (veluWFull W S)).a₆) = _
    have hb₂ : Φ W.b₂ = W.b₂ := by
      rw [WeierstrassCurve.b₂, map_add, map_pow, map_mul, map_ofNat, h.a₁, h.a₂]
    rw [veluCurve_a₆, map_sub, map_sub, map_mul, map_mul, map_ofNat, hb₂,
      veluVFull_semi W h S hS, veluWFull_semi W h S hS, h.a₆]
    rfl

/-! ## ★出典の紐付け(`.src`) -/

def veluGx_semi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluGx は半線型同変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluGy_semi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluGy は半線型同変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def image_semiPair_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ で保たれる有限集合は Φ の像で自分自身。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluVFull_semi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の集合が Φ で保たれるなら veluVFull は不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluWFull_semi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標の集合が Φ で保たれるなら veluWFull は不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def fixesCoeffs_veluQuotientFull.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Gal-安定な部分群による Vélu の商の係数は Gal で固定される。★無条件)",
    sectionId := "genell-lemma-3-5" }

def fixesCoeffs_veluQuotientFull.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-01（第 1153）**——`Skeleton/GenEll/LCyclicReading.lean` の" ++
       "節点 1 の**核**である。☆残る 1 歩は `L̄^{Gal(L̄/L)} = L` を使って" ++
       "「すべての `σ` で固定 ⟹ `L` の元」と言う段であり、" ++
       "mathlib の無限 Galois 理論の navigation である。") 6 ]

end ABC3.Found.GenEll
