import ABC3.Found.GenEll.Prop16

/-!
# [GenEll] Proposition 1.4, (iii) —— **計量を取り替えても高さは BD-同値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★何を取るか

> **因子が同じで計量が連続なら、`|ht_D(x) − ht_E(x)| ≤ C`(`C` は `F` にも `x` にも依らない)**

★これが `Proposition 1.4, (iii)`(生成ファイバーが同じなら高さは `≈`)の
**アルキメデス側の中身**である。因子表示では「生成ファイバーが同じ」は
「因子が同じ(垂直成分の差を除く)」であり、計量の違いだけが残る。

## ★★機構は `Proposition 1.6` と同じ

1. `X^arc` はコンパクト(`ArcModel.compactSpace`)
2. 連続な Green 関数は**上にも下にも**有界(`ArcModel.exists_bound`)
3. 正規化された和は `[F:ℚ]` で割るので、評価が `F` に依らない

★`ArchBound.lean` は下界(`archADiv_sum_ge`)しか持っていなかったので、
**上界を鏡像として足す**。

## ★本ファイルで取れるもの

| 定理 | 内容 |
|---|---|
| `archADiv_sum_le` | ★アルキメデス側の寄与の上界 |
| `archADiv_sum_div_finrank_le` | ★正規化した上界(`F` に依らない) |
| `htArith_eq_add` | ★高さ = 有限素点の寄与 + 正規化したアルキメデス側の寄与 |
| `htArith_sub_abs_le` | ★★★★**因子が同じなら高さの差は一様に有界** |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★アルキメデス側の寄与の上界 -/

/-- ★**Green 関数が `C` 以下なら、寄与は `C · [F : ℚ]` 以下**(`archADiv_sum_ge` の鏡像)。 -/
theorem archADiv_sum_le {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (C : ℝ) (hg : ∀ p, g p ≤ C) :
    (archADiv F g xF).sum (fun _ r => r) ≤ C * (Module.finrank ℚ F : ℝ) := by
  classical
  rw [Finsupp.sum]
  have hsub : (archADiv F g xF).support ⊆ (Finset.univ : Finset (InfinitePlace F)) :=
    fun v _ => Finset.mem_univ v
  have hle : ∀ v ∈ (Finset.univ : Finset (InfinitePlace F)),
      (archADiv F g xF) v ≤ C * (InfinitePlace.mult v : ℝ) := by
    intro v _
    rw [archADiv_apply]
    have := mul_le_mul_of_nonneg_left (hg (archPoint xF v))
      (Nat.cast_nonneg (InfinitePlace.mult v) : (0:ℝ) ≤ _)
    linarith
  have heq : ∑ v ∈ (archADiv F g xF).support, (archADiv F g xF) v
      = ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v :=
    Finset.sum_subset hsub (fun v _ hv => Finsupp.notMem_support_iff.1 hv)
  have hmult : ∑ v : InfinitePlace F, ((InfinitePlace.mult v : ℝ))
      = (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) (InfinitePlace.sum_mult_eq (K := F))
  rw [heq]
  calc ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v
      ≤ ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (C * (InfinitePlace.mult v : ℝ)) :=
        Finset.sum_le_sum hle
    _ = C * (Module.finrank ℚ F : ℝ) := by rw [← Finset.mul_sum, hmult]

/-- ★**正規化するとアルキメデス側の上界も `F` に依らない**。 -/
theorem archADiv_sum_div_finrank_le {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (C : ℝ) (hg : ∀ p, g p ≤ C) :
    (archADiv F g xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) ≤ C := by
  have hpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [div_le_iff₀ hpos]
  exact archADiv_sum_le F g xF C hg

/-! ## ★高さの分解 -/

/-- ★**高さ = 有限素点の寄与 + 正規化したアルキメデス側の寄与**。 -/
theorem htArith_eq_add {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    htArith F D xF
      = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
        + (archADiv F D.green xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) := by
  rw [htArith, degNormalized, degNormalized, deg_pullbackADiv]
  ring

/-! ## ★★★★因子が同じなら高さの差は一様に有界 -/

variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-- ★★★★**[GenEll] Proposition 1.4, (iii) のアルキメデス側**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

因子が同じで計量がどちらも連続なら、`|ht_D − ht_E|` は**一様に**有界である
——定数は `F` にも点にも依らない。

★★機構は `Proposition 1.6` と同じ 3 段(コンパクト性・連続関数の有界性・正規化)。 -/
theorem htArith_sub_abs_le (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D E : ArithCartier X) (hdiv : D.divisor = E.divisor)
    (hD : @Continuous _ _ M.topology _ D.green) (hE : @Continuous _ _ M.topology _ E.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        |htArith F D xF - htArith F E xF| ≤ C := by
  obtain ⟨C₁, hC₁, hlo₁, hhi₁⟩ := M.exists_bound D.green hD
  obtain ⟨C₂, hC₂, hlo₂, hhi₂⟩ := M.exists_bound E.green hE
  refine ⟨C₁ + C₂, by linarith, ?_⟩
  intro F _ _ xF
  have h1 := htArith_eq_add F D xF
  have h2 := htArith_eq_add F E xF
  rw [hdiv] at h1
  have hd1lo : -C₁ ≤ (archADiv F D.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) := archADiv_sum_div_finrank_ge F D.green xF C₁ hC₁ hlo₁
  have hd1hi : (archADiv F D.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) ≤ C₁ := archADiv_sum_div_finrank_le F D.green xF C₁ hhi₁
  have hd2lo : -C₂ ≤ (archADiv F E.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) := archADiv_sum_div_finrank_ge F E.green xF C₂ hC₂ hlo₂
  have hd2hi : (archADiv F E.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) ≤ C₂ := archADiv_sum_div_finrank_le F E.green xF C₂ hhi₂
  rw [abs_le]
  constructor
  · rw [h1, h2]; linarith
  · rw [h1, h2]; linarith

/-! ## ★出典の紐付け(`.src`) -/

def archADiv_sum_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(アルキメデス側の寄与の上界)",
    sectionId := "genell-prop-1-4" }

def htArith_sub_abs_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(計量を取り替えても高さは BD-同値——射影モデルを与えられたものとして)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
