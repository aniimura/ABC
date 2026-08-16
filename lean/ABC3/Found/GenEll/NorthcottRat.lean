import ABC3.Meta.Claim
import Mathlib.NumberTheory.Height.NumberField
import Mathlib.NumberTheory.Height.Northcott

/-!
# [GenEll] Proposition 1.4, (iv) の基底 —— `ℚ` 上の Northcott 性(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★これは mathlib が **TODO** と書いている場所である

`Mathlib/NumberTheory/Height/Northcott.lean` の冒頭:

> We provide instances showing that `K` also satisfies the Northcott property
> * for `logHeight₁`,
> * **(TODO)** for `Projectivization.mulHeight`,
> * **(TODO)** for `Projectivization.logHeight`.

★そして **`Northcott (mulHeight₁ (K := K))` の基底インスタンスが 1 つも無い**
(2026-08-17 実測——mathlib 全体を `instance.*Northcott` で検索して、
あるのは `Ideal.cardQuot` のものと、上の条件つき変換のみ)。

★★**ゆえに `Proposition 1.4, (iv)`(= Northcott の有限性)は、
「mathlib に高さがある」だけでは届かない。**
`genell-goal.md` §3-2 は「高さ | ★**ある**。5 ファイル 2006 行」と記録していたが、
**要る定理そのものは無い**——これは記述を精密化する測定である。

## ★本ファイルが取るもの

> **`ℚ` 上の `mulHeight₁` は Northcott 性を満たす。**

★`mathlib` の `mulHeight₁_eq_max : mulHeight₁ q = max |q.num| q.den` から、
高さが `B` 以下の有理数は分子・分母が有界、ゆえに有限個である。

★★**これは基底の場合である。** `Proposition 1.4, (iv)` が要求するのは
`X(ℚ̄)^{≤d}`(**次数有界**の代数点)上の有限性で、そちらは数え上げが要る。
**基底が取れたことと、全体が取れたことを混同しない。**

## ★★追記(2026-08-17 夜)—— 任意の数体へ一般化された

`Found/GenEll/DenominatorBound.lean` の `northcott_mulHeight₁` が
**任意の数体 `K`** について同じ instance を与えたので、
★**本ファイルの `northcott_mulHeight₁_rat` はその特別な場合になった**。

★本ファイルを残すのは、`Rat.mulHeight₁_eq_max`(高さ = `max |num| den`)による
**明示的な** 評価 `num_den_le_of_mulHeight₁_le` が別途有用だからである。
一般の場合の証明は分母イデアルを経由しており、`ℚ` のこの明示形は出てこない。
-/

namespace ABC3.Found.GenEll

open Height

/-- ★高さ `B` 以下の有理数は、分子と分母がともに `⌊B⌋₊` で抑えられる。 -/
theorem num_den_le_of_mulHeight₁_le {q : ℚ} {B : ℝ} (h : mulHeight₁ q ≤ B) :
    q.num.natAbs ≤ ⌊B⌋₊ ∧ q.den ≤ ⌊B⌋₊ := by
  have hq : ((max q.num.natAbs q.den : ℕ) : ℝ) ≤ B := by
    rw [← Rat.mulHeight₁_eq_max q]
    exact h
  have hB : (0 : ℝ) ≤ B := le_trans (by positivity) hq
  have hmax : max q.num.natAbs q.den ≤ ⌊B⌋₊ := Nat.le_floor hq
  exact ⟨le_trans (le_max_left _ _) hmax, le_trans (le_max_right _ _) hmax⟩

/-- ★★**`ℚ` 上の Northcott 性**。

★mathlib は `Northcott (mulHeight₁ (K := K))` を**仮定として要求する**だけで、
基底インスタンスを持たない(2026-08-17 実測)。ここで `ℚ` の場合を与える。 -/
instance northcott_mulHeight₁_rat : Northcott (mulHeight₁ (K := ℚ)) where
  finite_le B := by
    classical
    -- 分子・分母の対を取る写像で有限集合の像に埋める
    set N : ℕ := ⌊B⌋₊ with hN
    have hsub : {q : ℚ | mulHeight₁ q ≤ B}
        ⊆ (fun p : ℤ × ℕ => (p.1 : ℚ) / (p.2 : ℚ)) ''
            (↑((Finset.Icc (-(N : ℤ)) (N : ℤ)) ×ˢ (Finset.range (N + 1)))) := by
      intro q hq
      obtain ⟨hnum, hden⟩ := num_den_le_of_mulHeight₁_le hq
      refine ⟨(q.num, q.den), ?_, ?_⟩
      · simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_Icc,
          Finset.mem_range]
        refine ⟨⟨?_, ?_⟩, by omega⟩
        · have : q.num.natAbs ≤ N := hnum
          omega
        · have : q.num.natAbs ≤ N := hnum
          omega
      · exact_mod_cast Rat.num_div_den q
    exact Set.Finite.subset (Set.Finite.image _ (Finset.finite_toSet _)) hsub

/-! ## ★出典の紐付け(`.src`) -/

def northcott_mulHeight₁_rat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6, item := "Proposition 1.4, (iv)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
