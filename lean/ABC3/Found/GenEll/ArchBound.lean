import ABC3.Found.GenEll.ConductorHeight

/-!
# [GenEll] Proposition 1.6 のアルキメデス側 —— **正規化が定数を一様にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

## ★★★原文の証明はアルキメデス側をこう言う

> for the contributions at the archimedean primes, from the fact that the
> continuous function `|s|_L` on the compact topological space `X^arc` is bounded.

★つまり「有界だから `≲` に吸収される」である。
★★**しかし `≲` の定数が `F` に依ってはいけない**——
`log-cond_D ≲ ht_D` は `X(ℚ̄)` 全体、すなわち**すべての数体**の上での主張だからである。

## ★★★正規化がその一様性を与える

`ht` は `deg_F` を `[F : ℚ]` で割ったもの(`degNormalized`)である。
アルキメデス側の寄与は無限素点の個数だけ足されるが、

    `#{無限素点} = r₁ + r₂ ≤ r₁ + 2r₂ = [F : ℚ]`

なので、`[F : ℚ]` で割れば **`|g| ≤ C` から `|寄与| ≤ C` が出る**。
★★★**定数は `F` に依らない。**

★★これが `Definition 1.1` が「正規化した次数」を使う理由の 1 つである
——原文は理由を書いていないが、`Proposition 1.6` の証明が
それを**要求している**。

## ★仮定が discharge できることも示す

「連続関数は下に有界」は `IsCompact.exists_isMinOn` で出る。
★★ただし `X^arc = X(ℂ)` に位相を入れる所は**まだ繋がっていない**——
`ℙⁿ(ℂ)` のコンパクト性は `ProjTopology.lean` で取ったが、
`X ⊆ ℙⁿ` の埋め込みを経由する段が残っている。
★したがって本ファイルの下界補題は**任意のコンパクト空間について**述べる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★無限素点の個数は次数以下 -/

/-- ★**`#{無限素点} ≤ [F : ℚ]`**。

`r₁ + r₂ ≤ r₁ + 2r₂ = [F : ℚ]` である。
★★**これが正規化の効き所**——アルキメデス側の寄与を `[F:ℚ]` で割ると
`F` に依らない定数で抑えられる。 -/
theorem card_infinitePlace_le_finrank :
    Fintype.card (InfinitePlace F) ≤ Module.finrank ℚ F := by
  rw [InfinitePlace.card_eq_nrRealPlaces_add_nrComplexPlaces,
    ← InfinitePlace.card_add_two_mul_card_eq_rank]
  omega

/-! ## ★★アルキメデス側の寄与の下界 -/

/-- ★**Green 関数が `-C` 以上なら、寄与は `-C · #{無限素点}` 以上**。 -/
theorem archADiv_sum_ge {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (C : ℝ) (hg : ∀ p, -C ≤ g p) :
    -(C * (Fintype.card (InfinitePlace F) : ℝ))
      ≤ (archADiv F g xF).sum (fun _ r => r) := by
  classical
  rw [Finsupp.sum]
  have hsub : (archADiv F g xF).support ⊆ (Finset.univ : Finset (InfinitePlace F)) :=
    fun v _ => Finset.mem_univ v
  have hle : ∀ v ∈ (Finset.univ : Finset (InfinitePlace F)),
      -C ≤ (archADiv F g xF) v := by
    intro v _
    rw [archADiv_apply]
    exact hg _
  -- 台の外では値は `0` なので、`univ` 上の和に置き換えてよい
  have heq : ∑ v ∈ (archADiv F g xF).support, (archADiv F g xF) v
      = ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v :=
    Finset.sum_subset hsub (fun v _ hv => Finsupp.notMem_support_iff.1 hv)
  rw [heq]
  calc -(C * (Fintype.card (InfinitePlace F) : ℝ))
      = ∑ _v ∈ (Finset.univ : Finset (InfinitePlace F)), (-C) := by
        rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]; ring
    _ ≤ ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v :=
        Finset.sum_le_sum hle

/-- ★★★**正規化するとアルキメデス側の下界は `F` に依らない**。

★★`#{無限素点} ≤ [F : ℚ]` なので、`[F : ℚ]` で割れば `-C` で抑えられる。
★★★**これが `Definition 1.1` が正規化した次数を使う理由の 1 つ**である。 -/
theorem archADiv_sum_div_finrank_ge {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (C : ℝ) (hC : 0 ≤ C) (hg : ∀ p, -C ≤ g p) :
    -C ≤ (archADiv F g xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) := by
  have hpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [le_div_iff₀ hpos]
  refine le_trans ?_ (archADiv_sum_ge F g xF C hg)
  have hcard : (Fintype.card (InfinitePlace F) : ℝ) ≤ (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast card_infinitePlace_le_finrank F
  nlinarith

/-! ## ★★★`Proposition 1.6` —— 原文の `≲` の形 -/

/-- ★★★**導手は高さで抑えられる(定数つき)** —— 原文の `≲` の形。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

★★★**定数 `C` は `F` に依らない。** これが正規化の効き所である。

★`HeightNonneg.lean` の `logCond_le_htArith` は `g ≥ 0` を仮定していたが、
本定理は **`g ≥ -C`(下に有界)**でよい——それが原文の仮定である。 -/
theorem logCond_le_htArith_add {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (C : ℝ) (hC : 0 ≤ C)
    (h : pullbackIdeal F D.divisor xF ≠ 0) (hg : ∀ p, -C ≤ D.green p) :
    logCond F D.divisor xF ≤ htArith F D xF + C := by
  have h1 : logCond F D.divisor xF
      ≤ degNormalized (idealADiv F (pullbackIdeal F D.divisor xF)) :=
    logCond_le_degNormalized_pullback F D.divisor xF h
  have h2 : -C ≤ (archADiv F D.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) :=
    archADiv_sum_div_finrank_ge F D.green xF C hC hg
  have h3 : htArith F D xF
      = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
        + (archADiv F D.green xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) := by
    rw [htArith, degNormalized, degNormalized, deg_pullbackADiv]
    ring
  rw [h3]
  linarith

/-- ★★**BD-class の水準** —— 原文の `log-cond_D ≲ ht_D` そのもの。

★★★**定数は `F` に依らない**ので、`X(ℚ̄)` 全体(= すべての数体)の上で成り立つ。 -/
theorem logCond_bdge_htArith_of_bddBelow {X : Scheme.{0}} (D : ArithCartier X)
    (C : ℝ) (hC : 0 ≤ C) (hg : ∀ p, -C ≤ D.green p)
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0) :
    BDge (fun xF : specRingOfIntegers F ⟶ X => logCond F D.divisor xF)
      (fun xF => htArith F D xF) :=
  ⟨C, fun xF => by
    show logCond F D.divisor xF - htArith F D xF ≤ C
    linarith [logCond_le_htArith_add F D xF C hC (h xF) hg]⟩

/-! ## ★仮定が discharge できること -/

/-- ★**コンパクト空間上の連続関数は下に有界**。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

★★これが原文の「`|s|_L` は有界」に当たる。
★★★ただし `X^arc = X(ℂ)` に位相を入れる段は**まだ繋がっていない**——
`ℙⁿ(ℂ)` のコンパクト性は `ProjTopology.lean` にあるが、
`X ⊆ ℙⁿ` の埋め込みを経由する段が残っている。**混同しない。** -/
theorem exists_lower_bound_of_continuous {T : Type*} [TopologicalSpace T]
    [CompactSpace T] [Nonempty T] (g : T → ℝ) (hg : Continuous g) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p, -C ≤ g p := by
  obtain ⟨x₀, -, hmin⟩ :=
    isCompact_univ.exists_isMinOn Set.univ_nonempty hg.continuousOn
  refine ⟨|g x₀|, abs_nonneg _, fun p => ?_⟩
  have h1 : g x₀ ≤ g p := hmin (Set.mem_univ p)
  have h2 : -|g x₀| ≤ g x₀ := neg_abs_le _
  linarith

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `Proposition 1.6` は `L = O_X(D)` との対応を含み、
また `X^arc` の位相との接続が残っている。 -/

def logCond_bdge_htArith_of_bddBelow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(Green 関数が下に有界な場合——X^arc の位相との接続は含まない)",
    sectionId := "genell-prop-1-6" }

def card_infinitePlace_le_finrank.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(正規化が定数を一様にすることのみ)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
