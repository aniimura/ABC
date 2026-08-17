import ABC3.Found.GenEll.ArchBound

/-!
# [GenEll] Example 1.3, (ii) —— **compactly bounded は posit しなくてよい**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

## ★★★原文の定義は `archPoint` があれば書ける

原文は `compactly bounded subset` をこう定める:

> for the subset of points `x ∈ X(F) ⊆ X(Q̄)`, where `[F : Q] < ∞`, such that for each
> `v ∈ V^arc`(resp. `v ∈ V^non`), the set of `[F : Q]` points of `X^arc`
> (resp. `X(Q_v)`) determined by `x` is contained in `K_v`.

★★**「`x` が定める `X^arc` の `[F:ℚ]` 個の点」は `archPoint`(`ArchPoint.lean`、本日実装)
そのものである。** 無限素点 `v` ごとに `𝓞_F → ℂ` を取って `x_F` に合成したものであり、
その個数は `#{無限素点} ≤ [F:ℚ]` である(`ArchBound.lean`)。

★★★**したがって `Interface` の `CompactlyBounded` を posit する必要は無い**——
アルキメデス側は**定義できる**。
★`Example 1.3, (i)`(Galois-finite)が既に定義で書けていた
(`Skeleton/GenEll/Section1.lean`)のと同じ事情である。

## ★★これが何のためにあるか

`Theorem 2.1`(Effective Mordell/ABC/Vojta)は
**compactly bounded subset の上で**述べられる。
★その理由は、アルキメデス側の寄与を**一様に**抑えるためである。

★★★**本ファイルはその「一様に抑えられる」ことを実際に証明する**
(`archADiv_sum_le_of_boundedBy`)。

## ★posit を外せない部分は明示する

非アルキメデス側(`X(ℚ_v)` の compact domain)は、`ℚ_v` 上の点と
その位相を要する。★本ファイルは**アルキメデス側だけ**を扱う。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★`x` が定める `X^arc` の点の集合 -/

/-- ★★**原文「the set of `[F : Q]` points of `X^arc` determined by `x`」**。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

★無限素点ごとの ℂ-点の像である。 -/
def archPointSet {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    Set (complexPoints X) :=
  Set.range (archPoint xF)

/-- ★**`x` が `K` に囲われている**(原文「is contained in `K_v`」のアルキメデス側)。

★★`Example 1.3, (ii)` の定義本体である——**posit ではない**。 -/
def BoundedByArch {X : Scheme.{0}} (K : Set (complexPoints X))
    (xF : specRingOfIntegers F ⟶ X) : Prop :=
  archPointSet F xF ⊆ K

theorem boundedByArch_iff {X : Scheme.{0}} (K : Set (complexPoints X))
    (xF : specRingOfIntegers F ⟶ X) :
    BoundedByArch F K xF ↔ ∀ v : InfinitePlace F, archPoint xF v ∈ K := by
  constructor
  · intro h v; exact h ⟨v, rfl⟩
  · rintro h _ ⟨v, rfl⟩; exact h v

/-- ★**全体は誰でも囲う**(非空虚性の witness)。 -/
theorem boundedByArch_univ {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    BoundedByArch F (Set.univ : Set (complexPoints X)) xF :=
  fun _ _ => trivial

/-! ## ★★★囲われた点ではアルキメデス側の寄与が一様に抑えられる -/

/-- ★★★**囲われた点の上では、アルキメデス側の寄与は `[C, ...]` で一様に抑えられる**。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

★★★**定数は `F` に依らない。** 正規化(`[F:ℚ]` で割る)がそれを与える
——`#{無限素点} ≤ [F:ℚ]` だからである(`ArchBound.lean`)。

★★これが `Theorem 2.1` が compactly bounded subset の上で述べられる理由である。 -/
theorem archADiv_sum_le_of_boundedBy {X : Scheme.{0}} (g : GreenFn X)
    (K : Set (complexPoints X)) (xF : specRingOfIntegers F ⟶ X)
    (hb : BoundedByArch F K xF) (C : ℝ) (hC : 0 ≤ C) (hg : ∀ p ∈ K, g p ≤ C) :
    (archADiv F g xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) ≤ C := by
  classical
  have hpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [div_le_iff₀ hpos, Finsupp.sum]
  have hsub : (archADiv F g xF).support ⊆ (Finset.univ : Finset (InfinitePlace F)) :=
    fun v _ => Finset.mem_univ v
  have heq : ∑ v ∈ (archADiv F g xF).support, (archADiv F g xF) v
      = ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v :=
    Finset.sum_subset hsub (fun v _ hv => Finsupp.notMem_support_iff.1 hv)
  rw [heq]
  -- ★`mult` の重みが入ったので、`Σ_v mult v = [F:ℚ]` がそのまま効く
  have hle : ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)), (archADiv F g xF) v
      ≤ ∑ v ∈ (Finset.univ : Finset (InfinitePlace F)),
          C * (InfinitePlace.mult v : ℝ) := by
    refine Finset.sum_le_sum fun v _ => ?_
    rw [archADiv_apply]
    have := mul_le_mul_of_nonneg_left (hg _ (hb ⟨v, rfl⟩))
      (Nat.cast_nonneg (InfinitePlace.mult v) : (0:ℝ) ≤ _)
    linarith
  have hmult : ∑ v : InfinitePlace F, ((InfinitePlace.mult v : ℝ))
      = (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast congrArg (Nat.cast : ℕ → ℝ) (InfinitePlace.sum_mult_eq (K := F))
  rw [← Finset.mul_sum, hmult] at hle
  linarith

/-! ## ★仮定が discharge できること -/

/-- ★**コンパクト集合上の連続関数は上に有界**。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

★★これが原文の「compact domain」の効き所である。
★★★`X^arc` の位相を `ℙⁿ(ℂ)` から引くところは
`ProjClosed.lean` にあるが、`complexPoints X` との接続はまだである。**混同しない。** -/
theorem exists_upper_bound_on_compact {T : Type*} [TopologicalSpace T]
    {K : Set T} (hK : IsCompact K) {g : T → ℝ} (hg : Continuous g) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p ∈ K, g p ≤ C := by
  rcases K.eq_empty_or_nonempty with rfl | hne
  · exact ⟨0, le_refl 0, fun p hp => absurd hp (Set.notMem_empty p)⟩
  obtain ⟨x₀, hx₀, hmax⟩ := hK.exists_isMaxOn hne hg.continuousOn
  refine ⟨|g x₀|, abs_nonneg _, fun p hp => ?_⟩
  exact le_trans (hmax hp) (le_abs_self _)

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `Example 1.3, (ii)` は非アルキメデス側
(`X(ℚ_v)` の compact domain)も含む。 -/

def BoundedByArch.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Example 1.3, (ii)(アルキメデス側のみ——非アルキメデス側は未着手)",
    sectionId := "genell-ex-1-3" }

def archADiv_sum_le_of_boundedBy.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Example 1.3, (ii)(囲われた点でアルキメデス側の寄与が一様に抑えられること)",
    sectionId := "genell-ex-1-3" }

end ABC3.Found.GenEll
