import Mathlib.GroupTheory.MonoidLocalization.Basic
import Mathlib.GroupTheory.Subgroup.Saturated
import Mathlib.Data.ENat.Basic
import Mathlib.Algebra.Order.Group.Nat
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Algebra.Group.Submonoid.Membership

/-!
# [FrdI] §0 のモノイド語彙 —— `sharp` / `integral` / `saturated`

原典: S. Mochizuki, *The Geometry of Frobenioids I: The General Theory* [FrdI]、
物理 p.11(全 126 ページ。**400 dpi 目視確認 2026-08-15**)。

★**原文は加法で書く**:

原文 (FrdI p.11):
> of M will be written additively. We shall denote by

原文 (FrdI p.11):
> the submonoid [which, in fact, forms a group] of invertible elements of M, by

原文 (FrdI p.11):
> the natural homomorphism from M to its groupification M gp. Thus, M gp is the

そして3語の定義は1文にまとまっている:

原文 (FrdI p.11):
> that M is torsion-free if M has no torsion elements; we shall say that M is sharp if

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;

原文 (FrdI p.11):
> we shall say that M is saturated if every a ∈M gp for which n · a lies in the image

原文 (FrdI p.11):
> of M for some n ∈N≥1 lies in the image of M.

★`M ±` は原文では **`M^±`**(上付き)であり、`pdftotext` が上付きを空白区切りの別トークンに
するために `M ±` と出る。`papers.json` の `FrdI.notationNotes` に実測を記録した。

## ★なぜここを作るか(層0の仕分けの結論)

2026-08-15 の層0(129節点)の仕分けで、**いま着手できるもの**として6件を名指しした。
その筆頭がこの3語である。理由は2つ:

1. **mathlib に土台がある。** `IsAffineMonoid = IsCancelMul + Monoid.FG + IsMulTorsionFree`
   (`Mathlib/Algebra/AffineMonoid/Basic.lean`)。
2. ★**これは Frobenioid の入口である。** `Frobenioid` は
   `[IUTchIII] Definition 3.8` の解消条件(中層待ち2件のうち1件)であり、
   [FrdI] は機械概算で最も密(到達45項目・辺184本・平均出次数 4.1)。
   **層0の最も安い部分が、中層の最大の壁の入口になっている。**

## ★mathlib の何を使ったか(車輪の再発明をしない)

- `M^gp`(groupification): mathlib の **`AddLocalization (⊤ : AddSubmonoid M)`**。
  ★原文の同値関係「`a₁ + b₂ + c = b₁ + a₂ + c` for some `c ∈ M`」は
  `AddLocalization.r` と**同じ形**である(原文 p.11 目視)。
  mathlib には `GrothendieckAddGroup M` という `abbrev` もあるが、
  2026-08-15 時点の本ビルドでは `import` しても識別子が見えなかった(原因未特定)。
  `abbrev` の右辺そのものを使っているので、実質は同じものである。
- `saturated`: **`AddSubmonoid.NSMulSaturated`**(`Mathlib/GroupTheory/Subgroup/Saturated.lean`)。
  下の `frdI_saturated_iff` で、原文の定義がこれと一致することを示す。
- `sharp` に相当する名前は mathlib に無い(2026-08-15 実測: `\bsharp\b` は 17 件あるが
  **全部「sharp bound」**で別物)。1行なので自前で書く。
-/

namespace ABC3.Found.FrdI

open AddLocalization

variable (M : Type*) [AddCommMonoid M]

/-! ### `M^gp` —— groupification -/

/-- **`M^gp`** —— 原文の groupification。mathlib の `AddLocalization` そのもの。

原文 (FrdI p.11):
> the natural homomorphism from M to its groupification M gp. Thus, M gp is the -/
abbrev Gp : Type _ := AddLocalization (⊤ : AddSubmonoid M)

/-- 正準写像 `M → M^gp`。 -/
def toGp (m : M) : Gp M := AddLocalization.mk m 0

variable {M}

/-- `toGp` が等しいことの判定 —— 原文の同値関係そのもの
(`a₁ + b₂ + c = b₁ + a₂ + c` for some `c ∈ M`)。 -/
theorem toGp_eq_iff {m m' : M} : toGp M m = toGp M m' ↔ ∃ c : M, m + c = m' + c := by
  constructor
  · intro h
    rw [toGp, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at h
    obtain ⟨c, hc⟩ := h
    refine ⟨(c : M), ?_⟩
    simpa [add_comm] using hc
  · rintro ⟨c, hc⟩
    rw [toGp, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    refine ⟨⟨c, trivial⟩, ?_⟩
    simpa [add_comm] using hc

variable (M)

/-! ### 3語の定義 -/

/-- **[FrdI] §0 `sharp`** —— `M^± = 0`、すなわち可逆元が `0` だけ。

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective; -/
def IsSharp : Prop := ∀ a : M, IsAddUnit a → a = 0

/-- **[FrdI] §0 `integral`** —— `M → M^gp` が単射。

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective; -/
def IsIntegralMonoid : Prop := Function.Injective (toGp M)

/-- **[FrdI] §0 `saturated`** —— `n · a` が `M` の像に入るような `a ∈ M^gp`
(`n ≥ 1`)は、それ自身 `M` の像に入る。

原文 (FrdI p.11):
> we shall say that M is saturated if every a ∈M gp for which n · a lies in the image

原文 (FrdI p.11):
> of M for some n ∈N≥1 lies in the image of M. -/
def IsSaturatedMonoid : Prop :=
  ∀ (a : Gp M) (n : ℕ), 0 < n → n • a ∈ Set.range (toGp M) → a ∈ Set.range (toGp M)

/-! ### mathlib との接続 -/

/-- ★**原文の `saturated` の条件は、mathlib の `AddSubmonoid.NSMulSaturated` そのものである。**

原文は「`n ≥ 1` かつ `n·a` が入るなら `a` も入る」と述べ、mathlib は
「`n • g ∈ H` なら `n = 0` または `g ∈ H`」と述べる。**同じ主張**である。

★これを任意の部分モノイドについて示しておくのは、`M^gp` の中の `M` の像に
そのまま当てられるようにするためである。 -/
theorem frdI_saturated_iff {G : Type*} [AddMonoid G] (H : AddSubmonoid G) :
    (∀ (a : G) (n : ℕ), 0 < n → n • a ∈ H → a ∈ H) ↔ H.NSMulSaturated := by
  constructor
  · intro h n g hng
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · exact Or.inl rfl
    · exact Or.inr (h g n hn hng)
  · intro h a n hn hmem
    rcases h hmem with h0 | h1
    · omega
    · exact h1

/-- **cancellative なら integral**(mathlib の局所化の単射性が要求するのと同じ条件)。 -/
theorem isIntegralMonoid_of_isCancelAdd [IsCancelAdd M] : IsIntegralMonoid M := by
  intro m m' h
  obtain ⟨c, hc⟩ := toGp_eq_iff.mp h
  exact add_right_cancel hc

/-! ### ★非退化 —— 3語それぞれに「満たす例」と「満たさない例」

★これが無いと「型は付くが中身が空」を排除できない。 -/

/-- `ℕ` は **sharp** —— 可逆元は `0` だけ。 -/
theorem isSharp_nat : IsSharp ℕ := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h := u.val_neg
  omega

/-- ★`ℤ` は **sharp でない** —— `1` が可逆で `1 ≠ 0`。 -/
theorem not_isSharp_int : ¬ IsSharp ℤ := by
  intro h
  have h1 : IsAddUnit (1 : ℤ) := ⟨⟨1, -1, by omega, by omega⟩, rfl⟩
  have h2 := h 1 h1
  omega

/-- `ℕ` は **integral** —— `ℕ` は cancellative。 -/
theorem isIntegralMonoid_nat : IsIntegralMonoid ℕ :=
  isIntegralMonoid_of_isCancelAdd ℕ

/-- ★`ℕ∞` は **integral でない** —— `1 + ⊤ = 2 + ⊤` なので `toGp 1 = toGp 2`、
しかし `1 ≠ 2`。 -/
theorem not_isIntegralMonoid_enat : ¬ IsIntegralMonoid ℕ∞ := by
  intro h
  have heq : toGp ℕ∞ 1 = toGp ℕ∞ 2 := toGp_eq_iff.mpr ⟨⊤, by simp⟩
  have := h heq
  simp at this

/-- ★`ℕ` の中の**偶数**のなす部分モノイドは **saturated でない** ——
`2 • 1 = 2` は入るが `1` は入らない。

`frdI_saturated_iff` により、これは原文の saturated の条件が破れている実例である。 -/
theorem not_saturated_evens :
    ¬ (∀ (a : ℕ) (n : ℕ), 0 < n → n • a ∈ AddSubmonoid.closure ({2} : Set ℕ) →
        a ∈ AddSubmonoid.closure ({2} : Set ℕ)) := by
  intro h
  have h2 : (2 : ℕ) ∈ AddSubmonoid.closure ({2} : Set ℕ) :=
    AddSubmonoid.subset_closure rfl
  have h1 : (1 : ℕ) ∈ AddSubmonoid.closure ({2} : Set ℕ) := by
    refine h 1 2 (by omega) ?_
    simpa using h2
  rw [AddSubmonoid.mem_closure_singleton] at h1
  obtain ⟨n, hn⟩ := h1
  simp only [smul_eq_mul] at hn
  omega

/-- `⊤` は **saturated**(自明に)。「全部落とす」定義ではないことの対照。 -/
theorem saturated_top :
    ∀ (a : ℕ) (n : ℕ), 0 < n → n • a ∈ (⊤ : AddSubmonoid ℕ) → a ∈ (⊤ : AddSubmonoid ℕ) :=
  fun _ _ _ _ => trivial

end ABC3.Found.FrdI
