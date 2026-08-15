import Mathlib.GroupTheory.MonoidLocalization.Basic
import Mathlib.GroupTheory.Subgroup.Saturated
import Mathlib.Data.ENat.Basic
import Mathlib.Algebra.Order.Group.Nat
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Algebra.Group.Submonoid.Membership
import Mathlib.Data.ZMod.Basic

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
> the natural homomorphism from M to its groupification M gp. Thus, M gp is the
-/
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
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;
-/
def IsSharp : Prop := ∀ a : M, IsAddUnit a → a = 0

/-- **[FrdI] §0 `integral`** —— `M → M^gp` が単射。

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;
-/
def IsIntegralMonoid : Prop := Function.Injective (toGp M)

/-- **[FrdI] §0 `saturated`** —— `n · a` が `M` の像に入るような `a ∈ M^gp`
(`n ≥ 1`)は、それ自身 `M` の像に入る。

原文 (FrdI p.11):
> we shall say that M is saturated if every a ∈M gp for which n · a lies in the image

原文 (FrdI p.11):
> of M for some n ∈N≥1 lies in the image of M.
-/
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

/-! ### ★独立性 —— どの1語も他の2語から従わない -/

/-- `M → M^gp` が全射なら saturated(自明に)。 -/
theorem isSaturatedMonoid_of_range_eq_univ (h : Set.range (toGp M) = Set.univ) :
    IsSaturatedMonoid M := fun a _ _ _ => h ▸ Set.mem_univ a

/-- `ℤ` では `M → M^gp` は全射。 -/
theorem range_toGp_int : Set.range (toGp ℤ) = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  induction x using AddLocalization.ind with
  | _ p =>
    obtain ⟨a, b⟩ := p
    refine ⟨a - (b : ℤ), ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simp⟩

/-- `ℕ∞` では `M → M^gp` は全射(`⊤` がすべてを潰すので `M^gp` は1点)。 -/
theorem range_toGp_enat : Set.range (toGp ℕ∞) = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  induction x using AddLocalization.ind with
  | _ p =>
    obtain ⟨a, b⟩ := p
    refine ⟨0, ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    refine ⟨⟨⊤, trivial⟩, ?_⟩
    simp

/-- ★**独立性 (1)**: `ℤ` は **integral かつ saturated だが sharp でない**。 -/
theorem independence_sharp :
    IsIntegralMonoid ℤ ∧ IsSaturatedMonoid ℤ ∧ ¬ IsSharp ℤ :=
  ⟨isIntegralMonoid_of_isCancelAdd ℤ,
   isSaturatedMonoid_of_range_eq_univ ℤ range_toGp_int,
   not_isSharp_int⟩

/-- `ℕ∞` は sharp —— 正準順序つき加法モノイドでは `a + b = 0 → a = 0`。 -/
theorem isSharp_enat : IsSharp ℕ∞ := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h := u.val_neg
  exact (add_eq_zero.mp h).1

/-- ★**独立性 (2)**: `ℕ∞` は **sharp かつ saturated だが integral でない**。 -/
theorem independence_integral :
    IsSharp ℕ∞ ∧ IsSaturatedMonoid ℕ∞ ∧ ¬ IsIntegralMonoid ℕ∞ :=
  ⟨isSharp_enat,
   isSaturatedMonoid_of_range_eq_univ ℕ∞ range_toGp_enat,
   not_isIntegralMonoid_enat⟩

/-! ### 独立性 (3) —— `saturated` が他の2語から従わないこと -/

/-- `⟨2,3⟩ ⊆ ℕ` が生成する部分モノイド。`1` を含まないことが効く。 -/
abbrev TwoThree : AddSubmonoid ℕ := AddSubmonoid.closure ({2, 3} : Set ℕ)

theorem two_mem_twoThree : (2 : ℕ) ∈ TwoThree :=
  AddSubmonoid.subset_closure (by simp)

theorem three_mem_twoThree : (3 : ℕ) ∈ TwoThree :=
  AddSubmonoid.subset_closure (by simp)

theorem one_not_mem_twoThree : (1 : ℕ) ∉ TwoThree := by
  intro h
  rw [AddSubmonoid.mem_closure_pair] at h
  obtain ⟨m, n, hmn⟩ := h
  simp only [smul_eq_mul] at hmn
  omega

/-- ★**独立性 (3)**: `⟨2,3⟩` は **sharp かつ integral だが saturated でない**。

`a := mk 3 2`(= `1`)は `2 • a = mk 6 4`(= `2`)が像に入るのに、
それ自身は像に入らない —— `1 ∉ ⟨2,3⟩` だから。 -/
theorem independence_saturated :
    IsSharp (TwoThree) ∧ IsIntegralMonoid (TwoThree) ∧ ¬ IsSaturatedMonoid (TwoThree) := by
  refine ⟨?_, isIntegralMonoid_of_isCancelAdd _, ?_⟩
  · intro a ha
    obtain ⟨u, rfl⟩ := ha
    have h := u.val_neg
    have hsum : ((u : TwoThree) : ℕ) + ((u.neg : TwoThree) : ℕ) = 0 := congrArg Subtype.val h
    refine Subtype.ext ?_
    show ((u : TwoThree) : ℕ) = 0
    omega
  · intro hsat
    set a : Gp TwoThree :=
      AddLocalization.mk ⟨3, three_mem_twoThree⟩ ⟨⟨2, two_mem_twoThree⟩, trivial⟩ with ha
    have hmem : (2 : ℕ) • a ∈ Set.range (toGp TwoThree) := by
      refine ⟨⟨2, two_mem_twoThree⟩, ?_⟩
      rw [ha, AddLocalization.mk_nsmul, toGp, AddLocalization.mk_eq_mk_iff,
        AddLocalization.r_iff_exists]
      refine ⟨0, ?_⟩
      ext
      simp
    obtain ⟨m, hm⟩ := hsat a 2 (by omega) hmem
    rw [ha, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at hm
    obtain ⟨c, hc⟩ := hm
    have : ((m : TwoThree) : ℕ) + 2 = 3 := by
      have := congrArg Subtype.val hc
      simp at this
      omega
    refine one_not_mem_twoThree ?_
    have hm1 : ((m : TwoThree) : ℕ) = 1 := by omega
    simp [← hm1]

/-! ### `torsion-free` と、原文が主張する含意

原文 (FrdI p.11):
> inductive system I∗determines a natural morphism

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if

★**「which」の正体**: 直前の `M → M^pf`(**perfection** への自然な射)である。
`M^gp` ではない。原文は `M^pf` を「乗法 `n·` による帰納極限」として定める:

原文 (FrdI p.11):
> the perfection of M, that is to say, the inductive limit of the inductive system I∗of

★★**原文は「torsion element」を定義していない**(2026-08-15 実測: [FrdI] 全体で
`torsion` は 3 箇所——上の 2 箇所と、無関係な `torsion subgroup` の言及 1 箇所のみ)。
したがって `torsion-free` の意味は**原文からは決まらない**。標準的な読みを採り、
その旨をここに明記する。

★mathlib の `IsAddTorsionFree` は **「`n • ·` が単射」**であり、素朴な読み
(「`n • a = 0` なら `a = 0`」)とは**別の条件**である。群では一致するが、
モノイドでは一般に別物なので、両方を置いて関係を証明する。 -/

/-- **[FrdI] §0 `torsion-free`**(素朴な読み) —— torsion 元を持たない。

★原文はこの語を定義していない。ここは標準的な読み
「`n ≥ 1` かつ `n • a = 0` なら `a = 0`」を採る。 -/
def IsTorsionFreeNaive : Prop := ∀ (a : M) (n : ℕ), 0 < n → n • a = 0 → a = 0

variable {M}

/-- ★★**原文の「hence, in particular」の中身** —— **sharp なら torsion-free**。

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if

原文 (FrdI p.11):
> M is sharp, integral, and saturated. We shall say that M is perfect if multiplication

★証明は1行で、**`sharp` だけで足りる**(`integral` も `saturated` も要らない):
`n • a = 0` かつ `n ≥ 1` なら `a + (n-1) • a = 0` なので `a` は可逆、
`sharp` より `a = 0`。 -/
theorem isTorsionFreeNaive_of_isSharp (h : IsSharp M) : IsTorsionFreeNaive M := by
  intro a n hn hna
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  refine h a ⟨⟨a, m • a, ?_, ?_⟩, rfl⟩
  · calc a + m • a = m • a + a := add_comm _ _
      _ = (m + 1) • a := (succ_nsmul a m).symm
      _ = 0 := hna
  · calc m • a + a = (m + 1) • a := (succ_nsmul a m).symm
      _ = 0 := hna

/-- mathlib の `IsAddTorsionFree`(`n • ·` が単射)なら、素朴な読みも成り立つ。

★逆は一般には成り立たない。**同じ語で別の条件**である。 -/
theorem isTorsionFreeNaive_of_isAddTorsionFree [IsAddTorsionFree M] :
    IsTorsionFreeNaive M := by
  intro a n hn hna
  have h0 : n • a = n • (0 : M) := by simpa using hna
  exact IsAddTorsionFree.nsmul_right_injective (by omega) h0

/-! ### `M^±` と `of characteristic type`

原文 (FrdI p.11):
> the submonoid [which, in fact, forms a group] of invertible elements of M, by

原文 (FrdI p.11):
> the quotient monoid of M by M ±, which we shall refer to as the characteristic of

原文 (FrdI p.11):
> characteristic type if the fibers of the natural map M →M char are torsors over M ±.

★**`M^char` そのもの(商モノイド)は作っていない。** 作るにはモノイドの合同関係
(`AddCon`)が要る。ここで要るのは「**同じファイバーに属する**」という関係だけなので、
それを直接書く——`M^± ` で生成される合同関係は
「`∃ u v ∈ M^±, a + u = b + v`」である。**この選択は原文の言い換えであって、
原文が書いていないことを足してはいない。** -/

variable (M)

/-- **`M^±`** —— 可逆元のなす部分モノイド。 -/
def unitsSubmonoid : AddSubmonoid M where
  carrier := {a : M | IsAddUnit a}
  zero_mem' := isAddUnit_zero
  add_mem' := IsAddUnit.add

variable {M}

/-- `M → M^char` の**同じファイバー**に属する関係
(= `M^±` が生成する合同関係)。 -/
def CharRel (a b : M) : Prop := ∃ u v : M, IsAddUnit u ∧ IsAddUnit v ∧ a + u = b + v

theorem charRel_refl (a : M) : CharRel a a := ⟨0, 0, isAddUnit_zero, isAddUnit_zero, rfl⟩

variable (M)

/-- **[FrdI] §0 `of characteristic type`** —— `M → M^char` のファイバーが
`M^±` 上の**トーサー**(単純推移的)であること。

原文 (FrdI p.11):
> characteristic type if the fibers of the natural map M →M char are torsors over M ±.

★「トーサー」を「推移的かつ自由」に開いて書く。 -/
def IsOfCharacteristicType : Prop :=
  ∀ a b : M, CharRel a b → ∃! u : M, IsAddUnit u ∧ b = a + u

/-! ### 非退化 -/

/-- `ℕ` は **of characteristic type** —— 可逆元が `0` だけなのでファイバーは1点。 -/
theorem isOfCharacteristicType_nat : IsOfCharacteristicType ℕ := by
  intro a b hab
  obtain ⟨u, v, hu, hv, huv⟩ := hab
  obtain ⟨w, rfl⟩ := hu
  obtain ⟨w', rfl⟩ := hv
  have hw : ((w : ℕ)) = 0 := by have := w.val_neg; omega
  have hw' : ((w' : ℕ)) = 0 := by have := w'.val_neg; omega
  refine ⟨0, ⟨isAddUnit_zero, by omega⟩, ?_⟩
  rintro y ⟨hy, hy2⟩
  obtain ⟨z, rfl⟩ := hy
  have := z.val_neg
  omega

/-- ★`WithTop (ZMod 2)` は **of characteristic type でない** ——
`⊤` のファイバー上で `M^±` の作用が**自由でない**(`⊤ + 0 = ⊤ + 1 = ⊤`)。 -/
theorem not_isOfCharacteristicType_withTop_zmod2 :
    ¬ IsOfCharacteristicType (WithTop (ZMod 2)) := by
  intro h
  have h1 : IsAddUnit ((1 : ZMod 2) : WithTop (ZMod 2)) :=
    ⟨⟨((1 : ZMod 2) : WithTop (ZMod 2)), ((1 : ZMod 2) : WithTop (ZMod 2)), by decide, by decide⟩, rfl⟩
  obtain ⟨u, hu, huniq⟩ := h ⊤ ⊤ (charRel_refl ⊤)
  have e0 : (0 : WithTop (ZMod 2)) = u := huniq 0 ⟨isAddUnit_zero, by simp⟩
  have e1 : ((1 : ZMod 2) : WithTop (ZMod 2)) = u := huniq _ ⟨h1, by simp⟩
  rw [← e0] at e1
  exact absurd e1 (by decide)

/-! ### ★[FrdI] Definition 1.1, (i) —— `pre-divisorial` / `divisorial`

原典: [FrdI] 物理 p.19(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.19):
> (i) We shall say that M ∈Ob(Mon) is pre-divisorial if it is integral [cf. §0],

原文 (FrdI p.19):
> saturated [cf. §0], and of characteristic type [cf. §0]. Suppose that M is pre-

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is

★**Definition 1.1 (i) が §0 から引く語は `integral` / `saturated` /
`of characteristic type` の3つで、いま全部揃った。**

★ただし同じ (i) が定める `group-like`(「`M^char` が零」)は**書けていない**——
`M^char` を商モノイドとして作っていないため。 -/

/-- **[FrdI] Definition 1.1, (i)** —— `pre-divisorial`。

原文 (FrdI p.19):
> (i) We shall say that M ∈Ob(Mon) is pre-divisorial if it is integral [cf. §0],
-/
def IsPreDivisorial : Prop :=
  IsIntegralMonoid M ∧ IsSaturatedMonoid M ∧ IsOfCharacteristicType M

/-- **[FrdI] Definition 1.1, (i)** —— `divisorial` = pre-divisorial かつ sharp。

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is
-/
def IsDivisorial : Prop := IsPreDivisorial M ∧ IsSharp M

variable {M}

/-- 群は **of characteristic type** —— `M^± = M` なのでファイバーは全体、
`u = b - a` が一意。 -/
theorem isOfCharacteristicType_int : IsOfCharacteristicType ℤ := by
  intro a b _
  refine ⟨b - a, ⟨⟨⟨b - a, a - b, by omega, by omega⟩, rfl⟩, by omega⟩, ?_⟩
  rintro y ⟨-, hy⟩
  omega

/-- ★`ℤ` は **pre-divisorial**(Definition 1.1 (i) を満たす実例)。 -/
theorem isPreDivisorial_int : IsPreDivisorial ℤ :=
  ⟨isIntegralMonoid_of_isCancelAdd ℤ,
   isSaturatedMonoid_of_range_eq_univ ℤ range_toGp_int,
   isOfCharacteristicType_int⟩

/-- ★`ℤ` は **divisorial ではない** —— sharp でないから。
`pre-divisorial` と `divisorial` が別物であることの確認。 -/
theorem not_isDivisorial_int : ¬ IsDivisorial ℤ := fun h => not_isSharp_int h.2

end ABC3.Found.FrdI
