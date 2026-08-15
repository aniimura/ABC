import ABC3.Found.FrdI.MonoidVocabulary
import ABC3.Found.FrdI.CategoryVocabulary
import Mathlib.Data.PNat.Basic
import Mathlib.Algebra.Category.MonCat.Basic

/-!
# [FrdI] Definition 1.1, (ii)(iii)(iv) —— elementary Frobenioid と pre-Frobenioid

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.19–p.20(**400 dpi 目視確認 2026-08-15**)。

## ★何が既にあって、何が無かったか

前回 `MonoidVocabulary.lean` で `Definition 1.1, (i)` を閉じた。
それは **モノイドの側**(`pre-divisorial` / `divisorial` / `group-like`)である。
`Definition 1.2` が依拠するのは **(ii)(iii)(iv)、すなわち圏の側**であり、
そこはまだ何も無かった。このファイルがそこを作る。

## この節の構成

* `(iii)` **`𝔽_M`** = モノイド `M` に付随する elementary Frobenioid(モノイド。台は `M × ℕ≥1`)
* `(iii)` **`𝔽 = 𝔽_{ℤ≥0}`** = standard Frobenioid
* `(ii)` **`monoid on 𝒟`** = 反変関手 `Φ : 𝒟 → 𝔐𝔬𝔫` で2条件を満たすもの
* `(iii)` **`𝔽_Φ`** = `Φ` に付随する elementary Frobenioid(圏)
* `(iv)` **pre-Frobenioid structure** `𝒞 → 𝔽_Φ`

## ★mathlib の実測(2026-08-15)

`𝔽_M` の積 `(a₁,n₁)·(a₂,n₂) = (a₁ + n₁·a₂, n₁·n₂)` は**半直積の形**をしているが、
mathlib の `SemidirectProduct (φ : G →* MulAut N)` は

* `G`、`N` が**群**であることを要求し
* 作用が**自己同型** (`MulAut N`) であることを要求する

のに対し、ここでの `ℕ≥1` の `M` への作用 `a ↦ n · a` は
**自己同型ではない**(`M = ℕ` で `2 · (−)` は全射でない)。★したがって
mathlib の `SemidirectProduct` は使えない。自前で書く理由はこれである。

使うもの:

* `ℕ≥1` = mathlib の **`ℕ+`**(`PNat`、`CancelCommMonoid ℕ+`)
* `𝔐𝔬𝔫`(可換モノイドの圏) = mathlib の **`AddCommMonCat`**
  (原文 p.11: 「the category of commutative monoids」、演算は加法で書く)
* 反変関手 = `Dᵒᵖ ⥤ AddCommMonCat`
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

/-! ### ★`𝔽_M` —— モノイドに付随する elementary Frobenioid

原文 (FrdI p.20):
> set of FM is the product

原文 (FrdI p.20):
> equipped with the monoid structure is given as follows: if a1, a2 ∈M, n1, n2 ∈N≥1,

原文 (FrdI p.20):
> then (a1, n1) · (a2, n2) = (a1 + n1 · a2, n1 · n2). Also, we shall write F def = FZ≥0 and
-/

/-- **[FrdI] Definition 1.1, (iii)** —— `𝔽_M`。台は `M × ℕ≥1`。

`div` は零因子 `Div(−)`、`deg` は Frobenius 次数 `deg_Fr(−)` に対応する。 -/
@[ext]
structure ElemFrob (M : Type u) [AddCommMonoid M] where
  /-- 零因子(原文の `Div(−)`) -/
  div : M
  /-- Frobenius 次数(原文の `deg_Fr(−)`) -/
  deg : ℕ+

namespace ElemFrob

variable {M : Type u} [AddCommMonoid M]

instance : Mul (ElemFrob M) where
  mul x y := ⟨x.div + (x.deg : ℕ) • y.div, x.deg * y.deg⟩

instance : One (ElemFrob M) where
  one := ⟨0, 1⟩

@[simp] theorem mul_div (x y : ElemFrob M) :
    (x * y).div = x.div + (x.deg : ℕ) • y.div := rfl
@[simp] theorem mul_deg (x y : ElemFrob M) : (x * y).deg = x.deg * y.deg := rfl
@[simp] theorem one_div : (1 : ElemFrob M).div = 0 := rfl
@[simp] theorem one_deg : (1 : ElemFrob M).deg = 1 := rfl

instance : Monoid (ElemFrob M) where
  mul_assoc x y z := by
    ext
    · simp [smul_add, mul_smul, add_assoc]
    · simp [mul_assoc]
  one_mul x := by ext <;> simp
  mul_one x := by ext <;> simp

/-- **Frobenius 次数は乗法的** —— 原文 Remark 1.1.1 の
`deg_Fr(φ ∘ ψ) = deg_Fr(φ) · deg_Fr(ψ)`。

原文 (FrdI p.21):
> Indeed, this follows immediately from the definition of an elementary Frobenioid in
-/
def degHom : ElemFrob M →* ℕ+ where
  toFun := deg
  map_one' := rfl
  map_mul' _ _ := rfl

/-! ### 非退化 -/

/-- ★`𝔽_M` は**可換でない** —— `M = ℕ` での反例。

`(0,2)·(1,1) = (2,2)` だが `(1,1)·(0,2) = (1,2)`。

★これが無いと「`M × ℕ≥1` に何か積を入れた」だけで済んでしまう。
原文の積は**半直積であって直積ではない**。 -/
theorem not_commutative_elemFrob :
    ¬ ∀ x y : ElemFrob ℕ, x * y = y * x := by
  intro h
  have := congrArg div (h ⟨0, 2⟩ ⟨1, 1⟩)
  simp at this

/-- 可逆元の Frobenius 次数は `1` —— 単位群は「次数 1 の層」に閉じ込められる。

★「`𝔽_M` の単位群が全体ではない」ことの具体形。 -/
theorem deg_eq_one_of_isUnit {x : ElemFrob M} (hx : IsUnit x) : x.deg = 1 := by
  obtain ⟨u, rfl⟩ := hx
  have h1 : ((u : ElemFrob M) * (↑u⁻¹ : ElemFrob M)).deg = 1 := by
    rw [u.mul_inv]; rfl
  rw [mul_deg] at h1
  have h2 : ((u : ElemFrob M).deg : ℕ) * ((↑u⁻¹ : ElemFrob M).deg : ℕ) = 1 := by
    exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) h1
  exact PNat.coe_eq_one_iff.mp (Nat.eq_one_of_mul_eq_one_right h2)

/-- ★`𝔽_ℕ` には**可逆でない元がある**(例: `(0, 2)`)。 -/
theorem not_isUnit_deg_two : ¬ IsUnit (⟨0, 2⟩ : ElemFrob ℕ) := by
  intro h
  have := deg_eq_one_of_isUnit h
  exact absurd this (by decide)

/-- **[FrdI] Definition 1.1, (iii)** —— **standard Frobenioid** `𝔽 = 𝔽_{ℤ≥0}`。

原文 (FrdI p.20):
> then (a1, n1) · (a2, n2) = (a1 + n1 · a2, n1 · n2). Also, we shall write F def = FZ≥0 and

原文 (FrdI p.20):
> refer to F as the standard Frobenioid.

★`ℤ_{≥0}` は `ℕ` である。 -/
abbrev Standard : Type := ElemFrob ℕ

/-- standard Frobenioid も非可換(上の反例そのもの)。 -/
theorem standard_not_commutative : ¬ ∀ x y : Standard, x * y = y * x :=
  not_commutative_elemFrob

end ElemFrob

/-! ### ★`monoid on 𝒟` —— [FrdI] Definition 1.1, (ii)

原文 (FrdI p.19):
> as a monoid on D if the following conditions are satisfied: (a) every morphism of

原文 (FrdI p.19):
> monoids α∗ : Φ(A) →Φ(B) induced by a morphism α : B →A of D is char-

原文 (FrdI p.19):
> acteristically injective [cf. §0]; (b) if α is an FSM-morphism [cf. §0] of D, then

原文 (FrdI p.19):
> α∗ : Φ(A) →Φ(B) is an isomorphism of monoids.
-/

variable (D : Type u) [Category.{v} D]

/-- **[FrdI] Definition 1.1, (ii)** —— `monoid on 𝒟`。

`𝒟` 上の**反変**関手 `Φ : 𝒟 → 𝔐𝔬𝔫`(= `Dᵒᵖ ⥤ AddCommMonCat`)であって、
原文の (a)(b) を満たすもの。

★(a) は `MonoidVocabulary.lean` の `IsCharacteristicallyInjective`、
(b) は `CategoryVocabulary.lean` の `IsFSMMorphism` を使う——
**§0 の2語がここで初めて合流する。** -/
structure MonoidOn where
  /-- 反変関手 `Φ : 𝒟 → 𝔐𝔬𝔫` -/
  functor : Dᵒᵖ ⥤ AddCommMonCat.{w}
  /-- **(a)** `α*` は characteristically injective -/
  charInj : ∀ {A B : D} (α : B ⟶ A),
    IsCharacteristicallyInjective (functor.map α.op).hom
  /-- **(b)** `α` が FSM-morphism なら `α*` は同型 -/
  fsmIso : ∀ {A B : D} (α : B ⟶ A), IsFSMMorphism α →
    Function.Bijective (functor.map α.op).hom

variable {D}

namespace MonoidOn

variable (Φ : MonoidOn.{v, u, w} D)

/-- `Φ(A)` の台。 -/
abbrev val (A : D) : Type w := (Φ.functor.obj (op A) : AddCommMonCat.{w})

/-- **`α*`** —— 原文の誘導射 `Φ(A) → Φ(B)`(`α : B ⟶ A` に対し**逆向き**)。 -/
def map {A B : D} (α : B ⟶ A) : Φ.val A →+ Φ.val B := (Φ.functor.map α.op).hom

@[simp] theorem map_id (A : D) (x : Φ.val A) : Φ.map (𝟙 A) x = x := by
  simp [map]

@[simp] theorem map_comp {A B E : D} (α : B ⟶ A) (β : E ⟶ B) (x : Φ.val A) :
    Φ.map (β ≫ α) x = Φ.map β (Φ.map α x) := by
  simp [map, op_comp, Φ.functor.map_comp]

/-- **(a)** は特に `α*` が単射であることを含む。 -/
theorem map_injective {A B : D} (α : B ⟶ A) : Function.Injective (Φ.map α) :=
  (Φ.charInj α).1

/-- `Φ` が **divisorial** であること。

原文 (FrdI p.19):
> Φ(A) [as A ranges over the objects of D] (respectively, some monoid Φ(A) [where
-/
def IsDivisorialOn : Prop := ∀ A : D, IsDivisorial (Φ.val A)

end MonoidOn

/-! ### ★`𝔽_Φ` —— [FrdI] Definition 1.1, (iii) の一般形

原文 (FrdI p.20):
> defined as follows: The objects of FΦ are the objects of D. If A, B ∈Ob(FΦ), whose

原文 (FrdI p.20):
> of FΦ is defined to be a collection of data

原文 (FrdI p.20):
> where φD : AD →BD is a morphism of D; Zφ ∈Φ(AD); nφ ∈N≥1. Here,

原文 (FrdI p.20):
> is given as follows:
-/

/-- `𝔽_Φ` の対象 = `𝒟` の対象。★型としては `𝒟` そのものだが、
圏の構造が違うので新しい型として包む。 -/
structure ElemFrobCat (Φ : MonoidOn.{v, u, w} D) where
  /-- `𝒟` の対象 -/
  base : D

namespace ElemFrobCat

variable {Φ : MonoidOn.{v, u, w} D}

/-- `𝔽_Φ` の射 `A → B` = 3つ組 `(φ_𝒟, Z_φ, n_φ)`。

原文 (FrdI p.20):
> where φD : AD →BD is a morphism of D; Zφ ∈Φ(AD); nφ ∈N≥1. Here,
-/
@[ext]
structure Hom (A B : ElemFrobCat Φ) where
  /-- `Base(φ)` -/
  base : A.base ⟶ B.base
  /-- `Div(φ) ∈ Φ(A_𝒟)` -/
  div : Φ.val A.base
  /-- `deg_Fr(φ) ∈ ℕ≥1` -/
  deg : ℕ+

/-- 合成。

原文 (FrdI p.20):
> ψ ◦φ = (ψD ◦φD, φ∗ D(Zψ) + nψ · Zφ, nψ · nφ) : A →C
-/
def Hom.comp {A B E : ElemFrobCat Φ} (φ : Hom A B) (ψ : Hom B E) : Hom A E where
  base := φ.base ≫ ψ.base
  div := Φ.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div
  deg := ψ.deg * φ.deg

instance : Category (ElemFrobCat Φ) where
  Hom := Hom
  id A := ⟨𝟙 A.base, 0, 1⟩
  comp := Hom.comp
  id_comp φ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.id_comp _
    · show Φ.map (𝟙 _) φ.div + (φ.deg : ℕ) • (0 : Φ.val _) = φ.div
      simp
    · show φ.deg * 1 = φ.deg
      simp
  comp_id φ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.comp_id _
    · show Φ.map φ.base (0 : Φ.val _) + ((1 : ℕ+) : ℕ) • φ.div = φ.div
      simp
    · show (1 : ℕ+) * φ.deg = φ.deg
      simp
  assoc φ ψ χ := by
    refine Hom.ext ?_ ?_ ?_
    · exact Category.assoc _ _ _
    · show Φ.map (φ.base ≫ ψ.base) χ.div
          + (χ.deg : ℕ) • (Φ.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div)
        = Φ.map φ.base (Φ.map ψ.base χ.div + (χ.deg : ℕ) • ψ.div)
          + ((χ.deg * ψ.deg : ℕ+) : ℕ) • φ.div
      rw [Φ.map_comp, map_add, map_nsmul, smul_add, ← add_assoc, ← mul_smul]
      norm_cast
    · simp [Hom.comp, mul_assoc]

@[simp] theorem comp_base {A B E : ElemFrobCat Φ} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).base = φ.base ≫ ψ.base := rfl
@[simp] theorem comp_deg {A B E : ElemFrobCat Φ} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).deg = ψ.deg * φ.deg := rfl
@[simp] theorem id_base (A : ElemFrobCat Φ) : Hom.base (𝟙 A) = 𝟙 A.base := rfl
@[simp] theorem id_deg (A : ElemFrobCat Φ) : Hom.deg (𝟙 A) = 1 := rfl
@[simp] theorem id_div (A : ElemFrobCat Φ) : Hom.div (𝟙 A) = 0 := rfl

/-- **自然な射影関手 `𝔽_Φ → 𝒟`**。

原文 (FrdI p.20):
> functor:
-/
def proj : ElemFrobCat Φ ⥤ D where
  obj A := A.base
  map φ := φ.base

/-! ### ★Remark 1.1.1 の3式

原文 (FrdI p.21):
> Base(φ ◦ψ) = Base(φ) ◦Base(ψ)

原文 (FrdI p.21):
> Div(φ ◦ψ) = (Base(ψ))∗(Div(φ)) + degFr(φ) · Div(ψ)

原文 (FrdI p.21):
> degFr(φ ◦ψ) = degFr(φ) · degFr(ψ)

原文 (FrdI p.21):
> Indeed, this follows immediately from the definition of an elementary Frobenioid in

★原文は `φ ◦ ψ`(先に `ψ`)で書く。Lean の `ψ ≫ φ` が同じものである。
★「follows immediately」と書かれた3式を**実際に示す**。
-/

/-- `Base(φ ∘ ψ) = Base(φ) ∘ Base(ψ)`。 -/
theorem base_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).base = ψ.base ≫ φ.base := rfl

/-- `Div(φ ∘ ψ) = (Base(ψ))*(Div(φ)) + deg_Fr(φ) · Div(ψ)`。 -/
theorem div_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).div = Φ.map ψ.base φ.div + (φ.deg : ℕ) • ψ.div := rfl

/-- `deg_Fr(φ ∘ ψ) = deg_Fr(φ) · deg_Fr(ψ)`。 -/
theorem degFr_comp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (ψ ≫ φ).deg = φ.deg * ψ.deg := rfl

/-- ★**divisor monoid が sharp なら、`𝔽_Φ` の同型は零因子 0 を持つ**。

`φ ≫ φ⁻¹ = 𝟙` の零因子成分は
`Φ(Base φ)(Div φ⁻¹) + deg_Fr(φ⁻¹) · Div φ = 0`。
sharp なモノイドでは `x + y = 0` から `y` が可逆、よって `y = 0`。
さらに `deg_Fr(φ⁻¹) ≥ 1` なので、`MonoidVocabulary.lean` の
`isTorsionFreeNaive_of_isSharp` により `Div φ = 0`。

★これは `Definition 1.2` の `isometric` と `isomorphism` を結ぶ含意であり、
`Definition 1.2` の非退化 witness(`Witness.lean`)で使う。 -/
theorem div_eq_zero_of_isIso (hsharp : ∀ A : D, IsSharp (Φ.val A))
    {A B : ElemFrobCat Φ} (φ : A ⟶ B) [IsIso φ] : Hom.div φ = 0 := by
  have hd := congrArg Hom.div (IsIso.hom_inv_id φ)
  rw [div_comp, id_div] at hd
  have hunit : IsAddUnit (((Hom.deg (inv φ) : ℕ+) : ℕ) • Hom.div φ) :=
    ⟨⟨_, Φ.map (Hom.base φ) (Hom.div (inv φ)), by rw [add_comm]; exact hd, hd⟩, rfl⟩
  have hz := hsharp A.base _ hunit
  exact isTorsionFreeNaive_of_isSharp (hsharp A.base) _ _ (Hom.deg (inv φ)).pos hz

end ElemFrobCat

/-- ★**`𝔽_Φ` は totally epimorphic である**。

原文 (FrdI p.27):
> Proof. Since D is a connected, totally epimorphic category, the fact that FΦ is

原文 (FrdI p.27):
> tion 1.1, (iii); the fact that a pre-divisorial monoid is integral [cf. Definition 1.1,

★これは `Proposition 1.5` の証明の**第1段**である。原文が挙げる入力は3つ:

1. **`𝒟` が totally epimorphic**(`hD`)
2. **`Φ(A)` が integral**(`hint`)—— pre-divisorial なら integral
3. **`Definition 1.1, (ii), (a)` の単射性** —— `MonoidOn.charInj` として構造に入っている

★**3つがそのまま3成分に対応する**: `base` は 1 で、`div` は 2 と 3 で、`deg` は `ℕ+` の簡約で消える。

★**2026-08-15 訂正の記録**: 一度この補題を
「`𝒟` の hom が subsingleton」+「`Φ(A)` が `IsCancelAdd`」という**原文より強い仮定**で書いていた。
模型(`Vee` は前順序)では通るが**一般には使えない**形で、`Proposition 1.5` を阻む。
★「模型で通る」と「一般に使える」は違う、という実例。 -/
theorem isTotallyEpimorphic_elemFrobCat {Φ : MonoidOn.{v, u, w} D}
    (hD : IsTotallyEpimorphic D) (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) :
    IsTotallyEpimorphic (ElemFrobCat Φ) := by
  intro X Y f
  refine ⟨fun {Z} g h hgh => ?_⟩
  -- `deg`: `ℕ+` は簡約的
  have hdeg : ElemFrobCat.Hom.deg g = ElemFrobCat.Hom.deg h := by
    have hh := congrArg ElemFrobCat.Hom.deg hgh
    simp only [ElemFrobCat.comp_deg] at hh
    exact mul_right_cancel hh
  -- `base`: ★入力1 —— `𝒟` が totally epimorphic なので `Base f` は epi
  have hbase : ElemFrobCat.Hom.base g = ElemFrobCat.Hom.base h := by
    haveI : Epi (ElemFrobCat.Hom.base f) := hD _ _ _
    have hh := congrArg ElemFrobCat.Hom.base hgh
    simp only [ElemFrobCat.comp_base] at hh
    exact (cancel_epi (ElemFrobCat.Hom.base f)).mp hh
  -- `div`: ★入力2 と 3
  have hdiv : ElemFrobCat.Hom.div g = ElemFrobCat.Hom.div h := by
    have hh := congrArg ElemFrobCat.Hom.div hgh
    simp only [ElemFrobCat.div_comp, hdeg] at hh
    letI := isCancelAdd_of_isIntegralMonoid _ (hint X.base)
    exact Φ.map_injective (ElemFrobCat.Hom.base f) (add_right_cancel hh)
  exact ElemFrobCat.Hom.ext hbase hdiv hdeg

/-! ### ★`pre-Frobenioid structure` —— [FrdI] Definition 1.1, (iv)

原文 (FrdI p.20):
> (iv) Let D, Φ, FΦ be as in (iii); C a category. Assume further that Φ is divisorial,

原文 (FrdI p.20):
> and that C, D are connected, totally epimorphic categories [cf. §0]. Then we shall

原文 (FrdI p.20):
> as a pre-Frobenioid structure on C. The natural projection functor FΦ →D thus
-/

/-- **[FrdI] Definition 1.1, (iv)** —— `𝒞` 上の pre-Frobenioid structure。

★原文が課す仮定を**すべて**フィールドとして持つ:
`Φ` が divisorial、`𝒞` と `𝒟` が totally epimorphic、そして関手 `𝒞 → 𝔽_Φ`。

★`connected`(**圏**について)は §0 が**対象**についてしか定めていないため、
ここには入れていない。原文の曖昧さであり、測定結果として報告に残す。 -/
structure PreFrobenioid (C : Type u) [Category.{v} C] (Φ : MonoidOn.{v, u, w} D) where
  /-- 原文の `𝒞 → 𝔽_Φ` -/
  toElem : C ⥤ ElemFrobCat Φ
  /-- `Φ` は divisorial -/
  divisorial : Φ.IsDivisorialOn
  /-- `𝒞` は totally epimorphic -/
  totEpiC : IsTotallyEpimorphic C
  /-- `𝒟` は totally epimorphic -/
  totEpiD : IsTotallyEpimorphic D

namespace PreFrobenioid

variable {C : Type u} [Category.{v} C] {Φ : MonoidOn.{v, u, w} D}

/-- **自然な射影関手 `𝒞 → 𝒟`**。

原文 (FrdI p.20):
> restricts to a natural projection functor
-/
def proj (P : PreFrobenioid C Φ) : C ⥤ D := P.toElem ⋙ ElemFrobCat.proj

/-- **`Base(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev Base (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) :
    (P.toElem.obj A).base ⟶ (P.toElem.obj B).base := (P.toElem.map φ).base

/-- **`Div(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev Div (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) : Φ.val (P.toElem.obj A).base :=
  (P.toElem.map φ).div

/-- **`deg_Fr(−)`** —— `𝒞` の射に制限した操作。 -/
abbrev degFr (P : PreFrobenioid C Φ) {A B : C} (φ : A ⟶ B) : ℕ+ := (P.toElem.map φ).deg

@[simp] theorem Base_id (P : PreFrobenioid C Φ) (A : C) : P.Base (𝟙 A) = 𝟙 _ := by
  simp [PreFrobenioid.Base]

@[simp] theorem degFr_id (P : PreFrobenioid C Φ) (A : C) : P.degFr (𝟙 A) = 1 := by
  simp [PreFrobenioid.degFr]

@[simp] theorem Div_id (P : PreFrobenioid C Φ) (A : C) : P.Div (𝟙 A) = 0 := by
  simp [PreFrobenioid.Div]

/-- ★**Remark 1.1.1 の3式は `𝒞` の側でも成り立つ**(関手性から)。 -/
theorem degFr_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.degFr (ψ ≫ φ) = P.degFr φ * P.degFr ψ := by
  simp [degFr, P.toElem.map_comp]

/-- `Base(−)` の合成則(`𝒞` の側)。 -/
theorem Base_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.Base (ψ ≫ φ) = P.Base ψ ≫ P.Base φ := by
  simp [PreFrobenioid.Base, P.toElem.map_comp]

/-- `Div(−)` の合成則(`𝒞` の側)。 -/
theorem Div_comp (P : PreFrobenioid C Φ) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    P.Div (ψ ≫ φ) = Φ.map (P.Base ψ) (P.Div φ) + (P.degFr φ : ℕ) • P.Div ψ := by
  simp [PreFrobenioid.Div, PreFrobenioid.Base, PreFrobenioid.degFr, P.toElem.map_comp,
    ElemFrobCat.div_comp]

end PreFrobenioid

/-! ### ★出典の紐付け(`.src`) -/

def ElemFrob.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iii) — 𝔽_M = M × ℕ≥1",
    sectionId := "frdi-def-1-1-iii" }

def MonoidOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (ii) — monoid on 𝒟",
    sectionId := "frdi-def-1-1-ii" }

def ElemFrobCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iii) — 𝔽_Φ",
    sectionId := "frdi-def-1-1-iii" }

def PreFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iv) — pre-Frobenioid structure",
    sectionId := "frdi-def-1-1-iv" }

end ABC3.Found.FrdI
