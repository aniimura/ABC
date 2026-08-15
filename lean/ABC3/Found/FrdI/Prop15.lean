import ABC3.Found.FrdI.Prop14

/-!
# [FrdI] Proposition 1.5 —— Elementary Frobenioids are Frobenioids

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.27(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.27):
> Proposition 1.5. (Elementary Frobenioids are Frobenioids) Let Φ be a

原文 (FrdI p.27):
> pre-divisorial monoid on a connected, totally epimorphic category D. Then:

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-

## ★★`wP` との差 —— 何が模型固有だったか

`WitnessFrobenioid.lean` の `wIsFrobenioid` は `Definition 1.3` の22条を
`wP`(`𝒟 = Vee`、`Φ = ℕ`、`𝒞 = 𝔽_Φ`、`toElem = 𝟭`)について埋めた。
★**一般化にあたって、そのうち何が模型固有だったかを先に洗い出す。**

**模型固有だったもの(一般には使えない)**:

1. ★★**`wHom_ext`** —— 「`Vee` の hom は subsingleton なので射は `(wd, wn)` で決まる」。
   `wIsFrobenioid` の**ほぼ全条**がこれを使っている。一般の `𝒟` では
   底の射の等式を**別途示す**必要がある。
2. ★**`toElem = 𝟭`** —— `wP` は恒等関手だったので `Div φ = φ.div` だった。
   ★**原文の構造は `𝔽_Φ → 𝔽_{Φ^char}`** なので `Div φ = toChar φ.div` であり、
   **isometric は「零因子が 0」ではなく「零因子が可逆」を意味する**
   (`toChar_eq_zero_iff`)。`wP` では `Φ = ℕ` が sharp だったので両者が一致していた。

**模型固有でなかったもの(一般に成り立つ)**:

3. ★**`𝔽_Φ` は isotropic type** —— isometric な pre-step は
   「零因子が可逆・次数 1・底が同型」なので、**逆射を明示的に作れる**。
   下の `elemFrob_isotropic` で一般に証明する。
4. ★したがって **`𝔽_Φ` のすべての射が co-angular** ——
   これは `Proposition 1.4, (i)` の直接の帰結であり、
   ★**原文自身が同じ道を通る**:

原文 (FrdI p.27):
> sition 1.4, (i) [which is applicable to all morphisms of FΦ since FΦ is of isotropic

## ★このファイルの到達点

* `elemFrobToChar` —— 関手 `𝔽_Φ → 𝔽_{Φ^char}`
* `elemPreFrobenioid` —— ★**`𝔽_Φ` が `Φ^char` 上の pre-Frobenioid であること**
* `elemFrob_isotropic` / `elemFrob_coAngular` —— 上の 3, 4

★**`Definition 1.3` の22条は、まだ埋めていない。** どこで止まっているかは
ファイル末尾の記録を見よ。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] (Φ : MonoidOn.{v, u, w} D)

/-! ### ★関手 `𝔽_Φ → 𝔽_{Φ^char}`

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-
-/

/-- ★**`𝔽_Φ → 𝔽_{Φ^char}`** —— 対象はそのまま、零因子を `Φ(A)^char` へ送る。

関手性は `charMap_toChar`(`charMap g ∘ toChar = toChar ∘ g`)から出る。 -/
def elemFrobToChar : ElemFrobCat Φ ⥤ ElemFrobCat Φ.charOn where
  obj A := ⟨A.base⟩
  map {A B} f := ⟨f.base, toChar f.div, f.deg⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show toChar (0 : Φ.val A.base) = 0
    exact map_zero _
  map_comp {A B E} f g := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show toChar (Φ.map f.base g.div + ((g.deg : ℕ+) : ℕ) • f.div)
      = charMap (Φ.map f.base) (toChar g.div) + ((g.deg : ℕ+) : ℕ) • toChar f.div
    rw [map_add, map_nsmul, charMap_toChar]

@[simp] theorem elemFrobToChar_div {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).div = toChar f.div := rfl
@[simp] theorem elemFrobToChar_deg {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).deg = f.deg := rfl
@[simp] theorem elemFrobToChar_base {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).base = f.base := rfl

/-! ### ★`𝔽_Φ` は `Φ^char` 上の pre-Frobenioid

原文 (FrdI p.27):
> (i)]; and the injectivity condition of Definition 1.1, (ii), (a). Thus, FΦ is a pre-
-/

variable (hD : IsTotallyEpimorphic D) (hpd : ∀ A : D, IsPreDivisorial (Φ.val A))

/-- ★**`Proposition 1.5, (i)` の前半** —— `𝔽_Φ` は `Φ^char` 上の pre-Frobenioid。

原文が挙げる3つの入力がそのままフィールドを埋める:

* `divisorial` ← `charOn_isDivisorialOn`(`Φ` が pre-divisorial なら `Φ^char` は divisorial)
* `totEpiC` ← `isTotallyEpimorphic_elemFrobCat`(`𝒟` が totally epimorphic + `Φ` が integral)
* `totEpiD` ← 仮定 -/
def elemPreFrobenioid : PreFrobenioid (ElemFrobCat Φ) Φ.charOn where
  toElem := elemFrobToChar Φ
  divisorial := Φ.charOn_isDivisorialOn hpd
  totEpiC := isTotallyEpimorphic_elemFrobCat hD (fun A => (hpd A).1)
  totEpiD := hD

@[simp] theorem elemPreFrobenioid_Div {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).Div f = toChar f.div := rfl
@[simp] theorem elemPreFrobenioid_degFr {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).degFr f = f.deg := rfl
@[simp] theorem elemPreFrobenioid_Base {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).Base f = f.base := rfl

/-! ### ★`𝔽_Φ` は isotropic type

原文 (FrdI p.27):
> type]. This completes the proof of assertion (i). Assertion (iii) is immediate from
-/

/-- ★★**`𝔽_Φ` のすべての対象は isotropic**(一般)。

isometric な pre-step `φ` は
「零因子 `φ.div` が**可逆**(`toChar_eq_zero_iff`)・次数 1・底が同型」なので、
**逆射を明示的に構成できる**: `⟨inv φ.base, Φ.map (inv φ.base) U.neg, 1⟩`。

★`wP` では「零因子が **0**」だったが、一般には「**可逆**」である。
`Φ = ℕ` が sharp なので両者が一致していた——**これが模型固有の性質だった**。 -/
theorem elemFrob_isotropic (A : ElemFrobCat Φ) :
    IsIsotropic (elemPreFrobenioid Φ hD hpd) A := by
  intro Dd φ hiso hstep
  have hdeg : φ.deg = 1 := hstep.1
  haveI hbase : IsIso φ.base := hstep.2
  have hu : IsAddUnit φ.div := toChar_eq_zero_iff.mp hiso
  obtain ⟨U, hU⟩ := hu
  refine ⟨⟨inv φ.base, Φ.map (inv φ.base) U.neg, 1⟩, ?_, ?_⟩
  · refine ElemFrobCat.Hom.ext (IsIso.hom_inv_id _) ?_ ?_
    · show Φ.map φ.base (Φ.map (inv φ.base) U.neg) + ((1 : ℕ+) : ℕ) • φ.div
        = (0 : Φ.val A.base)
      rw [← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
      simp only [PNat.one_coe, one_smul, ← hU]
      rw [add_comm]
      exact U.val_neg
    · show (1 : ℕ+) * φ.deg = 1
      simp [hdeg]
  · refine ElemFrobCat.Hom.ext (IsIso.inv_hom_id _) ?_ ?_
    · show Φ.map (inv φ.base) φ.div + (φ.deg : ℕ) • Φ.map (inv φ.base) U.neg
        = (0 : Φ.val Dd.base)
      rw [hdeg]
      simp only [PNat.one_coe, one_smul, ← map_add, ← hU]
      rw [U.val_neg, map_zero]
    · show φ.deg * 1 = 1
      simp [hdeg]

/-- ★**`𝔽_Φ` のすべての射は co-angular**。

★これは `Proposition 1.4, (i)` の直接の帰結であり、**原文自身が同じ道を通る**。

原文 (FrdI p.27):
> sition 1.4, (i) [which is applicable to all morphisms of FΦ since FΦ is of isotropic
-/
theorem elemFrob_coAngular {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsCoAngular (elemPreFrobenioid Φ hD hpd) φ :=
  prop_1_4_i _ φ (fun X _ => elemFrob_isotropic Φ hD hpd X)

/-- ★**`Proposition 1.4, (i)` の「In particular」を `𝔽_Φ` に当てたもの** ——
Frobenius 型 ⟺ isometric な base-isomorphism。 -/
theorem elemFrob_frobeniusType_iff {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ ↔
      (IsIsometric (elemPreFrobenioid Φ hD hpd) φ ∧
        IsBaseIsomorphism (elemPreFrobenioid Φ hD hpd) φ) :=
  prop_1_4_i_frobeniusType _ φ (fun X _ => elemFrob_isotropic Φ hD hpd X)

/-! ## ★第1段 —— `Proposition 1.5, (ii)`: `𝒪^▷(A) ≅ Φ(A)`

原文 (FrdI p.27):
> (ii) There is a natural, functorial isomorphism

★原文の証明は「It is immediate from the definitions that assertion (ii) holds」の
一言である。**本当に immediate かを測る**のがここの目的。
-/

/-- `Φ(A)` の元 `a` から作る `𝒪^▷(A)` の元 —— 底は恒等、Frobenius 次数は 1。 -/
def otriOf (A : ElemFrobCat Φ) (a : Φ.val A.base) : End A :=
  (⟨𝟙 A.base, a, 1⟩ : A ⟶ A)

@[simp] theorem otriOf_div (A : ElemFrobCat Φ) (a : Φ.val A.base) :
    ElemFrobCat.Hom.div (otriOf Φ A a) = a := rfl
@[simp] theorem otriOf_base (A : ElemFrobCat Φ) (a : Φ.val A.base) :
    ElemFrobCat.Hom.base (otriOf Φ A a) = 𝟙 A.base := rfl
@[simp] theorem otriOf_deg (A : ElemFrobCat Φ) (a : Φ.val A.base) :
    ElemFrobCat.Hom.deg (otriOf Φ A a) = 1 := rfl


theorem otriOf_mem (A : ElemFrobCat Φ) (a : Φ.val A.base) :
    otriOf Φ A a ∈ OTri (elemPreFrobenioid Φ hD hpd) A :=
  ⟨((elemPreFrobenioid Φ hD hpd).Base_id A).symm, rfl⟩

/-- ★**`𝒪^▷(A)` は `otriOf` の像そのもの** —— base-identity かつ linear な自己射は
「底が `𝟙`、次数が 1」に他ならない。 -/
theorem mem_otri_iff {A : ElemFrobCat Φ} (x : End A) :
    x ∈ OTri (elemPreFrobenioid Φ hD hpd) A ↔ x = otriOf Φ A (ElemFrobCat.Hom.div x) := by
  constructor
  · rintro ⟨hb, hl⟩
    refine ElemFrobCat.Hom.ext ?_ rfl hl
    show ElemFrobCat.Hom.base x = 𝟙 A.base
    exact hb.trans ((elemPreFrobenioid Φ hD hpd).Base_id A)
  · intro h
    have hm := otriOf_mem Φ hD hpd A (ElemFrobCat.Hom.div x)
    rwa [← h] at hm

/-- ★★**`Proposition 1.5, (ii)`** —— `Φ(A) ≃ 𝒪^▷(A)` はモノイドの同型。

`Φ(A)` は加法モノイドなので `Multiplicative` を挟む(mathlib の標準の道具)。 -/
def otriEquiv (A : ElemFrobCat Φ) :
    Multiplicative (Φ.val A.base) ≃* OTri (elemPreFrobenioid Φ hD hpd) A where
  toFun a := ⟨otriOf Φ A a.toAdd, otriOf_mem Φ hD hpd A _⟩
  invFun x := Multiplicative.ofAdd (ElemFrobCat.Hom.div (x : End A))
  left_inv _ := rfl
  right_inv x := Subtype.ext ((mem_otri_iff Φ hD hpd _).mp x.2).symm
  map_mul' a b := by
    refine Subtype.ext (ElemFrobCat.Hom.ext (by simp) ?_ (by simp))
    show a.toAdd + b.toAdd
      = Φ.map (𝟙 A.base) a.toAdd + ((1 : ℕ+) : ℕ) • b.toAdd
    rw [Φ.map_id]
    simp

/-- ★原文の `[so 𝒪^×(A) ≅ Φ(A)^±]` —— 同型のもとで `𝒪^×(A)` は
**`Φ(A)` の可逆元**にちょうど対応する。

原文 (FrdI p.27):
> [so O×(A)
-/
theorem otriOf_mem_otimes (A : ElemFrobCat Φ) (a : Φ.val A.base) :
    otriOf Φ A a ∈ OTimes (elemPreFrobenioid Φ hD hpd) A ↔ IsAddUnit a := by
  constructor
  · rintro ⟨-, hu⟩
    have hiso : IsIso (otriOf Φ A a : A ⟶ A) := (isUnit_iff_isIso _).mp hu
    exact ((ElemFrobCat.isIso_iff _).mp hiso).2.1
  · intro hu
    refine ⟨otriOf_mem Φ hD hpd A a, ?_⟩
    exact (isUnit_iff_isIso _).mpr
      ((ElemFrobCat.isIso_iff _).mpr ⟨by simp only [otriOf_base]; infer_instance, hu, rfl⟩)

/-! ### ★関手性を測る

★`𝒪^▷` の推移写像は、原文が `Definition 1.3, (iii), (c)` で与える
「`φ ≫ β = α ≫ φ`」という**圏論的な条件**で決まる。それが `Φ.map` と
一致するかを測る。 -/

include hpd in
/-- ★★**推移写像の完全な形** —— `φ ≫ β = α ≫ φ` は零因子の等式
`Φ(Base φ)(b) = deg_Fr(φ) · a` と**同値**である。

★**Frobenius 次数が入る**のが要点。`deg_Fr(φ) = 1`(linear)のときだけ
推移写像は `Φ.map (Base φ)` そのものになる。 -/
theorem otriOf_comm_iff {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (a : Φ.val A.base) (b : Φ.val B.base) :
    (φ ≫ otriOf Φ B b : A ⟶ B) = otriOf Φ A a ≫ φ ↔
      Φ.map (ElemFrobCat.Hom.base φ) b
        = ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ) • a := by
  letI := isCancelAdd_of_isIntegralMonoid _ ((hpd A.base).1)
  have hkey : (Φ.map (ElemFrobCat.Hom.base φ) b + ElemFrobCat.Hom.div φ
        = ElemFrobCat.Hom.div φ + ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ) • a) ↔
      Φ.map (ElemFrobCat.Hom.base φ) b = ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ) • a := by
    rw [add_comm (ElemFrobCat.Hom.div φ)]
    exact ⟨fun h => add_right_cancel h, fun h => by rw [h]⟩
  constructor
  · intro h
    have hd := congrArg ElemFrobCat.Hom.div h
    rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp] at hd
    simp only [otriOf_div, otriOf_deg, otriOf_base, PNat.one_coe, one_smul,
      MonoidOn.map_id] at hd
    exact hkey.mp hd
  · intro h
    refine ElemFrobCat.Hom.ext (by simp) ?_ (by simp)
    rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp]
    simp only [otriOf_div, otriOf_deg, otriOf_base, PNat.one_coe, one_smul,
      MonoidOn.map_id]
    exact hkey.mpr h

/-- ★**線型射に沿った推移写像は `Φ.map (Base φ)` そのもの**。

★これが `Proposition 1.5, (ii)` の「natural, functorial」の中身である ——
`Φ.map` は `MonoidOn.map_id` / `MonoidOn.map_comp` で関手則を満たすから、
`𝒪^▷` の推移も関手則を満たす。

★**測定結果**: 関手性が言えるのは **linear な射に沿ってだけ**である。
一般の射では `otriOf_comm_iff` のとおり Frobenius 次数 `deg_Fr(φ)` が入るので、
`(𝔽_Φ)ᵒᵖ` 全体の上の関手にはならない。原文の「functorial」はこの範囲の主張である。 -/
theorem otriOf_natural {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hφ : IsLinear (elemPreFrobenioid Φ hD hpd) φ) (b : Φ.val B.base) :
    (φ ≫ otriOf Φ B b : A ⟶ B)
      = otriOf Φ A (Φ.map (ElemFrobCat.Hom.base φ) b) ≫ φ := by
  refine (otriOf_comm_iff Φ hpd φ _ b).mpr ?_
  show _ = ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ) • _
  rw [show ElemFrobCat.Hom.deg φ = 1 from hφ]
  simp

/-! ## ★第2段 —— 射が pull-back ⟺ **linear isometry**

原文 (FrdI p.27):
> definition of the category FΦ in Definition 1.1, (iii)] that a morphism of FΦ is a

原文 (FrdI p.27):
> pull-back morphism if and only if it is a linear isometry. The fact that FΦ satisfies

★原文は「one verifies immediately」と書く。**測る。**
-/

/-- ★★**`𝔽_Φ` の pull-back morphism の完全な特徴づけ** ——
`φ` が pull-back ⟺ `deg_Fr(φ) = 1` かつ `Div(φ) = 0`(= 零因子が**可逆**)。

★逆向き(`⇐`)では逆写像を**明示的に作る**: `(g, h) ↦ ⟨h, Div(g) + Φ(h)(U⁻), deg(g)⟩`。
`U` は `φ` の零因子の加法逆元で、**可逆性がここで効く**。

★順向き(`⇒`)は目標 `((φ_base, 0, 1), 𝟙)` **1点への全射性だけ**で両方が出る ——
次数成分から `deg_Fr(φ) = 1`、零因子成分から可逆性。 -/
theorem elemFrob_isPullBack_iff {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsPullBack (elemPreFrobenioid Φ hD hpd) φ ↔
      (IsLinear (elemPreFrobenioid Φ hD hpd) φ ∧
        IsIsometric (elemPreFrobenioid Φ hD hpd) φ) := by
  constructor
  · intro h
    obtain ⟨f, hf⟩ := (h A).2
      ⟨((⟨ElemFrobCat.Hom.base φ, 0, 1⟩ : A ⟶ B), 𝟙 A.base), (Category.id_comp _).symm⟩
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f ≫ φ : A ⟶ B) = (⟨ElemFrobCat.Hom.base φ, 0, 1⟩ : A ⟶ B) :=
      congrArg Prod.fst hp
    have h2 : ElemFrobCat.Hom.base f = 𝟙 A.base := congrArg Prod.snd hp
    have hdeg : ElemFrobCat.Hom.deg φ = 1 := by
      have hd := congrArg ElemFrobCat.Hom.deg h1
      rw [ElemFrobCat.comp_deg] at hd
      exact pnat_left_eq_one hd
    refine ⟨hdeg, ?_⟩
    have hd := congrArg ElemFrobCat.Hom.div h1
    rw [ElemFrobCat.div_comp, h2, hdeg] at hd
    simp only [PNat.one_coe, one_smul, MonoidOn.map_id] at hd
    exact toChar_eq_zero_iff.mpr
      ⟨⟨ElemFrobCat.Hom.div φ, ElemFrobCat.Hom.div f, hd, by rw [add_comm]; exact hd⟩, rfl⟩
  · rintro ⟨hlin, hiso⟩
    have hdeg : ElemFrobCat.Hom.deg φ = 1 := hlin
    obtain ⟨U, hU⟩ := toChar_eq_zero_iff.mp hiso
    intro X
    constructor
    · intro f₁ f₂ hf
      have hp := Subtype.ext_iff.mp hf
      have hb : ElemFrobCat.Hom.base f₁ = ElemFrobCat.Hom.base f₂ := congrArg Prod.snd hp
      have hc : (f₁ ≫ φ : X ⟶ B) = f₂ ≫ φ := congrArg Prod.fst hp
      letI := isCancelAdd_of_isIntegralMonoid _ ((hpd X.base).1)
      refine ElemFrobCat.Hom.ext hb ?_ ?_
      · have hd := congrArg ElemFrobCat.Hom.div hc
        rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hb, hdeg] at hd
        simp only [PNat.one_coe, one_smul] at hd
        exact add_left_cancel hd
      · have hd := congrArg ElemFrobCat.Hom.deg hc
        rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg] at hd
        exact mul_left_cancel hd
    · rintro ⟨⟨g, hb⟩, hgh⟩
      obtain ⟨k, hk⟩ : ∃ k : X.base ⟶ A.base, k = hb := ⟨hb, rfl⟩
      subst hk
      refine ⟨(⟨k, ElemFrobCat.Hom.div g + Φ.map k U.neg,
        ElemFrobCat.Hom.deg g⟩ : X ⟶ A), ?_⟩
      refine Subtype.ext (Prod.ext ?_ rfl)
      refine ElemFrobCat.Hom.ext hgh.symm ?_ ?_
      · show Φ.map k (ElemFrobCat.Hom.div φ)
            + ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ)
                • (ElemFrobCat.Hom.div g + Φ.map k U.neg)
          = ElemFrobCat.Hom.div g
        rw [hdeg]
        simp only [PNat.one_coe, one_smul, ← hU]
        rw [← add_assoc, add_comm (Φ.map k (U : Φ.val A.base)), add_assoc, ← map_add,
          U.val_neg, map_zero, add_zero]
      · show ElemFrobCat.Hom.deg φ * ElemFrobCat.Hom.deg g = ElemFrobCat.Hom.deg g
        rw [hdeg, one_mul]

/-- ★`Proposition 1.4, (ii)` は Frobenioid でしか言えないが、
`𝔽_Φ` では **pre-Frobenioid の段階で** 同じ結論が出る ——
`𝔽_Φ` のすべての射が co-angular だから、LB-invertible = isometric。 -/
theorem elemFrob_isPullBack_iff_lbInvertible {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsPullBack (elemPreFrobenioid Φ hD hpd) φ ↔
      (IsLBInvertible (elemPreFrobenioid Φ hD hpd) φ ∧
        IsLinear (elemPreFrobenioid Φ hD hpd) φ) := by
  rw [elemFrob_isPullBack_iff]
  exact ⟨fun h => ⟨⟨elemFrob_coAngular Φ hD hpd φ, h.2⟩, h.1⟩, fun h => ⟨h.2, h.1.2⟩⟩

/-! ## ★第3段 —— `Proposition 1.5, (i)` の 7 つの type

原文 (FrdI p.27):
> all objects of FΦ are Aut-ample, Autsub-ample, End-ample, base-trivial, Frobenius-

原文 (FrdI p.27):
> trivial, Frobenius-normalized, and isotropic. Also, one verifies immediately [cf. the

★`isotropic` は既に `elemFrob_isotropic` で示した。残り 6 つをここで示す。
-/

/-- ★**End-ample** —— `⟨g, 0, 1⟩` を取ればよい。 -/
theorem elemFrob_endAmple (A : ElemFrobCat Φ) : IsEndAmple (elemPreFrobenioid Φ hD hpd) A := by
  intro g
  obtain ⟨k, hk⟩ : ∃ k : A.base ⟶ A.base, k = g := ⟨g, rfl⟩
  subst hk
  exact ⟨(⟨k, 0, 1⟩ : A ⟶ A), rfl⟩

/-- ★**Aut-ample** —— `⟨g, 0, 1⟩` は `g` が同型なら同型(`ElemFrobCat.isIso_iff`)。 -/
theorem elemFrob_autAmple (A : ElemFrobCat Φ) : IsAutAmple (elemPreFrobenioid Φ hD hpd) A := by
  intro g hg
  obtain ⟨k, hk⟩ : ∃ k : A.base ⟶ A.base, k = g := ⟨g, rfl⟩
  subst hk
  exact ⟨(⟨k, 0, 1⟩ : A ⟶ A),
    (ElemFrobCat.isIso_iff _).mpr ⟨hg, isAddUnit_zero, rfl⟩, rfl⟩

/-- ★**Aut^sub-ample** —— `𝒟` 側の証人 `(Y, ψ, β)` を
`(⟨Y⟩, ⟨ψ, 0, 1⟩, ⟨β, 0, 1⟩)` へ**そのまま持ち上げる**。

★零因子をすべて `0` に取ると、`β ≫ ψ = ψ ≫ g` の零因子成分が
両辺 `0` になるので、`𝒟` の等式だけで済む。 -/
theorem elemFrob_autSubAmple (A : ElemFrobCat Φ) :
    IsAutSubAmple (elemPreFrobenioid Φ hD hpd) A := by
  intro g hg
  obtain ⟨k, hk⟩ : ∃ k : A.base ⟶ A.base, k = g := ⟨g, rfl⟩
  subst hk
  obtain ⟨Y, ψ, β, hβ, hcomm⟩ := hg
  refine ⟨(⟨k, 0, 1⟩ : A ⟶ A), ⟨(⟨Y⟩ : ElemFrobCat Φ), (⟨ψ, 0, 1⟩ : (⟨Y⟩ : ElemFrobCat Φ) ⟶ A),
    (⟨β, 0, 1⟩ : (⟨Y⟩ : ElemFrobCat Φ) ⟶ (⟨Y⟩ : ElemFrobCat Φ)), ?_, ?_⟩, rfl⟩
  · exact (ElemFrobCat.isIso_iff _).mpr ⟨hβ, isAddUnit_zero, rfl⟩
  · refine ElemFrobCat.Hom.ext hcomm ?_ rfl
    rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp]
    simp

/-- ★**base-trivial** —— 底の同型 `e` から `⟨e.inv, 0, 1⟩` を作る。 -/
theorem elemFrob_baseTrivial (A : ElemFrobCat Φ) :
    IsBaseTrivial (elemPreFrobenioid Φ hD hpd) A := by
  rintro Dd ⟨e⟩
  obtain ⟨k, hk⟩ : ∃ k : Dd.base ⟶ A.base, k = e.inv := ⟨e.inv, rfl⟩
  haveI : IsIso k := hk ▸ inferInstanceAs (IsIso e.inv)
  haveI : IsIso (⟨k, 0, 1⟩ : Dd ⟶ A) :=
    (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, isAddUnit_zero, rfl⟩
  exact ⟨asIso (⟨k, 0, 1⟩ : Dd ⟶ A)⟩

/-- ★**Frobenius-trivial** —— `ζ n := ⟨𝟙, 0, n⟩` がモノイド準同型になる。 -/
theorem elemFrob_frobeniusTrivial (A : ElemFrobCat Φ) :
    IsFrobeniusTrivial (elemPreFrobenioid Φ hD hpd) A := by
  refine ⟨⟨⟨fun n => (⟨𝟙 A.base, 0, n⟩ : A ⟶ A), ?_⟩, ?_⟩, fun n => rfl, fun n => ⟨?_, ?_, ?_⟩⟩
  · exact (ElemFrobCat.Hom.ext (ElemFrobCat.id_base A).symm
      (ElemFrobCat.id_div A).symm (ElemFrobCat.id_deg A).symm)
  · intro m n
    show (⟨𝟙 A.base, 0, m * n⟩ : A ⟶ A)
      = ((⟨𝟙 A.base, 0, n⟩ : A ⟶ A) ≫ (⟨𝟙 A.base, 0, m⟩ : A ⟶ A) : A ⟶ A)
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base]; simp
    · rw [ElemFrobCat.div_comp]; simp
    · rw [ElemFrobCat.comp_deg]
  · exact ((elemPreFrobenioid Φ hD hpd).Base_id A).symm
  · exact ⟨elemFrob_coAngular Φ hD hpd _, map_zero _⟩
  · show IsIso (𝟙 A.base)
    infer_instance

/-- `otriOf` の冪 —— モノイド同型 `otriEquiv` の帰結だが、直接示すほうが短い。 -/
theorem otriOf_pow (A : ElemFrobCat Φ) (a : Φ.val A.base) (n : ℕ) :
    (otriOf Φ A a) ^ n = otriOf Φ A (n • a) := by
  induction n with
  | zero =>
    rw [pow_zero, zero_nsmul]
    exact ElemFrobCat.Hom.ext (ElemFrobCat.id_base A) (ElemFrobCat.id_div A)
      (ElemFrobCat.id_deg A)
  | succ n ih =>
    rw [pow_succ, ih]
    refine ElemFrobCat.Hom.ext (by simp) ?_ (by simp)
    show Φ.map (𝟙 A.base) (n • a) + ((1 : ℕ+) : ℕ) • a = (n + 1) • a
    simp [MonoidOn.map_id, succ_nsmul]

/-- ★**Frobenius-normalized** —— `φ ≫ α^d = α ≫ φ`。

零因子成分は左辺が `d·a + Div(φ)`、右辺が `Div(φ) + d·a` で、
**加法の可換性だけ**で一致する。 -/
theorem elemFrob_frobeniusNormalized (A : ElemFrobCat Φ) :
    IsFrobeniusNormalized (elemPreFrobenioid Φ hD hpd) A := by
  intro φ hφ α hα
  have hb : ElemFrobCat.Hom.base φ = 𝟙 A.base :=
    hφ.trans ((elemPreFrobenioid Φ hD hpd).Base_id A)
  have hαe : α = otriOf Φ A (ElemFrobCat.Hom.div α) := (mem_otri_iff Φ hD hpd α).mp hα
  rw [hαe, otriOf_pow]
  refine ElemFrobCat.Hom.ext (by rw [ElemFrobCat.comp_base, ElemFrobCat.comp_base, hb]; simp) ?_
    (by rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg]; simp)
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hb]
  simp only [otriOf_div, otriOf_deg, otriOf_base, MonoidOn.map_id, PNat.one_coe, one_smul]
  exact add_comm _ _

/-! ## ★第5段 —— `Proposition 1.5, (iii)`: perfect / group-like type

原文 (FrdI p.27):
> (iii) If all of the monoids in the image of Φ are perfect (respectively, group-

原文 (FrdI p.27):
> like), then FΦ is of perfect (respectively, group-like) type.
-/

/-- ★**group-like** —— `Φ(X)` がすべて群なら、`𝔽_Φ` は group-like type。

★`IsGroupLikeObj` は **`Φ^char`** についての条件である(pre-Frobenioid の
divisor monoid が `Φ^char` だから)。`Φ(X)` が群なら `Φ(X)^char` は零で、
その `char` も零。 -/
theorem elemFrob_groupLike (hgl : ∀ X : D, IsGroupLike (Φ.val X)) (A : ElemFrobCat Φ) :
    IsGroupLikeObj (elemPreFrobenioid Φ hD hpd) A := by
  have key : ∀ a : MChar (Φ.val A.base), IsAddUnit a := by
    intro a
    obtain ⟨m, hm⟩ := toChar_surjective (Φ.val A.base) a
    rw [← hm, toChar_eq_zero_iff.mpr ((isGroupLike_iff _).mp (hgl A.base) m)]
    exact isAddUnit_zero
  exact (isGroupLike_iff _).mpr key

/-- ★**perfect** —— `Φ(X)` がすべて perfect なら、`𝔽_Φ` は perfect type。

★★**perfect 性が効くのは1箇所だけ**である —— pre-step `ψ'` を `B₁ ⟶ B₂` へ
降ろすとき、零因子について `n · Div(ψ) = t` を解く必要があり、
**`n` 倍の全単射性がちょうどそれを与える**(存在が全射性、一意性が単射性)。

★条件 (a)(「各次数の Frobenius 型射の終域として現れる」)には perfect 性は要らない。
★底成分と次数成分も perfect 性を使わない —— 底は `Base(φ₂)` が同型だから、
次数は `ℕ+` の簡約だから決まる。 -/
theorem elemFrob_perfect (hpf : ∀ X : D, IsPerfectMonoid (Φ.val X)) (A : ElemFrobCat Φ) :
    IsPerfectObj (elemPreFrobenioid Φ hD hpd) A := by
  intro n
  refine ⟨fun B _ => ⟨B, (⟨𝟙 B.base, 0, n⟩ : B ⟶ B),
    ⟨⟨elemFrob_coAngular Φ hD hpd _, map_zero _⟩, ?_⟩, rfl⟩, ?_⟩
  · show IsIso (𝟙 B.base)
    infer_instance
  intro B₁ B₁' B₂ B₂' φ₁ φ₂ hf₁ hd₁ hf₂ hd₂ _ _ ψ' hψ'
  haveI : IsIso (ElemFrobCat.Hom.base φ₁) := hf₁.2
  haveI : IsIso (ElemFrobCat.Hom.base φ₂) := hf₂.2
  haveI : IsIso (ElemFrobCat.Hom.base ψ') := hψ'.2
  obtain ⟨V, hV⟩ := toChar_eq_zero_iff.mp hf₂.1.2
  set b : B₁.base ⟶ B₂.base :=
    ElemFrobCat.Hom.base φ₁ ≫ ElemFrobCat.Hom.base ψ' ≫ inv (ElemFrobCat.Hom.base φ₂) with hb
  haveI : IsIso b := by rw [hb]; infer_instance
  set t : Φ.val B₁.base :=
    Φ.map (ElemFrobCat.Hom.base φ₁) (ElemFrobCat.Hom.div ψ')
      + ElemFrobCat.Hom.div φ₁ + Φ.map b V.neg with ht
  obtain ⟨d, hdt⟩ := (hpf B₁.base n).2 t
  -- 零因子の等式は、`Div(ψ) = d` のとき両辺一致する
  have hdiv : ∀ e : Φ.val B₁.base, ((n : ℕ+) : ℕ) • e = t →
      Φ.map (ElemFrobCat.Hom.base φ₁) (ElemFrobCat.Hom.div ψ') + ElemFrobCat.Hom.div φ₁
        = Φ.map b (ElemFrobCat.Hom.div φ₂) + ((n : ℕ+) : ℕ) • e := by
    intro e he
    rw [he, ht, ← hV, ← add_assoc, add_comm (Φ.map b (V : Φ.val B₂.base)), add_assoc,
      add_assoc, ← map_add, V.val_neg, map_zero, add_zero]
  have hstep : IsPreStep (elemPreFrobenioid Φ hD hpd) (⟨b, d, 1⟩ : B₁ ⟶ B₂) :=
    ⟨rfl, (inferInstance : IsIso b)⟩
  have heq : (φ₁ ≫ ψ' : B₁ ⟶ B₂') = (⟨b, d, 1⟩ : B₁ ⟶ B₂) ≫ φ₂ := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base, ElemFrobCat.comp_base, hb]
      simp
    · rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp]
      show Φ.map (ElemFrobCat.Hom.base φ₁) (ElemFrobCat.Hom.div ψ')
          + ((ElemFrobCat.Hom.deg ψ' : ℕ+) : ℕ) • ElemFrobCat.Hom.div φ₁
        = Φ.map b (ElemFrobCat.Hom.div φ₂)
          + ((ElemFrobCat.Hom.deg φ₂ : ℕ+) : ℕ) • d
      rw [show ElemFrobCat.Hom.deg ψ' = 1 from hψ'.1,
        show ElemFrobCat.Hom.deg φ₂ = n from hd₂]
      simp only [PNat.one_coe, one_smul]
      exact hdiv d hdt
    · rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg,
        show ElemFrobCat.Hom.deg ψ' = 1 from hψ'.1,
        show ElemFrobCat.Hom.deg φ₁ = n from hd₁,
        show ElemFrobCat.Hom.deg φ₂ = n from hd₂]
      rw [one_mul, mul_one]
  refine ⟨(⟨b, d, 1⟩ : B₁ ⟶ B₂), ⟨hstep, heq⟩, ?_⟩
  rintro ψ₀ ⟨hstep₀, heq₀⟩
  have hbase₀ : ElemFrobCat.Hom.base ψ₀ = b := by
    have hh := congrArg ElemFrobCat.Hom.base heq₀
    rw [ElemFrobCat.comp_base, ElemFrobCat.comp_base] at hh
    calc ElemFrobCat.Hom.base ψ₀
        = (ElemFrobCat.Hom.base ψ₀ ≫ ElemFrobCat.Hom.base φ₂)
            ≫ inv (ElemFrobCat.Hom.base φ₂) := by simp
      _ = (ElemFrobCat.Hom.base φ₁ ≫ ElemFrobCat.Hom.base ψ')
            ≫ inv (ElemFrobCat.Hom.base φ₂) := by rw [hh]
      _ = b := by rw [hb, Category.assoc]
  refine ElemFrobCat.Hom.ext hbase₀ ?_ hstep₀.1
  -- ★零因子: `n` 倍の**単射性**でここが決まる
  refine (hpf B₁.base n).1 ?_
  have hh := congrArg ElemFrobCat.Hom.div heq₀
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hbase₀,
    show ElemFrobCat.Hom.deg ψ' = 1 from hψ'.1,
    show ElemFrobCat.Hom.deg φ₂ = n from hd₂] at hh
  simp only [PNat.one_coe, one_smul] at hh
  letI := isCancelAdd_of_isIntegralMonoid _ ((hpd B₁.base).1)
  show ((n : ℕ+) : ℕ) • ElemFrobCat.Hom.div ψ₀ = ((n : ℕ+) : ℕ) • d
  have := hh.symm.trans (hdiv d hdt)
  exact add_left_cancel this

/-! ## ★第4段 —— `Definition 1.3` の 22 条

原文 (FrdI p.27):
> the conditions of Definition 1.3 now follows immediately from the definition of the

原文 (FrdI p.27):
> category FΦ in Definition 1.1, (iii), together with assertion (ii) and the “explicit

★原文が挙げる入力は3つ: **`Definition 1.1, (iii)` の定義**、**(ii)**、
**`Proposition 1.4, (i)` の explicit description**。
★**1条ずつ独立の定理として示す**(構造体を途中まで埋めることはできないため)。
-/

/-- **(i)(a)** —— `⟨Y⟩` を取れば底はちょうど `Y`、しかも Frobenius-trivial。 -/
theorem elemFrob_baseSurj (Y : D) :
    ∃ A : ElemFrobCat Φ, IsFrobeniusTrivial (elemPreFrobenioid Φ hD hpd) A ∧
      Nonempty (((elemPreFrobenioid Φ hD hpd).toElem.obj A).base ≅ Y) :=
  ⟨(⟨Y⟩ : ElemFrobCat Φ), elemFrob_frobeniusTrivial Φ hD hpd _, ⟨Iso.refl Y⟩⟩

/-- **(i)(b)** —— span は `(𝟙_A, ⟨α, 0, 1⟩)` で取れる。★`𝒟` 側の同型が
そのまま `𝔽_Φ` の pre-step に持ち上がる。 -/
theorem elemFrob_preStepSpan (A B : ElemFrobCat Φ)
    (α : ((elemPreFrobenioid Φ hD hpd).toElem.obj A).base ⟶
      ((elemPreFrobenioid Φ hD hpd).toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : ElemFrobCat Φ) (φ : X ⟶ A) (ψ : X ⟶ B)
      (hφ : IsPreStep (elemPreFrobenioid Φ hD hpd) φ),
      IsPreStep (elemPreFrobenioid Φ hD hpd) ψ ∧
        α = @inv _ _ _ _ ((elemPreFrobenioid Φ hD hpd).Base φ) hφ.2
              ≫ (elemPreFrobenioid Φ hD hpd).Base ψ := by
  refine ⟨A, 𝟙 A, (⟨α, 0, 1⟩ : A ⟶ B), isPreStep_id _ A, ⟨rfl, hα⟩, ?_⟩
  refine (@IsIso.eq_inv_comp _ _ _ _ _ ((elemPreFrobenioid Φ hD hpd).Base (𝟙 A))
    (isPreStep_id (elemPreFrobenioid Φ hD hpd) A).2 _ _).mpr ?_
  rw [(elemPreFrobenioid Φ hD hpd).Base_id A, Category.id_comp]
  rfl

/-- **(ii)** —— `⟨𝟙, 0, n⟩` が次数 `n` の Frobenius 型射。 -/
theorem elemFrob_frobDegSurj (A : ElemFrobCat Φ) (n : ℕ+) :
    ∃ (B : ElemFrobCat Φ) (φ : A ⟶ B),
      IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ ∧
        (elemPreFrobenioid Φ hD hpd).degFr φ = n := by
  refine ⟨A, (⟨𝟙 A.base, 0, n⟩ : A ⟶ A),
    ⟨⟨elemFrob_coAngular Φ hD hpd _, map_zero _⟩, ?_⟩, rfl⟩
  show IsIso (𝟙 A.base)
  infer_instance

/-- **(iii)(a)** —— `𝔽_Φ` ではすべての射が co-angular なので自明。 -/
theorem elemFrob_coAngularComp {A B E : ElemFrobCat Φ} (ψ : A ⟶ B) (φ : B ⟶ E) :
    IsCoAngular (elemPreFrobenioid Φ hD hpd) ψ → IsCoAngular (elemPreFrobenioid Φ hD hpd) φ →
      IsCoAngular (elemPreFrobenioid Φ hD hpd) (ψ ≫ φ) :=
  fun _ _ => elemFrob_coAngular Φ hD hpd _

/-- **(iii)(b)** —— 同上。 -/
theorem elemFrob_coAngularOfPreStep {A' A : ElemFrobCat Φ} (α : A' ⟶ A) :
    IsCoAngular (elemPreFrobenioid Φ hD hpd) α →
      IsPreStep (elemPreFrobenioid Φ hD hpd) α →
        ∀ φ : A' ⟶ A, IsCoAngular (elemPreFrobenioid Φ hD hpd) φ :=
  fun _ _ φ => elemFrob_coAngular Φ hD hpd φ

/-- **(iv)(b)** —— 第2段の帰結。 -/
theorem elemFrob_pullBackLB {A B : ElemFrobCat Φ} (α : A ⟶ B)
    (h : IsPullBack (elemPreFrobenioid Φ hD hpd) α) :
    IsLBInvertible (elemPreFrobenioid Φ hD hpd) α ∧
      IsLinear (elemPreFrobenioid Φ hD hpd) α :=
  (elemFrob_isPullBack_iff_lbInvertible Φ hD hpd α).mp h

/-- **(iv)(a)** —— 分解は `(⟨𝟙,0,deg⟩, ⟨𝟙,div,1⟩, ⟨base,0,1⟩)` で取れる。

★**3成分がそのまま3因子に分かれる**のが `𝔽_Φ` の特徴である ——
Frobenius 型が次数を、pre-step が零因子を、pull-back が底を担う。 -/
theorem elemFrob_arbFactor {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    ∃ (X Y : ElemFrobCat Φ) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (elemPreFrobenioid Φ hD hpd) γ ∧
        IsPreStep (elemPreFrobenioid Φ hD hpd) β ∧
        IsPullBack (elemPreFrobenioid Φ hD hpd) α := by
  refine ⟨A, A, (⟨𝟙 A.base, 0, ElemFrobCat.Hom.deg φ⟩ : A ⟶ A),
    (⟨𝟙 A.base, ElemFrobCat.Hom.div φ, 1⟩ : A ⟶ A),
    (⟨ElemFrobCat.Hom.base φ, 0, 1⟩ : A ⟶ B), ?_,
    ⟨⟨elemFrob_coAngular Φ hD hpd _, map_zero _⟩, ?_⟩, ⟨rfl, ?_⟩,
    (elemFrob_isPullBack_iff Φ hD hpd _).mpr ⟨rfl, map_zero _⟩⟩
  · refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base, ElemFrobCat.comp_base]
      simp
    · rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp]
      simp
    · rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg]
      simp
  · show IsIso (𝟙 A.base)
    infer_instance
  · show IsIso (𝟙 A.base)
    infer_instance

/-- **(v)(b)** —— `φ = φ ≫ 𝟙`。★すべての射が co-angular なので、
「co-angular pre-step のあと isometric pre-step」は `𝟙` を足すだけでよい。 -/
theorem elemFrob_preStepFactor {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hφ : IsPreStep (elemPreFrobenioid Φ hD hpd) φ) :
    ∃ (X : ElemFrobCat Φ) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular (elemPreFrobenioid Φ hD hpd) β ∧
        IsPreStep (elemPreFrobenioid Φ hD hpd) β ∧
        IsIsometric (elemPreFrobenioid Φ hD hpd) α ∧
        IsPreStep (elemPreFrobenioid Φ hD hpd) α :=
  ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, elemFrob_coAngular Φ hD hpd φ, hφ,
    (elemPreFrobenioid Φ hD hpd).Div_id B, isPreStep_id _ B⟩

/-- **(v)(c)** —— `φ = 𝟙 ≫ φ`。 -/
theorem elemFrob_preStepFactor' {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hφ : IsPreStep (elemPreFrobenioid Φ hD hpd) φ) :
    ∃ (X : ElemFrobCat Φ) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric (elemPreFrobenioid Φ hD hpd) β ∧
        IsPreStep (elemPreFrobenioid Φ hD hpd) β ∧
        IsCoAngular (elemPreFrobenioid Φ hD hpd) α ∧
        IsPreStep (elemPreFrobenioid Φ hD hpd) α :=
  ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm, (elemPreFrobenioid Φ hD hpd).Div_id A,
    isPreStep_id _ A, elemFrob_coAngular Φ hD hpd φ, hφ⟩

/-- **(vii)(a)** —— すべての対象が isotropic なので、`𝟙` が isotropic hull。 -/
theorem elemFrob_isotropicHullExists (A : ElemFrobCat Φ) :
    ∃ (B : ElemFrobCat Φ) (φ : A ⟶ B), IsIsotropicHull (elemPreFrobenioid Φ hD hpd) φ :=
  ⟨A, 𝟙 A, (elemPreFrobenioid Φ hD hpd).Div_id A, isPreStep_id _ A,
    elemFrob_isotropic Φ hD hpd A,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

/-- **(vii)(b)** —— 同上。 -/
theorem elemFrob_isotropicClosed {A B : ElemFrobCat Φ} (_φ : A ⟶ B) :
    IsIsotropic (elemPreFrobenioid Φ hD hpd) A → IsIsotropic (elemPreFrobenioid Φ hD hpd) B :=
  fun _ => elemFrob_isotropic Φ hD hpd B

/-- **(ii)** の本質的一意性 —— `β := ⟨inv(Base φ) ≫ Base ψ, Φ(inv Base φ)(Div ψ − Div φ), 1⟩`。

★零因子の「引き算」は、`Frobenius 型 ⟹ 零因子が可逆`(第2段と同じ事実)で意味を持つ。 -/
theorem elemFrob_frobDegUniq (A B E : ElemFrobCat Φ) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ)
    (hψ : IsFrobeniusType (elemPreFrobenioid Φ hD hpd) ψ)
    (hd : (elemPreFrobenioid Φ hD hpd).degFr φ = (elemPreFrobenioid Φ hD hpd).degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI : IsIso (ElemFrobCat.Hom.base φ) := hφ.2
  haveI : IsIso (ElemFrobCat.Hom.base ψ) := hψ.2
  obtain ⟨U, hU⟩ := toChar_eq_zero_iff.mp hφ.1.2
  have hψu : IsAddUnit (ElemFrobCat.Hom.div ψ) := toChar_eq_zero_iff.mp hψ.1.2
  have hcu : IsAddUnit (ElemFrobCat.Hom.div ψ + U.neg) := hψu.add ⟨-U, rfl⟩
  refine ⟨(⟨inv (ElemFrobCat.Hom.base φ) ≫ ElemFrobCat.Hom.base ψ,
    Φ.map (inv (ElemFrobCat.Hom.base φ)) (ElemFrobCat.Hom.div ψ + U.neg), 1⟩ : B ⟶ E), ?_, ?_⟩
  · exact (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, hcu.map _, rfl⟩
  · refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base]
      simp
    · rw [ElemFrobCat.div_comp]
      show Φ.map (ElemFrobCat.Hom.base φ)
            (Φ.map (inv (ElemFrobCat.Hom.base φ)) (ElemFrobCat.Hom.div ψ + U.neg))
          + ((1 : ℕ+) : ℕ) • ElemFrobCat.Hom.div φ = ElemFrobCat.Hom.div ψ
      rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
      simp only [PNat.one_coe, one_smul, ← hU, add_assoc, U.neg_val, add_zero]
    · rw [ElemFrobCat.comp_deg, one_mul]
      exact hd

/-- **(v)(a)** —— pre-step は monomorphism。

★3成分すべてが個別に決まる: 底は `Base φ` が同型(epi)、次数は `ℕ+` の簡約、
零因子は `Φ(A)` が integral(= cancel)であること。 -/
theorem elemFrob_preStepMono {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hφ : IsPreStep (elemPreFrobenioid Φ hD hpd) φ) : Mono φ := by
  haveI : IsIso (ElemFrobCat.Hom.base φ) := hφ.2
  refine ⟨fun {Z} f g hfg => ?_⟩
  letI := isCancelAdd_of_isIntegralMonoid _ ((hpd Z.base).1)
  have hb : ElemFrobCat.Hom.base f = ElemFrobCat.Hom.base g := by
    have hh := congrArg ElemFrobCat.Hom.base hfg
    rw [ElemFrobCat.comp_base, ElemFrobCat.comp_base] at hh
    exact (cancel_mono (ElemFrobCat.Hom.base φ)).mp hh
  refine ElemFrobCat.Hom.ext hb ?_ ?_
  · have hh := congrArg ElemFrobCat.Hom.div hfg
    rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hb,
      show ElemFrobCat.Hom.deg φ = 1 from hφ.1] at hh
    simp only [PNat.one_coe, one_smul] at hh
    exact add_left_cancel hh
  · have hh := congrArg ElemFrobCat.Hom.deg hfg
    rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg] at hh
    exact mul_left_cancel hh

/-- **(iii)(c)** 順方向 —— `𝒪^▷(A) → 𝒪^▷(B)` の全単射。

★`Φ.map (Base φ)` が **`Base φ` が同型なので全単射**であることに帰着する。 -/
theorem elemFrob_otriFwd {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hst : IsPreStep (elemPreFrobenioid Φ hD hpd) φ)
    (α : End A) (hα : α ∈ OTri (elemPreFrobenioid Φ hD hpd) A) :
    ∃! β : End B, β ∈ OTri (elemPreFrobenioid Φ hD hpd) B ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  haveI : IsIso (ElemFrobCat.Hom.base φ) := hst.2
  have hαe := (mem_otri_iff Φ hD hpd α).mp hα
  have hkey : ∀ b : Φ.val B.base,
      Φ.map (ElemFrobCat.Hom.base φ) b
          = ((ElemFrobCat.Hom.deg φ : ℕ+) : ℕ) • ElemFrobCat.Hom.div α ↔
        b = Φ.map (inv (ElemFrobCat.Hom.base φ)) (ElemFrobCat.Hom.div α) := by
    intro b
    rw [show ElemFrobCat.Hom.deg φ = 1 from hst.1]
    simp only [PNat.one_coe, one_smul]
    constructor
    · intro h
      rw [← h, ← MonoidOn.map_comp, IsIso.inv_hom_id, MonoidOn.map_id]
    · intro h
      rw [h, ← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
  refine ⟨otriOf Φ B (Φ.map (inv (ElemFrobCat.Hom.base φ)) (ElemFrobCat.Hom.div α)),
    ⟨otriOf_mem Φ hD hpd B _, ?_⟩, ?_⟩
  · rw [hαe]
    exact (otriOf_comm_iff Φ hpd φ _ _).mpr ((hkey _).mpr rfl)
  · rintro β ⟨hβ, hβe⟩
    have hβ2 := (mem_otri_iff Φ hD hpd β).mp hβ
    rw [hαe, hβ2] at hβe
    rw [hβ2]
    exact congrArg (otriOf Φ B) ((hkey _).mp ((otriOf_comm_iff Φ hpd φ _ _).mp hβe))

/-- **(iii)(c)** 逆方向。★こちらは `Φ.map (Base φ)` を**当てるだけ**で決まる。 -/
theorem elemFrob_otriBwd {A B : ElemFrobCat Φ} (φ : A ⟶ B)
    (hst : IsPreStep (elemPreFrobenioid Φ hD hpd) φ)
    (β : End B) (hβ : β ∈ OTri (elemPreFrobenioid Φ hD hpd) B) :
    ∃! α : End A, α ∈ OTri (elemPreFrobenioid Φ hD hpd) A ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  haveI : IsIso (ElemFrobCat.Hom.base φ) := hst.2
  have hβe := (mem_otri_iff Φ hD hpd β).mp hβ
  have hdeg : ElemFrobCat.Hom.deg φ = 1 := hst.1
  refine ⟨otriOf Φ A (Φ.map (ElemFrobCat.Hom.base φ) (ElemFrobCat.Hom.div β)),
    ⟨otriOf_mem Φ hD hpd A _, ?_⟩, ?_⟩
  · rw [hβe]
    refine (otriOf_comm_iff Φ hpd φ _ _).mpr ?_
    rw [hdeg]
    simp
  · rintro α ⟨hα, hαe⟩
    have hα2 := (mem_otri_iff Φ hD hpd α).mp hα
    rw [hβe, hα2] at hαe
    rw [hα2]
    refine congrArg (otriOf Φ A) ?_
    have := (otriOf_comm_iff Φ hpd φ _ _).mp hαe
    rw [hdeg] at this
    simpa using this.symm

/-- **(iii)(c)** 全単射は `Base(φ)` にしか依らない。

★`otriOf_comm_iff` の右辺が `Base φ` と `deg_Fr φ` しか含まないので、
**そのまま従う**。 -/
theorem elemFrob_otriBase {A B : ElemFrobCat Φ} (φ φ' : A ⟶ B)
    (hst : IsPreStep (elemPreFrobenioid Φ hD hpd) φ)
    (hst' : IsPreStep (elemPreFrobenioid Φ hD hpd) φ')
    (hbase : (elemPreFrobenioid Φ hD hpd).Base φ = (elemPreFrobenioid Φ hD hpd).Base φ')
    (α : End A) (hα : α ∈ OTri (elemPreFrobenioid Φ hD hpd) A)
    (β : End B) (hβ : β ∈ OTri (elemPreFrobenioid Φ hD hpd) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' := by
  have hαe := (mem_otri_iff Φ hD hpd α).mp hα
  have hβe := (mem_otri_iff Φ hD hpd β).mp hβ
  rw [hαe, hβe] at h ⊢
  refine (otriOf_comm_iff Φ hpd φ' _ _).mpr ?_
  have hh := (otriOf_comm_iff Φ hpd φ _ _).mp h
  rw [show ElemFrobCat.Hom.deg φ = 1 from hst.1] at hh
  rw [show ElemFrobCat.Hom.deg φ' = 1 from hst'.1]
  rw [show ElemFrobCat.Hom.base φ' = ElemFrobCat.Hom.base φ from hbase.symm]
  exact hh

/-- **(vi)** 単元を除く忠実性 —— `Div` が `Φ^char` で一致するので、
零因子は **`Φ(A)` の可逆元の差**しかない。それを `𝒪^×(B)` の元にする。 -/
theorem elemFrob_faithfulUpToUnits {A B : ElemFrobCat Φ} (φ ψ : A ⟶ B)
    (hb : BaseEquivalent (elemPreFrobenioid Φ hD hpd) φ ψ)
    (hm : MetricallyEquivalent (elemPreFrobenioid Φ hD hpd) φ ψ)
    (hφ : IsPreStep (elemPreFrobenioid Φ hD hpd) φ)
    (hψ : IsPreStep (elemPreFrobenioid Φ hD hpd) ψ) :
    ∃ α : End B, α ∈ OTimes (elemPreFrobenioid Φ hD hpd) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  haveI : IsIso (ElemFrobCat.Hom.base ψ) := hψ.2
  obtain ⟨u, v, ⟨U, hU⟩, ⟨V, hV⟩, huv⟩ := toChar_eq_iff.mp hm
  have hv : IsAddUnit v := ⟨V, hV⟩
  have hcu : IsAddUnit (v + U.neg) := hv.add ⟨-U, rfl⟩
  refine ⟨otriOf Φ B (Φ.map (inv (ElemFrobCat.Hom.base ψ)) (v + U.neg)),
    (otriOf_mem_otimes Φ hD hpd B _).mpr (hcu.map _), ?_⟩
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · rw [ElemFrobCat.comp_base]
    show ElemFrobCat.Hom.base φ = ElemFrobCat.Hom.base ψ ≫ 𝟙 B.base
    rw [Category.comp_id]
    exact hb
  · rw [ElemFrobCat.div_comp]
    show ElemFrobCat.Hom.div φ
      = Φ.map (ElemFrobCat.Hom.base ψ)
          (Φ.map (inv (ElemFrobCat.Hom.base ψ)) (v + U.neg))
        + ((1 : ℕ+) : ℕ) • ElemFrobCat.Hom.div ψ
    rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
    simp only [PNat.one_coe, one_smul]
    have hz : u + U.neg = 0 := by rw [← hU]; exact U.val_neg
    calc ElemFrobCat.Hom.div φ
        = ElemFrobCat.Hom.div φ + (u + U.neg) := by rw [hz, add_zero]
      _ = ElemFrobCat.Hom.div φ + u + U.neg := by rw [add_assoc]
      _ = ElemFrobCat.Hom.div ψ + v + U.neg := by rw [huv]
      _ = v + U.neg + ElemFrobCat.Hom.div ψ := by
          rw [add_comm (ElemFrobCat.Hom.div ψ) v, add_assoc,
            add_comm (ElemFrobCat.Hom.div ψ) U.neg, ← add_assoc]
  · rw [ElemFrobCat.comp_deg]
    show ElemFrobCat.Hom.deg φ = 1 * ElemFrobCat.Hom.deg ψ
    rw [one_mul, show ElemFrobCat.Hom.deg φ = 1 from hφ.1,
      show ElemFrobCat.Hom.deg ψ = 1 from hψ.1]

/-- ★★**(v)(b) と (v)(c) に共通する構成**。

原文は (b) と (c) で「isometric ∘ co-angular」「co-angular ∘ isometric」と
順序を入れ替えた2条を並べるが、★**`𝔽_Φ` では同じ構成が両方を与える** ——
必要なのは「`Div β'` が `Div β` と**可逆元だけ違う**」ことだけで、
その可逆元がどちらの因子から来るかは効かない。

同型は `γ := ⟨inv(Base β) ≫ Base β', Φ(inv Base β)(u), 1⟩`。 -/
theorem elemFrob_preStepIso {A B X X' : ElemFrobCat Φ}
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβ : IsPreStep (elemPreFrobenioid Φ hD hpd) β)
    (hα : IsLinear (elemPreFrobenioid Φ hD hpd) α)
    (hβ' : IsPreStep (elemPreFrobenioid Φ hD hpd) β')
    (hα' : IsLinear (elemPreFrobenioid Φ hD hpd) α')
    (u : Φ.val A.base) (huu : IsAddUnit u)
    (hu : ElemFrobCat.Hom.div β' = ElemFrobCat.Hom.div β + u) :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI : IsIso (ElemFrobCat.Hom.base β) := hβ.2
  haveI : IsIso (ElemFrobCat.Hom.base β') := hβ'.2
  letI := isCancelAdd_of_isIntegralMonoid _ ((hpd A.base).1)
  have hbase : ElemFrobCat.Hom.base β ≫ ElemFrobCat.Hom.base α
      = ElemFrobCat.Hom.base β' ≫ ElemFrobCat.Hom.base α' := by
    have hh := congrArg ElemFrobCat.Hom.base heq
    rwa [ElemFrobCat.comp_base, ElemFrobCat.comp_base] at hh
  have hdivq : Φ.map (ElemFrobCat.Hom.base β) (ElemFrobCat.Hom.div α)
      = Φ.map (ElemFrobCat.Hom.base β') (ElemFrobCat.Hom.div α') + u := by
    have hh := congrArg ElemFrobCat.Hom.div heq
    rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp,
      show ElemFrobCat.Hom.deg α = 1 from hα,
      show ElemFrobCat.Hom.deg α' = 1 from hα', hu] at hh
    simp only [PNat.one_coe, one_smul] at hh
    refine add_right_cancel (b := ElemFrobCat.Hom.div β) ?_
    rw [hh, add_assoc, add_comm u (ElemFrobCat.Hom.div β)]
  have hgu : IsAddUnit (Φ.map (inv (ElemFrobCat.Hom.base β)) u) := huu.map _
  haveI hiso : IsIso ((⟨inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β',
      Φ.map (inv (ElemFrobCat.Hom.base β)) u, 1⟩ : X ⟶ X')) :=
    (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, hgu, rfl⟩
  have hβγ : (β ≫ (⟨inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β',
      Φ.map (inv (ElemFrobCat.Hom.base β)) u, 1⟩ : X ⟶ X') : A ⟶ X') = β' := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base]
      show ElemFrobCat.Hom.base β ≫ inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β'
        = ElemFrobCat.Hom.base β'
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    · rw [ElemFrobCat.div_comp]
      show Φ.map (ElemFrobCat.Hom.base β) (Φ.map (inv (ElemFrobCat.Hom.base β)) u)
          + ((1 : ℕ+) : ℕ) • ElemFrobCat.Hom.div β = ElemFrobCat.Hom.div β'
      rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
      simp only [PNat.one_coe, one_smul]
      rw [hu, add_comm]
    · rw [ElemFrobCat.comp_deg]
      show (1 : ℕ+) * ElemFrobCat.Hom.deg β = ElemFrobCat.Hom.deg β'
      rw [one_mul, show ElemFrobCat.Hom.deg β = 1 from hβ.1,
        show ElemFrobCat.Hom.deg β' = 1 from hβ'.1]
  have hγα : ((⟨inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β',
      Φ.map (inv (ElemFrobCat.Hom.base β)) u, 1⟩ : X ⟶ X') ≫ α' : X ⟶ B) = α := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · rw [ElemFrobCat.comp_base]
      show (inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β') ≫ ElemFrobCat.Hom.base α'
        = ElemFrobCat.Hom.base α
      rw [Category.assoc, ← hbase, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    · rw [ElemFrobCat.div_comp]
      show Φ.map (inv (ElemFrobCat.Hom.base β) ≫ ElemFrobCat.Hom.base β')
            (ElemFrobCat.Hom.div α')
          + ((ElemFrobCat.Hom.deg α' : ℕ+) : ℕ) • Φ.map (inv (ElemFrobCat.Hom.base β)) u
        = ElemFrobCat.Hom.div α
      rw [show ElemFrobCat.Hom.deg α' = 1 from hα']
      simp only [PNat.one_coe, one_smul]
      rw [MonoidOn.map_comp, ← map_add, ← hdivq, ← MonoidOn.map_comp,
        IsIso.inv_hom_id, MonoidOn.map_id]
    · rw [ElemFrobCat.comp_deg]
      show ElemFrobCat.Hom.deg α' * 1 = ElemFrobCat.Hom.deg α
      rw [mul_one, show ElemFrobCat.Hom.deg α' = 1 from hα',
        show ElemFrobCat.Hom.deg α = 1 from hα]
  refine ⟨asIso _, ?_, hβγ.symm⟩
  rw [← hγα]
  simp

/-- **(v)(b)** の一意性 —— 可逆元は `Div α`, `Div α'`(isometric)から作る。 -/
theorem elemFrob_preStepFactorUniq {A B : ElemFrobCat Φ} (X X' : ElemFrobCat Φ)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβ : IsPreStep (elemPreFrobenioid Φ hD hpd) β)
    (hαi : IsIsometric (elemPreFrobenioid Φ hD hpd) α)
    (hα : IsLinear (elemPreFrobenioid Φ hD hpd) α)
    (hβ' : IsPreStep (elemPreFrobenioid Φ hD hpd) β')
    (hαi' : IsIsometric (elemPreFrobenioid Φ hD hpd) α')
    (hα' : IsLinear (elemPreFrobenioid Φ hD hpd) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  letI := isCancelAdd_of_isIntegralMonoid _ ((hpd A.base).1)
  obtain ⟨P', hP'⟩ := toChar_eq_zero_iff.mp hαi'
  have hau : IsAddUnit (ElemFrobCat.Hom.div α) := toChar_eq_zero_iff.mp hαi
  have hpn : IsAddUnit P'.neg := ⟨-P', rfl⟩
  have hh := congrArg ElemFrobCat.Hom.div heq
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp,
    show ElemFrobCat.Hom.deg α = 1 from hα,
    show ElemFrobCat.Hom.deg α' = 1 from hα'] at hh
  simp only [PNat.one_coe, one_smul] at hh
  refine elemFrob_preStepIso Φ hD hpd β α β' α' heq hβ hα hβ' hα'
    (Φ.map (ElemFrobCat.Hom.base β) (ElemFrobCat.Hom.div α)
      + Φ.map (ElemFrobCat.Hom.base β') P'.neg)
    ((hau.map _).add (hpn.map _)) ?_
  calc ElemFrobCat.Hom.div β'
      = Φ.map (ElemFrobCat.Hom.base β') (ElemFrobCat.Hom.div α')
          + ElemFrobCat.Hom.div β' + Φ.map (ElemFrobCat.Hom.base β') P'.neg := by
        rw [add_assoc, add_comm (ElemFrobCat.Hom.div β'), ← add_assoc, ← map_add, ← hP',
          P'.val_neg, map_zero, zero_add]
    _ = Φ.map (ElemFrobCat.Hom.base β) (ElemFrobCat.Hom.div α)
          + ElemFrobCat.Hom.div β + Φ.map (ElemFrobCat.Hom.base β') P'.neg := by rw [hh]
    _ = ElemFrobCat.Hom.div β
          + (Φ.map (ElemFrobCat.Hom.base β) (ElemFrobCat.Hom.div α)
            + Φ.map (ElemFrobCat.Hom.base β') P'.neg) := by
        rw [add_assoc, ← add_assoc, add_comm (Φ.map (ElemFrobCat.Hom.base β)
          (ElemFrobCat.Hom.div α)) (ElemFrobCat.Hom.div β), add_assoc]

/-- **(v)(c)** の一意性 —— 可逆元は `Div β`, `Div β'`(isometric)から作る。

★同じ `elemFrob_preStepIso` に落ちる。**原文の2条の差は `u` の作り方だけ**である。 -/
theorem elemFrob_preStepFactorUniq' {A B : ElemFrobCat Φ} (X X' : ElemFrobCat Φ)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβi : IsIsometric (elemPreFrobenioid Φ hD hpd) β)
    (hβ : IsPreStep (elemPreFrobenioid Φ hD hpd) β)
    (hα : IsPreStep (elemPreFrobenioid Φ hD hpd) α)
    (hβi' : IsIsometric (elemPreFrobenioid Φ hD hpd) β')
    (hβ' : IsPreStep (elemPreFrobenioid Φ hD hpd) β')
    (hα' : IsPreStep (elemPreFrobenioid Φ hD hpd) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨Q, hQ⟩ := toChar_eq_zero_iff.mp hβi
  have hbu' : IsAddUnit (ElemFrobCat.Hom.div β') := toChar_eq_zero_iff.mp hβi'
  refine elemFrob_preStepIso Φ hD hpd β α β' α' heq hβ hα.1 hβ' hα'.1
    (ElemFrobCat.Hom.div β' + Q.neg) (hbu'.add ⟨-Q, rfl⟩) ?_
  rw [← add_assoc, add_comm (ElemFrobCat.Hom.div β), add_assoc, ← hQ, Q.val_neg, add_zero]

/-- **(iv)(a)** の分解の一意性。

★★**新しい構成は要らない** —— 原文の3因子分解の一意性は、
既に示した2つに**分解して**出る:

1. `ε` は **(ii) の本質的一意性**(`elemFrob_frobDegUniq`)そのもの。
2. `𝔽_Φ` は totally epimorphic なので `γ` は epi。**消去**して
   `β ≫ α = (ε ≫ β') ≫ α'` を得る。
3. `δ` は **(v)(b) の一意性**(`elemFrob_preStepFactorUniq`)そのもの。
   ★ここで**第2段**(pull-back ⟺ linear isometry)が効く —— `α`, `α'` が
   isometric でなければ (v)(b) の補題に渡せない。 -/
theorem elemFrob_arbFactorUniq {A B : ElemFrobCat Φ} (X Y X' Y' : ElemFrobCat Φ)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : (γ ≫ β ≫ α : A ⟶ B) = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (elemPreFrobenioid Φ hD hpd) γ)
    (hβ : IsPreStep (elemPreFrobenioid Φ hD hpd) β)
    (hα : IsPullBack (elemPreFrobenioid Φ hD hpd) α)
    (hγ' : IsFrobeniusType (elemPreFrobenioid Φ hD hpd) γ')
    (hβ' : IsPreStep (elemPreFrobenioid Φ hD hpd) β')
    (hα' : IsPullBack (elemPreFrobenioid Φ hD hpd) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  obtain ⟨hαl, hαi⟩ := (elemFrob_isPullBack_iff Φ hD hpd α).mp hα
  obtain ⟨hαl', hαi'⟩ := (elemFrob_isPullBack_iff Φ hD hpd α').mp hα'
  have hdegγ : (elemPreFrobenioid Φ hD hpd).degFr γ
      = (elemPreFrobenioid Φ hD hpd).degFr γ' := by
    have hh := congrArg ElemFrobCat.Hom.deg heq
    rw [ElemFrobCat.comp_deg, ElemFrobCat.comp_deg, ElemFrobCat.comp_deg,
      ElemFrobCat.comp_deg,
      show ElemFrobCat.Hom.deg β = 1 from hβ.1,
      show ElemFrobCat.Hom.deg α = 1 from hαl,
      show ElemFrobCat.Hom.deg β' = 1 from hβ'.1,
      show ElemFrobCat.Hom.deg α' = 1 from hαl'] at hh
    simpa using hh
  obtain ⟨εm, hεiso, hγε⟩ := elemFrob_frobDegUniq Φ hD hpd A X X' γ γ' hγ hγ' hdegγ
  obtain ⟨hεb, -, hεn⟩ := (ElemFrobCat.isIso_iff εm).mp hεiso
  haveI := hεiso
  haveI : Epi γ := isTotallyEpimorphic_elemFrobCat hD (fun Z => (hpd Z).1) _ _ γ
  have hcancel : (β ≫ α : X ⟶ B) = (εm ≫ β') ≫ α' := by
    refine (cancel_epi γ).mp ?_
    rw [heq, ← hγε]
    simp
  have hβ'' : IsPreStep (elemPreFrobenioid Φ hD hpd) (εm ≫ β') := by
    refine ⟨?_, ?_⟩
    · show ElemFrobCat.Hom.deg (εm ≫ β') = 1
      rw [ElemFrobCat.comp_deg, hεn, show ElemFrobCat.Hom.deg β' = 1 from hβ'.1, mul_one]
    · haveI hb1 : IsIso (ElemFrobCat.Hom.base εm) := hεb
      haveI hb2 : IsIso (ElemFrobCat.Hom.base β') := hβ'.2
      show IsIso (ElemFrobCat.Hom.base (εm ≫ β'))
      rw [ElemFrobCat.comp_base]
      infer_instance
  obtain ⟨δ, hδ1, hδ2⟩ := elemFrob_preStepFactorUniq Φ hD hpd Y Y' β α (εm ≫ β') α'
    hcancel hβ hαi hαl hβ'' hαi' hαl'
  refine ⟨δ, asIso εm, hδ1, ?_, hγε.symm⟩
  rw [← hδ2]
  simp

/-- **(i)(c)** —— `(𝔽_Φ^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値。

★★**第2段(pull-back ⟺ linear isometry)が3つの成分すべてで効く**:

* **忠実**: pull-back は次数 1 なので `deg` は自動、`div` は `Φ(A)` の簡約で決まる。
* **充満**: `IsPullBack` の**全射性そのもの**が原像を与える。
  持ち上げた射がまた pull-back であることは、次数と零因子の可逆性から出る。
* **本質的全射**: `⟨q, 0, 1⟩` が pull-back(次数 1・零因子 0)。 -/
theorem elemFrob_plBkEquiv (A : ElemFrobCat Φ) :
    (plBkOverFunctor (elemPreFrobenioid Φ hD hpd) A).IsEquivalence := by
  haveI hfaith : (plBkOverFunctor (elemPreFrobenioid Φ hD hpd) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : ElemFrobCat.Hom.base f.left.hom = ElemFrobCat.Hom.base g.left.hom :=
      congrArg CommaMorphism.left hfg
    obtain ⟨hWl, -⟩ := (elemFrob_isPullBack_iff Φ hD hpd W.hom.hom).mp W.hom.property
    obtain ⟨hfl, -⟩ := (elemFrob_isPullBack_iff Φ hD hpd f.left.hom).mp f.left.property
    obtain ⟨hgl, -⟩ := (elemFrob_isPullBack_iff Φ hD hpd g.left.hom).mp g.left.property
    have hwf : (f.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have hwg : (g.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    letI := isCancelAdd_of_isIntegralMonoid _ ((hpd Z.left.obj.base).1)
    have hdiv : ElemFrobCat.Hom.div f.left.hom = ElemFrobCat.Hom.div g.left.hom := by
      have h1 := congrArg ElemFrobCat.Hom.div (hwf.trans hwg.symm)
      rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp, hb,
        show ElemFrobCat.Hom.deg W.hom.hom = 1 from hWl] at h1
      simp only [PNat.one_coe, one_smul] at h1
      exact add_left_cancel h1
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      (ElemFrobCat.Hom.ext hb hdiv (hfl.trans hgl.symm)))
  haveI hfull : (plBkOverFunctor (elemPreFrobenioid Φ hD hpd) A).Full := by
    constructor
    intro Z W h
    obtain ⟨hWl, hWi⟩ := (elemFrob_isPullBack_iff Φ hD hpd W.hom.hom).mp W.hom.property
    obtain ⟨hZl, hZi⟩ := (elemFrob_isPullBack_iff Φ hD hpd Z.hom.hom).mp Z.hom.property
    obtain ⟨S, hS⟩ := toChar_eq_zero_iff.mp hZi
    obtain ⟨f₀, hf₀⟩ := (W.hom.property Z.left.obj).2
      ⟨(Z.hom.hom, h.left), (Over.w h).symm⟩
    have hp := Subtype.ext_iff.mp hf₀
    have h1 : (f₀ ≫ W.hom.hom) = Z.hom.hom := congrArg Prod.fst hp
    have h2 : ElemFrobCat.Hom.base f₀ = h.left := congrArg Prod.snd hp
    have hdeg : ElemFrobCat.Hom.deg f₀ = 1 := by
      have hh := congrArg ElemFrobCat.Hom.deg h1
      rw [ElemFrobCat.comp_deg, show ElemFrobCat.Hom.deg W.hom.hom = 1 from hWl,
        one_mul, show ElemFrobCat.Hom.deg Z.hom.hom = 1 from hZl] at hh
      exact hh
    have hd : Φ.map (ElemFrobCat.Hom.base f₀) (ElemFrobCat.Hom.div W.hom.hom)
        + ElemFrobCat.Hom.div f₀ = ElemFrobCat.Hom.div Z.hom.hom := by
      have hh := congrArg ElemFrobCat.Hom.div h1
      rw [ElemFrobCat.div_comp, show ElemFrobCat.Hom.deg W.hom.hom = 1 from hWl] at hh
      simpa using hh
    have hdivu : IsAddUnit (ElemFrobCat.Hom.div f₀) := by
      refine ⟨⟨ElemFrobCat.Hom.div f₀,
        Φ.map (ElemFrobCat.Hom.base f₀) (ElemFrobCat.Hom.div W.hom.hom) + S.neg, ?_, ?_⟩, rfl⟩
      · rw [← add_assoc, add_comm (ElemFrobCat.Hom.div f₀), hd, ← hS, S.val_neg]
      · rw [add_assoc, add_comm S.neg, ← add_assoc, hd, ← hS, S.val_neg]
    refine ⟨Over.homMk (⟨f₀, (elemFrob_isPullBack_iff Φ hD hpd f₀).mpr
      ⟨hdeg, toChar_eq_zero_iff.mpr hdivu⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext h1), Over.OverMorphism.ext h2⟩
  haveI hess : (plBkOverFunctor (elemPreFrobenioid Φ hD hpd) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨q, hq⟩ : ∃ q : Y.left ⟶ A.base, q = Y.hom := ⟨Y.hom, rfl⟩
    refine ⟨Over.mk (show (⟨(⟨Y.left⟩ : ElemFrobCat Φ)⟩ :
        PlBk (elemPreFrobenioid Φ hD hpd)) ⟶ (⟨A⟩ : PlBk (elemPreFrobenioid Φ hD hpd)) from
      ⟨(⟨q, 0, 1⟩ : (⟨Y.left⟩ : ElemFrobCat Φ) ⟶ A),
        (elemFrob_isPullBack_iff Φ hD hpd _).mpr ⟨rfl, map_zero _⟩⟩), ⟨?_⟩⟩
    refine Over.isoMk (Iso.refl _) ?_
    show 𝟙 Y.left ≫ Y.hom = q
    rw [Category.id_comp, hq]
  exact ⟨hfaith, hfull, hess⟩

/-- ★★★**`Definition 1.3` の core 21 条をすべて満たす**。

原文 (FrdI p.27):
> the conditions of Definition 1.3 now follows immediately from the definition of the
-/
theorem elemFrob_frobenioidCore : FrobenioidCore (elemPreFrobenioid Φ hD hpd) where
  baseSurj := elemFrob_baseSurj Φ hD hpd
  preStepSpan := elemFrob_preStepSpan Φ hD hpd
  plBkEquiv := elemFrob_plBkEquiv Φ hD hpd
  frobDegSurj := elemFrob_frobDegSurj Φ hD hpd
  frobDegUniq := elemFrob_frobDegUniq Φ hD hpd
  coAngularComp := elemFrob_coAngularComp Φ hD hpd
  coAngularOfPreStep := elemFrob_coAngularOfPreStep Φ hD hpd
  otriFwd := fun φ _ hst α hα => elemFrob_otriFwd Φ hD hpd φ hst α hα
  otriBwd := fun φ _ hst β hβ => elemFrob_otriBwd Φ hD hpd φ hst β hβ
  otriBase := fun φ φ' _ hst _ hst' hbase α hα β hβ h =>
    elemFrob_otriBase Φ hD hpd φ φ' hst hst' hbase α hα β hβ h
  arbFactor := elemFrob_arbFactor Φ hD hpd
  arbFactorUniq := fun X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα' =>
    elemFrob_arbFactorUniq Φ hD hpd X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα'
  pullBackLB := elemFrob_pullBackLB Φ hD hpd
  preStepMono := elemFrob_preStepMono Φ hD hpd
  preStepFactor := elemFrob_preStepFactor Φ hD hpd
  preStepFactorUniq := fun X X' β α β' α' heq _ hβ hαi hα _ hβ' hαi' hα' =>
    elemFrob_preStepFactorUniq Φ hD hpd X X' β α β' α' heq hβ hαi hα.1 hβ' hαi' hα'.1
  preStepFactor' := elemFrob_preStepFactor' Φ hD hpd
  preStepFactorUniq' := fun X X' β α β' α' heq hβi hβ _ hα hβi' hβ' _ hα' =>
    elemFrob_preStepFactorUniq' Φ hD hpd X X' β α β' α' heq hβi hβ hα hβi' hβ' hα'
  faithfulUpToUnits := fun φ ψ hb hm _ hφ _ hψ =>
    elemFrob_faithfulUpToUnits Φ hD hpd φ ψ hb hm hφ hψ
  isotropicHullExists := elemFrob_isotropicHullExists Φ hD hpd
  isotropicClosed := elemFrob_isotropicClosed Φ hD hpd

/-! ### ★`Definition 1.3, (iii), (d)` —— 残る 2 本の圏同値

原文 (FrdI p.25):
> categories.

★**忠実性は `𝔽_Φ` が totally epimorphic であること(コスライス側)と
pre-step が mono であること(スライス側)から出る** —— 目標が前順序なので、
射が一意であることを示せばよい。
★**充満性は `Φ^char` の `MLe` を `Φ` へ持ち上げる**作業で、
そこで `CharRel`(可逆元の差)が現れる。 -/

section CoaPre

variable [MorphismProperty.IsMultiplicative (coaPreProp (elemPreFrobenioid Φ hD hpd))]

/-- **(iii)(d)** コスライス側 `_A(𝔽_Φ^coa-pre) → Order(Φ^char(A))` は圏同値。 -/
theorem elemFrob_coaPreUnderEquiv (A : ElemFrobCat Φ) :
    (coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).IsEquivalence := by
  haveI hfaith : (coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    have h2 : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w g)
    haveI : Epi Z.hom.hom := isTotallyEpimorphic_elemFrobCat hD (fun X => (hpd X).1) _ _ _
    exact Under.UnderMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_epi Z.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).Full := by
    constructor
    intro Z W h
    haveI hz : IsIso (ElemFrobCat.Hom.base Z.hom.hom) := Z.hom.property.2.2
    obtain ⟨x, hx⟩ := leOfHom h
    obtain ⟨x₀, hx₀⟩ : ∃ x₀ : Φ.val A.base, toChar x₀ = x := by
      obtain ⟨y, hy⟩ := toChar_surjective _ x
      exact ⟨y, hy⟩
    obtain ⟨u, v, ⟨U, hU⟩, ⟨V, hV⟩, huv⟩ :
        CharRel (ElemFrobCat.Hom.div Z.hom.hom + x₀) (ElemFrobCat.Hom.div W.hom.hom) :=
      toChar_eq_iff.mp (by rw [map_add, hx₀]; exact hx)
    have hkey : Φ.map (ElemFrobCat.Hom.base Z.hom.hom)
          (Φ.map (inv (ElemFrobCat.Hom.base Z.hom.hom)) (x₀ + u + V.neg))
        + ((1 : ℕ+) : ℕ) • ElemFrobCat.Hom.div Z.hom.hom
        = ElemFrobCat.Hom.div W.hom.hom := by
      rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
      simp only [PNat.one_coe, one_smul]
      calc x₀ + u + V.neg + ElemFrobCat.Hom.div Z.hom.hom
          = ElemFrobCat.Hom.div Z.hom.hom + x₀ + u + V.neg := by
            rw [add_comm (x₀ + u + V.neg), ← add_assoc, ← add_assoc]
        _ = ElemFrobCat.Hom.div W.hom.hom + v + V.neg := by rw [huv]
        _ = ElemFrobCat.Hom.div W.hom.hom := by rw [add_assoc, ← hV, V.val_neg, add_zero]
    haveI := Preorder.subsingleton_hom
      ((coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).obj Z)
      ((coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).obj W)
    refine ⟨Under.homMk (⟨(⟨inv (ElemFrobCat.Hom.base Z.hom.hom)
        ≫ ElemFrobCat.Hom.base W.hom.hom,
        Φ.map (inv (ElemFrobCat.Hom.base Z.hom.hom)) (x₀ + u + V.neg), 1⟩ :
          Z.right.obj ⟶ W.right.obj),
      ⟨elemFrob_coAngular Φ hD hpd _, rfl, ?_⟩⟩ : Z.right ⟶ W.right) ?_,
      Subsingleton.elim _ _⟩
    · show IsIso (inv (ElemFrobCat.Hom.base Z.hom.hom) ≫ ElemFrobCat.Hom.base W.hom.hom)
      haveI : IsIso (ElemFrobCat.Hom.base W.hom.hom) := W.hom.property.2.2
      infer_instance
    · refine InducedWideCategory.Hom.ext ?_
      simp only [WideSubcategory.comp_def]
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · rw [ElemFrobCat.comp_base, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
      · rw [ElemFrobCat.div_comp]
        exact hkey
      · rw [ElemFrobCat.comp_deg, one_mul,
          show ElemFrobCat.Hom.deg Z.hom.hom = 1 from Z.hom.property.2.1,
          show ElemFrobCat.Hom.deg W.hom.hom = 1 from W.hom.property.2.1]
  haveI hess : (coaPreUnderFunctor (elemPreFrobenioid Φ hD hpd) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a, ha⟩ : ∃ a : Φ.val A.base, toChar a = c := by
      obtain ⟨y, hy⟩ := toChar_surjective _ c
      exact ⟨y, hy⟩
    refine ⟨Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
        ⟶ (⟨A⟩ : WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
      ⟨otriOf Φ A a, elemFrob_coAngular Φ hD hpd _, rfl, ?_⟩), ⟨eqToIso ?_⟩⟩
    · show IsIso (𝟙 A.base)
      infer_instance
    · exact congrArg toOrderCat ha
  exact ⟨hfaith, hfull, hess⟩

/-- `otriOf` の底は `𝟙` なので、その逆射での引き戻しは何もしない。

★`inv` の instance 引数は `IsIso` が `Prop` なので**証明無関係**であり、
標準の instance で述べた `IsIso.inv_id` がそのまま使える。 -/
theorem elemFrob_map_inv_otriOf {A : ElemFrobCat Φ} (a : Φ.val A.base)
    (hi : IsIso ((elemPreFrobenioid Φ hD hpd).Base (otriOf Φ A a))) :
    Φ.charOn.map (@inv _ _ _ _ _ hi)
        ((elemPreFrobenioid Φ hD hpd).Div (otriOf Φ A a)) = toChar a := by
  rw [show (@inv _ _ _ _ _ hi) = 𝟙 A.base from IsIso.inv_id]
  exact Φ.charOn.map_id A.base _

/-- **(iii)(d)** スライス側 `(𝔽_Φ^coa-pre)_A → Order(Φ^char(A))^opp` は圏同値。 -/
theorem elemFrob_coaPreOverEquiv (A : ElemFrobCat Φ) :
    (coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).IsEquivalence := by
  haveI hfaith : (coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have h2 : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : Mono W.hom.hom := elemFrob_preStepMono Φ hD hpd _ W.hom.property.2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_mono W.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).Full := by
    constructor
    intro Z W h
    haveI hz : IsIso (ElemFrobCat.Hom.base Z.hom.hom) := Z.hom.property.2.2
    haveI hw : IsIso (ElemFrobCat.Hom.base W.hom.hom) := W.hom.property.2.2
    obtain ⟨x, hx⟩ := leOfHom h.unop
    obtain ⟨x₀, hx₀⟩ : ∃ x₀ : Φ.val A.base, toChar x₀ = x := by
      obtain ⟨y, hy⟩ := toChar_surjective _ x
      exact ⟨y, hy⟩
    obtain ⟨u, v, ⟨U, hU⟩, ⟨V, hV⟩, huv⟩ :
        CharRel (Φ.map (inv (ElemFrobCat.Hom.base W.hom.hom))
            (ElemFrobCat.Hom.div W.hom.hom) + x₀)
          (Φ.map (inv (ElemFrobCat.Hom.base Z.hom.hom))
            (ElemFrobCat.Hom.div Z.hom.hom)) :=
      toChar_eq_iff.mp (by
        rw [map_add, hx₀]
        exact hx)
    haveI := Preorder.subsingleton_hom
      ((coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).obj W).unop
      ((coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).obj Z).unop
    refine ⟨Over.homMk (⟨(⟨ElemFrobCat.Hom.base Z.hom.hom
        ≫ inv (ElemFrobCat.Hom.base W.hom.hom),
        Φ.map (ElemFrobCat.Hom.base Z.hom.hom) (x₀ + u + V.neg), 1⟩ :
          Z.left.obj ⟶ W.left.obj),
      ⟨elemFrob_coAngular Φ hD hpd _, rfl, ?_⟩⟩ : Z.left ⟶ W.left) ?_,
      Quiver.Hom.unop_inj (Subsingleton.elim _ _)⟩
    · show IsIso (ElemFrobCat.Hom.base Z.hom.hom ≫ inv (ElemFrobCat.Hom.base W.hom.hom))
      infer_instance
    · refine InducedWideCategory.Hom.ext ?_
      simp only [WideSubcategory.comp_def]
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · rw [ElemFrobCat.comp_base, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
      · rw [ElemFrobCat.div_comp,
          show ElemFrobCat.Hom.deg W.hom.hom = 1 from W.hom.property.2.1]
        simp only [PNat.one_coe, one_smul]
        have hinner : Φ.map (inv (ElemFrobCat.Hom.base W.hom.hom))
              (ElemFrobCat.Hom.div W.hom.hom) + (x₀ + u + V.neg)
            = Φ.map (inv (ElemFrobCat.Hom.base Z.hom.hom))
              (ElemFrobCat.Hom.div Z.hom.hom) := by
          rw [← add_assoc, ← add_assoc, huv, add_assoc, ← hV, V.val_neg, add_zero]
        rw [MonoidOn.map_comp, ← map_add, hinner, ← MonoidOn.map_comp,
          IsIso.hom_inv_id, MonoidOn.map_id]
      · rw [ElemFrobCat.comp_deg, mul_one,
          show ElemFrobCat.Hom.deg W.hom.hom = 1 from W.hom.property.2.1,
          show ElemFrobCat.Hom.deg Z.hom.hom = 1 from Z.hom.property.2.1]
  haveI hess : (coaPreOverFunctor (elemPreFrobenioid Φ hD hpd) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a, ha⟩ : ∃ a : Φ.val A.base, toChar a = c.unop := by
      obtain ⟨y, hy⟩ := toChar_surjective _ c.unop
      exact ⟨y, hy⟩
    refine ⟨Over.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd)))
        ⟶ (⟨A⟩ : WideSubcategory (coaPreProp (elemPreFrobenioid Φ hD hpd))) from
      ⟨otriOf Φ A a, elemFrob_coAngular Φ hD hpd _, rfl, ?_⟩), ⟨eqToIso ?_⟩⟩
    · show IsIso (𝟙 A.base)
      infer_instance
    · refine Opposite.unop_injective ?_
      rw [← ha]
      exact elemFrob_map_inv_otriOf Φ hD hpd a _
  exact ⟨hfaith, hfull, hess⟩

end CoaPre

/-- ★★★**[FrdI] Proposition 1.5, (i) の Frobenioid 部分** ——
`𝔽_Φ` は関手 `𝔽_Φ → 𝔽_{Φ^char}` を備えて **Frobenioid** である。 -/
theorem elemFrob_isFrobenioid : Frobenioid (elemPreFrobenioid Φ hD hpd) := by
  haveI := coaPreProp_isMultiplicative (elemPreFrobenioid Φ hD hpd)
    (elemFrob_frobenioidCore Φ hD hpd).coAngularComp
  exact ⟨elemFrob_frobenioidCore Φ hD hpd,
    elemFrob_coaPreUnderEquiv Φ hD hpd, elemFrob_coaPreOverEquiv Φ hD hpd⟩

/-! ### ★出典の紐付け(`.src`)

★`.src` は「その原典項目を**完全に**実装した」という主張である。
`Proposition 1.5` は (i)(ii)(iii) がすべて揃ったので、ここで初めて付ける。

* **(i)**: `elemFrob_isFrobenioid`(Frobenioid 本体)＋ 7 つの type
  (`elemFrob_autAmple` / `elemFrob_autSubAmple` / `elemFrob_endAmple` /
  `elemFrob_baseTrivial` / `elemFrob_frobeniusTrivial` /
  `elemFrob_frobeniusNormalized` / `elemFrob_isotropic`)
* **(ii)**: `otriEquiv`(モノイド同型)＋ `otriOf_mem_otimes`(`𝒪^×` ↔ `Φ(A)^±`)
* **(iii)**: `elemFrob_perfect` / `elemFrob_groupLike`
-/

def elemFrob_isFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 27, item := "Proposition 1.5, (i)",
    sectionId := "frdi-prop-1-5-i" }

def otriEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 27, item := "Proposition 1.5, (ii)",
    sectionId := "frdi-prop-1-5-ii" }

def elemFrob_perfect.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 27, item := "Proposition 1.5, (iii)",
    sectionId := "frdi-prop-1-5-iii" }

end ABC3.Found.FrdI
