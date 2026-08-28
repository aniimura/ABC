/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.APicToSheaf
import Mathlib.RingTheory.Flat.TorsionFree
import Mathlib.RingTheory.PicardGroup
import Mathlib.LinearAlgebra.TensorProduct.RightExactness
import ABC3.Meta.Claim

/-!
# 段 E の核 —— **指数はテンソル積で掛け算になる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

`deg_F` の**加法性**（台帳 `arakelov-degF-finite-places` の段 E）の、
有限部分側の中身は次の**純粋な可換環論**である:

    `#( (A ⊗ B) / R·(a ⊗ b) ) = #( A / R·a ) · #( B / R·b )`

（`A`, `B` は可逆 `R`-加群、`R` は整域、`a ≠ 0`。）

★アルキメデス側の加法性は在庫の `archDeg_mul` である。

## ★★★★★機構（4 段）

    (1) card_tensor_invertible : `#(Q ⊗ B) = #(Q)`（`Q` 有限、`B` 可逆）
    (2) mk_injective_of_ne_zero: `y ↦ a ⊗ y` は単射
    (3) range_rTensor_span     : その像は `(R·a) ⊗ B` の像に一致する
    (4) card_quotient_tmul     : 塔の公式で組む

★★★(1) が鍵である。証明は次のとおり:
`Q` は有限なので `S := R / Ann(Q)` は**有限環**である（`S` は `Q` に忠実に作用する）。
有限可換環は Artin 環なので極大スペクトルが有限、
したがって `Pic S` は自明（mathlib の `Finite (MaximalSpectrum S) → Subsingleton (Pic S)`）。
★★よって `S ⊗_R B` は `S` 上**自由階数 1**、すなわち `≅ S` である。
`Q ⊗_R B ≅ Q ⊗_S (S ⊗_R B) ≅ Q ⊗_S S ≅ Q`。

★(4) は `A ⊗ B` の中で `R·(a⊗b) ≤ range(y ↦ a⊗y) ≤ A ⊗ B` という塔を作り、
在庫の `card_quotient_mul_card_map` / `card_map_mkQ_eq` /
`comap_subtype_span_singleton` に渡す。
★★上の段は `rTensor.equiv`（テンソルの右完全性）で `(A/R·a) ⊗ B` になり、(1) で `A/R·a` になる。
★★★下の段は (2) の単射性から `B / R·b` になる。
-/

namespace ABC3.Found.Arakelov

open scoped TensorProduct

/-! ## ★★★★★★(1) 有限加群と可逆加群のテンソルは位数を変えない -/

/-- ★**有限加群の零化イデアルによる商は有限環である**。

★`R / Ann(Q)` は `Q` に忠実に作用するので `Q → Q` の中に単射で入る。 -/
theorem finite_quotient_annihilator (R : Type) [CommRing R] (Q : Type) [AddCommGroup Q]
    [Module R Q] [Finite Q] : Finite (R ⧸ Module.annihilator R Q) := by
  letI := Module.quotientAnnihilator (R := R) (M := Q)
  have hinj : Function.Injective
      (fun (a : R ⧸ Module.annihilator R Q) => (fun q : Q => a • q)) := by
    intro a b hab
    induction a using Quotient.ind with
    | _ a =>
      induction b using Quotient.ind with
      | _ b =>
        have h : ∀ q : Q, a • q = b • q := fun q => congrFun hab q
        have hm : (a - b) ∈ Module.annihilator R Q := by
          rw [Module.mem_annihilator]
          intro q
          rw [sub_smul, h q, sub_self]
        exact (Submodule.Quotient.eq _).mpr hm
  exact Finite.of_injective _ hinj

/-- ★★★★★★★★**`#(Q ⊗ B) = #(Q)`**（`Q` 有限、`B` 可逆）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★機構: `S := R / Ann(Q)` は有限環、ゆえに Artin 環、ゆえに極大スペクトルが有限、
ゆえに `Pic S` が自明。★★★したがって可逆 `S`-加群 `S ⊗_R B` は自由階数 1、
すなわち `≅ S` である。あとは基底変換の打ち消し
（`AlgebraTensorModule.cancelBaseChange`）で `Q ⊗_R B ≅ Q` となる。

★これは「可逆加群でテンソルしても長さも位数も変わらない」という古典的事実の、
**有限環に落として mathlib の `Pic` の自明性を使う**形の証明である。 -/
theorem card_tensor_invertible (R : Type) [CommRing R] (Q : Type) [AddCommGroup Q] [Module R Q]
    [Finite Q] (B : Type) [AddCommGroup B] [Module R B] [Module.Invertible R B] :
    Nat.card (Q ⊗[R] B) = Nat.card Q := by
  haveI : Finite (R ⧸ Module.annihilator R Q) := finite_quotient_annihilator R Q
  letI := Module.quotientAnnihilator (R := R) (M := Q)
  haveI : IsScalarTower R (R ⧸ Module.annihilator R Q) Q :=
    (Module.isTorsionBySet_annihilator R Q).isScalarTower
  obtain ⟨e⟩ := Module.Invertible.free_iff_linearEquiv.mp
    (inferInstance : Module.Free (R ⧸ Module.annihilator R Q)
      ((R ⧸ Module.annihilator R Q) ⊗[R] B))
  have e2 : Q ⊗[R ⧸ Module.annihilator R Q] ((R ⧸ Module.annihilator R Q) ⊗[R] B)
      ≃ₗ[R ⧸ Module.annihilator R Q] Q ⊗[R] B :=
    TensorProduct.AlgebraTensorModule.cancelBaseChange R (R ⧸ Module.annihilator R Q)
      (R ⧸ Module.annihilator R Q) Q B
  have e3 : Q ⊗[R ⧸ Module.annihilator R Q] ((R ⧸ Module.annihilator R Q) ⊗[R] B)
      ≃ₗ[R ⧸ Module.annihilator R Q]
        Q ⊗[R ⧸ Module.annihilator R Q] (R ⧸ Module.annihilator R Q) :=
    TensorProduct.congr (LinearEquiv.refl _ Q) e
  have e4 : Q ⊗[R ⧸ Module.annihilator R Q] (R ⧸ Module.annihilator R Q)
      ≃ₗ[R ⧸ Module.annihilator R Q] Q :=
    TensorProduct.rid _ Q
  have hn : Nonempty ((Q ⊗[R] B) ≃ Q) := ⟨(e2.symm.trans (e3.trans e4)).toEquiv⟩
  exact Nat.card_congr hn.some

/-! ## ★★★★★★★(2)(3) `y ↦ a ⊗ y` の像 -/

/-- ★★★★**`y ↦ a ⊗ y` は単射である**（`a ≠ 0`、`A`・`B` 可逆、`R` 整域）。

★可逆加群は平坦、平坦加群は捻れを持たないので `1 ↦ a` は単射。
★★`B` が平坦なので `rTensor B` がその単射性を保つ。 -/
theorem mk_injective_of_ne_zero (R : Type) [CommRing R] [IsDomain R] (A B : Type)
    [AddCommGroup A] [Module R A] [Module.Invertible R A]
    [AddCommGroup B] [Module R B] [Module.Invertible R B]
    (a : A) (ha : a ≠ 0) : Function.Injective (TensorProduct.mk R A B a) := by
  have h1 : Function.Injective (LinearMap.toSpanSingleton R A a) := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro r hr
    by_contra hne
    have hreg := Module.Flat.isSMulRegular_of_isRegular (R := R) (M := A)
      (IsRegular.of_ne_zero hne)
    exact ha (hreg (by simpa [LinearMap.toSpanSingleton] using hr))
  have h2 := Module.Flat.rTensor_preserves_injective_linearMap
    (M := B) (LinearMap.toSpanSingleton R A a) h1
  have h3 : (TensorProduct.mk R A B a)
      = (LinearMap.rTensor B (LinearMap.toSpanSingleton R A a)) ∘ₗ
        (TensorProduct.lid R B).symm.toLinearMap := by
    ext y; simp [LinearMap.toSpanSingleton]
  rw [h3]
  exact h2.comp (TensorProduct.lid R B).symm.injective

/-- ★★★**`(R·a) ⊗ B` の像は `{a ⊗ y}` の張る部分加群である**。

★`(r·a) ⊗ y = a ⊗ (r·y)` という一行が中身である。 -/
theorem range_rTensor_span (R : Type) [CommRing R] (A B : Type)
    [AddCommGroup A] [Module R A] [AddCommGroup B] [Module R B] (a : A) :
    LinearMap.range (LinearMap.rTensor B (Submodule.span R {a}).subtype)
      = LinearMap.range (TensorProduct.mk R A B a) := by
  apply le_antisymm
  · rintro z ⟨w, rfl⟩
    induction w using TensorProduct.induction_on with
    | zero => simp
    | tmul x y =>
        obtain ⟨r, hr⟩ := Submodule.mem_span_singleton.mp x.2
        refine ⟨r • y, ?_⟩
        show a ⊗ₜ[R] (r • y) = (LinearMap.rTensor B (Submodule.span R {a}).subtype) (x ⊗ₜ y)
        rw [LinearMap.rTensor_tmul]
        show a ⊗ₜ[R] (r • y) = ((x : A)) ⊗ₜ[R] y
        rw [← hr, TensorProduct.smul_tmul]
    | add w1 w2 h1 h2 =>
        rw [map_add]
        exact Submodule.add_mem _ h1 h2
  · rintro z ⟨y, rfl⟩
    exact ⟨(⟨a, Submodule.mem_span_singleton_self a⟩ : Submodule.span R {a}) ⊗ₜ[R] y, rfl⟩

/-! ## ★★★★★★★★★★(4) 指数の乗法性 -/

/-- ★★★★★★★★★★**指数はテンソル積で掛け算になる**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `#( (A ⊗ B) / R·(a ⊗ b) ) = #( A / R·a ) · #( B / R·b )`

★★これが台帳 `arakelov-degF-finite-places` の**段 E の有限部分側**である。
★★★塔は `R·(a⊗b) ≤ range(y ↦ a⊗y) ≤ A ⊗ B`。
上の段はテンソルの右完全性（`rTensor.equiv`）で `(A/R·a) ⊗ B` になり、
`card_tensor_invertible` で `A/R·a` に落ちる。
下の段は `y ↦ a⊗y` の単射性で `B/R·b` になる。 -/
theorem card_quotient_tmul (R : Type) [CommRing R] [IsDomain R] (A B : Type)
    [AddCommGroup A] [Module R A] [Module.Invertible R A]
    [AddCommGroup B] [Module R B] [Module.Invertible R B]
    (a : A) (ha : a ≠ 0) (b : B) [Finite (A ⧸ (R ∙ a))] :
    Nat.card ((A ⊗[R] B) ⧸ (R ∙ (a ⊗ₜ[R] b)))
      = Nat.card (A ⧸ (R ∙ a)) * Nat.card (B ⧸ (R ∙ b)) := by
  have hinj := mk_injective_of_ne_zero R A B a ha
  have hmem : (a ⊗ₜ[R] b : A ⊗[R] B) ∈ LinearMap.range (TensorProduct.mk R A B a) := ⟨b, rfl⟩
  have hle : (R ∙ (a ⊗ₜ[R] b)) ≤ LinearMap.range (TensorProduct.mk R A B a) :=
    (Submodule.span_singleton_le_iff_mem _ _).mpr hmem
  rw [card_quotient_mul_card_map R (A ⊗[R] B) _ _ hle, card_map_mkQ_eq,
    comap_subtype_span_singleton R (A ⊗[R] B) _ (a ⊗ₜ[R] b) hmem]
  congr 1
  · have h1 : ((A ⊗[R] B) ⧸ LinearMap.range (TensorProduct.mk R A B a))
        ≃ₗ[R] ((A ⊗[R] B) ⧸ LinearMap.range
          (LinearMap.rTensor B (Submodule.span R {a}).subtype)) :=
      Submodule.quotEquivOfEq _ _ (range_rTensor_span R A B a).symm
    have h2 := rTensor.equiv (f := (Submodule.span R {a}).subtype)
      (g := (Submodule.span R {a}).mkQ) B
      (LinearMap.exact_subtype_mkQ (Submodule.span R {a}))
      (Submodule.mkQ_surjective (Submodule.span R {a}))
    have h3 : Nat.card ((A ⊗[R] B) ⧸ LinearMap.range (TensorProduct.mk R A B a))
        = Nat.card ((A ⧸ (R ∙ a)) ⊗[R] B) :=
      Nat.card_congr ((h1.trans h2).toEquiv)
    rw [h3, card_tensor_invertible R (A ⧸ (R ∙ a)) B]
  · have hmap : Submodule.map
        ((LinearEquiv.ofInjective (TensorProduct.mk R A B a) hinj) : B →ₗ[R] _) (R ∙ b)
        = R ∙ (⟨a ⊗ₜ[R] b, hmem⟩ : LinearMap.range (TensorProduct.mk R A B a)) := by
      rw [Submodule.map_span]
      simp only [Set.image_singleton, LinearEquiv.coe_coe]
      congr 1
    exact (Nat.card_congr (Submodule.Quotient.equiv (R ∙ b) _
      (LinearEquiv.ofInjective (TensorProduct.mk R A B a) hinj) hmap).toEquiv).symm

/-! ### ★出典の紐付け(`.src`) -/

def card_tensor_invertible.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(有限加群と可逆加群のテンソルは位数を変えない)",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_tmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 E の有限部分側——指数はテンソル積で掛け算になる)",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_tmul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Finite (MaximalSpectrum S) → Subsingleton (Pic S)(半局所環の Pic は自明)"
      (.inMathlib "CommRing.Pic") 4,
    .citation "[mathlib]" "Module.Invertible.free_iff_linearEquiv(可逆かつ自由なら R に同型)"
      (.inMathlib "Module.Invertible.free_iff_linearEquiv") 4,
    .citation "[mathlib]" "rTensor.equiv(テンソル積の右完全性)"
      (.inMathlib "rTensor.equiv") 4,
    .citation "[mathlib]" "Module.Flat.isSMulRegular_of_isRegular(平坦加群は捻れを持たない)"
      (.inMathlib "Module.Flat.isSMulRegular_of_isRegular") 4,
    .citation "[ABC3]" "card_quotient_mul_card_map / card_map_mkQ_eq(塔の公式、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.card_quotient_mul_card_map") 4,
    .implicitStep
      ("★原文は deg_F が準同型であることを [Szp] Prop 1.1 の同型を引いて済ませている。" ++
       "★★本ブロックはその同型を経由せず、古典的な定義から直接示す道の" ++
       "**有限部分側**である。アルキメデス側は在庫の archDeg_mul") 4 ]

end ABC3.Found.Arakelov
