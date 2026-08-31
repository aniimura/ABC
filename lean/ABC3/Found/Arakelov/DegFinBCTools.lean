/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.BaseChangeIso
import ABC3.Found.Arakelov.CardQuotientBC
import Mathlib.LinearAlgebra.TensorProduct.RightExactness
import ABC3.Meta.Claim

/-!
# 底変換の**有限側**の道具 —— 係数環の取り替えと商の位数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L)

## ★★★★★★★★これは何か

底変換の有限側

    `degFin_K(f^*L̄) = [K:F] · degFin_F(L̄)`

を出すのに要る 2 つの橋を入れる:

| もの | 内容 |
|---|---|
| `card_quotient_span_gammaModPre` | ★`Γ(X,⊤)`-span と `R`-span は**同じ部分集合**である |
| `range_lTensor_span_singleton` | ★★`T ⊗ (R·a)` の像は `T·(1 ⊗ a)` である |
| `card_quotient_tensor_span` | ★★★`#((T ⊗ A)/T·(1⊗a)) = #(T ⊗ (A/R·a))` |

★★`degFinPre`（`§9-780`）は `Γ(X,⊤)`-span で定義されているが、
`§9-791` の同型と `§9-794` の数え上げは `R = 𝓞_F`／`T = 𝓞_K` の側にある。
一つ目の補題がその**摩擦を渡る**。

★★★三つ目はテンソルの右完全性（`lTensor.equiv`）である
——`(T ⊗ A)/T·(1⊗a) ≅ T ⊗ (A/R·a)`。

## ★配管の記録

部分加群の商と加法部分群の商は**定義的に等しい**（`§9-794` の測定）ので、
係数環が違っても**台が同じなら位数は等しい**——`congrArg` 一発で渡れる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open scoped TensorProduct

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`Γ(X,⊤)`-span と `R`-span は同じ部分集合である**。

★`Γ-Spec` 同型で係数が移るだけなので、生成する部分加群の台は変わらない。
★★部分加群の商と加法部分群の商は定義的に等しいので、位数もそのまま等しい。 -/
theorem card_quotient_span_gammaModPre (R : CommRingCat.{0}) (M : (Spec R).PresheafOfModules)
    (x : (gammaModPre R M : Type)) :
    Nat.card ((gammaModPre R M : Type) ⧸ ((Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) ∙ x))
      = Nat.card ((gammaModPre R M : Type) ⧸ ((R : Type) ∙ x)) := by
  have hsub : (((Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) ∙ x).toAddSubgroup)
      = (((R : Type) ∙ x).toAddSubgroup) := by
    ext y
    constructor
    · intro hy
      obtain ⟨c, rfl⟩ := Submodule.mem_span_singleton.mp hy
      refine Submodule.mem_span_singleton.mpr ⟨(Scheme.ΓSpecIso R).hom.hom c, ?_⟩
      rw [gammaModPre_smul]
      congr 1
      exact congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m c)
        (Scheme.ΓSpecIso R).hom_inv_id
    · intro hy
      obtain ⟨t, rfl⟩ := Submodule.mem_span_singleton.mp hy
      exact Submodule.mem_span_singleton.mpr ⟨(Scheme.ΓSpecIso R).inv.hom t, rfl⟩
  exact congrArg (fun (H : AddSubgroup (gammaModPre R M : Type)) =>
    Nat.card ((gammaModPre R M : Type) ⧸ H)) hsub

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`T ⊗ (R·a)` の像は `T·(1 ⊗ a)` である**。

★`t ⊗ (r·a) = (r·t) ⊗ a = (r·t)·(1 ⊗ a)` という一行が中身である。 -/
theorem range_lTensor_span_singleton (R T A : Type) [CommRing R] [CommRing T] [Algebra R T]
    [AddCommGroup A] [Module R A] (a : A) :
    ((LinearMap.range (LinearMap.lTensor T (Submodule.span R {a}).subtype)) : Set (T ⊗[R] A))
      = ((Submodule.span T {(1 : T) ⊗ₜ[R] a} : Submodule T (T ⊗[R] A)) : Set (T ⊗[R] A)) := by
  have hkey : ∀ (s : T), s • ((1 : T) ⊗ₜ[R] a) = s ⊗ₜ[R] a := by
    intro s
    rw [TensorProduct.smul_tmul', smul_eq_mul, mul_one]
  ext z
  constructor
  · rintro ⟨w, rfl⟩
    induction w using TensorProduct.induction_on with
    | zero => simp
    | tmul t y =>
        obtain ⟨r, hr⟩ := Submodule.mem_span_singleton.mp y.2
        refine Submodule.mem_span_singleton.mpr ⟨r • t, ?_⟩
        rw [hkey, LinearMap.lTensor_tmul]
        show (r • t) ⊗ₜ[R] a = t ⊗ₜ[R] ((y : A))
        rw [← hr, ← TensorProduct.smul_tmul]
    | add u v hu hv =>
        rw [map_add]
        exact Submodule.add_mem _ hu hv
  · intro hz
    obtain ⟨t, rfl⟩ := Submodule.mem_span_singleton.mp hz
    refine ⟨t ⊗ₜ[R] (⟨a, Submodule.mem_span_singleton_self a⟩ : Submodule.span R {a}), ?_⟩
    rw [hkey, LinearMap.lTensor_tmul]
    rfl

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**`#((T ⊗ A)/T·(1⊗a)) = #(T ⊗ (A/R·a))`**（`T` は平坦）。

★テンソルの右完全性（`lTensor.equiv`）と、上の像の同定の合成である。 -/
theorem card_quotient_tensor_span (R T A : Type) [CommRing R] [CommRing T] [Algebra R T]
    [Module.Flat R T] [AddCommGroup A] [Module R A] (a : A) :
    Nat.card ((T ⊗[R] A) ⧸ (Submodule.span T {(1 : T) ⊗ₜ[R] a}))
      = Nat.card (T ⊗[R] (A ⧸ (Submodule.span R {a}))) := by
  have hsub : (Submodule.span T {(1 : T) ⊗ₜ[R] a}).toAddSubgroup
      = (LinearMap.range (LinearMap.lTensor T (Submodule.span R {a}).subtype)).toAddSubgroup := by
    ext y
    exact (Set.ext_iff.mp (range_lTensor_span_singleton R T A a) y).symm
  have e := lTensor.equiv (f := (Submodule.span R {a}).subtype)
    (g := (Submodule.span R {a}).mkQ) T
    (LinearMap.exact_subtype_mkQ (Submodule.span R {a}))
    (Submodule.mkQ_surjective (Submodule.span R {a}))
  have h1 : Nat.card ((T ⊗[R] A) ⧸ (Submodule.span T {(1 : T) ⊗ₜ[R] a}))
      = Nat.card ((T ⊗[R] A) ⧸ (LinearMap.range (LinearMap.lTensor T
          (Submodule.span R {a}).subtype))) :=
    congrArg (fun (H : AddSubgroup (T ⊗[R] A)) => Nat.card ((T ⊗[R] A) ⧸ H)) hsub
  rw [h1, Nat.card_congr e.toEquiv]

/-! ### ★出典の紐付け(`.src`) -/

def card_quotient_span_gammaModPre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Γ(X,⊤)-span と R-span は同じ部分集合である)",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_tensor_span.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(#((T ⊗ A)/T·(1⊗a)) = #(T ⊗ (A/R·a)))",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_tensor_span.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "lTensor.equiv(テンソル積の右完全性)"
      (.inMathlib "lTensor.equiv") 4,
    .implicitStep
      ("★これらは底変換の有限側 degFin_K(f^*L̄) = [K:F]·degFin_F(L̄) の道具である。" ++
       "★★§9-791 の同型 Γ(f^*L) ≅ 𝓞_K ⊗ Γ(L) と §9-794 の数え上げを繋ぐのに使う") 4 ]

end ABC3.Found.Arakelov
