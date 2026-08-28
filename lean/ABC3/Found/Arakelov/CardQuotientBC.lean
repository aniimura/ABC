/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.Ideal.Norm.RelNorm
import Mathlib.NumberTheory.NumberField.Basic
import ABC3.Meta.Claim

/-!
# 底変換の**数え上げ**（巡回加群の場合）——`#(𝓞_K/𝔞𝓞_K) = #(𝓞_F/𝔞)^{[K:F]}`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★これは何か

底変換 `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の**有限側**は、`§9-791` の同型

    `Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L)`

と合わせると

    `#(𝓞_K ⊗_{𝓞_F} Q) = #(Q)^{[K:F]}`   （`Q` は有限 `𝓞_F`-加群）

に帰着する。★本ファイルはその**巡回加群の場合**（`Q = 𝓞_F/𝔞`）を入れる。

## ★★★測定の記録（2026-08-28）

mathlib に **`Ideal.absNorm_algebraMap`** がある:

    `absNorm (I.map (algebraMap R S)) = (absNorm I) ^ finrank (FractionRing R) (FractionRing S)`

★`Ideal.absNorm I = Nat.card (R ⧸ I)` は定義そのものなので、そのまま数え上げになる。
★★使うには `attribute [local instance] FractionRing.liftAlgebra` が要る
——mathlib の当該節がそうしている。

## ★残っている段（明示）

★★`finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = [K:F]` の同一視。
★★★一般の有限 `Q` への持ち上げ（巡回商による濾過での dévissage）。
-/

namespace ABC3.Found.Arakelov

open NumberField

attribute [local instance] FractionRing.liftAlgebra in
/-- ★★★★★★**イデアルを延長したときの商の位数**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `#(𝓞_K / 𝔞·𝓞_K) = #(𝓞_F / 𝔞) ^ [FractionRing 𝓞_K : FractionRing 𝓞_F]`

★mathlib の `Ideal.absNorm_algebraMap` そのままである
（`Ideal.absNorm I = Nat.card (R ⧸ I)` は定義）。 -/
theorem card_quotient_map_algebraMap (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (a : Ideal (𝓞 F)) :
    Nat.card ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K))))
      = (Nat.card ((𝓞 F) ⧸ a))
        ^ (Module.finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K))) :=
  Ideal.absNorm_algebraMap (R := 𝓞 F) (S := 𝓞 K) a

attribute [local instance] FractionRing.liftAlgebra FractionRing.isScalarTower_liftAlgebra in
/-- ★★★★★**`FractionRing` の次数は体の次数である**。

★`FractionRing (𝓞 F) ≃ₐ F`（`FractionRing.algEquiv`）を 2 つ使って
`Algebra.finrank_eq_of_equiv_equiv` に渡す。
★★互換性は局所化の一意性（`IsLocalization.ringHom_ext`）で出る。 -/
theorem finrank_fractionRing_ringOfIntegers (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] :
    Module.finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = Module.finrank F K := by
  refine Algebra.finrank_eq_of_equiv_equiv
    (FractionRing.algEquiv (𝓞 F) F).toRingEquiv
    (FractionRing.algEquiv (𝓞 K) K).toRingEquiv ?_
  apply IsLocalization.ringHom_ext (nonZeroDivisors (𝓞 F))
  ext x
  simp only [RingHom.coe_comp, Function.comp_apply, RingEquiv.toRingHom_eq_coe,
    RingEquiv.coe_toRingHom, AlgEquiv.coe_ringEquiv, AlgEquiv.commutes]
  rw [← IsScalarTower.algebraMap_apply (𝓞 F) (FractionRing (𝓞 F)) (FractionRing (𝓞 K)),
    IsScalarTower.algebraMap_apply (𝓞 F) (𝓞 K) (FractionRing (𝓞 K)),
    AlgEquiv.commutes,
    ← IsScalarTower.algebraMap_apply (𝓞 F) F K,
    IsScalarTower.algebraMap_apply (𝓞 F) (𝓞 K) K]

attribute [local instance] FractionRing.liftAlgebra in
/-- ★★★★★★★**イデアルを延長したときの商の位数**（`[K:F]` で書いた形）。

    `#(𝓞_K / 𝓞·𝓞_K) = #(𝓞_F / 𝓞) ^ [K:F]` -/
theorem card_quotient_map_algebraMap' (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (a : Ideal (𝓞 F)) :
    Nat.card ((𝓞 K) ⧸ (a.map (algebraMap (𝓞 F) (𝓞 K))))
      = (Nat.card ((𝓞 F) ⧸ a)) ^ (Module.finrank F K) := by
  rw [card_quotient_map_algebraMap F K a, finrank_fractionRing_ringOfIntegers F K]

/-! ### ★出典の紐付け(`.src`) -/

def card_quotient_map_algebraMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(イデアルを延長したときの商の位数——底変換の数え上げの巡回の場合)",
    sectionId := "genell-def-1-1-ii" }

def card_quotient_map_algebraMap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.absNorm_algebraMap(absNorm (I.map) = absNorm I ^ finrank)"
      (.inMathlib "Ideal.absNorm_algebraMap") 4,
    .implicitStep
      ("★残っている段: (a) finrank (FractionRing (𝓞 F)) (FractionRing (𝓞 K)) = [K:F] の同一視、" ++
       "(b) 一般の有限 Q への持ち上げ(巡回商による濾過での dévissage)。" ++
       "★★これと §9-791 の同型 Γ(f^*L) ≅ 𝓞_K ⊗_{𝓞_F} Γ(L) を合わせると" ++
       "底変換の有限側が出る") 4 ]

end ABC3.Found.Arakelov
