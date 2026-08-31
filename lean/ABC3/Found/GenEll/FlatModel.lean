/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ComapMul
import Mathlib.RingTheory.Flat.TorsionFree

/-!
# `ℚ`-スキームの閉包は **`ℤ`-平坦**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★★原文第 2 文の括弧の中身

`Remark 1.4.1` の第 2 文はこう書く:

> In particular, this theory may be applied to any normal, projective scheme `Y` over `ℚ`
> [i.e., by regarding `Y` as the "`X_ℚ`" determined by some **`ℤ`-flat, `ℤ`-projective model**
> "`X`" of `Y` **that arises from a projective embedding of `Y`**].

★★その「`ℤ`-平坦なモデル」は**`ℚ`-スキームの `ℤ`-スキームの中での閉包**である。
本ファイルはその**平坦性の側**を取る。

## ★★★★★★機構は 3 行

`X ≝` スキーム論的像 `= P / ker(𝒪_P → f_*𝒪_Y)` なので、アフィン開 `U` の上では

    Γ(X, U) = Γ(P, U) / ker(f.app U)  ↪  Γ(Y, f⁻¹U)

★右辺は `ℚ`-代数なので `ℤ` 上平坦（`ℚ` が `ℤ` 上平坦・体上の加群は自由）。
★★単射に埋め込めば**ねじれ無し**が遺伝し、`ℤ` は Bezout なので**平坦**に戻る
（`Module.Flat.flat_iff_torsion_eq_bot_of_isBezout`）。

★★★**吹き上げも `Proj` も要らない**——`Scheme.Hom.ker_apply`（mathlib）と
`Module.Flat.trans` だけである。

## ★残っている段（明示）

★`ℤ`-射影性（`ℙⁿ_ℤ` への閉埋め込み）は本ファイルには無い
——`ample-and-projective-embedding` の段 C（`ℙⁿ_S` と `O(1)`）である。
★★`X_ℚ ≅ Y`（生成ファイバーが元に戻ること）も本ファイルには無い。
-/

namespace ABC3.Found.GenEll

open Module AlgebraicGeometry CategoryTheory

/-! ## ★ねじれ無しは単射で遺伝する -/

/-- ★★**単射な線型写像に沿ってねじれ無しは遺伝する**。

★`IsTorsionFree R M` は「正則元の作用が単射」なので、
`f` が単射なら合成の単射性から出る。 -/
theorem isTorsionFree_of_injective {R M N : Type} [CommRing R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] (f : M →ₗ[R] N) (hf : Function.Injective f)
    [Module.IsTorsionFree R N] : Module.IsTorsionFree R M := by
  refine ⟨fun r hr => ?_⟩
  intro m₁ m₂ h
  refine hf ?_
  have h2 : r • f m₁ = r • f m₂ := by
    rw [← f.map_smul, ← f.map_smul]
    exact congrArg f h
  exact Module.IsTorsionFree.isSMulRegular hr h2

/-- ★★★**`ℤ` 上では平坦性が単射で遺伝する**。

★`ℤ` は Bezout 整域なので「平坦 ⟺ ねじれ無し」
（`Module.Flat.flat_iff_torsion_eq_bot_of_isBezout`）。
★★一般の環では成り立たない——**Bezout 性が本質**である。 -/
theorem flat_int_of_injective {M N : Type} [AddCommGroup M] [Module ℤ M]
    [AddCommGroup N] [Module ℤ N] (f : M →ₗ[ℤ] N) (hf : Function.Injective f)
    [Module.Flat ℤ N] : Module.Flat ℤ M := by
  haveI : Module.IsTorsionFree ℤ N := inferInstance
  haveI : Module.IsTorsionFree ℤ M := isTorsionFree_of_injective f hf
  rw [Module.Flat.flat_iff_torsion_eq_bot_of_isBezout,
    ← Submodule.isTorsionFree_iff_torsion_eq_bot]
  infer_instance

/-! ## ★★★★核による商の平坦性 -/

/-- ★★**行き先が `ℤ`-平坦なら、核による商も `ℤ`-平坦**。 -/
theorem flat_int_quotient_ker {A B : Type} [CommRing A] [CommRing B] (φ : A →+* B)
    [Module.Flat ℤ B] : Module.Flat ℤ (A ⧸ RingHom.ker φ) := by
  refine flat_int_of_injective
    (AddMonoidHom.toIntLinearMap (RingHom.kerLift φ).toAddMonoidHom) ?_
  exact RingHom.kerLift_injective φ

/-- ★★★★**`ℚ`-代数への環準同型の核による商は `ℤ`-平坦**。

★`ℚ` は `ℤ` 上平坦、体上の加群は自由（したがって平坦）、
`Module.Flat.trans` で合成する。 -/
theorem flat_int_quotient_ker_of_ratAlgebra {A B : Type} [CommRing A] [CommRing B]
    [Algebra ℚ B] (φ : A →+* B) : Module.Flat ℤ (A ⧸ RingHom.ker φ) := by
  haveI : Module.Flat ℚ B := Module.Flat.of_free
  haveI : Module.Flat ℤ B := Module.Flat.trans ℤ ℚ B
  exact flat_int_quotient_ker φ

/-- ★**イデアルが `ℚ`-代数への写像の核として書けるなら、商は `ℤ`-平坦**。 -/
theorem flat_int_quotient_of_ker_eq {A B : Type} [CommRing A] [CommRing B] [Algebra ℚ B]
    (φ : A →+* B) (I : Ideal A) (hI : I = RingHom.ker φ) : Module.Flat ℤ (A ⧸ I) := by
  subst hI
  exact flat_int_quotient_ker_of_ratAlgebra φ

/-! ## ★★★★★★★★★スキーム論的像へ -/

/-- ★★★★★★★★**スキーム論的像のイデアル商は `ℤ`-平坦**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

★機構は `Scheme.Hom.ker_apply`（mathlib、準コンパクトな射の核をアフィンで読む）
1 本である。 -/
theorem flat_int_ker_ideal {P Y : Scheme.{0}} (f : Y ⟶ P) [QuasiCompact f] (U : P.affineOpens)
    [Algebra ℚ Γ(Y, f ⁻¹ᵁ (U : P.Opens))] :
    Module.Flat ℤ (Γ(P, (U : P.Opens)) ⧸ (Scheme.Hom.ker f).ideal U) :=
  flat_int_quotient_of_ker_eq (Scheme.Hom.app f (U : P.Opens)).hom _
    (Scheme.Hom.ker_apply f U)

/-- ★★★★★★★★★**スキーム論的像そのものの切断が `ℤ`-平坦**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

> **`ℚ`-スキーム `Y` の `ℤ`-スキーム `P` の中での閉包は `ℤ`-平坦である。**

★★これが原文第 2 文の括弧「`ℤ`-flat, `ℤ`-projective model」の**平坦性の側**である。
★★★`Scheme.IdealSheafData.subschemeObjIso`（mathlib）で
`Γ(像, ι⁻¹U) ≅ Γ(P,U)/核` に直すだけ。 -/
theorem flat_int_subschemeObj {P Y : Scheme.{0}} (f : Y ⟶ P) [QuasiCompact f]
    (U : P.affineOpens) [Algebra ℚ Γ(Y, f ⁻¹ᵁ (U : P.Opens))] :
    Module.Flat ℤ Γ((Scheme.Hom.ker f).subscheme,
      (Scheme.Hom.ker f).subschemeι ⁻¹ᵁ (U : P.Opens)) := by
  haveI := flat_int_ker_ideal f U
  let e := Scheme.IdealSheafData.subschemeObjIso (Scheme.Hom.ker f) U
  have hinj : Function.Injective (e.hom.hom) :=
    ConcreteCategory.injective_of_mono_of_preservesPullback e.hom
  exact flat_int_of_injective (AddMonoidHom.toIntLinearMap (e.hom.hom).toAddMonoidHom) hinj

/-! ### ★出典の紐付け(`.src`)

★★`Remark 1.4.1` の第 2 文の括弧のうち**平坦性の側だけ**である。
射影性（`ℙⁿ_ℤ` への閉埋め込み）と `X_ℚ ≅ Y` は含まない。 -/

def flat_int_of_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ 上は平坦性が単射で遺伝すること——Bezout 性)",
    sectionId := "genell-rem-1-4-1" }

def flat_int_subschemeObj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℚ-スキームの閉包が ℤ-平坦であること。射影性は含まない)",
    sectionId := "genell-rem-1-4-1" }

def flat_int_subschemeObj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.ker_apply(準コンパクトな射の核をアフィンで読む)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.ker_apply") 8,
    .citation "[mathlib]"
      "Scheme.IdealSheafData.subschemeObjIso(閉部分スキームのアフィン記述)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.subschemeObjIso") 8,
    .citation "[mathlib]"
      "Module.Flat.flat_iff_torsion_eq_bot_of_isBezout(Bezout 整域上で平坦 ⟺ ねじれ無し)"
      (.inMathlib "Module.Flat.flat_iff_torsion_eq_bot_of_isBezout") 8,
    .implicitStep
      ("★★★★原文は「ℤ-flat, ℤ-projective model that arises from a projective " ++
       "embedding of Y」と括弧 1 つで済ませている。" ++
       "★本ファイルはそのうち**平坦性の側**を取った——" ++
       "吹き上げも Proj も要らず、ker_apply と Module.Flat.trans だけである") 8,
    .implicitStep
      ("★★残る段: ℤ-射影性(ℙⁿ_ℤ への閉埋め込み)は " ++
       "ample-and-projective-embedding の段 C(ℙⁿ_S と O(1))である。" ++
       "★X_ℚ ≅ Y(生成ファイバーが元に戻ること)も本ファイルには無い") 8 ]

end ABC3.Found.GenEll
