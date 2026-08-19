/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ElementaryFrobenioid
import Mathlib.LinearAlgebra.TensorProduct.Tower
import Mathlib.Data.NNReal.Defs
import Mathlib.Data.NNRat.Defs

/-!
# [FrdI] 係数拡大 `S ⊗_ℕ Φ` —— 実化 `Φ^rlf` と完全化 `Φ^pf` の**関手的**模型

`Def24.lean` は `Definition 2.4, (i)` の実化 `M^rlf` を、原文どおり
**素点分解** `M^rlf ⊆ ∏_{p ∈ Prime(M)} ℝ≥0` として作った。
★そこには**素点ごとの同一視 `ι`(選択)**が入るので、`𝒟` 上の**関手**にするのが難しい。

★★`Theorem 5.2` は因子の単系・有理関数の単系をどちらも「`𝒟` 上の単系」として要求し、
`Proposition 5.3` は `𝒞^rlf` を `(Φ^rlf, ℝ·Φ^birat)` の、
`(𝒞^un-tr)^pf` を `(Φ^pf, ℚ·Φ^birat)` の model Frobenioid と**定義する**。
したがって `Φ^rlf` と `Φ^pf` を `MonoidOn 𝒟` として立てる必要がある。

## ★設計 —— `ℕ` 上のテンソルで関手性を無料にする

原文自身が `M^rlf` を「`M^pf ⊗ ℝ`」の一言で述べている。そこで**係数を引数に出して**

  `M ⊗_S := S ⊗_ℕ M`   (`S` は可換半環)

と置く。★`ℕ` 上のテンソルなので**任意の可換単系に対して定義でき**、
* `Module S (S ⊗_ℕ M)` が**自動で付く**(`Algebra ℕ S` からの基底変換)、
* 関手性 `scMap` が `TensorProduct.map` から**無料で出る**、
* 全射性の保存(`scMap_surjective`)も直ちに出る。

★`S := ℝ≥0` が実化(`RlfT`)、`S := ℚ≥0` が完全化(`PfT`)である。

## ★★逸脱(記録)

★これは `Def24.lean` の `Rlf`(素点分解版)と**別の定義**である。
perf-factorial な `M` では `M^pf ≅ ⊕_p ℚ≥0` なので両者は一致するはずだが、
**その一致はまだ証明していない**(依存グラフの鎖 `rlf` の節点 `rlf-agree`)。

★もう 1 つ、`MonoidOn` の条件 (a)(characteristic injectivity)は
**テンソルの平坦性**にあたり、一般の可換単系では成り立たない。
本ファイルは**それを仮定として受け取る**形(`phiScOn` の `hcharInj`)にしてある
—— すなわち `Φ^rlf : MonoidOn 𝒟` は**ただ 1 つの言明**に還元されている
(鎖 `rlf` の節点 `rlf-flat`)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w

/-! ## ★1. `M ⊗_S = S ⊗_ℕ M` -/

variable (S : Type) [CommSemiring S]

/-- ★★**係数拡大の関手的模型** —— `S ⊗_ℕ M`。 -/
abbrev ScT (M : Type w) [AddCommMonoid M] : Type w := S ⊗[ℕ] M

/-- ★★**実化** `M^rlf = ℝ≥0 ⊗_ℕ M`。 -/
abbrev RlfT (M : Type w) [AddCommMonoid M] : Type w := ScT NNReal M

/-- ★★**完全化** `M^pf = ℚ≥0 ⊗_ℕ M`。 -/
abbrev PfT (M : Type w) [AddCommMonoid M] : Type w := ScT NNRat M

variable {S}

/-- ★`M → M ⊗_S`(`m ↦ 1 ⊗ m`)。 -/
noncomputable def toSc {M : Type w} [AddCommMonoid M] : M →+ ScT S M where
  toFun m := (1 : S) ⊗ₜ m
  map_zero' := by simp
  map_add' _ _ := by rw [TensorProduct.tmul_add]

/-- ★★**関手性** —— `TensorProduct.map` から無料で出る。 -/
noncomputable def scMap {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N) :
    ScT S M →+ ScT S N :=
  (TensorProduct.map (LinearMap.id (R := ℕ) (M := S)) f.toNatLinearMap).toAddMonoidHom

@[simp] theorem scMap_tmul {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    (r : S) (m : M) : scMap f (r ⊗ₜ m) = r ⊗ₜ f m := rfl

theorem scMap_id {M : Type w} [AddCommMonoid M] :
    scMap (S := S) (AddMonoidHom.id M) = AddMonoidHom.id (ScT S M) := by
  ext x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul r m => rfl
  | add x y hx hy => simp [map_add, hx, hy]

theorem scMap_comp {M N O : Type w} [AddCommMonoid M] [AddCommMonoid N] [AddCommMonoid O]
    (f : M →+ N) (g : N →+ O) :
    scMap (S := S) (g.comp f) = (scMap g).comp (scMap f) := by
  ext x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul r m => rfl
  | add x y hx hy => simp [map_add, hx, hy]

/-- ★`toSc` の自然性。 -/
@[simp] theorem scMap_toSc {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    (m : M) : scMap (S := S) f (toSc m) = toSc (f m) := rfl

/-- ★★**全射は保たれる**(テンソルの右完全性のうち易しい方)。 -/
theorem scMap_surjective {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] {f : M →+ N}
    (hf : Function.Surjective f) : Function.Surjective (scMap (S := S) f) := by
  intro y
  induction y using TensorProduct.induction_on with
  | zero => exact ⟨0, map_zero _⟩
  | tmul r n =>
      obtain ⟨m, rfl⟩ := hf n
      exact ⟨r ⊗ₜ m, rfl⟩
  | add x y hx hy =>
      obtain ⟨a, rfl⟩ := hx
      obtain ⟨b, rfl⟩ := hy
      exact ⟨a + b, map_add _ _ _⟩

/-! ## ★2. `𝒟` 上の反変関手 -/

variable {D : Type u} [Category.{v} D]

variable (S) in
/-- ★`Φ ⊗_S` の台となる反変関手。 -/
noncomputable def scFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (ScT S (Φ.val X.unop))
  map f := AddCommMonCat.ofHom (scMap (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (𝟙 X.unop) = AddMonoidHom.id _ := by ext a; exact Φ.map_id _ a
    show scMap (Φ.map (𝟙 X.unop)) x = _
    rw [h, scMap_id]
    rfl
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (g.unop ≫ f.unop) = (Φ.map g.unop).comp (Φ.map f.unop) := by
      ext a; exact Φ.map_comp _ _ a
    show scMap (Φ.map (g.unop ≫ f.unop)) x = _
    rw [h, scMap_comp]
    rfl

variable (S) in
/-- ★★★★**`Φ ⊗_S` を `𝒟` 上の単系として立てる**。

★`hcharInj`(条件 (a))だけを仮定に置く —— これが**テンソルの平坦性**にあたる
唯一の残りである(鎖 `rlf` の節点 `rlf-flat`)。
★条件 (b)(FSM 射なら同型)は `Φ.fsmIso` から**証明できる**:
単射性は `hcharInj` の第 1 成分、全射性は `scMap_surjective`。 -/
noncomputable def phiScOn (Φ : MonoidOn.{v, u, w} D)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α))) :
    MonoidOn.{v, u, w} D where
  functor := scFunctor S Φ
  charInj α := hcharInj α
  fsmIso α hα := ⟨(hcharInj α).1, scMap_surjective (Φ.fsmIso α hα).2⟩

/-! ### ★出典の紐付け -/

/-- ★locator —— `Definition 2.4, (i)` の実化(★**関手的模型**の側。
`Def24.lean` の素点分解版とは別の模型であり、一致は未証明)。 -/
def RlfT.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48,
    item := "Definition 2.4, (i) — realification(関手的模型 ℝ≥0 ⊗_ℕ M)",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
