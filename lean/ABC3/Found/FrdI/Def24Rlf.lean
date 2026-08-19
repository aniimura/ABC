/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ElementaryFrobenioid
import Mathlib.LinearAlgebra.TensorProduct.Tower
import Mathlib.Data.NNReal.Defs

/-!
# [FrdI] 実化 `Φ^rlf` の**関手的な**模型

`Def24.lean` は `Definition 2.4, (i)` の実化 `M^rlf` を、原文どおり
**素点分解** `M^rlf ⊆ ∏_{p ∈ Prime(M)} ℝ≥0` として作った。
★そこには**素点ごとの同一視 `ι`(選択)**が入るので、`𝒟` 上の**関手**にするのが難しい。

★★`Theorem 5.2` は有理関数の単系 `B` を「`𝒟` 上の単系」として要求し、
`Proposition 5.3` の第 1 文は `𝒞^rlf` を `(Φ^rlf, ℝ·Φ^birat)` の model Frobenioid と
**定義する**。したがって `Φ^rlf` を `MonoidOn 𝒟` として立てる必要がある。

## ★設計 —— テンソルで関手性を無料にする

原文自身が `M^rlf` を「`M^pf ⊗ ℝ`」の一言で述べている。そこで

  `M^rlf := ℝ≥0 ⊗_ℕ M`

と置く。★`ℕ` 上のテンソルなので**任意の可換単系に対して定義でき**、
* `Module ℝ≥0 (ℝ≥0 ⊗_ℕ M)` が**自動で付く**(`Algebra ℕ ℝ≥0` からの基底変換)、
* 関手性 `rlfMap` が `TensorProduct.map` から**無料で出る**、
* 全射性の保存(`rlfMap_surjective`)も直ちに出る。

★`M^pf` も同じ形(`ℚ≥0 ⊗_ℕ M`)で書ける —— 本ファイルでは実化だけを扱う。

## ★★逸脱(記録)

★これは `Def24.lean` の `Rlf`(素点分解版)と**別の定義**である。
perf-factorial な `M` では `M^pf ≅ ⊕_p ℚ≥0` なので両者は一致するはずだが、
**その一致はまだ証明していない**(依存グラフの鎖 `rlf` の節点 `rlf-agree`)。
★したがって現時点では「実化の 2 つの模型」が並立している。

★もう 1 つ、`MonoidOn` の条件 (a)(characteristic injectivity)は
**テンソルの平坦性**にあたり、一般の可換単系では成り立たない。
本ファイルは**それを仮定として受け取る**形(`phiRlfOn` の `hcharInj`)にしてある
—— すなわち `Φ^rlf : MonoidOn 𝒟` は**ただ 1 つの言明**に還元されている
(鎖 `rlf` の節点 `rlf-flat`)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w

/-! ## ★1. `M^rlf = ℝ≥0 ⊗_ℕ M` -/

/-- ★★**実化の関手的模型** —— `M^rlf := ℝ≥0 ⊗_ℕ M`。 -/
abbrev RlfT (M : Type w) [AddCommMonoid M] : Type w := NNReal ⊗[ℕ] M

/-- ★`M → M^rlf`(`m ↦ 1 ⊗ m`)。 -/
noncomputable def toRlf {M : Type w} [AddCommMonoid M] : M →+ RlfT M where
  toFun m := (1 : NNReal) ⊗ₜ m
  map_zero' := by simp
  map_add' _ _ := by rw [TensorProduct.tmul_add]

/-- ★★**関手性** —— `TensorProduct.map` から無料で出る。 -/
noncomputable def rlfMap {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N) :
    RlfT M →+ RlfT N :=
  (TensorProduct.map (LinearMap.id (R := ℕ) (M := NNReal)) f.toNatLinearMap).toAddMonoidHom

@[simp] theorem rlfMap_tmul {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    (r : NNReal) (m : M) : rlfMap f (r ⊗ₜ m) = r ⊗ₜ f m := rfl

theorem rlfMap_id {M : Type w} [AddCommMonoid M] :
    rlfMap (AddMonoidHom.id M) = AddMonoidHom.id (RlfT M) := by
  ext x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul r m => rfl
  | add x y hx hy => simp [map_add, hx, hy]

theorem rlfMap_comp {M N O : Type w} [AddCommMonoid M] [AddCommMonoid N] [AddCommMonoid O]
    (f : M →+ N) (g : N →+ O) : rlfMap (g.comp f) = (rlfMap g).comp (rlfMap f) := by
  ext x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul r m => rfl
  | add x y hx hy => simp [map_add, hx, hy]

/-- ★`toRlf` の自然性。 -/
@[simp] theorem rlfMap_toRlf {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    (m : M) : rlfMap f (toRlf m) = toRlf (f m) := rfl

/-- ★★**全射は保たれる**(テンソルの右完全性のうち易しい方)。 -/
theorem rlfMap_surjective {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] {f : M →+ N}
    (hf : Function.Surjective f) : Function.Surjective (rlfMap f) := by
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

/-- ★`Φ^rlf` の台となる反変関手。 -/
noncomputable def rlfFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (RlfT (Φ.val X.unop))
  map f := AddCommMonCat.ofHom (rlfMap (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (𝟙 X.unop) = AddMonoidHom.id _ := by ext a; exact Φ.map_id _ a
    show rlfMap (Φ.map (𝟙 X.unop)) x = _
    rw [h, rlfMap_id]
    rfl
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (g.unop ≫ f.unop) = (Φ.map g.unop).comp (Φ.map f.unop) := by
      ext a; exact Φ.map_comp _ _ a
    show rlfMap (Φ.map (g.unop ≫ f.unop)) x = _
    rw [h, rlfMap_comp]
    rfl

/-- ★★★★**`Φ^rlf` を `𝒟` 上の単系として立てる**。

★`hcharInj`(条件 (a))だけを仮定に置く —— これが**テンソルの平坦性**にあたる
唯一の残りである(鎖 `rlf` の節点 `rlf-flat`)。
★条件 (b)(FSM 射なら同型)は `Φ.fsmIso` から**証明できる**:
単射性は `hcharInj` の第 1 成分、全射性は `rlfMap_surjective`。 -/
noncomputable def phiRlfOn (Φ : MonoidOn.{v, u, w} D)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (rlfMap (Φ.map α))) :
    MonoidOn.{v, u, w} D where
  functor := rlfFunctor Φ
  charInj α := hcharInj α
  fsmIso α hα := ⟨(hcharInj α).1, rlfMap_surjective (Φ.fsmIso α hα).2⟩

/-! ### ★出典の紐付け -/

/-- ★locator —— `Definition 2.4, (i)` の実化(★**関手的模型**の側。
`Def24.lean` の素点分解版とは別の模型であり、一致は未証明)。 -/
def RlfT.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48,
    item := "Definition 2.4, (i) — realification(関手的模型 ℝ≥0 ⊗_ℕ M)",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
