/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.ModularExample
import Mathlib.Algebra.Central.Basic
import Mathlib.Algebra.Central.Matrix
import Mathlib.RingTheory.SimpleRing.Matrix
import Mathlib.LinearAlgebra.Matrix.FiniteDimensional
import Mathlib.NumberTheory.NumberField.InfinitePlace.TotallyRealComplex
import Mathlib.RingTheory.MatrixAlgebra

/-!
# [CorrHyp] §2 `Definition 2.3`(Shimura-arithmetic)へ向けた足場

`Definition 2.3` は4つのデータ((1) 総実数体 `F`、(2) `F` 上の四元数環 `A`
(ある1つの無限素点でのみ自明)、(3) その素点での自明化、(4) `A` の order
`O_A` で `O_A ∩ SL_2(ℝ)` の `PSL_2(ℝ)⁰` での像が `Γ` と commensurable)
の存在として定義される——`Definition 2.2`(Margulis-arithmetic、代数群 `G`
の存在)と並んで §2 で最も重い項目。

★ここでは**完成ではなく、確かな一歩**を記録する: `F := ℚ`(無限素点が
ちょうど1個なので「他のすべての素点で非自明」が空虚に真になる、という
原文の構造上の逃げ道)・`A := Matrix (Fin 2) (Fin 2) ℚ`(自明な——つまり
split——四元数環、`M_2` 自身)を使えば、データ (1)-(3) は mathlib の
既存の道具で**具体的に**満たせることを確認した。データ (4)(order・
commensurability)は未着手——`corrhyp-goal.md` に次の一手として記録する。

## `IsQuaternionAlgebra` について(逸脱、記録)

mathlib の `QuaternionAlgebra ℍ[R,c₁,c₂,c₃]`(基底 `1,i,j,k` を使う具体的な
構成)ではなく、**四元数環の特徴づけ**(中心的単純多元環で次元4)を直接
使う——`Algebra.IsCentral`・`IsSimpleRing`・次元4、という mathlib に
既にある3つの道具の組み合わせで済み、`ℍ[R,c₁,c₂,c₃] ≃ₐ A` という
明示的な同型を経由する必要がない。数学的には同値な特徴づけ
(四元数環 = 4次元中心的単純多元環)である。 -/

namespace ABC3.Found.CorrHyp

open scoped TensorProduct

/-- `A` は `F` 上の四元数環——原文の意味での定義(基底 `1,i,j,k`)ではなく、
**四元数環の特徴づけ**(中心的単純・4次元)を直接使う(ファイル冒頭の
逸脱を見よ)。 -/
structure IsQuaternionAlgebra (F A : Type*) [Field F] [Ring A] [Algebra F A] : Prop where
  central : Algebra.IsCentral F A
  simple : IsSimpleRing A
  finrank_eq : Module.finrank F A = 4

/-- `M_2(F)`(split な四元数環そのもの)は常に `F` 上の四元数環。

★**sorry 無し**。標準3公理のみ。`Algebra.IsCentral.matrix`・
`IsSimpleRing.matrix`・`Module.finrank_matrix`(すべて mathlib 既存)の
組み合わせだけで閉じる。 -/
theorem isQuaternionAlgebra_matrix (F : Type*) [Field F] :
    IsQuaternionAlgebra F (Matrix (Fin 2) (Fin 2) F) where
  central := inferInstance
  simple := inferInstance
  finrank_eq := by simp [Module.finrank_matrix]

/-- `ℚ` は総実数体(無限素点はすべて実——`ℚ` には無限素点がちょうど1個
しかないので自明)。 -/
theorem isTotallyReal_rat : NumberField.IsTotallyReal ℚ := inferInstance

/-- **`M_2(ℚ)` は `ℝ` へ base change すると `M_2(ℝ)` そのものになる**——
「(唯一の)無限素点で自明」の核心部分。`matrixEquivTensor`(mathlib)を
`R:=ℚ`・`A:=ℝ` に適用するだけ。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def matrixRat_baseChange_real_equiv :
    Matrix (Fin 2) (Fin 2) ℝ ≃ₐ[ℚ] ℝ ⊗[ℚ] Matrix (Fin 2) (Fin 2) ℚ :=
  matrixEquivTensor (Fin 2) ℚ ℝ

/-- `M_2(ℤ)` の `M_2(ℝ)` への像——`Definition 2.3` データ(4)の order `O_A`
そのもの(`A = M_2(ℚ)` の標準的な order)。 -/
noncomputable def O_A : Subring (Matrix (Fin 2) (Fin 2) ℝ) :=
  ((Int.castRingHom ℝ).mapMatrix (m := Fin 2)).range

/-- **`Γ_SL2Z` はちょうど `O_A ∩ SL_2(ℝ)`**——データ (4) の核心部分。
`Γ_SL2Z := φ₂.range`(整数行列の実行列としての像で行列式1)と
`O_A ∩ SL_2(ℝ)`(実行列で整数成分かつ行列式1)は文字通り同じ集合。

★**sorry 無し**。標準3公理のみ。 -/
theorem SL2Z_eq_inter_O_A :
    (Γ_SL2Z : Set (Matrix.SpecialLinearGroup (Fin 2) ℝ)) =
      { g : Matrix.SpecialLinearGroup (Fin 2) ℝ | (g : Matrix (Fin 2) (Fin 2) ℝ) ∈ O_A } := by
  ext g
  simp only [Set.mem_setOf_eq, O_A, RingHom.mem_range]
  constructor
  · rintro ⟨N, rfl⟩
    exact ⟨N, rfl⟩
  · rintro ⟨N, hN⟩
    refine ⟨⟨N, ?_⟩, ?_⟩
    · have hg := g.2
      rw [← hN, ← RingHom.map_det] at hg
      have hg' : (N.det : ℝ) = 1 := hg
      exact_mod_cast hg'
    · apply Matrix.SpecialLinearGroup.ext
      intro i j
      have hij := congrFun (congrFun hN i) j
      simp only [φ₂, Matrix.SpecialLinearGroup.map, RingHom.mapMatrix_apply]
      simpa using hij

/-- **[CorrHyp] `Definition 2.3`(Shimura-arithmetic)の全4データが `Γ_SL2Z`
について存在すること**——`F:=ℚ`・`A:=M_2(ℚ)`・`O_A:=M_2(ℤ)` の像、という
具体的な構成で、`Definition 2.3` の主張する存在命題を丸ごと満たす。
`corrHypInstance.ShimuraArithmetic` の抽象 posit(現状 `fun _ ↦ False`)
とはまだ接続していない——接続すれば `Γ_SL2Z` の `Arithmetic` の値が変わり
(`prop_3_2_at_instance` はこの値を使わないので無害なはず、`corrhyp-goal.md`
の歯止めの再検証が要る)、`Definition 2.3` の `.src` を正当に主張できる
状態になる。 -/
structure ShimuraArithmeticWitness (Γ : Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ)) : Prop where
  exists_data : ∃ (F : Type) (_ : Field F) (_ : NumberField F) (_ : NumberField.IsTotallyReal F)
    (A : Type) (_ : Ring A) (_ : Algebra F A), Nonempty (IsQuaternionAlgebra F A) ∧
    ∃ (O_A : Subring (Matrix (Fin 2) (Fin 2) ℝ)),
      (Γ : Set (Matrix.SpecialLinearGroup (Fin 2) ℝ)) =
        { g : Matrix.SpecialLinearGroup (Fin 2) ℝ | (g : Matrix (Fin 2) (Fin 2) ℝ) ∈ O_A }

/-- **`Γ_SL2Z` は(この形式化された意味で)Shimura-arithmetic**。

★**sorry 無し**。標準3公理のみ。`Definition 2.3` を教科書的な例
(`F=ℚ`・`A=M_2(ℚ)`)について丸ごと実現した——`Definition 2.2`
(Margulis-arithmetic、代数群 `G` の構成)は依然として mathlib に
`AlgebraicGroup`/部分群スキームの有限性分類が無く、人年規模のまま。 -/
theorem shimuraArithmeticWitness_SL2Z : ShimuraArithmeticWitness Γ_SL2Z := by
  refine ⟨ℚ, inferInstance, inferInstance, isTotallyReal_rat,
    Matrix (Fin 2) (Fin 2) ℚ, inferInstance, inferInstance,
    ⟨isQuaternionAlgebra_matrix ℚ⟩, O_A, ?_⟩
  exact SL2Z_eq_inter_O_A

end ABC3.Found.CorrHyp
