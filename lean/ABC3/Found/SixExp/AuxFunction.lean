/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.RingTheory.IntegralClosure.IsIntegralClosure.Basic

/-!
# 補助関数の構成 —— チェーン `sixexp` の葉 `aux-function`

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーンの葉。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を作るか

超越性証明の補助関数

  ★★★**`F(z) = ∑_{(p,q) ∈ s} c(p,q) · exp((p x₀ + q x₁) z)`**

と、その 3 つの性質:

| | 宣言 |
|---|---|
| 整関数である | `auxFun_differentiable` |
| 増大度 `‖F(z)‖ ≤ (∑‖c‖)·exp(B‖z‖)` | `norm_auxFun_le` |
| ★★**格子点での値は 6 個の `exp(x_j y_i)` の単項式の線形結合** | `auxFun_latticePt` |
| Siegel の補題との接続(線形系が解ければ `F` は消える) | `auxFun_eq_zero_of_system` |
| 係数が代数的整数なら単項式も代数的整数 | `isIntegral_auxMatrix` |

## ★★★ここが「解析」と「算術」の接点である

`F(m₀y₀ + m₁y₁ + m₂y₂)` を展開すると

  `= ∑_{(p,q)} c(p,q) · ∏_i (E(0,i)^p · E(1,i)^q)^{m_i}`,  `E(j,i) = exp(x_j y_i)`

となる(`auxFun_latticePt`)。★**右辺は `c` について線形**なので、
「多くの格子点で `F` が消える」は `c` についての**線形系**である
—— そこに **数体上の Siegel の補題**
(`NumberField.house.exists_ne_zero_int_vec_house_le`、mathlib にある)を当てて
`c` を選ぶ。★★そして 6 個の `E(j,i)` が代数的なら値は代数的整数なので、
`Found/SixExp/Liouville.lean` の下界と突き合わせられる。

★これが six exponentials theorem の**枠**であり、残るのは
`extrapolation`(解析側の小ささと算術側の下界を回す帰納)である。

## ★測定(2026-08-18)

★mathlib に補助関数の構成は無い(超越性証明そのものが `Liouville` と
`Lindemann/AnalyticalPart` しか無い)。★一方で部品——`Complex.exp` の
加法性・`exp_nat_mul`・`exp_sum`、`IsIntegral.pow` / `.mul` / `.prod`——はすべてある。
-/

namespace ABC3.Found.SixExp

open Complex

/-! ## ★1. 補助関数 -/

/-- 補助関数の指数 `p·x₀ + q·x₁`。 -/
noncomputable def auxExp (x : Fin 2 → ℂ) (pq : ℕ × ℕ) : ℂ := pq.1 * x 0 + pq.2 * x 1

/-- ★★**補助関数** `F(z) = ∑_{(p,q) ∈ s} c(p,q) · exp((p x₀ + q x₁) z)`。 -/
noncomputable def auxFun (x : Fin 2 → ℂ) (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℂ) (z : ℂ) : ℂ :=
  ∑ pq ∈ s, c pq * Complex.exp (auxExp x pq * z)

/-- ★`F` は整関数。 -/
theorem auxFun_differentiable (x : Fin 2 → ℂ) (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℂ) :
    Differentiable ℂ (auxFun x s c) := by
  unfold auxFun
  refine Differentiable.fun_sum ?_
  intro pq _
  exact ((Complex.differentiable_exp.comp
    (differentiable_id.const_mul (auxExp x pq))).const_mul (c pq))

/-- ★**増大度** —— 指数和は有限指数の整関数。

★`Schwarz` の多零点版(`SchwarzZeros.lean`)に渡す `M` はこの評価から作る。 -/
theorem norm_auxFun_le (x : Fin 2 → ℂ) (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℂ)
    {B : ℝ} (hB : ∀ pq ∈ s, ‖auxExp x pq‖ ≤ B) (z : ℂ) :
    ‖auxFun x s c z‖ ≤ (∑ pq ∈ s, ‖c pq‖) * Real.exp (B * ‖z‖) := by
  rw [Finset.sum_mul]
  refine le_trans (norm_sum_le _ _) (Finset.sum_le_sum ?_)
  intro pq hpq
  rw [norm_mul]
  refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
  rw [Complex.norm_exp]
  refine Real.exp_le_exp.mpr (le_trans (Complex.re_le_norm _) ?_)
  rw [norm_mul]
  exact mul_le_mul_of_nonneg_right (hB pq hpq) (norm_nonneg _)

/-! ## ★2. 格子点での値 —— 解析と算術の接点 -/

/-- ★**六指数の 6 個の値** `E j i = exp(x_j · y_i)`。 -/
noncomputable def sixExpVals (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (j : Fin 2) (i : Fin 3) : ℂ :=
  Complex.exp (x j * y i)

/-- ★格子点 `∑ m_i y_i`。 -/
noncomputable def latticePt (y : Fin 3 → ℂ) (m : Fin 3 → ℕ) : ℂ := ∑ i : Fin 3, (m i : ℂ) * y i

/-- ★線形系の係数 —— `A[m,(p,q)] = ∏_i (E(0,i)^p · E(1,i)^q)^{m_i}`。 -/
noncomputable def auxMatrix (E : Fin 2 → Fin 3 → ℂ) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) : ℂ :=
  ∏ i : Fin 3, (E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i)

/-- ★★★**補助関数の格子点での値は、6 個の `exp(x_j y_i)` の単項式の線形結合**。

★これが「解析側の関数」と「算術側の代数的数」を結ぶ等式である。
証明は `exp` の加法性(`exp_sum` / `exp_add` / `exp_nat_mul`)だけ。 -/
theorem auxFun_latticePt (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → ℂ) (m : Fin 3 → ℕ) :
    auxFun x s c (latticePt y m) = ∑ pq ∈ s, auxMatrix (sixExpVals x y) m pq * c pq := by
  unfold auxFun latticePt auxMatrix sixExpVals
  refine Finset.sum_congr rfl ?_
  intro pq _
  rw [mul_comm (c pq)]
  congr 1
  rw [Finset.mul_sum, Complex.exp_sum]
  refine Finset.prod_congr rfl ?_
  intro i _
  rw [show auxExp x pq * ((m i : ℂ) * y i) = (m i : ℂ) * (auxExp x pq * y i) by ring,
    Complex.exp_nat_mul]
  congr 1
  rw [auxExp, add_mul, Complex.exp_add]
  congr 1
  · rw [show (pq.1 : ℂ) * x 0 * y i = (pq.1 : ℂ) * (x 0 * y i) by ring, Complex.exp_nat_mul]
  · rw [show (pq.2 : ℂ) * x 1 * y i = (pq.2 : ℂ) * (x 1 * y i) by ring, Complex.exp_nat_mul]

/-- ★★**Siegel の補題との接続** —— 線形系が解ければ、`F` は対応する格子点で消える。

★`T` は「消させたい格子点の添字」の有限集合。`c` をこの線形系の
**0 でない解**に取るのが Siegel の補題
(`NumberField.house.exists_ne_zero_int_vec_house_le`、mathlib にある)の役目である。 -/
theorem auxFun_eq_zero_of_system (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → ℂ) (T : Finset (Fin 3 → ℕ))
    (hsys : ∀ m ∈ T, ∑ pq ∈ s, auxMatrix (sixExpVals x y) m pq * c pq = 0) :
    ∀ m ∈ T, auxFun x s c (latticePt y m) = 0 := by
  intro m hm
  rw [auxFun_latticePt]
  exact hsys m hm

/-- ★**6 個の値が代数的整数なら、線形系の係数も代数的整数**。

★これで `Found/SixExp/Liouville.lean` の下界(`eq_zero_of_norm_embedding_lt`)に渡せる。 -/
theorem isIntegral_auxMatrix {E : Fin 2 → Fin 3 → ℂ} (hE : ∀ j i, IsIntegral ℤ (E j i))
    (m : Fin 3 → ℕ) (pq : ℕ × ℕ) : IsIntegral ℤ (auxMatrix E m pq) := by
  refine IsIntegral.prod _ (fun i _ => ?_)
  exact (((hE 0 i).pow pq.1).mul ((hE 1 i).pow pq.2)).pow (m i)

/-! ## ★3. 非空虚性の確認 -/

/-- ★**構成は空虚でない** —— `s = {(0,0)}`・`c ≡ 1` なら `F ≡ 1`。 -/
theorem auxFun_singleton_zero (x : Fin 2 → ℂ) (z : ℂ) :
    auxFun x {(0, 0)} (fun _ => 1) z = 1 := by
  simp [auxFun, auxExp]

end ABC3.Found.SixExp
