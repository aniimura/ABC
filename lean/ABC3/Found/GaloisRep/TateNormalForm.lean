import ABC3.Found.GaloisRep.TateNodeMove

/-!
# Galois (G6) 第 307 ブロック —— **★★★★★★★★★Tate 標準形への正規化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> 分裂乗法還元をもつ曲線は、変数変換で
> **`y² + xy = x³ + a₄x + a₆`(`a₄, a₆ ∈ 𝔪`)**に直せる
> (`exists_tate_normal_form_of_split`)

★★★これが **Tate 標準形**である。G6 の第 2 段がこれで閉じた。

## ★★★★★★4 つの変数変換

| 段 | 変換 | 効き目 |
|---|---|---|
| 1 | `(1, r₀, 0, t₀)` | 節点を原点へ(第 306) |
| 2 | `(1, 0, α, 0)` | 接線で剪断——`a₂ ∈ 𝔪`(第 305 の `α`) |
| 3 | `(a₁ + 2α, 0, 0, 0)` | 尺度を合わせて **`a₁ = 1`** |
| 4 | `(1 + 2s, −a₃, s, 0)` | **`a₂ = a₃ = 0`**(Artin–Schreier の Hensel) |

★★★★第 3 段が効くのは **`a₁ + 2α` が単元**だから。これは
`(c₄(a₁ + 2α))² = −c₄c₆`(第 305 の判別式)から**根の式だけで**出る——
剰余体の議論に戻らなくてよい。

## ★★★★★★★第 4 段は 1 本の Artin–Schreier で足りる

`a₂ = a₃ = 0` を同時に達成するには、素朴には 2 本の方程式が要るように見える。
★しかし `u = 1 + 2s` と結ぶと `a₁ = 1` が保たれ、`r = −a₃` で `a₃ = 0` が出て、
残るのは **`s² + s = a₂ − 3a₃`** の 1 本だけになる。
★★★★`s² + s = c`(`c ∈ 𝔪`)は **`X² + X − c` の微分が `1`** なので
標数によらず Hensel が効く(`exists_artin_schreier`)。
★★★★★**剰余標数 2 でも平方根を取らずに済む**——ここが分岐点だった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isUnit_a1_add_two_root` | ★★★★★`a₁ + 2α` は単元 |
| `shear_a2_mem` | ★★★★★剪断で `a₂ ∈ 𝔪` |
| `exists_artin_schreier` | ★★★★`s² + s = c` は解ける |
| `exists_normal_a1` | ★★★★★★`a₁ = 1` まで |
| `exists_tate_normal_form_of_split` | ★★★★★★★★★**Tate 標準形** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

section Coeff

variable {R : Type} [CommRing R]

/-- ★★接線の 2 次式の定数項(符号込み)。 -/
def tangentQ (E : WeierstrassCurve R) : R := 54 * E.b₆ - 3 * E.b₂ * E.b₄ + E.a₂ * E.c₄

/-- ★★★接線の 2 次式は**平行移動で変わらない**(節点の `x` 座標がずれる分と相殺する)。 -/
theorem tangentQ_trans (E : WeierstrassCurve R) (r t : R) :
    tangentQ ((⟨1, r, 0, t⟩ : VariableChange R) • E) = tangentQ E := by
  simp only [tangentQ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.c₄, WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂,
    WeierstrassCurve.variableChange_a₃, WeierstrassCurve.variableChange_a₄,
    WeierstrassCurve.variableChange_a₆, inv_one, Units.val_one, one_pow, one_mul]
  ring

theorem c4_trans (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).c₄ = E.c₄ := by
  rw [WeierstrassCurve.variableChange_c₄]
  simp only [inv_one, Units.val_one, one_pow, one_mul]

theorem c6_trans (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).c₆ = E.c₆ := by
  rw [WeierstrassCurve.variableChange_c₆]
  simp only [inv_one, Units.val_one, one_pow, one_mul]

theorem shear_a1 (E : WeierstrassCurve R) (s : R) :
    ((⟨1, 0, s, 0⟩ : VariableChange R) • E).a₁ = E.a₁ + 2 * s := by
  rw [WeierstrassCurve.variableChange_a₁]
  simp only [inv_one, Units.val_one, one_mul]

theorem shear_a3 (E : WeierstrassCurve R) (s : R) :
    ((⟨1, 0, s, 0⟩ : VariableChange R) • E).a₃ = E.a₃ := by
  rw [WeierstrassCurve.variableChange_a₃]
  simp only [inv_one, Units.val_one, one_pow, one_mul, zero_mul, mul_zero, add_zero]

theorem shear_a4 (E : WeierstrassCurve R) (s : R) :
    ((⟨1, 0, s, 0⟩ : VariableChange R) • E).a₄ = E.a₄ - s * E.a₃ := by
  rw [WeierstrassCurve.variableChange_a₄]
  simp only [inv_one, Units.val_one, one_pow, one_mul]
  ring

theorem shear_a6 (E : WeierstrassCurve R) (s : R) :
    ((⟨1, 0, s, 0⟩ : VariableChange R) • E).a₆ = E.a₆ := by
  rw [WeierstrassCurve.variableChange_a₆]
  simp only [inv_one, Units.val_one, one_pow, one_mul]
  ring

theorem scale_a1 (E : WeierstrassCurve R) (u : Rˣ) :
    ((⟨u, 0, 0, 0⟩ : VariableChange R) • E).a₁ = ((u⁻¹ : Rˣ) : R) * E.a₁ := by
  rw [WeierstrassCurve.variableChange_a₁]
  ring

theorem scale_a2 (E : WeierstrassCurve R) (u : Rˣ) :
    ((⟨u, 0, 0, 0⟩ : VariableChange R) • E).a₂ = ((u⁻¹ : Rˣ) : R) ^ 2 * E.a₂ := by
  rw [WeierstrassCurve.variableChange_a₂]
  ring

theorem scale_a3 (E : WeierstrassCurve R) (u : Rˣ) :
    ((⟨u, 0, 0, 0⟩ : VariableChange R) • E).a₃ = ((u⁻¹ : Rˣ) : R) ^ 3 * E.a₃ := by
  rw [WeierstrassCurve.variableChange_a₃]
  ring

theorem scale_a4 (E : WeierstrassCurve R) (u : Rˣ) :
    ((⟨u, 0, 0, 0⟩ : VariableChange R) • E).a₄ = ((u⁻¹ : Rˣ) : R) ^ 4 * E.a₄ := by
  rw [WeierstrassCurve.variableChange_a₄]
  ring

theorem scale_a6 (E : WeierstrassCurve R) (u : Rˣ) :
    ((⟨u, 0, 0, 0⟩ : VariableChange R) • E).a₆ = ((u⁻¹ : Rˣ) : R) ^ 6 * E.a₆ := by
  rw [WeierstrassCurve.variableChange_a₆]
  ring

theorem full_a1 (E : WeierstrassCurve R) (u : Rˣ) (r s t : R) :
    ((⟨u, r, s, t⟩ : VariableChange R) • E).a₁ = ((u⁻¹ : Rˣ) : R) * (E.a₁ + 2 * s) := by
  rw [WeierstrassCurve.variableChange_a₁]

theorem full_a2 (E : WeierstrassCurve R) (u : Rˣ) (r s t : R) :
    ((⟨u, r, s, t⟩ : VariableChange R) • E).a₂
      = ((u⁻¹ : Rˣ) : R) ^ 2 * (E.a₂ - s * E.a₁ + 3 * r - s ^ 2) := by
  rw [WeierstrassCurve.variableChange_a₂]

theorem full_a3 (E : WeierstrassCurve R) (u : Rˣ) (r s t : R) :
    ((⟨u, r, s, t⟩ : VariableChange R) • E).a₃
      = ((u⁻¹ : Rˣ) : R) ^ 3 * (E.a₃ + r * E.a₁ + 2 * t) := by
  rw [WeierstrassCurve.variableChange_a₃]

theorem full_a4 (E : WeierstrassCurve R) (u : Rˣ) (r s t : R) :
    ((⟨u, r, s, t⟩ : VariableChange R) • E).a₄
      = ((u⁻¹ : Rˣ) : R) ^ 4 * (E.a₄ - s * E.a₃ + 2 * r * E.a₂ - (t + r * s) * E.a₁
        + 3 * r ^ 2 - 2 * s * t) := by
  rw [WeierstrassCurve.variableChange_a₄]

theorem full_a6 (E : WeierstrassCurve R) (u : Rˣ) (r s t : R) :
    ((⟨u, r, s, t⟩ : VariableChange R) • E).a₆
      = ((u⁻¹ : Rˣ) : R) ^ 6 * (E.a₆ + r * E.a₄ + r ^ 2 * E.a₂ + r ^ 3 - t * E.a₃ - t ^ 2
        - r * t * E.a₁) := by
  rw [WeierstrassCurve.variableChange_a₆]

/-- ★★★★★**`a₁ + 2α` は単元**——判別式 `(c₄(a₁+2α))² = −c₄c₆` から。 -/
theorem isUnit_a1_add_two_root (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄) (hc6 : IsUnit E.c₆)
    (α : R) (hα : E.c₄ * α ^ 2 + E.a₁ * E.c₄ * α - tangentQ E = 0) :
    IsUnit (E.a₁ + 2 * α) := by
  have hd := tangent_disc_eq E
  have hkey : (E.c₄ * (E.a₁ + 2 * α)) ^ 2 = -(E.c₄ * E.c₆) := by
    rw [tangentQ] at hα
    linear_combination 4 * E.c₄ * hα + hd
  have h1 : IsUnit ((E.c₄ * (E.a₁ + 2 * α)) ^ 2) := by
    rw [hkey]
    exact (hc4.mul hc6).neg
  exact isUnit_of_mul_isUnit_right ((isUnit_pow_iff two_ne_zero).1 h1)

end Coeff

section LocalRing

variable {R : Type} [CommRing R] [IsLocalRing R]

/-- ★★★★★**剪断すれば `a₂ ∈ 𝔪`**——接線の根の効き目。 -/
theorem shear_a2_mem (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄)
    (h3 : E.a₃ ∈ IsLocalRing.maximalIdeal R) (h4 : E.a₄ ∈ IsLocalRing.maximalIdeal R)
    (h6 : E.a₆ ∈ IsLocalRing.maximalIdeal R) (α : R)
    (hα : E.c₄ * α ^ 2 + E.a₁ * E.c₄ * α - tangentQ E = 0) :
    ((⟨1, 0, α, 0⟩ : VariableChange R) • E).a₂ ∈ IsLocalRing.maximalIdeal R := by
  have hb4 : E.b₄ ∈ IsLocalRing.maximalIdeal R := by
    have hb : E.b₄ = 2 * E.a₄ + E.a₁ * E.a₃ := rfl
    rw [hb]
    exact Ideal.add_mem _ (Ideal.mul_mem_left _ _ h4) (Ideal.mul_mem_left _ _ h3)
  have hb6 : E.b₆ ∈ IsLocalRing.maximalIdeal R := by
    have hb : E.b₆ = E.a₃ ^ 2 + 4 * E.a₆ := rfl
    rw [hb]
    exact Ideal.add_mem _ (by rw [pow_two]; exact Ideal.mul_mem_left _ _ h3)
      (Ideal.mul_mem_left _ _ h6)
  rw [WeierstrassCurve.variableChange_a₂]
  simp only [inv_one, Units.val_one, one_pow, one_mul]
  refine (Ideal.unit_mul_mem_iff_mem _ hc4).1 ?_
  have hval : E.c₄ * (E.a₂ - α * E.a₁ + 3 * 0 - α ^ 2) = -(54 * E.b₆ - 3 * E.b₂ * E.b₄) := by
    rw [tangentQ] at hα
    linear_combination -hα
  rw [hval]
  exact neg_mem (Ideal.sub_mem _ (Ideal.mul_mem_left _ _ hb6) (Ideal.mul_mem_left _ _ hb4))

/-- ★★局所環では `1 + 𝔪` は単元。 -/
theorem isUnit_one_add_mem {x : R} (hx : x ∈ IsLocalRing.maximalIdeal R) : IsUnit (1 + x) := by
  by_contra h
  have hm : (1 + x) ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal _).2 h
  have h1 : (1 : R) ∈ IsLocalRing.maximalIdeal R := by
    have hx1 : (1 : R) = (1 + x) - x := by ring
    rw [hx1]
    exact Ideal.sub_mem _ hm hx
  exact (IsLocalRing.maximalIdeal R).ne_top_iff_one.1 (Ideal.IsMaximal.ne_top inferInstance) h1

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**`a₁ = 1` まで正規化できる**(節点 → 剪断 → 尺度)。 -/
theorem exists_normal_a1 (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄) (hc6 : IsUnit E.c₆)
    (hΔ : E.Δ ∈ IsLocalRing.maximalIdeal R) (α : R)
    (hα : E.c₄ * α ^ 2 + E.a₁ * E.c₄ * α - tangentQ E = 0) :
    ∃ C : VariableChange R, (C • E).a₁ = 1
      ∧ (C • E).a₂ ∈ IsLocalRing.maximalIdeal R ∧ (C • E).a₃ ∈ IsLocalRing.maximalIdeal R
      ∧ (C • E).a₄ ∈ IsLocalRing.maximalIdeal R ∧ (C • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨r₀, t₀, h3, h4, h6⟩ := exists_node_translation E hc4 hΔ
  have hc4₁ : IsUnit ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E).c₄ := by
    rw [c4_trans]; exact hc4
  have hc6₁ : IsUnit ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E).c₆ := by
    rw [c6_trans]; exact hc6
  have hα₁ : ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E).c₄ * α ^ 2
      + ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E).a₁
        * ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E).c₄ * α
      - tangentQ ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E) = 0 := by
    rw [c4_trans, trans_a₁, tangentQ_trans]
    exact hα
  have hu₂ : IsUnit (((⟨1, 0, α, 0⟩ : VariableChange R)
      • ((⟨1, r₀, 0, t₀⟩ : VariableChange R) • E)).a₁) := by
    rw [shear_a1]
    exact isUnit_a1_add_two_root _ hc4₁ hc6₁ α hα₁
  refine ⟨(⟨hu₂.unit, 0, 0, 0⟩ : VariableChange R)
    * ((⟨1, 0, α, 0⟩ : VariableChange R) * (⟨1, r₀, 0, t₀⟩ : VariableChange R)), ?_, ?_, ?_, ?_, ?_⟩
  all_goals rw [mul_smul, mul_smul]
  · rw [scale_a1]
    have h2 : ((hu₂.unit⁻¹ : Rˣ) : R) * ((hu₂.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [hu₂.unit_spec] at h2
    exact h2
  · rw [scale_a2]
    exact Ideal.mul_mem_left _ _ (shear_a2_mem _ hc4₁ h3 h4 h6 α hα₁)
  · rw [scale_a3, shear_a3]
    exact Ideal.mul_mem_left _ _ h3
  · rw [scale_a4, shear_a4]
    exact Ideal.mul_mem_left _ _ (Ideal.sub_mem _ h4 (Ideal.mul_mem_left _ _ h3))
  · rw [scale_a6, shear_a6]
    exact Ideal.mul_mem_left _ _ h6

end LocalRing

section Complete

variable {R : Type} [CommRing R] [IsLocalRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]

/-- ★★★★**`s² + s = c`(`c ∈ 𝔪`)は `𝔪` の中で解ける**——微分が `1` なので標数によらない。 -/
theorem exists_artin_schreier (c : R) (hc : c ∈ IsLocalRing.maximalIdeal R) :
    ∃ s : R, s ∈ IsLocalRing.maximalIdeal R ∧ s ^ 2 + s = c := by
  set f : Polynomial R := Polynomial.X ^ 2 + Polynomial.X - Polynomial.C c with hf
  have hmonic : f.Monic := by
    rw [hf]
    monicity!
  have heval : Polynomial.eval 0 f ∈ IsLocalRing.maximalIdeal R := by
    rw [hf]
    simp only [Polynomial.eval_sub, Polynomial.eval_add, Polynomial.eval_pow, Polynomial.eval_X,
      Polynomial.eval_C]
    simpa using neg_mem hc
  have hderiv : Polynomial.derivative f = 2 * Polynomial.X + 1 := by
    rw [hf]
    simp [Polynomial.derivative_sub, Polynomial.derivative_add]
    ring
  have hunit : IsUnit ((Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))
      (Polynomial.eval 0 (Polynomial.derivative f))) := by
    rw [hderiv]
    simp
  obtain ⟨s, hs, hsm⟩ := (IsAdicComplete.henselianRing R (IsLocalRing.maximalIdeal R)).is_henselian
    f hmonic 0 heval hunit
  refine ⟨s, by simpa using hsm, ?_⟩
  have h := hs
  rw [Polynomial.IsRoot, hf] at h
  simp only [Polynomial.eval_sub, Polynomial.eval_add, Polynomial.eval_pow, Polynomial.eval_X,
    Polynomial.eval_C] at h
  linear_combination h

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**`a₁ = 1` からは `a₂ = a₃ = 0` まで行ける**——Artin–Schreier 1 本で。 -/
theorem exists_normal_form_of_a1 (E : WeierstrassCurve R) (h1 : E.a₁ = 1)
    (h2 : E.a₂ ∈ IsLocalRing.maximalIdeal R) (h3 : E.a₃ ∈ IsLocalRing.maximalIdeal R)
    (h4 : E.a₄ ∈ IsLocalRing.maximalIdeal R) (h6 : E.a₆ ∈ IsLocalRing.maximalIdeal R) :
    ∃ C : VariableChange R, (C • E).a₁ = 1 ∧ (C • E).a₂ = 0 ∧ (C • E).a₃ = 0
      ∧ (C • E).a₄ ∈ IsLocalRing.maximalIdeal R
      ∧ (C • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨s, hs, hseq⟩ := exists_artin_schreier (E.a₂ - 3 * E.a₃)
    (Ideal.sub_mem _ h2 (Ideal.mul_mem_left _ _ h3))
  have hu : IsUnit (1 + 2 * s) := isUnit_one_add_mem (Ideal.mul_mem_left _ _ hs)
  have huv : ((hu.unit⁻¹ : Rˣ) : R) * (1 + 2 * s) = 1 := by
    have hx : ((hu.unit⁻¹ : Rˣ) : R) * ((hu.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [hu.unit_spec] at hx
    exact hx
  refine ⟨⟨hu.unit, -E.a₃, s, 0⟩, ?_, ?_, ?_, ?_, ?_⟩
  · rw [full_a1, h1]
    exact huv
  · rw [full_a2, h1]
    have hz : E.a₂ - s * 1 + 3 * (-E.a₃) - s ^ 2 = 0 := by
      linear_combination -hseq
    rw [hz, mul_zero]
  · rw [full_a3, h1]
    have hz : E.a₃ + (-E.a₃) * 1 + 2 * 0 = 0 := by ring
    rw [hz, mul_zero]
  · rw [full_a4, h1]
    refine Ideal.mul_mem_left _ _ ?_
    have hz : E.a₄ - s * E.a₃ + 2 * (-E.a₃) * E.a₂ - (0 + (-E.a₃) * s) * 1
        + 3 * (-E.a₃) ^ 2 - 2 * s * 0
        = E.a₄ - s * E.a₃ + (-2 * E.a₂) * E.a₃ + s * E.a₃ + (3 * E.a₃) * E.a₃ := by ring
    rw [hz]
    exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _
      (Ideal.sub_mem _ h4 (Ideal.mul_mem_left _ _ h3)) (Ideal.mul_mem_left _ _ h3))
      (Ideal.mul_mem_left _ _ h3)) (Ideal.mul_mem_left _ _ h3)
  · rw [full_a6, h1]
    refine Ideal.mul_mem_left _ _ ?_
    have hz : E.a₆ + (-E.a₃) * E.a₄ + (-E.a₃) ^ 2 * E.a₂ + (-E.a₃) ^ 3 - 0 * E.a₃ - 0 ^ 2
        - (-E.a₃) * 0 * 1
        = E.a₆ + (-E.a₄) * E.a₃ + (E.a₃ * E.a₂) * E.a₃ + (-E.a₃ * E.a₃) * E.a₃ := by ring
    rw [hz]
    exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ h6 (Ideal.mul_mem_left _ _ h3))
      (Ideal.mul_mem_left _ _ h3)) (Ideal.mul_mem_left _ _ h3)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**Tate 標準形**——`y² + xy = x³ + a₄x + a₆`(`a₄, a₆ ∈ 𝔪`)。 -/
theorem exists_tate_normal_form (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄) (hc6 : IsUnit E.c₆)
    (hΔ : E.Δ ∈ IsLocalRing.maximalIdeal R) (α : R)
    (hα : E.c₄ * α ^ 2 + E.a₁ * E.c₄ * α - tangentQ E = 0) :
    ∃ C : VariableChange R, (C • E).a₁ = 1 ∧ (C • E).a₂ = 0 ∧ (C • E).a₃ = 0
      ∧ (C • E).a₄ ∈ IsLocalRing.maximalIdeal R
      ∧ (C • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨C₁, k1, k2, k3, k4, k6⟩ := exists_normal_a1 E hc4 hc6 hΔ α hα
  obtain ⟨C₂, m1, m2, m3, m4, m6⟩ := exists_normal_form_of_a1 (C₁ • E) k1 k2 k3 k4 k6
  refine ⟨C₂ * C₁, ?_, ?_, ?_, ?_, ?_⟩
  all_goals rw [mul_smul]
  · exact m1
  · exact m2
  · exact m3
  · exact m4
  · exact m6

end Complete

/-! ## ★★★★★★★★★分裂乗法還元からの Tate 標準形 -/

section Split

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★**分裂乗法還元をもつ曲線は Tate 標準形に直せる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tate_normal_form_of_split (W : WeierstrassCurve K)
    [hmul : W.HasMultiplicativeReduction R] (hsplit : W.HasSplitMultiplicativeReduction R)
    (hΔ0 : W.Δ ≠ 0) :
    ∃ C : VariableChange R, (C • integralModel R W).a₁ = 1
      ∧ (C • integralModel R W).a₂ = 0 ∧ (C • integralModel R W).a₃ = 0
      ∧ (C • integralModel R W).a₄ ∈ IsLocalRing.maximalIdeal R
      ∧ (C • integralModel R W).a₆ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨α, hα⟩ := exists_tangent_root W hsplit hΔ0
  exact exists_tate_normal_form (integralModel R W) (integralModel_c4_isUnit W)
    (integralModel_c6_isUnit W hΔ0) (integralModel_Delta_mem W hΔ0) α hα

end Split

/-! ## ★出典の紐付け(`.src`) -/

def exists_artin_schreier.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Artin-Schreier 型の Hensel)",
    sectionId := "genell-def-3-3" }

def exists_tate_normal_form.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 標準形への正規化)",
    sectionId := "genell-def-3-3" }

def exists_tate_normal_form_of_split.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(分裂乗法還元から Tate 標準形へ)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
