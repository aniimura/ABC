import ABC3.Found.GaloisRep.TateAnalyticTrunc

/-!
# Galois (G6) 第 232 ブロック —— **★★★★★★★`(z,τ)` と `(u,w)` を繋ぐ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★二つの座標系を繋ぐ

第 220 の解析側の Tate 方程式は `(z, τ)`(上半平面の 2 点)で書かれている。
第 231 の評価は `(u, w)`(ℂ の 2 元、`q = uw`)で書かれている。本ブロックで繋ぐ。

    u := exp(2πi z),   w := exp(2πi (τ − z)),   uw = exp(2πi τ) = q

★★**`im z < im τ` なら `τ − z` も上半平面の点である**(`subUH`)。
これで `‖w‖ < 1` が `UpperHalfPlane.norm_exp_two_pi_I_lt_one` から自動的に出る。
★`him` は第 216 から持ち歩いている仮定で、ここでも同じものが効いている
(第 225 では `z ∉ 格子` を、ここでは `w` が上半平面の点であることを与える)。

## ★★★★★★繋がった形

第 226 の分解と第 227 の `s₁` の ℤ 和(`∑_{n∈ℤ}f(q^{−n}) = 2∑_{m≥1}f(qᵐ)`)を使うと:

    tateXfun z τ = tateXanalytic u w
    tateYfun z τ = tateYanalytic u w

★★★したがって**第 220 の Tate 方程式がそのまま `(u,w)` の形で成り立つ**
(`tate_equation_uw`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `subUH` | ★★`τ − z` は上半平面の点 |
| `exp_mul_exp_sub` | ★`u·w = q` |
| `tateXfun_eq_analytic` | ★★★★★★`tateXfun z τ = tateXanalytic u w` |
| `tateYfun_eq_analytic` | ★★★★★★`tateYfun z τ = tateYanalytic u w` |
| `tate_equation_uw` | ★★★★★★★**解析側の Tate 方程式(`(u,w)` の形)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★`τ − z` も上半平面の点 -/

/-- ★`im z < im τ` なら `τ − z` も上半平面の点である。 -/
noncomputable def subUH (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    UpperHalfPlane :=
  ⟨(τ : ℂ) - (z : ℂ), by
    simp only [Complex.sub_im]
    linarith⟩

@[simp] theorem subUH_coe (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    ((subUH z τ him : UpperHalfPlane) : ℂ) = (τ : ℂ) - (z : ℂ) := rfl

/-- ★`u·w = q`。 -/
theorem exp_mul_exp_sub (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    Complex.exp (2 * ↑π * I * z) * Complex.exp (2 * ↑π * I * (subUH z τ him))
      = Complex.exp (2 * ↑π * I * τ) := by
  rw [subUH_coe, ← Complex.exp_add]
  ring_nf

/-! ## ★★★★★★二つの形が一致する -/

/-- ★★★★★★**`tateXfun` は `(u,w)` の形の `tateXanalytic` と同じもの**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXfun_eq_analytic (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    tateXfun z τ
      = tateXanalytic (Complex.exp (2 * ↑π * I * z))
        (Complex.exp (2 * ↑π * I * (subUH z τ him))) := by
  have hu : ‖Complex.exp (2 * ↑π * I * (z : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have hw : ‖Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one (subUH z τ him)
  have hu0 : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hw0 : Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ)) ≠ 0 :=
    Complex.exp_ne_zero _
  have hmul := exp_mul_exp_sub z τ him
  have hq : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  have hq0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  rw [tateXfun, tateXanalytic, ← hmul,
    tsum_int_tateXterm_split _ _ hu hw hu0 hw0]
  rw [hmul, tsum_int_tateXterm_one _ hq hq0, ← hmul]

/-- ★★★★★★**`tateYfun` は `(u,w)` の形の `tateYanalytic` と同じもの**。 -/
theorem tateYfun_eq_analytic (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    tateYfun z τ
      = tateYanalytic (Complex.exp (2 * ↑π * I * z))
        (Complex.exp (2 * ↑π * I * (subUH z τ him))) := by
  have hu : ‖Complex.exp (2 * ↑π * I * (z : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have hw : ‖Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one (subUH z τ him)
  have hu0 : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hw0 : Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ)) ≠ 0 :=
    Complex.exp_ne_zero _
  have hmul := exp_mul_exp_sub z τ him
  have hq : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  have hq0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  rw [tateYfun, tateYanalytic, ← hmul,
    tsum_int_tateYterm_split _ _ hu hw hu0 hw0]
  rw [hmul, tsum_int_tateXterm_one _ hq hq0, ← hmul]
  ring

/-! ## ★★★★★★★`(u,w)` の形の Tate 方程式 -/

/-- ★★★★★★★**解析側の Tate 方程式(`(u,w)` の形)**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_equation_uw (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    tateYanalytic (Complex.exp (2 * ↑π * I * z))
        (Complex.exp (2 * ↑π * I * (subUH z τ him))) ^ 2
      + tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
        * tateYanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
      = tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him))) ^ 3
        + (-5 * sigmaSum 3 τ) * tateXanalytic (Complex.exp (2 * ↑π * I * z))
          (Complex.exp (2 * ↑π * I * (subUH z τ him)))
        - (5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12 := by
  rw [← tateXfun_eq_analytic z τ him, ← tateYfun_eq_analytic z τ him]
  exact tate_equation_analytic' z τ him

/-! ## ★出典の紐付け(`.src`) -/

def subUH.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——tau - z も上半平面の点)",
    sectionId := "genell-def-3-3" }

def tateXfun_eq_analytic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——(z,tau) と (u,w) を繋ぐ)",
    sectionId := "genell-def-3-3" }

def tate_equation_uw.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の方程式の (u,w) 形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
