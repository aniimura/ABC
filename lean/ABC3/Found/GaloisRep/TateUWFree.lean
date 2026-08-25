import ABC3.Found.GaloisRep.PolyVanish

/-!
# Galois (G6) 第 237 ブロック —— **★★★★★★★`w` を自由にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★`τ` を消して `w` を自由にする

第 235 の評価は `z, τ : ℍ` の言葉で書かれている。しかし第 236 の
「`‖f(w)‖ ≤ C‖w‖ᵐ` なら `Xᵐ ∣ f`」を当てるには、**`w` を小さい円板の中で自由に動かせる**
必要がある。

★★**`0 < ‖w‖ < 1` なら `w = exp(2πi w')` なる上半平面の点 `w'` がある**
(`exists_uh_exp_eq`)——`w' = log w / (2πi)` とおけば `im w' = −log‖w‖/(2π) > 0`。
`τ := z + w'` とすれば `im z < im τ` で `subUH z τ = w'`、よって
`exp(2πi·subUH) = w`、`exp(2πiτ) = u·w` となる。

★★★これで第 235 が `(u, w)` の言葉になる(`norm_tateDefectTrunc_uw_le`)。

## ★★★★分子を `w` の多項式と見る

`P ∈ ℤ[A,W]` を `A ↦ u`(定数)、`W ↦ X` で `ℂ[w]` に送る(`evalW`)。
★`w` を代入すれば `tateEval u w P` に戻る(`eval_evalW`)——
`MvPolynomial.ringHom_ext`(`C` と `X i` での一致で環準同型が決まる)で示す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalW` / `eval_evalW` | ★★★★分子を `w` の多項式と見る |
| `exists_uh_exp_eq` | ★★★`0 < ‖w‖ < 1` は上半平面から来る |
| `addUH` / `subUH_addUH` | ★`τ := z + w'` の配線 |
| `norm_tateDefectTrunc_uw_le` | ★★★★★★★**評価の `(u,w)` 形** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★分子を `w` の多項式と見る -/

/-- ★★`A ↦ u`(定数)、`W ↦ X` として `ℤ[A,W]` を `ℂ[w]` に送る。 -/
noncomputable def evalW (u : ℂ) : TateBase →+* Polynomial ℂ :=
  MvPolynomial.eval₂Hom ((Polynomial.C : ℂ →+* Polynomial ℂ).comp (Int.castRingHom ℂ))
    ![Polynomial.C u, Polynomial.X]

theorem evalW_univA (u : ℂ) : evalW u univA = Polynomial.C u := by
  rw [evalW, univA, MvPolynomial.eval₂Hom_X']
  rfl

theorem evalW_univW (u : ℂ) : evalW u univW = Polynomial.X := by
  rw [evalW, univW, MvPolynomial.eval₂Hom_X']
  rfl

/-- ★★★★**`w` を代入すると `tateEval u w` に戻る**。 -/
theorem eval_evalW (u w : ℂ) (P : TateBase) :
    (evalW u P).eval w = tateEval u w P := by
  have hhom : (Polynomial.evalRingHom w).comp (evalW u) = tateEval u w := by
    refine MvPolynomial.ringHom_ext ?_ ?_
    · intro r
      simp [evalW, tateEval]
    · intro i
      fin_cases i
      · simp [evalW, tateEval, MvPolynomial.eval₂Hom_X']
      · simp [evalW, tateEval, MvPolynomial.eval₂Hom_X']
  exact congrArg (fun f => f P) hhom

/-! ## ★★★`0 < ‖w‖ < 1` は上半平面から来る -/

/-- ★★★**`0 < ‖w‖ < 1` なら `w = exp(2πi w')` なる上半平面の点 `w'` がある**。

★`w' = log w / (2πi)` とおけば `im w' = −log‖w‖/(2π) > 0`。 -/
theorem exists_uh_exp_eq (w : ℂ) (hw0 : w ≠ 0) (hw : ‖w‖ < 1) :
    ∃ w' : UpperHalfPlane, Complex.exp (2 * ↑π * I * (w' : ℂ)) = w := by
  have hpi : (0 : ℝ) < π := Real.pi_pos
  set t : ℂ := Complex.log w / (2 * ↑π * I) with htdef
  have hne : (2 * (π : ℂ) * I) ≠ 0 := two_pi_I_ne_zero
  have hval : t = -Complex.I * Complex.log w / (2 * (π : ℝ)) := by
    rw [htdef]
    field_simp
    ring_nf
    rw [Complex.I_sq]
    ring
  have h2 : t.im = -(Complex.log w).re / (2 * π) := by
    rw [hval, Complex.div_im]
    simp [Complex.normSq_apply]
    field_simp
  have him : 0 < t.im := by
    rw [h2, Complex.log_re]
    have hpos : 0 < ‖w‖ := norm_pos_iff.2 hw0
    have hlogneg : Real.log ‖w‖ < 0 := Real.log_neg hpos hw
    have h2pi : 0 < 2 * π := by linarith
    exact div_pos (by linarith) h2pi
  refine ⟨⟨t, him⟩, ?_⟩
  have hmul : 2 * (π : ℂ) * I * t = Complex.log w := by
    rw [htdef]
    field_simp
  show Complex.exp (2 * ↑π * I * t) = w
  rw [hmul, Complex.exp_log hw0]

/-! ## ★`τ := z + w'` の配線 -/

/-- ★上半平面の 2 点の和も上半平面の点。 -/
noncomputable def addUH (z w' : UpperHalfPlane) : UpperHalfPlane :=
  ⟨(z : ℂ) + (w' : ℂ), by
    have h1 : 0 < (z : ℂ).im := z.2
    have h2 : 0 < (w' : ℂ).im := w'.2
    simp only [Complex.add_im]
    linarith⟩

@[simp] theorem addUH_coe (z w' : UpperHalfPlane) :
    ((addUH z w' : UpperHalfPlane) : ℂ) = (z : ℂ) + (w' : ℂ) := rfl

theorem im_lt_addUH (z w' : UpperHalfPlane) :
    (z : ℂ).im < ((addUH z w' : UpperHalfPlane) : ℂ).im := by
  have h2 : 0 < (w' : ℂ).im := w'.2
  rw [addUH_coe]
  simp only [Complex.add_im]
  linarith

theorem subUH_addUH (z w' : UpperHalfPlane) :
    ((subUH z (addUH z w') (im_lt_addUH z w') : UpperHalfPlane) : ℂ) = (w' : ℂ) := by
  rw [subUH_coe, addUH_coe]
  ring

/-! ## ★★★★★★★評価の `(u,w)` 形 -/

/-- ★★★★★★★**切り詰めた差の評価(`(u,w)` の形)**——`τ` を消して `w` を自由にした。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_tateDefectTrunc_uw_le (z : UpperHalfPlane) (w : ℂ) (hw0 : w ≠ 0)
    (hw : ‖w‖ ≤ 1 / 2)
    (hq : ‖Complex.exp (2 * ↑π * I * (z : ℂ)) * w‖ ≤ 1 / 128) (n : ℕ) :
    ‖tateDefectTrunc (n + 1) (Complex.exp (2 * ↑π * I * z)) w
        (Complex.exp (2 * ↑π * I * z) * w)‖
      ≤ (3 * (‖tateXterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖
          + ‖tateYterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖ + 30 * ((n : ℝ) + 1) + 80) ^ 2
        + 6 * (‖tateXterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖
          + ‖tateYterm (Complex.exp (2 * ↑π * I * (z : ℂ)))‖ + 30 * ((n : ℝ) + 1) + 80) + 1)
        * (50 * ((64 : ℝ) * ‖Complex.exp (2 * ↑π * I * (z : ℂ)) * w‖) ^ (n + 1)) := by
  obtain ⟨w', hw'⟩ := exists_uh_exp_eq w hw0 (by linarith)
  set τ := addUH z w' with hτdef
  have him : (z : ℂ).im < (τ : ℂ).im := im_lt_addUH z w'
  have hsub : Complex.exp (2 * ↑π * I * ((subUH z τ him : UpperHalfPlane) : ℂ)) = w := by
    rw [subUH_addUH z w']
    exact hw'
  have hτ : Complex.exp (2 * ↑π * I * (τ : ℂ))
      = Complex.exp (2 * ↑π * I * (z : ℂ)) * w := by
    rw [← hsub]
    exact (exp_mul_exp_sub z τ him).symm
  have h := norm_tateDefectTrunc_le z τ him (by rw [hτ]; exact hq) (by rw [hsub]; exact hw) n
  rw [hsub, hτ] at h
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def evalW.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子を w の多項式と見る)",
    sectionId := "genell-def-3-3" }

def exists_uh_exp_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——0 < ‖w‖ < 1 は上半平面から来る)",
    sectionId := "genell-def-3-3" }

def norm_tateDefectTrunc_uw_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——評価の (u,w) 形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
