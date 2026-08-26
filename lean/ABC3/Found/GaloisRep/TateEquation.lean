import ABC3.Found.GaloisRep.TateDescent

/-!
# Galois (G6) 第 240 ブロック —— **★★★★★★★★★★葉 (b) 完成**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★葉 (b) —— Tate 級数は Weierstrass 方程式を満たす

    Y(u,q)² + X(u,q)·Y(u,q) = X(u,q)³ + a₄(q)·X(u,q) + a₆(q)

を**任意の `I` 進完備環**で示した(`tate_equation`)。第 212 で `I` を法として取れていた
ものが、これで**厳密な等式**になった。

## ★★★もう一方の整除性は対称性で出る

第 239 で `W^{n+1} ∣ P` が出た。`A^{n+1} ∣ P` は同じ議論の**役割の入れ替え**である:

★★`tateDefectTrunc (n+1) u w (uw) = tateDefectTrunc (n+1) w u (wu)`(第 224 の
`tateDefectTrunc_symm`)なので、**`w` を固定して `u → 0` を見る**だけでよい。
評価(第 237)はそのまま `(w, u)` の順で当たる。
★局所化の側に「`A` と `W` を入れ替える自己同型」を作る必要は無かった——
`tateDefectTrunc` の対称性だけで足りる。

`ℤ[A,W] ≅ ℤ[w][a]`(`A` を外側に)を `toPP2` / `ofPP2` で作り、
第 239 と同じ流れで `A^{n+1} ∣ P` を出す。

## ★★仕上げ

第 224 の `tateDefect_eq_zero_of_base` は `∀ n` の形なので、`n = 0` は
`A⁰ = 1 ∣ P` で自明、`n = m+1` は本ブロックと第 239 で埋まる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalA` / `eval_evalA` | 分子を `u` の多項式と見る |
| `toPP2` / `univA_pow_dvd_iff` | ★★★★`ℤ[A,W] ≅ ℤ[w][a]` |
| `coeff_evalA_eq_zero` | ★★★★★★★★`A` 側の低次係数も消える |
| `univA_pow_dvd_numerator` | ★★★★★★★★分子は `A^{n+1}` でも割れる |
| `tateDefect_eq_zero` | ★★★★★★★★★**方程式の差は 0** |
| `tate_equation` | ★★★★★★★★★★**葉 (b) 完成** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★分子を `u` の多項式と見る -/

/-- ★★`W ↦ w`(定数)、`A ↦ X` として `ℤ[A,W]` を `ℂ[a]` に送る。 -/
noncomputable def evalA (w : ℂ) : TateBase →+* Polynomial ℂ :=
  MvPolynomial.eval₂Hom ((Polynomial.C : ℂ →+* Polynomial ℂ).comp (Int.castRingHom ℂ))
    ![Polynomial.X, Polynomial.C w]

theorem evalA_univA (w : ℂ) : evalA w univA = Polynomial.X := by
  rw [evalA, univA, MvPolynomial.eval₂Hom_X']
  rfl

theorem evalA_univW (w : ℂ) : evalA w univW = Polynomial.C w := by
  rw [evalA, univW, MvPolynomial.eval₂Hom_X']
  rfl

theorem eval_evalA (u w : ℂ) (P : TateBase) :
    (evalA w P).eval u = tateEval u w P := by
  have hhom : (Polynomial.evalRingHom u).comp (evalA w) = tateEval u w := by
    refine MvPolynomial.ringHom_ext ?_ ?_
    · intro r
      simp [evalA, tateEval]
    · intro i
      fin_cases i
      · simp [evalA, tateEval, MvPolynomial.eval₂Hom_X']
      · simp [evalA, tateEval, MvPolynomial.eval₂Hom_X']
  exact congrArg (fun f => f P) hhom

/-! ## ★★★★`ℤ[A,W] ≅ ℤ[w][a]`(`A` を外側に) -/

noncomputable def toPP2 : TateBase →+* Polynomial (Polynomial ℤ) :=
  MvPolynomial.eval₂Hom
    ((Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ)).comp
      (Polynomial.C : ℤ →+* Polynomial ℤ))
    ![Polynomial.X, Polynomial.C Polynomial.X]

noncomputable def ofPP2 : Polynomial (Polynomial ℤ) →+* TateBase :=
  Polynomial.eval₂RingHom
    (Polynomial.eval₂RingHom (MvPolynomial.C : ℤ →+* TateBase) univW) univA

theorem toPP2_univA : toPP2 univA = Polynomial.X := by
  rw [toPP2, univA, MvPolynomial.eval₂Hom_X']
  rfl

theorem toPP2_univW : toPP2 univW = Polynomial.C Polynomial.X := by
  rw [toPP2, univW, MvPolynomial.eval₂Hom_X']
  rfl

theorem ofPP2_toPP2 : ofPP2.comp toPP2 = RingHom.id TateBase := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [toPP2, ofPP2]
  · intro i
    fin_cases i
    · show ofPP2 (toPP2 univA) = univA
      rw [toPP2_univA, ofPP2]
      simp
    · show ofPP2 (toPP2 univW) = univW
      rw [toPP2_univW, ofPP2]
      simp

/-- ★★★★**`A^m ∣ P` は `X^m ∣ toPP2 P` と同じ**。 -/
theorem univA_pow_dvd_iff (m : ℕ) (P : TateBase) :
    univA ^ m ∣ P ↔ Polynomial.X ^ m ∣ toPP2 P := by
  constructor
  · rintro ⟨c, hc⟩
    exact ⟨toPP2 c, by rw [hc, map_mul, map_pow, toPP2_univA]⟩
  · rintro ⟨g, hg⟩
    refine ⟨ofPP2 g, ?_⟩
    have h := congrArg ofPP2 hg
    rw [map_mul, map_pow] at h
    have hA : ofPP2 Polynomial.X = univA := by
      rw [ofPP2]
      simp
    rw [hA] at h
    have hid : ofPP2 (toPP2 P) = P := congrArg (fun f => f P) ofPP2_toPP2
    rw [hid] at h
    exact h

theorem evalA_eq_map_toPP2 (w : ℂ) :
    evalA w = (Polynomial.mapRingHom (Polynomial.eval₂RingHom (Int.castRingHom ℂ) w)).comp
      toPP2 := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [evalA, toPP2]
  · intro i
    fin_cases i
    · show evalA w univA = _
      rw [evalA_univA]
      show _ = (Polynomial.mapRingHom _) (toPP2 univA)
      rw [toPP2_univA]
      simp
    · show evalA w univW = _
      rw [evalA_univW]
      show _ = (Polynomial.mapRingHom _) (toPP2 univW)
      rw [toPP2_univW]
      simp

theorem coeff_evalA_eq (w : ℂ) (P : TateBase) (i : ℕ) :
    (evalA w P).coeff i
      = Polynomial.eval₂ (Int.castRingHom ℂ) w ((toPP2 P).coeff i) := by
  have h : evalA w P
      = (Polynomial.mapRingHom (Polynomial.eval₂RingHom (Int.castRingHom ℂ) w)) (toPP2 P) :=
    congrArg (fun f : TateBase →+* Polynomial ℂ => f P) (evalA_eq_map_toPP2 w)
  rw [h]
  simp [Polynomial.coeff_map]

set_option maxHeartbeats 400000 in
theorem univA_pow_dvd_of_coeff (n : ℕ) (P : TateBase)
    (h : ∀ z : UpperHalfPlane, ∀ i, i < n + 1 →
      (evalA (Complex.exp (2 * ↑π * I * (z : ℂ))) P).coeff i = 0) :
    univA ^ (n + 1) ∣ P := by
  rw [univA_pow_dvd_iff, Polynomial.X_pow_dvd_iff]
  intro i hi
  have hgc : ((toPP2 P).coeff i).map (Int.castRingHom ℂ) = 0 := by
    refine map_eq_zero_of_exp_eval _ fun z => ?_
    rw [← coeff_evalA_eq]
    exact h z i hi
  have hcast : Function.Injective ((Int.castRingHom ℂ) : ℤ →+* ℂ) := by
    intro a b hab
    simpa using hab
  have hinj : Function.Injective (Polynomial.map (Int.castRingHom ℂ)) :=
    Polynomial.map_injective _ hcast
  have h0 : ((0 : Polynomial ℤ).map (Int.castRingHom ℂ)) = 0 := by simp
  exact hinj (by rw [hgc, h0])

/-! ## ★★★★★★★★`A` 側の係数の消滅 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**分子の `A` 側の低次係数も消える**——`w` を固定して `u → 0` を見る。 -/
theorem coeff_evalA_eq_zero (z : UpperHalfPlane) (n : ℕ) (P : TateBase) (d : tateDenoms)
    (hPd : tateDefectTrunc (n + 1) uA uW (uA * uW)
      * algebraMap TateBase TateUniv (d : TateBase) = algebraMap TateBase TateUniv P)
    (i : ℕ) (hi : i < n + 1) :
    (evalA (Complex.exp (2 * ↑π * I * (z : ℂ))) P).coeff i = 0 := by
  set w := Complex.exp (2 * ↑π * I * (z : ℂ)) with hwdef
  have hw : ‖w‖ < 1 := UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have hw1 : ‖w‖ ≤ 1 := hw.le
  obtain ⟨B, hB0, hB⟩ := exists_bound_tateEval (d : TateBase)
  set Cn : ℝ := 3 * (‖tateXterm w‖ + ‖tateYterm w‖ + 30 * ((n : ℝ) + 1) + 80) ^ 2
      + 6 * (‖tateXterm w‖ + ‖tateYterm w‖ + 30 * ((n : ℝ) + 1) + 80) + 1 with hCndef
  have hCn0 : (0 : ℝ) ≤ Cn := by
    simp only [hCndef]
    have h1 := norm_nonneg (tateXterm w)
    have h2 := norm_nonneg (tateYterm w)
    have h3 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
    positivity
  refine coeff_eq_zero_of_norm_le _ (Cn * 50 * 64 ^ (n + 1) * B) (1 / 128) (by norm_num)
    (n + 1) ?_ i hi
  intro u hu0 hur
  have hu2 : ‖u‖ ≤ 1 / 2 := by linarith
  have hu1 : ‖u‖ < 1 := by linarith
  have hwu : ‖w * u‖ ≤ 1 / 128 := by
    rw [norm_mul]
    nlinarith [norm_nonneg u, norm_nonneg w]
  have hnum := tateEval_numerator (n + 1) P d hPd u w hu1 hw
  have hsym : tateDefectTrunc (n + 1) u w (u * w) = tateDefectTrunc (n + 1) w u (w * u) := by
    rw [show w * u = u * w by ring]
    exact (tateDefectTrunc_symm (n + 1) u w (u * w)).symm
  have hbound := norm_tateDefectTrunc_uw_le z u hu0 hu2 hwu n
  have hd := hB u w hu1.le hw1
  rw [eval_evalA, hnum, hsym, norm_mul]
  have hmono : ((64 : ℝ) * ‖w * u‖) ^ (n + 1) ≤ 64 ^ (n + 1) * ‖u‖ ^ (n + 1) := by
    have hle : (64 : ℝ) * ‖w * u‖ ≤ 64 * ‖u‖ := by
      rw [norm_mul]
      nlinarith [norm_nonneg u, norm_nonneg w]
    calc ((64 : ℝ) * ‖w * u‖) ^ (n + 1) ≤ ((64 : ℝ) * ‖u‖) ^ (n + 1) :=
          pow_le_pow_left₀ (by positivity) hle _
      _ = 64 ^ (n + 1) * ‖u‖ ^ (n + 1) := by rw [mul_pow]
  have hupow : (0 : ℝ) ≤ ‖u‖ ^ (n + 1) := by positivity
  calc ‖tateDefectTrunc (n + 1) w u (w * u)‖ * ‖tateEval u w (d : TateBase)‖
      ≤ (Cn * (50 * ((64 : ℝ) * ‖w * u‖) ^ (n + 1))) * B :=
        mul_le_mul hbound hd (norm_nonneg _) (by positivity)
    _ ≤ (Cn * (50 * (64 ^ (n + 1) * ‖u‖ ^ (n + 1)))) * B := by
        have hstep : Cn * (50 * ((64 : ℝ) * ‖w * u‖) ^ (n + 1))
            ≤ Cn * (50 * (64 ^ (n + 1) * ‖u‖ ^ (n + 1))) := by
          nlinarith [hCn0, hmono]
        nlinarith [hB0, hstep]
    _ = Cn * 50 * 64 ^ (n + 1) * B * ‖u‖ ^ (n + 1) := by ring

/-- ★★★★★★★★**分子は `A^{n+1}` でも割れる**。 -/
theorem univA_pow_dvd_numerator (n : ℕ) (P : TateBase) (d : tateDenoms)
    (hPd : tateDefectTrunc (n + 1) uA uW (uA * uW)
      * algebraMap TateBase TateUniv (d : TateBase) = algebraMap TateBase TateUniv P) :
    univA ^ (n + 1) ∣ P :=
  univA_pow_dvd_of_coeff n P (fun z i hi => coeff_evalA_eq_zero z n P d hPd i hi)

/-! ## ★★★★★★★★★★葉 (b) 完成 -/

section Final

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★★**Weierstrass 方程式の差は 0 である**。 -/
theorem tateDefect_eq_zero [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDefect a w q hq = 0 := by
  refine tateDefect_eq_zero_of_base a w q hq haw ha hw ?_ ?_
  · intro n P d hPd
    cases n with
    | zero => simp
    | succ m => exact univA_pow_dvd_numerator m P d hPd
  · intro n P d hPd
    cases n with
    | zero => simp
    | succ m => exact univW_pow_dvd_numerator m P d hPd

/-- ★★★★★★★★★★**葉 (b) 完成 —— Tate 級数は Weierstrass 方程式を満たす**。

    Y(u,q)² + X(u,q)·Y(u,q) = X(u,q)³ + a₄(q)·X(u,q) + a₆(q)

★第 212 で `I` を法として取れていたものが、**厳密な等式**になった。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_equation [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateYpair a w q hq ^ 2 + tateXpair a w q hq * tateYpair a w q hq
      = tateXpair a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
        + (tateCurveAt q hq).a₆ :=
  tate_equation_of_defect a w q hq (tateDefect_eq_zero a w q hq haw ha hw)

end Final

/-! ## ★出典の紐付け(`.src`) -/

def univA_pow_dvd_numerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子は A^{n+1} でも割れる)",
    sectionId := "genell-def-3-3" }

def tateDefect_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——方程式の差は 0)",
    sectionId := "genell-def-3-3" }

def tate_equation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (b) Weierstrass 方程式)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
