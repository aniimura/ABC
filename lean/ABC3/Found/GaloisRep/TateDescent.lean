import ABC3.Found.GaloisRep.TateCoeffVanish

/-!
# Galois (G6) 第 239 ブロック —— **★★★★★★★★分子は `W^{n+1}` で割れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★ℂ から ℤ へ降ろす

第 238 で「各 `u = exp(2πiz)` について `(evalW u P)` の低次係数が 0」まで来た。
これを **`ℤ[A,W]` の中の `W^{n+1} ∣ P`** に降ろす。

★★★★★**`ℤ[A,W] ≅ ℤ[a][w]`(`W` を外側の変数に)を作る**のが鍵である
(`toPP` / `ofPP`)。すると

    W^m ∣ P  ⟺  X^m ∣ toPP P  ⟺  ∀ j < m, (toPP P).coeff j = 0

★`(toPP P).coeff j ∈ ℤ[a]` は `u` の多項式である。`evalW u` はちょうど
その係数を `u` で評価したものになる(`coeff_evalW_eq`)。

## ★★★`u` の動く先は無限集合

`z ∈ ℍ` を動かすと `exp(2πiz)` は無限個の値を取る——`z = i(k+1)` とすれば
`‖exp(2πiz)‖ = exp(−2π(k+1))` がすべて相異なる(`infinite_exp_range`)。
★これで「無限個の点で消える多項式は 0」(第 236)が使え、`ℂ[a]` の側で 0 になる。
★★`ℤ[a] → ℂ[a]` は単射なので `ℤ[a]` の側でも 0 ✓。

## ★★配線の落とし穴

`poly_eq_zero_of_infinite_zeros _ (Set.range fun z : ℍ => …) infinite_exp_range` と
**集合を明示すると `isDefEq` が時間切れになる**。`_ _ infinite_exp_range` と書いて
補題側から推論させると一瞬で通る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `toPP` / `ofPP` / `univW_pow_dvd_iff` | ★★★★★`ℤ[A,W] ≅ ℤ[a][w]` と整除性の言い換え |
| `coeff_evalW_eq` | ★★★★`evalW u` は係数を `u` で評価したもの |
| `infinite_exp_range` | ★★★`exp(2πi ℍ)` の像は無限 |
| `univW_pow_dvd_of_coeff` | ★★★★★★★`u` を動かして `W^{n+1} ∣ P` |
| `univW_pow_dvd_numerator` | ★★★★★★★★**分子は `W^{n+1}` で割れる** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★★`ℤ[A,W] ≅ ℤ[a][w]` -/

/-- ★★★★`ℤ[A,W] → ℤ[a][w]`(`W` を外側の変数に)。 -/
noncomputable def toPP : TateBase →+* Polynomial (Polynomial ℤ) :=
  MvPolynomial.eval₂Hom
    ((Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ)).comp
      (Polynomial.C : ℤ →+* Polynomial ℤ))
    ![Polynomial.C Polynomial.X, Polynomial.X]

/-- ★逆向き。 -/
noncomputable def ofPP : Polynomial (Polynomial ℤ) →+* TateBase :=
  Polynomial.eval₂RingHom
    (Polynomial.eval₂RingHom (MvPolynomial.C : ℤ →+* TateBase) univA) univW

theorem toPP_univA : toPP univA = Polynomial.C Polynomial.X := by
  rw [toPP, univA, MvPolynomial.eval₂Hom_X']
  rfl

theorem toPP_univW : toPP univW = Polynomial.X := by
  rw [toPP, univW, MvPolynomial.eval₂Hom_X']
  rfl

theorem ofPP_toPP : ofPP.comp toPP = RingHom.id TateBase := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [toPP, ofPP]
  · intro i
    fin_cases i
    · show ofPP (toPP univA) = univA
      rw [toPP_univA, ofPP]
      simp
    · show ofPP (toPP univW) = univW
      rw [toPP_univW, ofPP]
      simp

theorem toPP_ofPP : toPP.comp ofPP = RingHom.id (Polynomial (Polynomial ℤ)) := by
  refine Polynomial.ringHom_ext ?_ ?_
  · intro a
    have hinner : (toPP.comp ofPP).comp
        (Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ))
        = (Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ)) := by
      refine Polynomial.ringHom_ext ?_ ?_
      · intro r
        simp [toPP, ofPP]
      · show toPP (ofPP (Polynomial.C Polynomial.X)) = Polynomial.C Polynomial.X
        rw [ofPP]
        simp [toPP_univA]
    exact congrArg (fun f => f a) hinner
  · show toPP (ofPP Polynomial.X) = Polynomial.X
    rw [ofPP]
    simp [toPP_univW]

/-- ★★★★**`W^m ∣ P` は `X^m ∣ toPP P` と同じ**。 -/
theorem univW_pow_dvd_iff (m : ℕ) (P : TateBase) :
    univW ^ m ∣ P ↔ Polynomial.X ^ m ∣ toPP P := by
  constructor
  · rintro ⟨c, hc⟩
    exact ⟨toPP c, by rw [hc, map_mul, map_pow, toPP_univW]⟩
  · rintro ⟨g, hg⟩
    refine ⟨ofPP g, ?_⟩
    have h := congrArg ofPP hg
    rw [map_mul, map_pow] at h
    have hW : ofPP Polynomial.X = univW := by
      rw [ofPP]
      simp
    rw [hW] at h
    have hid : ofPP (toPP P) = P := congrArg (fun f => f P) ofPP_toPP
    rw [hid] at h
    exact h

/-! ## ★★★★`evalW` と `toPP` の関係 -/

theorem evalW_eq_map_toPP (u : ℂ) :
    evalW u = (Polynomial.mapRingHom (Polynomial.eval₂RingHom (Int.castRingHom ℂ) u)).comp
      toPP := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [evalW, toPP]
  · intro i
    fin_cases i
    · show evalW u univA = _
      rw [evalW_univA]
      show _ = (Polynomial.mapRingHom _) (toPP univA)
      rw [toPP_univA]
      simp
    · show evalW u univW = _
      rw [evalW_univW]
      show _ = (Polynomial.mapRingHom _) (toPP univW)
      rw [toPP_univW]
      simp

/-- ★★★★**`evalW u` は `toPP` の係数を `u` で評価したもの**。 -/
theorem coeff_evalW_eq (u : ℂ) (P : TateBase) (j : ℕ) :
    (evalW u P).coeff j
      = Polynomial.eval₂ (Int.castRingHom ℂ) u ((toPP P).coeff j) := by
  have h : evalW u P
      = (Polynomial.mapRingHom (Polynomial.eval₂RingHom (Int.castRingHom ℂ) u)) (toPP P) :=
    congrArg (fun f : TateBase →+* Polynomial ℂ => f P) (evalW_eq_map_toPP u)
  rw [h]
  simp [Polynomial.coeff_map]

/-! ## ★★★`u` の動く先は無限集合 -/

/-- ★★★**`exp(2πi ℍ)` の像は無限集合**——`z = i(k+1)` でノルムがすべて相異なる。 -/
theorem infinite_exp_range :
    (Set.range (fun z : UpperHalfPlane => Complex.exp (2 * ↑π * I * (z : ℂ)))).Infinite := by
  have hpi : (0 : ℝ) < π := Real.pi_pos
  refine Set.infinite_of_injective_forall_mem
    (f := fun k : ℕ => Complex.exp (2 * (π : ℂ) * I * (I * ((k : ℂ) + 1)))) ?_ ?_
  · intro a b hab
    have hnorm : ∀ k : ℕ, ‖Complex.exp (2 * (π : ℂ) * I * (I * ((k : ℂ) + 1)))‖
        = Real.exp (-(2 * π * ((k : ℝ) + 1))) := by
      intro k
      rw [Complex.norm_exp]
      congr 1
      have hz : 2 * (π : ℂ) * I * (I * ((k : ℂ) + 1))
          = -(2 * (π : ℂ) * ((k : ℂ) + 1)) := by
        have hI : (I : ℂ) * I = -1 := by
          rw [← sq, Complex.I_sq]
        linear_combination (2 * (π : ℂ) * ((k : ℂ) + 1)) * hI
      rw [hz]
      simp
    have h2 := congrArg (fun x : ℂ => ‖x‖) hab
    simp only [hnorm] at h2
    have h3 := Real.exp_injective h2
    have h5 : (2 * π) * ((a : ℝ) + 1) = (2 * π) * ((b : ℝ) + 1) := by linarith
    have h4 := mul_left_cancel₀ (by positivity : (2 * π : ℝ) ≠ 0) h5
    have h6 : (a : ℝ) = (b : ℝ) := by linarith
    exact_mod_cast h6
  · intro k
    have him : 0 < (Complex.I * ((k : ℂ) + 1)).im := by
      simp
      positivity
    exact ⟨⟨Complex.I * ((k : ℂ) + 1), him⟩, rfl⟩

set_option maxHeartbeats 400000 in
/-- ★★★★★`exp(2πiz)` のすべてで消える `ℤ` 係数多項式は ℂ に送って 0。

★`poly_eq_zero_of_infinite_zeros` に集合を**明示すると `isDefEq` が時間切れになる**。
`_ _ infinite_exp_range` と書いて補題側から推論させること。 -/
theorem map_eq_zero_of_exp_eval (g : Polynomial ℤ)
    (h : ∀ z : UpperHalfPlane,
      Polynomial.eval₂ (Int.castRingHom ℂ) (Complex.exp (2 * ↑π * I * (z : ℂ))) g = 0) :
    g.map (Int.castRingHom ℂ) = 0 := by
  refine poly_eq_zero_of_infinite_zeros _ _ infinite_exp_range ?_
  intro u hu
  obtain ⟨z, hz⟩ := hu
  rw [← hz, Polynomial.eval_map]
  exact h z

/-! ## ★★★★★★★★分子の整除性 -/

set_option maxHeartbeats 400000 in
/-- ★★★★★★★**`u` を動かして `W^{n+1} ∣ P` を出す**。 -/
theorem univW_pow_dvd_of_coeff (n : ℕ) (P : TateBase)
    (h : ∀ z : UpperHalfPlane, ∀ j, j < n + 1 →
      (evalW (Complex.exp (2 * ↑π * I * (z : ℂ))) P).coeff j = 0) :
    univW ^ (n + 1) ∣ P := by
  rw [univW_pow_dvd_iff, Polynomial.X_pow_dvd_iff]
  intro j hj
  have hgc : ((toPP P).coeff j).map (Int.castRingHom ℂ) = 0 := by
    refine map_eq_zero_of_exp_eval _ fun z => ?_
    rw [← coeff_evalW_eq]
    exact h z j hj
  have hcast : Function.Injective ((Int.castRingHom ℂ) : ℤ →+* ℂ) := by
    intro a b hab
    simpa using hab
  have hinj : Function.Injective (Polynomial.map (Int.castRingHom ℂ)) :=
    Polynomial.map_injective _ hcast
  have h0 : ((0 : Polynomial ℤ).map (Int.castRingHom ℂ)) = 0 := by simp
  exact hinj (by rw [hgc, h0])

/-- ★★★★★★★★**分子は `W^{n+1}` で割れる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem univW_pow_dvd_numerator (n : ℕ) (P : TateBase) (d : tateDenoms)
    (hPd : tateDefectTrunc (n + 1) uA uW (uA * uW)
      * algebraMap TateBase TateUniv (d : TateBase) = algebraMap TateBase TateUniv P) :
    univW ^ (n + 1) ∣ P :=
  univW_pow_dvd_of_coeff n P (fun z j hj => coeff_evalW_eq_zero z n P d hPd j hj)

/-! ## ★出典の紐付け(`.src`) -/

def univW_pow_dvd_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Z[A,W] ≅ Z[a][w] と整除性)",
    sectionId := "genell-def-3-3" }

def univW_pow_dvd_numerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子は W^{n+1} で割れる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
