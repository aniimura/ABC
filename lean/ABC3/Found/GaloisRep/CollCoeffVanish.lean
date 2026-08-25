import ABC3.Found.GaloisRep.CollBound

/-!
# Galois (G6) 第 255 ブロック —— **★★★★★★★★分子は `W^{n+1}` で割れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★ℂ への特殊化 → 分子の値

第 225 の `tateSpecializeC`・`tateEval_numerator` を 3 変数にしたもの。
`‖u‖, ‖v‖, ‖w‖ < 1` なら分母はすべて ℂ で単元なので、局所化から ℂ への射が作れる。
そこから

    collEval u v w P = collDefectTrunc (n+1) u v w (uvw) · collEval u v w d

が出る(`collEval_numerator`)。★**分子の値が「切り詰めた差 × 分母」**という形が、
第 254 の解析的評価と第 236 の `X_pow_dvd_of_norm_le` をつなぐ橋である。

## ★★★★★★★★係数の消滅

`u, v` を `‖·‖ ≤ 1/8` に固定して `w` を動かすと、第 254 の評価は

    ‖collDefectTrunc (n+1) u v w (uvw)‖ ≤ Cₙ · 4^{n+1} ‖w‖^{n+1}

となる(`‖uvw‖ ≤ ‖w‖` を使う)。★`coeff_eq_zero_of_norm_le`(第 236)より
`evalUV u v P` の `j < n+1` 次の係数は消える。

## ★動かす先は「小さい点」の無限集合でよい

第 239 では上半平面の指数像を使ったが、ここでは **`‖·‖ ≤ 1/8` である必要**がある
(3 変数とも動くので `tateXterm` を数値で潰した——§9-567)。
★`smallSet := {1/(k+9) : k ∈ ℕ}` を取れば、無限で、0 でなく、ノルムが `1/9` 以下。
**上半平面を経由しなくてよい**ぶん簡単になった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `collSpecializeC` | ★★★★★★万有な環から ℂ への特殊化 |
| `collEval_numerator` | ★★★★★★**分子の値 = 切り詰めた差 × 分母** |
| `exists_bound_collEval` | ★★★★分母の有界性 |
| `coeff_evalUV_eq_zero` | ★★★★★★★★**分子の低次係数は消える** |
| `smallSet` 他 | ★動かす先の無限集合 |
| `cW_pow_dvd_numerator` | ★★★★★★★★**分子は `W^{n+1}` で割れる** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★★★ℂ への特殊化 -/

theorem norm_collPts_lt_one (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1)
    (i : Fin 6) : ‖collPts u v w i‖ < 1 := by
  have hprod : ∀ a b : ℂ, ‖a‖ < 1 → ‖b‖ < 1 → ‖a * b‖ < 1 := fun a b ha hb =>
    norm_mul_lt_one ha hb
  fin_cases i
  · exact hu
  · exact hv
  · exact hw
  · exact hprod v w hv hw
  · exact hprod u w hu hw
  · exact hprod u v hu hv

theorem collEval_isUnit_denomSet_complex (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1)
    (hw : ‖w‖ < 1) {x : CollBase} (hx : x ∈ collDenomSet) : IsUnit (collEval u v w x) := by
  have hq : ‖u * v * w‖ < 1 := norm_mul_lt_one (norm_mul_lt_one hu hv) hw
  rcases hx with ⟨i, rfl⟩ | ⟨⟨m, i⟩, rfl⟩
  · rw [map_sub, map_one, collEval_denomBases]
    exact isUnit_one_sub_of_norm_lt_one (norm_collPts_lt_one u v w hu hv hw i)
  · rw [map_sub, map_one, map_mul, map_pow, collEval_Q, collEval_denomBases]
    refine isUnit_one_sub_of_norm_lt_one ?_
    rw [norm_mul, norm_pow]
    have h1 : ‖u * v * w‖ ^ (m + 1) ≤ 1 := pow_le_one₀ (norm_nonneg _) hq.le
    have h2 := norm_collPts_lt_one u v w hu hv hw i
    nlinarith [norm_nonneg (collPts u v w i), pow_nonneg (norm_nonneg (u * v * w)) (m + 1)]

theorem collEval_isUnit_denoms_complex (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1)
    (hw : ‖w‖ < 1) (y : collDenoms) : IsUnit (collEval u v w (y : CollBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := collDenomSet)
    (motive := fun z _ => IsUnit (collEval u v w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact collEval_isUnit_denomSet_complex u v w hu hv hw hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-- ★★★★★★**万有な環から ℂ への特殊化**——`‖u‖,‖v‖,‖w‖ < 1` だけで足りる。 -/
noncomputable def collSpecializeC (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    CollUniv →+* ℂ :=
  IsLocalization.lift (collEval_isUnit_denoms_complex u v w hu hv hw)

theorem collSpecializeC_base (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1)
    (x : CollBase) :
    collSpecializeC u v w hu hv hw (algebraMap CollBase CollUniv x) = collEval u v w x :=
  IsLocalization.lift_eq _ x

theorem collSpecializeC_kU (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    collSpecializeC u v w hu hv hw kU = u := by rw [kU, collSpecializeC_base, collEval_U]

theorem collSpecializeC_kV (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    collSpecializeC u v w hu hv hw kV = v := by rw [kV, collSpecializeC_base, collEval_V]

theorem collSpecializeC_kW (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    collSpecializeC u v w hu hv hw kW = w := by rw [kW, collSpecializeC_base, collEval_W]

/-- ★★★★★ℂ 側での切り詰めた差の値。 -/
theorem collSpecializeC_collDefectTrunc (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1)
    (hw : ‖w‖ < 1) (n : ℕ) :
    collSpecializeC u v w hu hv hw (collDefectTrunc n kU kV kW (kU * kV * kW))
      = collDefectTrunc n u v w (u * v * w) := by
  have hcu := collUnits_k
  have h := map_collDefectTrunc (collSpecializeC u v w hu hv hw) n kU kV kW (kU * kV * kW)
    hcu.hu hcu.hv hcu.hw hcu.hvw hcu.huw hcu.huv
    (fun m _ => hcu.hqu m) (fun m _ => hcu.hqv m) (fun m _ => hcu.hqw m)
    (fun m _ => hcu.hqvw m) (fun m _ => hcu.hquw m) (fun m _ => hcu.hquv m)
  rw [h, collSpecializeC_kU, collSpecializeC_kV, collSpecializeC_kW, map_mul, map_mul,
    collSpecializeC_kU, collSpecializeC_kV, collSpecializeC_kW]

/-- ★★★★★★**分子の値は「切り詰めた差 × 分母」である**——解析の段への橋。 -/
theorem collEval_numerator (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc n kU kV kW (kU * kV * kW) * algebraMap CollBase CollUniv (d : CollBase)
      = algebraMap CollBase CollUniv P)
    (u v w : ℂ) (hu : ‖u‖ < 1) (hv : ‖v‖ < 1) (hw : ‖w‖ < 1) :
    collEval u v w P
      = collDefectTrunc n u v w (u * v * w) * collEval u v w (d : CollBase) := by
  have h := congrArg (collSpecializeC u v w hu hv hw) hPd
  rw [map_mul, collSpecializeC_base, collSpecializeC_base,
    collSpecializeC_collDefectTrunc] at h
  exact h.symm

/-! ## ★★★★分母の有界性 -/

theorem collEval_C (u v w : ℂ) (r : ℤ) : collEval u v w (MvPolynomial.C r) = (r : ℂ) := by
  rw [collEval, MvPolynomial.eval₂Hom_C]
  rfl

/-- ★★★★**`‖u‖,‖v‖,‖w‖ ≤ 1` の上で `collEval u v w Q` は有界**。 -/
theorem exists_bound_collEval (Q : CollBase) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ u v w : ℂ, ‖u‖ ≤ 1 → ‖v‖ ≤ 1 → ‖w‖ ≤ 1 →
      ‖collEval u v w Q‖ ≤ B := by
  refine MvPolynomial.induction_on (motive := fun p =>
    ∃ B : ℝ, 0 ≤ B ∧ ∀ u v w : ℂ, ‖u‖ ≤ 1 → ‖v‖ ≤ 1 → ‖w‖ ≤ 1 →
      ‖collEval u v w p‖ ≤ B) Q ?_ ?_ ?_
  · intro r
    refine ⟨‖((r : ℤ) : ℂ)‖, norm_nonneg _, fun u v w _ _ _ => ?_⟩
    rw [collEval_C]
  · rintro p q ⟨Bp, hBp0, hBp⟩ ⟨Bq, hBq0, hBq⟩
    refine ⟨Bp + Bq, by linarith, fun u v w hu hv hw => ?_⟩
    rw [map_add]
    refine (norm_add_le _ _).trans ?_
    have h1 := hBp u v w hu hv hw
    have h2 := hBq u v w hu hv hw
    linarith
  · rintro p i ⟨Bp, hBp0, hBp⟩
    refine ⟨Bp, hBp0, fun u v w hu hv hw => ?_⟩
    rw [map_mul, norm_mul]
    have hxi : ‖collEval u v w (MvPolynomial.X i)‖ ≤ 1 := by
      fin_cases i
      · show ‖collEval u v w cU‖ ≤ 1
        rw [collEval_U]; exact hu
      · show ‖collEval u v w cV‖ ≤ 1
        rw [collEval_V]; exact hv
      · show ‖collEval u v w cW‖ ≤ 1
        rw [collEval_W]; exact hw
    have hp := hBp u v w hu hv hw
    nlinarith [norm_nonneg (collEval u v w p), norm_nonneg (collEval u v w (MvPolynomial.X i))]

/-! ## ★★★★★★★★分子の低次係数は消える -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**分子の低次係数は消える**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem coeff_evalUV_eq_zero (u v : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0)
    (hu : ‖u‖ ≤ 1 / 8) (hv : ‖v‖ ≤ 1 / 8) (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P)
    (j : ℕ) (hj : j < n + 1) :
    (evalUV u v P).coeff j = 0 := by
  obtain ⟨B, hB0, hB⟩ := exists_bound_collEval (d : CollBase)
  have hu1 : ‖u‖ ≤ 1 := by linarith
  have hv1 : ‖v‖ ≤ 1 := by linarith
  set Cn : ℝ := 12 * (25 * ((n : ℝ) + 1) + 8) * 50 with hCndef
  have hn0 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  have hCn0 : (0 : ℝ) ≤ Cn := by simp only [hCndef]; nlinarith
  refine coeff_eq_zero_of_norm_le _ (Cn * 4 ^ (n + 1) * B) (1 / 8) (by norm_num)
    (n + 1) ?_ j hj
  intro w hw0 hwr
  have hw : ‖w‖ ≤ 1 / 8 := hwr.le
  have hw1 : ‖w‖ ≤ 1 := by linarith
  have hnum := collEval_numerator (n + 1) P d hPd u v w (by linarith) (by linarith) (by linarith)
  have hbound := norm_collDefectTrunc_le u v w hu0 hv0 hw0 hu hv hw n
  have hd := hB u v w hu1 hv1 hw1
  rw [eval_evalUV, hnum, norm_mul]
  have hqle : ‖u * v * w‖ ≤ ‖w‖ := by
    rw [norm_mul, norm_mul]
    have hnu := norm_nonneg u
    have hnv := norm_nonneg v
    have hnw := norm_nonneg w
    have h1 : ‖u‖ * ‖v‖ ≤ 1 := by nlinarith
    nlinarith
  have hmono : (4 * ‖u * v * w‖) ^ (n + 1) ≤ 4 ^ (n + 1) * ‖w‖ ^ (n + 1) := by
    have hle : (4 : ℝ) * ‖u * v * w‖ ≤ 4 * ‖w‖ := by linarith
    calc (4 * ‖u * v * w‖) ^ (n + 1) ≤ ((4 : ℝ) * ‖w‖) ^ (n + 1) :=
          pow_le_pow_left₀ (by positivity) hle _
      _ = 4 ^ (n + 1) * ‖w‖ ^ (n + 1) := by rw [mul_pow]
  have hwpow : (0 : ℝ) ≤ ‖w‖ ^ (n + 1) := by positivity
  have hstep : ‖collDefectTrunc (n + 1) u v w (u * v * w)‖
      ≤ Cn * (4 ^ (n + 1) * ‖w‖ ^ (n + 1)) := by
    simp only [hCndef] at hbound ⊢
    nlinarith [hbound, hmono, hn0]
  calc ‖collDefectTrunc (n + 1) u v w (u * v * w)‖ * ‖collEval u v w (d : CollBase)‖
      ≤ (Cn * (4 ^ (n + 1) * ‖w‖ ^ (n + 1))) * B :=
        mul_le_mul hstep hd (norm_nonneg _) (by positivity)
    _ = Cn * 4 ^ (n + 1) * B * ‖w‖ ^ (n + 1) := by ring

/-! ## ★動かす先の無限集合 -/

/-- ★小さい点の無限集合——`1/(k+9)`。★上半平面を経由しなくてよい。 -/
noncomputable def smallSet : Set ℂ := Set.range (fun k : ℕ => (1 : ℂ) / ((k : ℂ) + 9))

theorem cast_add_nine_ne_zero (k : ℕ) : ((k : ℂ) + 9) ≠ 0 := by
  have h : ((k : ℂ) + 9) = ((k + 9 : ℕ) : ℂ) := by push_cast; ring
  rw [h]
  exact Nat.cast_ne_zero.2 (by omega)

theorem infinite_smallSet : smallSet.Infinite := by
  refine Set.infinite_range_of_injective ?_
  intro a b hab
  simp only at hab
  have ha := cast_add_nine_ne_zero a
  have hb := cast_add_nine_ne_zero b
  field_simp at hab
  have h2 : ((b : ℂ)) = ((a : ℂ)) := by linear_combination hab
  have h3 : b = a := by exact_mod_cast h2
  omega

theorem smallSet_ne_zero {t : ℂ} (ht : t ∈ smallSet) : t ≠ 0 := by
  obtain ⟨k, rfl⟩ := ht
  simp only [ne_eq, div_eq_zero_iff, one_ne_zero, false_or]
  exact cast_add_nine_ne_zero k

theorem smallSet_norm_le {t : ℂ} (ht : t ∈ smallSet) : ‖t‖ ≤ 1 / 8 := by
  obtain ⟨k, rfl⟩ := ht
  have h : ((k : ℂ) + 9) = ((k + 9 : ℕ) : ℂ) := by push_cast; ring
  simp only [h, norm_div, norm_one, Complex.norm_natCast]
  have h9 : (8 : ℝ) ≤ ((k + 9 : ℕ) : ℝ) := by
    have he : ((k + 9 : ℕ) : ℝ) = (k : ℝ) + 9 := by push_cast; ring
    rw [he]
    have : (0 : ℝ) ≤ (k : ℝ) := Nat.cast_nonneg k
    linarith
  exact one_div_le_one_div_of_le (by norm_num) h9

/-! ## ★★★★★★★★分子の整除性 -/

/-- ★★★★★★★★**分子は `W^{n+1}` で割れる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem cW_pow_dvd_numerator (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P) :
    cW ^ (n + 1) ∣ P :=
  cW_pow_dvd_of_coeff n P smallSet infinite_smallSet
    (fun u hu v hv j hj => coeff_evalUV_eq_zero u v (smallSet_ne_zero hu) (smallSet_ne_zero hv)
      (smallSet_norm_le hu) (smallSet_norm_le hv) n P d hPd j hj)

/-! ## ★出典の紐付け(`.src`) -/

def collEval_numerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子の値 = 切り詰めた差 × 分母)",
    sectionId := "genell-def-3-3" }

def cW_pow_dvd_numerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子は W^{n+1} で割れる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
