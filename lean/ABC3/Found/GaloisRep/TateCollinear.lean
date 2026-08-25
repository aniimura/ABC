import ABC3.Found.GaloisRep.CollCoeffVanish

/-!
# Galois (G6) 第 256 ブロック —— **★★★★★★★★★葉 (c) が終わった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> `u v w = q` なら Tate 級数の 3 点 `(X(u), Y(u))`, `(X(v), Y(v))`, `(X(w), Y(w))` は
> 一直線上にある(`tate_collinear`)。

これが**葉 (c)(一意化が準同型であること)の中身**である。

## ★★★★★変数の入れ替えで `U`・`V` 側を出す

第 255 で `W^{n+1} ∣ P` を示した。`U`・`V` 側は**変数の名前替え**で出る。

    swapUW := rename (Equiv.swap 0 2),   swapVW := rename (Equiv.swap 1 2)

これは `CollBase` の環同型(対合)で、`swapUW cW = cU` なので

    `W^m ∣ swapUW P`  ⟹  `U^m ∣ P`          (`cU_pow_dvd_of_swap`)

となる。★★そして

    (evalUV u v (swapUW P)).eval t = collEval t v u P

なので、**第 253 の `cW_pow_dvd_of_coeff` をそのまま `swapUW P` に当てられる**。
新しい多項式環の塔(`ℤ[v][w][u]` 等)を作る必要がない。

## ★★★★★★評価は 3 変数対称だった

第 254 の `norm_collDefectTrunc_le` は `(u,v,w)` について対称な形をしている。
そこで「動かす変数 `t` はどのスロットでもよい」形に一度だけ書き直しておく
(`norm_collEval_le`——仮定は `‖uvw‖ ≤ ‖t‖` だけ)。
★これで 3 つの側が**同じ 1 本の評価**から出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `norm_collEval_le` | ★★★★★★動かす変数はどのスロットでもよい |
| `swapUW`・`swapVW` | ★★★★★変数の入れ替え |
| `eval_evalUV_swapUW` 他 | ★★★★入れ替えたあとの評価 |
| `cU_pow_dvd_numerator`・`cV_pow_dvd_numerator` | ★★★★★★★★残る二つの整除性 |
| `tate_collinear` | ★★★★★★★★★**葉 (c) の到達点** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★★★動かす変数はどのスロットでもよい -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**分子の値の評価**——動かす変数 `t` はどのスロットでもよい。

★仮定は `‖uvw‖ ≤ ‖t‖` だけ。第 254 の評価が `(u,v,w)` 対称だから書ける形である。 -/
theorem norm_collEval_le (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P)
    (B : ℝ) (_hB0 : 0 ≤ B)
    (hB : ∀ u v w : ℂ, ‖u‖ ≤ 1 → ‖v‖ ≤ 1 → ‖w‖ ≤ 1 → ‖collEval u v w (d : CollBase)‖ ≤ B)
    (u v w t : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0) (hw0 : w ≠ 0)
    (hu : ‖u‖ ≤ 1 / 8) (hv : ‖v‖ ≤ 1 / 8) (hw : ‖w‖ ≤ 1 / 8)
    (ht : ‖u * v * w‖ ≤ ‖t‖) :
    ‖collEval u v w P‖
      ≤ 12 * (25 * ((n : ℝ) + 1) + 8) * 50 * 4 ^ (n + 1) * B * ‖t‖ ^ (n + 1) := by
  have hn0 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  set Cn : ℝ := 12 * (25 * ((n : ℝ) + 1) + 8) * 50 with hCndef
  have hCn0 : (0 : ℝ) ≤ Cn := by simp only [hCndef]; nlinarith
  have hnum := collEval_numerator (n + 1) P d hPd u v w (by linarith) (by linarith) (by linarith)
  have hbound := norm_collDefectTrunc_le u v w hu0 hv0 hw0 hu hv hw n
  have hd := hB u v w (by linarith) (by linarith) (by linarith)
  rw [hnum, norm_mul]
  have hmono : (4 * ‖u * v * w‖) ^ (n + 1) ≤ 4 ^ (n + 1) * ‖t‖ ^ (n + 1) := by
    have hle : (4 : ℝ) * ‖u * v * w‖ ≤ 4 * ‖t‖ := by linarith
    calc (4 * ‖u * v * w‖) ^ (n + 1) ≤ ((4 : ℝ) * ‖t‖) ^ (n + 1) :=
          pow_le_pow_left₀ (by positivity) hle _
      _ = 4 ^ (n + 1) * ‖t‖ ^ (n + 1) := by rw [mul_pow]
  have htpow : (0 : ℝ) ≤ ‖t‖ ^ (n + 1) := by positivity
  have hstep : ‖collDefectTrunc (n + 1) u v w (u * v * w)‖
      ≤ Cn * (4 ^ (n + 1) * ‖t‖ ^ (n + 1)) := by
    simp only [hCndef] at hbound ⊢
    nlinarith [hbound, hmono, hn0]
  calc ‖collDefectTrunc (n + 1) u v w (u * v * w)‖ * ‖collEval u v w (d : CollBase)‖
      ≤ (Cn * (4 ^ (n + 1) * ‖t‖ ^ (n + 1))) * B :=
        mul_le_mul hstep hd (norm_nonneg _) (by positivity)
    _ = Cn * 4 ^ (n + 1) * B * ‖t‖ ^ (n + 1) := by ring

/-! ## ★★★★★変数の入れ替え -/

/-- ★`U` と `W` を入れ替える環同型。 -/
noncomputable def swapUW : CollBase →+* CollBase :=
  (MvPolynomial.rename (Equiv.swap (0 : Fin 3) 2)).toRingHom

/-- ★`V` と `W` を入れ替える環同型。 -/
noncomputable def swapVW : CollBase →+* CollBase :=
  (MvPolynomial.rename (Equiv.swap (1 : Fin 3) 2)).toRingHom

theorem swapUW_X (i : Fin 3) : swapUW (MvPolynomial.X i) = MvPolynomial.X (Equiv.swap 0 2 i) :=
  MvPolynomial.rename_X _ i

theorem swapVW_X (i : Fin 3) : swapVW (MvPolynomial.X i) = MvPolynomial.X (Equiv.swap 1 2 i) :=
  MvPolynomial.rename_X _ i

theorem swapUW_comp : swapUW.comp swapUW = RingHom.id CollBase := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r; simp [swapUW]
  · intro i
    show swapUW (swapUW (MvPolynomial.X i)) = MvPolynomial.X i
    rw [swapUW_X, swapUW_X, Equiv.swap_apply_self]

theorem swapVW_comp : swapVW.comp swapVW = RingHom.id CollBase := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r; simp [swapVW]
  · intro i
    show swapVW (swapVW (MvPolynomial.X i)) = MvPolynomial.X i
    rw [swapVW_X, swapVW_X, Equiv.swap_apply_self]

theorem swapUW_cW : swapUW cW = cU := by
  show swapUW (MvPolynomial.X 2) = MvPolynomial.X 0
  rw [swapUW_X, Equiv.swap_apply_right]

theorem swapVW_cW : swapVW cW = cV := by
  show swapVW (MvPolynomial.X 2) = MvPolynomial.X 1
  rw [swapVW_X, Equiv.swap_apply_right]

/-- ★★★★`W^m ∣ swapUW P` なら `U^m ∣ P`。 -/
theorem cU_pow_dvd_of_swap (m : ℕ) (P : CollBase) (h : cW ^ m ∣ swapUW P) : cU ^ m ∣ P := by
  obtain ⟨c, hc⟩ := h
  refine ⟨swapUW c, ?_⟩
  have hid : swapUW (swapUW P) = P := congrArg (fun f => f P) swapUW_comp
  calc P = swapUW (swapUW P) := hid.symm
    _ = swapUW (cW ^ m * c) := by rw [hc]
    _ = cU ^ m * swapUW c := by rw [map_mul, map_pow, swapUW_cW]

/-- ★★★★`W^m ∣ swapVW P` なら `V^m ∣ P`。 -/
theorem cV_pow_dvd_of_swap (m : ℕ) (P : CollBase) (h : cW ^ m ∣ swapVW P) : cV ^ m ∣ P := by
  obtain ⟨c, hc⟩ := h
  refine ⟨swapVW c, ?_⟩
  have hid : swapVW (swapVW P) = P := congrArg (fun f => f P) swapVW_comp
  calc P = swapVW (swapVW P) := hid.symm
    _ = swapVW (cW ^ m * c) := by rw [hc]
    _ = cV ^ m * swapVW c := by rw [map_mul, map_pow, swapVW_cW]

/-! ## ★★★★入れ替えたあとの評価 -/

/-- ★★★★**入れ替えたあとの評価**——動かす変数が第 1 スロットに来る。 -/
theorem eval_evalUV_swapUW (u v t : ℂ) (P : CollBase) :
    (evalUV u v (swapUW P)).eval t = collEval t v u P := by
  have hhom : ((Polynomial.evalRingHom t).comp (evalUV u v)).comp swapUW = collEval t v u := by
    refine MvPolynomial.ringHom_ext ?_ ?_
    · intro r
      simp [swapUW, evalUV, collEval]
    · intro i
      fin_cases i
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapUW (MvPolynomial.X 0)) = _
        rw [swapUW_X, Equiv.swap_apply_left]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapUW (MvPolynomial.X 1)) = _
        rw [swapUW_X, show Equiv.swap (0 : Fin 3) 2 1 = 1 from by decide]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapUW (MvPolynomial.X 2)) = _
        rw [swapUW_X, Equiv.swap_apply_right]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
  exact congrArg (fun f => f P) hhom

/-- ★★★★**入れ替えたあとの評価**——動かす変数が第 2 スロットに来る。 -/
theorem eval_evalUV_swapVW (u v t : ℂ) (P : CollBase) :
    (evalUV u v (swapVW P)).eval t = collEval u t v P := by
  have hhom : ((Polynomial.evalRingHom t).comp (evalUV u v)).comp swapVW = collEval u t v := by
    refine MvPolynomial.ringHom_ext ?_ ?_
    · intro r
      simp [swapVW, evalUV, collEval]
    · intro i
      fin_cases i
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapVW (MvPolynomial.X 0)) = _
        rw [swapVW_X, show Equiv.swap (1 : Fin 3) 2 0 = 0 from by decide]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapVW (MvPolynomial.X 1)) = _
        rw [swapVW_X, Equiv.swap_apply_left]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · show ((Polynomial.evalRingHom t).comp (evalUV u v)) (swapVW (MvPolynomial.X 2)) = _
        rw [swapVW_X, Equiv.swap_apply_right]
        simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
  exact congrArg (fun f => f P) hhom

/-! ## ★★★★★★★★残る二つの整除性 -/

set_option maxHeartbeats 800000 in
theorem coeff_evalUV_swapUW_eq_zero (u v : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0)
    (hu : ‖u‖ ≤ 1 / 8) (hv : ‖v‖ ≤ 1 / 8) (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P)
    (j : ℕ) (hj : j < n + 1) :
    (evalUV u v (swapUW P)).coeff j = 0 := by
  obtain ⟨B, hB0, hB⟩ := exists_bound_collEval (d : CollBase)
  refine coeff_eq_zero_of_norm_le _
    (12 * (25 * ((n : ℝ) + 1) + 8) * 50 * 4 ^ (n + 1) * B) (1 / 8) (by norm_num) (n + 1) ?_ j hj
  intro t ht0 htr
  rw [eval_evalUV_swapUW]
  refine norm_collEval_le n P d hPd B hB0 hB t v u t ht0 hv0 hu0 htr.le hv hu ?_
  rw [norm_mul, norm_mul]
  have h3 : ‖v‖ * ‖u‖ ≤ 1 := by nlinarith [norm_nonneg u, norm_nonneg v]
  have h4 : ‖t‖ * ‖v‖ * ‖u‖ = (‖v‖ * ‖u‖) * ‖t‖ := by ring
  rw [h4]
  nlinarith [norm_nonneg t]

set_option maxHeartbeats 800000 in
theorem coeff_evalUV_swapVW_eq_zero (u v : ℂ) (hu0 : u ≠ 0) (hv0 : v ≠ 0)
    (hu : ‖u‖ ≤ 1 / 8) (hv : ‖v‖ ≤ 1 / 8) (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P)
    (j : ℕ) (hj : j < n + 1) :
    (evalUV u v (swapVW P)).coeff j = 0 := by
  obtain ⟨B, hB0, hB⟩ := exists_bound_collEval (d : CollBase)
  refine coeff_eq_zero_of_norm_le _
    (12 * (25 * ((n : ℝ) + 1) + 8) * 50 * 4 ^ (n + 1) * B) (1 / 8) (by norm_num) (n + 1) ?_ j hj
  intro t ht0 htr
  rw [eval_evalUV_swapVW]
  refine norm_collEval_le n P d hPd B hB0 hB u t v t hu0 ht0 hv0 hu htr.le hv ?_
  rw [norm_mul, norm_mul]
  have h3 : ‖u‖ * ‖v‖ ≤ 1 := by nlinarith [norm_nonneg u, norm_nonneg v]
  have h4 : ‖u‖ * ‖t‖ * ‖v‖ = (‖u‖ * ‖v‖) * ‖t‖ := by ring
  rw [h4]
  nlinarith [norm_nonneg t]

/-- ★★★★★★★★**分子は `U^{n+1}` で割れる**。 -/
theorem cU_pow_dvd_numerator (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P) :
    cU ^ (n + 1) ∣ P :=
  cU_pow_dvd_of_swap (n + 1) P
    (cW_pow_dvd_of_coeff n (swapUW P) smallSet infinite_smallSet
      (fun a ha b hb j hj => coeff_evalUV_swapUW_eq_zero a b (smallSet_ne_zero ha)
        (smallSet_ne_zero hb) (smallSet_norm_le ha) (smallSet_norm_le hb) n P d hPd j hj))

/-- ★★★★★★★★**分子は `V^{n+1}` で割れる**。 -/
theorem cV_pow_dvd_numerator (n : ℕ) (P : CollBase) (d : collDenoms)
    (hPd : collDefectTrunc (n + 1) kU kV kW (kU * kV * kW)
      * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P) :
    cV ^ (n + 1) ∣ P :=
  cV_pow_dvd_of_swap (n + 1) P
    (cW_pow_dvd_of_coeff n (swapVW P) smallSet infinite_smallSet
      (fun a ha b hb j hj => coeff_evalUV_swapVW_eq_zero a b (smallSet_ne_zero ha)
        (smallSet_ne_zero hb) (smallSet_norm_le ha) (smallSet_norm_le hb) n P d hPd j hj))

/-! ## ★★★★★★★★★葉 (c) の到達点 -/

namespace TateCollinearSection

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★★**Tate 級数の 3 点は共線**——葉 (c)(一意化が準同型)の中身。

`u v w = q` のとき、3 点 `(X(u),Y(u))`, `(X(v),Y(v))`, `(X(w),Y(w))` の
共線性の行列式が 0 になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_collinear [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    tateXpair u (v * w) q hq * (tateYpair v (u * w) q hq - tateYpair w (u * v) q hq)
      + tateXpair v (u * w) q hq * (tateYpair w (u * v) q hq - tateYpair u (v * w) q hq)
      + tateXpair w (u * v) q hq * (tateYpair u (v * w) q hq - tateYpair v (u * w) q hq) = 0 := by
  have h : collDefect u v w q hq = 0 := by
    refine collDefect_eq_zero_of_base u v w q hq huvw hcp ?_ ?_ ?_
    · intro n P d hPd
      cases n with
      | zero => simp
      | succ m => exact cU_pow_dvd_numerator m P d hPd
    · intro n P d hPd
      cases n with
      | zero => simp
      | succ m => exact cV_pow_dvd_numerator m P d hPd
    · intro n P d hPd
      cases n with
      | zero => simp
      | succ m => exact cW_pow_dvd_numerator m P d hPd
  exact h

end TateCollinearSection

/-! ## ★出典の紐付け(`.src`) -/

def norm_collEval_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——動かす変数はどのスロットでもよい)",
    sectionId := "genell-def-3-3" }

def TateCollinearSection.tate_collinear.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——3 点は共線、葉 (c) の到達点)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
