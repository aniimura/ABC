import ABC3.Found.GaloisRep.TateDivisibility
import ABC3.Found.GaloisRep.TateEquationAnalytic

/-!
# Galois (G6) 第 225 ブロック —— **★★★★★★万有な環から ℂ への特殊化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★解析の段への橋を架ける

第 224 で葉 (b) は「万有な環の分子 `P ∈ ℤ[A,W]` が `A^n`・`W^n` で割れる」に落ちた。
分子の係数の消滅を言うには、**`A ↦ u`、`W ↦ w`(複素数)と特殊化**して
`q = uw → 0` の様子を見る。そのための橋が本ブロックである。

★★**完備性は要らない**——`TateUniv` から ℂ への環準同型は、分母が ℂ で単元でありさえ
すれば作れる。`‖u‖ < 1` かつ `‖w‖ < 1` なら:

| 分母 | ℂ での単元性 |
|---|---|
| `1 − u` | ★`‖u‖ < 1` なので `u ≠ 1` |
| `1 − w` | ★同上 |
| `1 − (uw)^{m+1}u` | ★★`‖(uw)^{m+1}u‖ ≤ ‖u‖ < 1` |
| `1 − (uw)^{m+1}w` | ★★同上 |

これで `tateSpecializeC : TateUniv →+* ℂ` が得られる。

## ★★★★★★橋の形

`x·(分母) = (分子)` を ℂ に落とすと:

    tateEval u w P = tateDefectTrunc n u w (uw) · tateEval u w d       (`tateEval_numerator`)

★左辺は `P` の値、右辺は**切り詰めた差の複素数値**である。
★★あとは「`q = uw → 0` のとき右辺が `O(q^n)`」を解析側から出せばよい
——それが残る唯一の解析の段である。

## ★★`im z < im τ` は格子の外を意味する

`z = nτ + m` なら `im z = n·im τ`。`0 < im z < im τ` から `0 < n < 1` となり、
整数 `n` としては不可能である(`notMem_lattice_of_im_lt`)。
★これで第 220 の `tate_equation_analytic` から仮定 `hz` が落ちる
(`tate_equation_analytic'`)——解析側の恒等式は **`im z < im τ` だけで成り立つ**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `notMem_lattice_of_im_lt` | ★★★★`im z < im τ` なら `z` は格子の外 |
| `tate_equation_analytic'` | ★★★仮定を `him` だけにした解析側の方程式 |
| `isUnit_one_sub_of_norm_lt_one` | ★`‖t‖ < 1` なら `1 − t` は単元 |
| `tateSpecializeC` | ★★★★★★**万有な環から ℂ への特殊化** |
| `tateSpecializeC_tateDefectTrunc` | ★★★★★ℂ 側での切り詰めた差の値 |
| `tateEval_numerator` | ★★★★★★**分子の値 = 切り詰めた差 × 分母** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial PeriodPair

/-! ## ★★★★格子の外であること -/

/-- ★★★★**`im z < im τ` なら `z` は格子 `ℤτ + ℤ` の点でない**。

★`z = nτ + m` なら `im z = n·im τ`。`0 < im z < im τ` から `0 < n < 1` となり、
整数 `n` としては不可能である。 -/
theorem notMem_lattice_of_im_lt (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    (z : ℂ) ∉ (tauPair τ).lattice := by
  intro hz
  obtain ⟨x, hx⟩ := (tauPair τ).latticeEquivProd.symm.surjective ⟨(z : ℂ), hz⟩
  have hval : ((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ)) = (z : ℂ) := by
    rw [← tauPair_symm_apply τ x, hx]
  have him2 : (x.1 : ℝ) * (τ : ℂ).im = (z : ℂ).im := by
    have hc := congrArg Complex.im hval
    simpa using hc
  have hzpos : 0 < (z : ℂ).im := z.2
  have htpos : 0 < (τ : ℂ).im := τ.2
  have h1 : (0 : ℝ) < (x.1 : ℝ) := by nlinarith
  have h2 : (x.1 : ℝ) < 1 := by nlinarith
  have hi1 : (0 : ℤ) < x.1 := by exact_mod_cast h1
  have hi2 : x.1 < 1 := by exact_mod_cast h2
  omega

/-- ★★★**解析側の Tate 方程式(`im z < im τ` だけを仮定する形)**。 -/
theorem tate_equation_analytic' (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    tateYfun z τ ^ 2 + tateXfun z τ * tateYfun z τ
      = tateXfun z τ ^ 3 + (-5 * sigmaSum 3 τ) * tateXfun z τ
        - (5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12 :=
  tate_equation_analytic z τ him (notMem_lattice_of_im_lt z τ him)

/-! ## ★★分母は ℂ で単元になる -/

theorem isUnit_one_sub_of_norm_lt_one {t : ℂ} (h : ‖t‖ < 1) : IsUnit (1 - t) := by
  rw [isUnit_iff_ne_zero]
  intro hc
  have ht : t = 1 := by linear_combination -hc
  rw [ht] at h
  simp at h

theorem norm_lt_one_of_mul_pow {u w : ℂ} (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) (m : ℕ) :
    ‖(u * w) ^ (m + 1) * u‖ < 1 := by
  rw [norm_mul, norm_pow, norm_mul, mul_pow]
  have h1 : ‖u‖ ^ (m + 1) ≤ 1 := pow_le_one₀ (norm_nonneg u) hu.le
  have h2 : ‖w‖ ^ (m + 1) ≤ 1 := pow_le_one₀ (norm_nonneg w) hw.le
  have key : ‖u‖ ^ (m + 1) * ‖w‖ ^ (m + 1) * ‖u‖ ≤ 1 * 1 * ‖u‖ := by
    refine mul_le_mul_of_nonneg_right ?_ (norm_nonneg u)
    exact mul_le_mul h1 h2 (pow_nonneg (norm_nonneg w) _) zero_le_one
  linarith

theorem tateEval_isUnit_denomSet_complex (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1)
    {x : TateBase} (hx : x ∈ tateDenomSet) : IsUnit (tateEval u w x) := by
  rcases hx with (hx | ⟨m, rfl⟩) | ⟨m, rfl⟩
  · rcases hx with rfl | rfl
    · simpa using isUnit_one_sub_of_norm_lt_one hu
    · simpa using isUnit_one_sub_of_norm_lt_one hw
  · have hval : tateEval u w (1 - univQ ^ (m + 1) * univA) = 1 - (u * w) ^ (m + 1) * u := by
      simp
    rw [hval]
    exact isUnit_one_sub_of_norm_lt_one (norm_lt_one_of_mul_pow hu hw m)
  · have hval : tateEval u w (1 - univQ ^ (m + 1) * univW) = 1 - (u * w) ^ (m + 1) * w := by
      simp
    rw [hval]
    refine isUnit_one_sub_of_norm_lt_one ?_
    have hkey := norm_lt_one_of_mul_pow hw hu m
    rw [show w * u = u * w by ring] at hkey
    exact hkey

theorem tateEval_isUnit_denoms_complex (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1)
    (y : tateDenoms) : IsUnit (tateEval u w (y : TateBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := tateDenomSet)
    (motive := fun z _ => IsUnit (tateEval u w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact tateEval_isUnit_denomSet_complex u w hu hw hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-! ## ★★★★★★ℂ への特殊化 -/

/-- ★★★★★★**万有な環から ℂ への特殊化**——`‖u‖ < 1`、`‖w‖ < 1` だけで足りる。

★完備性は要らない(分母が単元でありさえすればよい)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tateSpecializeC (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) :
    TateUniv →+* ℂ :=
  IsLocalization.lift (tateEval_isUnit_denoms_complex u w hu hw)

theorem tateSpecializeC_base (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) (x : TateBase) :
    tateSpecializeC u w hu hw (algebraMap TateBase TateUniv x) = tateEval u w x :=
  IsLocalization.lift_eq _ x

theorem tateSpecializeC_uA (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) :
    tateSpecializeC u w hu hw uA = u := by
  rw [uA, tateSpecializeC_base, tateEval_A]

theorem tateSpecializeC_uW (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) :
    tateSpecializeC u w hu hw uW = w := by
  rw [uW, tateSpecializeC_base, tateEval_W]

/-- ★★★★★**ℂ 側での切り詰めた差の値**。 -/
theorem tateSpecializeC_tateDefectTrunc (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) (n : ℕ) :
    tateSpecializeC u w hu hw (tateDefectTrunc n uA uW (uA * uW))
      = tateDefectTrunc n u w (u * w) := by
  have h := map_tateDefectTrunc (tateSpecializeC u w hu hw) n uA uW (uA * uW)
    isUnit_one_sub_uA isUnit_one_sub_uW (fun m _ => isUnit_one_sub_uQA m)
    (fun m _ => isUnit_one_sub_uQW m)
  rw [h, tateSpecializeC_uA, tateSpecializeC_uW, map_mul, tateSpecializeC_uA,
    tateSpecializeC_uW]

/-- ★★★★★★**分子の値は「切り詰めた差 × 分母」である**——これが解析の段への橋である。

★あとは「`q = uw → 0` のとき右辺が `O(q^n)`」を解析側から出せばよい。 -/
theorem tateEval_numerator (n : ℕ) (P : TateBase) (d : tateDenoms)
    (hPd : tateDefectTrunc n uA uW (uA * uW) * algebraMap TateBase TateUniv (d : TateBase)
      = algebraMap TateBase TateUniv P)
    (u w : ℂ) (hu : ‖u‖ < 1) (hw : ‖w‖ < 1) :
    tateEval u w P = tateDefectTrunc n u w (u * w) * tateEval u w (d : TateBase) := by
  have h := congrArg (tateSpecializeC u w hu hw) hPd
  rw [map_mul, tateSpecializeC_base, tateSpecializeC_base,
    tateSpecializeC_tateDefectTrunc] at h
  exact h.symm

/-! ## ★出典の紐付け(`.src`) -/

def notMem_lattice_of_im_lt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——im z < im tau なら格子の外)",
    sectionId := "genell-def-3-3" }

def tateSpecializeC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——万有な環から ℂ への特殊化)",
    sectionId := "genell-def-3-3" }

def tateEval_numerator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分子の値と切り詰めた差)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
