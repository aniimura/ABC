import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GaloisRep.LatticeInvariant
import ABC3.Found.GaloisRep.TateEquationAnalytic
import Mathlib.NumberTheory.ModularForms.LevelOne.GradedRing

/-!
# GenEll 第 346 ブロック —— **★★★★★★★(i) 判別式の非消失が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★到達点

> **任意の周期対 `L` について `latticeDisc L = g₂³ − 27g₃² ≠ 0`**(`latticeDisc_ne_zero`)

★これで `latticeCurve L` はつねに楕円曲線になり(`isElliptic_latticeCurve'`)、
一意化 (ii)(群同型)・(iii)(`j` の全射性)の土台ができた。

★★おまけに **`j(ℤ + τℤ) = E₄(τ)³ / Δ(τ)`**(`latticeCurve_j_tauPair`)——
格子から作った曲線の `j` 不変量が**古典的なモジュラー `j` 関数**であることが出た。
★★★これは (iii) の入口そのものである。

## ★★★★★★★★2026-08-26 の在庫調査——**8 ブロックぶんの回り道をしていた**

★本ブロックの直前に、mathlib を「概念で」引き直したところ:

| 探していたもの | 実は | 場所 |
|---|---|---|
| `E₄³ − E₆² = 1728Δ` | ★**あった** | `ModularForms/LevelOne/GradedRing.lean`<br>`discriminant_eq_E₄_cube_sub_E₆_sq` |
| `E₄` の `q` 係数 `240` | ★**あった** | `LevelOne/DimensionFormula.lean` `E₄_qExpansion_coeff_one` |
| `E₆` の `q` 係数 `−504` | ★**あった** | 同上 `E₆_qExpansion_coeff_one` |
| 格子和 = `ζ(k)·Eisenstein` | ★**あった** | `tsum_eisSummand_eq_riemannZeta_mul_eisensteinSeries` |

★★★さらに悪いことに、4 つめは**自分の第 215 ブロック
(`Found/GaloisRep/LatticeInvariant.lean`、`G_eq_two_zeta_mul_E`)が既に使っていた**。

★★★★つまり第 337-345 の 9 ブロック(gcd 分解 → 全単射 → 積 → 絶対収束 → Fubini →
`ζ(4)`,`ζ(6)` → Eisenstein 同定 → 係数)は、
**すでに手元にある `G_eq_two_zeta_mul_E` の一行**で置き換えられた。
★本ブロックの証明が `g₂_tauPair` の 3 行で済んでいるのがその証拠である。

★★★★★**なぜ落ちたか**——`grep` を「**これから使うつもりの道具の名前**」
(`sturm_bound_levelOne`・`qExpansion_mul`)で打った。
`1728` という**概念**で打てば 1 秒で当たっていた。
`tools/lean-idioms.md` の「在庫は名前でなく概念で引く」の**2 度目の違反**である。

## ★★★★★★組み立て(残ったのは 5 行)

    g₂(Λ_τ) = 60·G₄ = 60·2ζ(4)·E₄ = (4π⁴/3)·E₄        (G_eq_two_zeta_mul_E + riemannZeta_four)
    g₃(Λ_τ) = 140·G₆ = 140·2ζ(6)·E₆ = (8π⁶/27)·E₆      (同 + riemannZeta_six、第 342 の在庫)
    Δ_lat   = (4π⁴/3)³E₄³ − 27(8π⁶/27)²E₆²
            = (64π¹²/27)(E₄³ − E₆²)
            = (64π¹²/27)·1728·Δ(τ) = 4096π¹²·Δ(τ)      (discriminant_eq_E₄_cube_sub_E₆_sq)
            ≠ 0                                        (discriminant_ne_zero)

★★一般の `L` へは**スケール則**(第 333 `latticeDisc_ne_zero_iff`)で降ろす:
`L = ω₂ · ⟨ω₁/ω₂, 1⟩`。★`im(ω₁/ω₂) < 0` のときは `ω₁ ↔ ω₂` を入れ替える
——**束 `span{ω₁,ω₂}` は入れ替えで変わらない**(`swap_lattice`)ので判別式も変わらない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `g₂_tauPair`・`g₃_tauPair` | ★★★★★`g₂ = (4π⁴/3)E₄`・`g₃ = (8π⁶/27)E₆` |
| `latticeDisc_tauPair` | ★★★★★★**`Δ_lat = 4096π¹²·Δ(τ)`** |
| `latticeDisc_tauPair_ne_zero` | ★★★★★正規化された束での非消失 |
| `omega₂_ne_zero`・`ratio_im_ne_zero` | ★周期比は実でない |
| `G_congr`・`latticeDisc_congr`・`swap_lattice` | ★★束が同じなら不変量も同じ |
| `latticeDisc_ne_zero` | ★★★★★★★**(i) 判別式の非消失** |
| `isElliptic_latticeCurve'` | ★★★★★★格子の曲線はつねに楕円曲線 |
| `latticeCurve_j_tauPair` | ★★★★★★★**`j = E₄³/Δ`**((iii) の入口) |
-/

namespace ABC3.Found.GenEll

open Complex Real PeriodPair ModularForm ABC3.Found.GaloisRep

/-! ## ★★★★★`g₂`・`g₃` を `E₄`・`E₆` で書く -/

/-- ★★★★★**`g₂(ℤ + τℤ) = (4π⁴/3)·E₄(τ)`**。 -/
theorem g₂_tauPair (τ : UpperHalfPlane) :
    (tauPair τ).g₂ = (4 * (π : ℂ) ^ 4 / 3) * ModularForm.E₄ τ := by
  rw [PeriodPair.g₂, G_eq_two_zeta_mul_E τ (k := 4) (by norm_num)]
  norm_num [riemannZeta_four]
  ring

/-- ★★★★★**`g₃(ℤ + τℤ) = (8π⁶/27)·E₆(τ)`**。 -/
theorem g₃_tauPair (τ : UpperHalfPlane) :
    (tauPair τ).g₃ = (8 * (π : ℂ) ^ 6 / 27) * ModularForm.E₆ τ := by
  rw [PeriodPair.g₃, G_eq_two_zeta_mul_E τ (k := 6) (by norm_num)]
  norm_num [ABC3.Found.GaloisRep.riemannZeta_six]
  ring

/-! ## ★★★★★★判別式 = `4096π¹²·Δ` -/

/-- ★★★★★★**`g₂³ − 27g₃² = 4096π¹²·Δ(τ)`**——`(2π)¹²Δ` である。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`E₄³ − E₆² = 1728Δ` は mathlib の `discriminant_eq_E₄_cube_sub_E₆_sq`。 -/
theorem latticeDisc_tauPair (τ : UpperHalfPlane) :
    latticeDisc (tauPair τ) = 4096 * (π : ℂ) ^ 12 * ModularForm.discriminant τ := by
  rw [latticeDisc, g₂_tauPair, g₃_tauPair, ModularForm.discriminant_eq_E₄_cube_sub_E₆_sq]
  ring

/-- ★★★★★正規化された束 `ℤ + τℤ` では判別式は消えない。 -/
theorem latticeDisc_tauPair_ne_zero (τ : UpperHalfPlane) : latticeDisc (tauPair τ) ≠ 0 := by
  rw [latticeDisc_tauPair]
  refine mul_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero 12 ?_))
    (ModularForm.discriminant_ne_zero τ)
  exact_mod_cast Real.pi_ne_zero

/-! ## ★一般の周期対への降ろし方 -/

/-- ★周期対は成分で決まる。 -/
theorem periodPair_ext {L L' : PeriodPair} (h1 : L.ω₁ = L'.ω₁) (h2 : L.ω₂ = L'.ω₂) : L = L' := by
  cases L; cases L'; simp_all

/-- ★`ω₂ ≠ 0`。 -/
theorem omega₂_ne_zero (L : PeriodPair) : L.ω₂ ≠ 0 := by
  have h := (linearIndependent_fin2 (K := ℝ)).1 L.indep
  simpa using h.1

/-- ★★**周期比 `ω₁/ω₂` は実数でない**——`ℝ` 上の一次独立から。 -/
theorem ratio_im_ne_zero (L : PeriodPair) : (L.ω₁ / L.ω₂).im ≠ 0 := by
  have h := (linearIndependent_fin2 (K := ℝ)).1 L.indep
  have h2 := omega₂_ne_zero L
  intro h0
  refine h.2 ((L.ω₁ / L.ω₂).re) ?_
  have hre : (L.ω₁ / L.ω₂) = ((L.ω₁ / L.ω₂).re : ℂ) := by
    apply Complex.ext <;> simp [h0]
  have hmul : L.ω₁ = ((L.ω₁ / L.ω₂).re : ℂ) * L.ω₂ := by
    rw [← hre]; field_simp
  simpa [Complex.real_smul] using hmul.symm

/-! ## ★★束が同じなら不変量も同じ -/

/-- ★★★束が同じなら `G` も同じ。 -/
theorem G_congr {L L' : PeriodPair} (h : L.lattice = L'.lattice) (n : ℕ) :
    L.G n = L'.G n := by
  rw [PeriodPair.G, PeriodPair.G]
  exact (Equiv.subtypeEquivRight (fun x => by rw [h])).tsum_eq
    (fun l : L'.lattice => ((l : ℂ) ^ n)⁻¹)

/-- ★★束が同じなら判別式も同じ。 -/
theorem latticeDisc_congr {L L' : PeriodPair} (h : L.lattice = L'.lattice) :
    latticeDisc L = latticeDisc L' := by
  rw [latticeDisc, latticeDisc, PeriodPair.g₂, PeriodPair.g₂, PeriodPair.g₃, PeriodPair.g₃,
    G_congr h 4, G_congr h 6]

/-- ★★**周期の入れ替えは束を変えない**——`span{ω₂,ω₁} = span{ω₁,ω₂}`。 -/
theorem swap_lattice (L : PeriodPair) (h : LinearIndependent ℝ ![L.ω₂, L.ω₁]) :
    (⟨L.ω₂, L.ω₁, h⟩ : PeriodPair).lattice = L.lattice := by
  rw [PeriodPair.lattice, PeriodPair.lattice]
  exact congrArg _ (Set.pair_comm _ _)

/-- ★入れ替えても一次独立。 -/
theorem indep_swap (L : PeriodPair) : LinearIndependent ℝ ![L.ω₂, L.ω₁] := by
  have h := L.indep.comp (Equiv.swap (0 : Fin 2) 1) (Equiv.injective _)
  have he : (![L.ω₁, L.ω₂] ∘ (Equiv.swap (0 : Fin 2) 1)) = ![L.ω₂, L.ω₁] := by
    funext i; fin_cases i <;> simp
  rwa [he] at h

/-- ★★★★`im(ω₁/ω₂) > 0` のとき、`L` は `⟨τ, 1⟩` の `ω₂` 倍である。 -/
theorem eq_scalePair_tauPair (L : PeriodPair) (h : 0 < (L.ω₁ / L.ω₂).im) :
    scalePair (tauPair ⟨L.ω₁ / L.ω₂, h⟩) L.ω₂ (omega₂_ne_zero L) = L := by
  have h2 := omega₂_ne_zero L
  refine periodPair_ext ?_ ?_
  · show L.ω₂ * (L.ω₁ / L.ω₂) = L.ω₁
    rw [mul_comm, div_mul_cancel₀ _ h2]
  · show L.ω₂ * 1 = L.ω₂
    ring

/-- ★★★★★向きが正の周期対での非消失。 -/
theorem latticeDisc_ne_zero_of_im_pos (L : PeriodPair) (h : 0 < (L.ω₁ / L.ω₂).im) :
    latticeDisc L ≠ 0 := by
  rw [← eq_scalePair_tauPair L h, latticeDisc_ne_zero_iff]
  exact latticeDisc_tauPair_ne_zero _

/-! ## ★★★★★★★(i) 判別式の非消失 -/

/-- ★★★★★★★**任意の周期対の判別式は消えない**——一意化 (i)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem latticeDisc_ne_zero (L : PeriodPair) : latticeDisc L ≠ 0 := by
  rcases lt_or_gt_of_ne (ratio_im_ne_zero L) with hneg | hpos
  · have hind := indep_swap L
    have him : 0 < ((⟨L.ω₂, L.ω₁, hind⟩ : PeriodPair).ω₁
        / (⟨L.ω₂, L.ω₁, hind⟩ : PeriodPair).ω₂).im := by
      show 0 < (L.ω₂ / L.ω₁).im
      rw [show L.ω₂ / L.ω₁ = (L.ω₁ / L.ω₂)⁻¹ from (inv_div _ _).symm]
      have hz : L.ω₁ / L.ω₂ ≠ 0 := fun h0 => by simp [h0] at hneg
      rw [Complex.inv_im]
      exact div_pos (by linarith) (Complex.normSq_pos.2 hz)
    rw [latticeDisc_congr (swap_lattice L hind).symm]
    exact latticeDisc_ne_zero_of_im_pos _ him
  · exact latticeDisc_ne_zero_of_im_pos L hpos

/-- ★★★★★★**格子から作った曲線はつねに楕円曲線である**。 -/
instance isElliptic_latticeCurve' (L : PeriodPair) : (latticeCurve L).IsElliptic :=
  isElliptic_latticeCurve L (latticeDisc_ne_zero L)

/-! ## ★★★★★★★`j` は古典的なモジュラー `j` 関数 -/

/-- ★★★★★★★**`j(ℤ + τℤ) = E₄(τ)³ / Δ(τ)`**——(iii) の入口。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem latticeCurve_j_tauPair (τ : UpperHalfPlane) :
    (latticeCurve (tauPair τ)).j = ModularForm.E₄ τ ^ 3 / ModularForm.discriminant τ := by
  have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hΔ := ModularForm.discriminant_ne_zero τ
  rw [latticeCurve_j, g₂_tauPair, latticeDisc_tauPair]
  field_simp
  ring

/-! ## ★出典の紐付け(`.src`) -/

def latticeDisc_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def latticeDisc_tauPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def latticeCurve_j_tauPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
