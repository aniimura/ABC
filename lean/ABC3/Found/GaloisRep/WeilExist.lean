import ABC3.Found.GaloisRep.WeilVal

/-!
# Galois (G5) 第 179 ブロック —— **★★★★★★★★Weil 対の存在と第 1 性質**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★`WeilSpec` を満たす値が存在する

第 178 で定義した `weilPairingVal` が「本物の」値になるには、`WeilSpec` を満たす `c` が
存在しなければならない。★材料はすべて揃っている:

| 部品 | 出どころ |
|---|---|
| `f_P`(`XYIdeal(P)^n` の生成元) | 第 126 `xyIdeal_pow_isPrincipal_integral` |
| `g`(`μ f_P` の `n` 乗根) | **D2**(第 149・142・175・139 の合流) |
| `τ`(`Q` の平行移動) | 第 124 `exists_translateAut_all` |
| `τ (μ f_P) = μ f_P` | 第 168 `aut_comp_mulByN` |
| `τ g / g` が定数 | 第 176 `exists_const_aut_div` |

★★`exists_nthRoot_mulByN` は `Skeleton` の `exists_nthRoot_comp_mulByN` を
**`Found` 側で組み直したもの**である(`Found` は `Skeleton` を import しない)。

## ★★★★★★第 1 性質——`e_n(P,Q)^n = 1`

`(τ g / g)^n = τ(g^n)/g^n = τ(μ f_P)/(μ f_P) = 1` から直ちに出る。
★`τ (μ f_P) = μ f_P` には `n·Q = 0` が要る(第 168)——スケルトンの述べ方と一致する。
★★データが無い場合は既定値 1 なので `1^n = 1` ✓。

## ★★残る仮定

`μ`(= `[n]^*`)の**単射性**は仮定に置いている。★これは既知の未解決事項であり、
`Skeleton` の `dvd_count_pullback` の `.needs` に「`x([n]·)` の超越性を別途示す必要がある
(3-8 ブロック)」として記録済みである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_nthRoot_mulByN` | ★★★★★★★★`Found` 側の `n` 乗根の存在 |
| `weilSpec_of_data` | ★★★★★★★★**`WeilSpec` の witness** |
| `weilPairingVal_pow_eq_one` | ★★★★★★**`e_n(P,Q)^n = 1`** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep
open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [IsAlgClosed F] [DecidableEq F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-! ## ★★★★★★★★`n` 乗根の存在(`Found` 側) -/

/-- ★★★★★★★★**`Found` 側で組み直した `n` 乗根の存在**。

★第 149(D1')・第 142(`n` 乗の構成)・第 175(D2)・第 139(定数の `n` 乗根)の合流。 -/
theorem exists_nthRoot_mulByN (h2 : IsUnit (2 : F))
    {x y : F} (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ g : W.FunctionField, g ^ n = μ fP := by
  haveI := isDedekindDomain_coordinateRing h2 W
  have hdvd := dvd_count_pullback W h2 μ hμinj hμF h n hP fP hfP
  have hμfP : μ fP ≠ 0 := fun h0 =>
    (fP_ne_zero W n fP hfP) (hμinj (by rw [h0, map_zero]))
  have hIne : FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) ≠ 0 := by
    rw [ne_eq, FractionalIdeal.spanSingleton_eq_zero_iff]; exact hμfP
  obtain ⟨J, hJ⟩ := exists_pow_of_dvd_count hIne hdvd
  obtain ⟨g, hg⟩ := fractionalIdeal_isPrincipal_proof h2 W h n hn hchar hP μ hμinj hμF
    hns hμP hμx hμy fP hfP J hJ
  exact exists_nthRoot_of_fractionalIdeal W hn hJ hg

/-! ## ★★★★★★★★`WeilSpec` の witness -/

/-- ★★★★★★★★**`WeilSpec` を満たす値が存在する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`f_P`・`g`・`τ` を順に取り、第 176 で `τ g / g` が定数であることを言うだけ。 -/
theorem weilSpec_of_data [Infinite F] (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    {x y : F} (hP : W.Nonsingular x y) (hPt : n • (Point.some x y hP) = 0)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (hQt : n • (Point.some x₀ y₀ hQ) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn) :
    ∃ c : F, WeilSpec W n (Point.some x y hP) (Point.some x₀ y₀ hQ) c := by
  obtain ⟨fP, hfPne, hfP⟩ := xyIdeal_pow_isPrincipal_integral hP n hPt
  obtain ⟨g, hg⟩ := exists_nthRoot_mulByN W h2 hP n hn hchar hPt μ hμinj hμF hns hμP hμx hμy
    fP hfP
  obtain ⟨τ, hτx, hτy⟩ := exists_translateAut_all W h4 hQ
  have hcomp := aut_comp_mulByN W τ hQ hτx hτy n hQt μ hμF hns hμx hμy hμP
  have hτz : τ (μ fP) = μ fP :=
    congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hcomp) fP
  have hz : μ fP ≠ 0 := fun h0 => hfPne (hμinj (by rw [h0, map_zero]))
  obtain ⟨c, hc0, hcdiv⟩ := exists_const_aut_div W h2 hn hg hz hτz
  exact ⟨c, x, y, hP, x₀, y₀, hQ, fP, μ, g, τ, xn, yn, hns, rfl, rfl, hfP, hμinj, hμF,
    hμx, hμy, hμP, hg, hτx, hτy, hcdiv⟩

/-! ## ★★★★★★第 1 性質 -/

/-- ★★★★★★**`e_n(P,Q)` は 1 の `n` 乗根である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`(τ g / g)^n = τ(g^n)/g^n = τ(μ f_P)/(μ f_P) = 1`。
★★データが無い場合は既定値 1 なので `1^n = 1`。 -/
theorem weilPairingVal_pow_eq_one {n : ℕ} (P Q : W.Point) (hQt : n • Q = 0) :
    (weilPairingVal W n P Q) ^ n = 1 := by
  by_cases hex : ∃ c : F, WeilSpec W n P Q c
  · obtain ⟨x, y, hP, x₀, y₀, hQn, fP, μ, g, τ, xn, yn, hns,
      hPe, hQe, hfP, hinj, hμF, hμx, hμy, hμP, hg, hτx, hτy, hdiv⟩ :=
      weilPairingVal_spec W n P Q hex
    have hQt' : n • Point.some x₀ y₀ hQn = 0 := by rw [← hQe]; exact hQt
    have hcomp := aut_comp_mulByN W τ hQn hτx hτy n hQt' μ hμF hns hμx hμy hμP
    have hτz : τ (μ fP) = μ fP :=
      congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hcomp) fP
    have hz : μ fP ≠ 0 := fun h0 => (fP_ne_zero W n fP hfP) (hinj (by rw [h0, map_zero]))
    have hone := pow_aut_div_eq_one W hg hz hτz
    rw [hdiv, ← map_pow] at hone
    exact (algebraMap F W.FunctionField).injective (hone.trans (map_one _).symm)
  · rw [weilPairingVal_of_not W n P Q hex, one_pow]

/-! ## ★出典の紐付け(`.src`) -/

def weilSpec_of_data.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の値が存在すること)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_pow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対が 1 の n 乗根であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
