import ABC3.Found.GaloisRep.RootOfUnity

/-!
# Galois (G5) 第 177 ブロック —— **★★★★★★Weil 対が選択に依らないこと**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★定義に入る選択は 3 つある

    e_n(P, Q) := τ_Q(g_P) / g_P

に入る選択は `f_P`(`XYIdeal(P)^n` の生成元)・`g_P`(`μ f_P` の `n` 乗根)・
`τ_Q`(平行移動)の 3 つである。★どれを取っても値が変わらないことを示す。

| 選択 | 依らない理由 | 本ブロック |
|---|---|---|
| `f_P` | 生成元は定数倍の差(第 128) | `aut_div_eq_of_const_mul` |
| `g_P` | `g'/g` は定数(下の `const_of_pow_eq_const`) | 同上 |
| `τ_Q` | 座標での値が決まれば一意 | `aut_ext` |

### ★★★★★定数の `n` 乗根も定数である

`h^n = c`(`c ∈ F^×`)なら、`F` が代数閉なので `c = e^n` と書けて
`(h/e)^n = 1`。★第 176 の `const_of_pow_eq_one` から `h/e` は定数、よって `h` も定数。

### ★★★★★`τ_Q` の一意性

`F(W)` の `F`-自己同型は `coordX`・`coordY` での値で決まる。★座標環への制限が
第 119 の `coordinateRing_hom_ext` で決まり、`F(W) = Frac(F[W])` なので
mathlib の `IsLocalization.ringHom_ext` で `F(W)` 全体に伸びる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `const_of_pow_eq_const` | ★★★★★**定数の `n` 乗根も定数** |
| `aut_div_eq_of_const_mul` | ★★★★★★**`τ(g)/g` は `g` の定数倍に依らない** |
| `aut_ext` | ★★★★★**自己同型は座標での値で決まる** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★★★★自己同型の一意性 -/

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- ★★★★★**関数体の `F`-自己同型は座標での値で決まる**。

★座標環への制限が第 119 の `coordinateRing_hom_ext` で決まり、
`F(W) = Frac(F[W])` なので `IsLocalization.ringHom_ext` で伸びる。 -/
theorem aut_ext {τ τ' : W.FunctionField ≃ₐ[F] W.FunctionField}
    (hx : τ (coordX W) = τ' (coordX W)) (hy : τ (coordY W) = τ' (coordY W)) : τ = τ' := by
  have hcomp : (τ.toAlgHom.toRingHom).comp (algebraMap W.CoordinateRing W.FunctionField)
      = (τ'.toAlgHom.toRingHom).comp (algebraMap W.CoordinateRing W.FunctionField) := by
    refine coordinateRing_hom_ext _ _ (fun d => ?_) hx hy
    simp only [RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe,
      ← IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField d]
    exact (τ.commutes d).trans (τ'.commutes d).symm
  have hring : (τ.toAlgHom.toRingHom : W.FunctionField →+* W.FunctionField)
      = τ'.toAlgHom.toRingHom :=
    IsLocalization.ringHom_ext W.CoordinateRing⁰ hcomp
  refine AlgEquiv.ext (fun z => ?_)
  exact congrFun (congrArg (fun f => (f : W.FunctionField →+* W.FunctionField).toFun) hring) z

/-! ## ★★★★★定数の `n` 乗根 -/

variable [DecidableEq F] [IsAlgClosed F] [W.IsElliptic]

/-- ★★★★★**定数の `n` 乗根も定数である**。 -/
theorem const_of_pow_eq_const (h2 : IsUnit (2 : F)) {h : W.FunctionField} {c : F} {n : ℕ}
    (hn : 1 ≤ n) (hc : c ≠ 0) (hpow : h ^ n = algebraMap F W.FunctionField c) :
    ∃ d : F, d ≠ 0 ∧ h = algebraMap F W.FunctionField d := by
  obtain ⟨e, he⟩ := IsAlgClosed.exists_pow_nat_eq (k := F) c (n := n) (by omega)
  have he0 : e ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at he
    exact hc he.symm
  have hem : algebraMap F W.FunctionField e ≠ 0 := by simpa using he0
  have hq : (h / algebraMap F W.FunctionField e) ^ n = 1 := by
    rw [div_pow, hpow, ← map_pow, he, div_self (by simpa using hc)]
  obtain ⟨k, hk0, hk⟩ := const_of_pow_eq_one W h2 hn hq
  refine ⟨e * k, mul_ne_zero he0 hk0, ?_⟩
  have hd := (div_eq_iff hem).1 hk
  rw [hd, map_mul, mul_comm]

/-- ★★★★★★**`τ(g)/g` は `g` の取り方(定数倍まで)に依らない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`f_P` の取り替え(定数倍)と `g_P` の取り替え(1 の `n` 乗根倍)を同時に吸収する。 -/
theorem aut_div_eq_of_const_mul (h2 : IsUnit (2 : F))
    {τ : W.FunctionField ≃ₐ[F] W.FunctionField}
    {g g' z : W.FunctionField} {c : F} {n : ℕ} (hn : 1 ≤ n) (hc : c ≠ 0)
    (hz : z ≠ 0) (hg : g ^ n = z)
    (hg' : g' ^ n = algebraMap F W.FunctionField c * z) :
    τ g' / g' = τ g / g := by
  have hgne : g ≠ 0 := by
    intro h0; rw [h0, zero_pow (by omega : n ≠ 0)] at hg; exact hz hg.symm
  have hq : (g' / g) ^ n = algebraMap F W.FunctionField c := by
    rw [div_pow, hg, hg', mul_div_assoc, div_self hz, mul_one]
  obtain ⟨d, hd0, hd⟩ := const_of_pow_eq_const W h2 hn hc hq
  have hdm : algebraMap F W.FunctionField d ≠ 0 := by simpa using hd0
  have hg'eq : g' = algebraMap F W.FunctionField d * g := by
    rw [← hd]; field_simp
  rw [hg'eq, map_mul, τ.commutes, mul_div_mul_left _ _ hdm]

/-! ## ★出典の紐付け(`.src`) -/

def aut_ext.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——自己同型が座標での値で決まること)",
    sectionId := "genell-thm-3-8" }

def aut_div_eq_of_const_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——値が g の選択に依らないこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
