import ABC3.Found.GaloisRep.WeilWellDef

/-!
# Galois (G5) 第 178 ブロック —— **★★★★★★★Weil 対の値が定義できた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★存在形の特徴づけで全域関数にする

スケルトンの `weilPairing` は `W.Point → W.Point → F` という**仮定なしの全域関数**である。
★そこで値を**特徴づける性質**を存在形で書く:

    WeilSpec W n P Q c :=
      ∃ (f_P) (μ) (g) (τ), [P, Q がアフィン点で、f_P が XYIdeal(P)^n の生成元、
                            μ = [n]^*、g^n = μ f_P、τ が Q の平行移動]
                           ∧ τ g / g = algebraMap c

そして

    weilPairingVal W n P Q := if h : ∃ c, WeilSpec W n P Q c then h.choose else 1

★★**存在形にしたのが要点**である。`P` や `Q` が `n` 等分点でないときは
データが無いので既定値 1 になり、`e_n(P,Q)^n = 1` のような性質が
**退化した場合でも成り立つ**。★★★∀ 形にすると空虚に真になって値が定まらない。

## ★★★★★★★値は一意である

第 177 の 3 つの補題がここで効く:

| 選択 | 一意性の根拠 |
|---|---|
| `x, y`(`P` の座標) | `Point.some` の単射性 |
| `x₀, y₀`(`Q` の座標) | 同上 |
| `xn, yn`(`n·生成点` の座標) | 同上 |
| `μ` | 第 119 `coordinateRing_hom_ext` |
| `τ` | 第 177 `aut_ext` |
| `f_P` | `span` が等しい ⟹ 同伴 ⟹ 定数倍(第 128) |
| `g` | 第 177 `aut_div_eq_of_const_mul` |

★したがって `weilPairingVal` は `WeilSpec` を満たす唯一の値を返す
(`weilPairingVal_eq`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `WeilSpec` | ★★★★★★値を特徴づける性質 |
| `weilPairingVal` | ★★★★★★**Weil 対の値**(全域関数) |
| `weilSpec_unique` | ★★★★★★★**値は一意** |
| `weilPairingVal_eq` | ★★★★★★★特徴づけを満たす値は Weil 対の値 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F]

/-! ## ★★★★★★特徴づけ -/

/-- ★★★★★★**Weil 対の値を特徴づける性質**。

★存在形にしてあるので、データが無いとき(`P` や `Q` が `n` 等分点でないとき)は
どの `c` についても偽になり、`weilPairingVal` は既定値 1 を返す。 -/
def WeilSpec (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (P Q : W.Point) (c : F) : Prop :=
  ∃ (x y : F) (hP : W.Nonsingular x y) (x₀ y₀ : F) (hQ : W.Nonsingular x₀ y₀)
    (fP : W.CoordinateRing) (μ : W.CoordinateRing →+* W.FunctionField)
    (g : W.FunctionField) (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    (xn yn : W.FunctionField)
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn),
    P = Point.some x y hP ∧ Q = Point.some x₀ y₀ hQ ∧
    (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP} ∧
    Function.Injective μ ∧
    (∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d) ∧
    μ (genX W) = xn ∧ μ (genY W) = yn ∧
    n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns ∧
    g ^ n = μ fP ∧
    τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀ ∧
    τ g / g = algebraMap F W.FunctionField c

/-- ★★★★★★**Weil 対の値**——`e_n(P,Q) = τ_Q(g_P)/g_P`。 -/
noncomputable def weilPairingVal (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (P Q : W.Point) : F :=
  open Classical in
  if h : ∃ c : F, WeilSpec W n P Q c then h.choose else 1

theorem weilPairingVal_spec (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (P Q : W.Point) (h : ∃ c : F, WeilSpec W n P Q c) :
    WeilSpec W n P Q (weilPairingVal W n P Q) := by
  classical
  rw [weilPairingVal, dif_pos h]
  exact h.choose_spec

theorem weilPairingVal_of_not (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (P Q : W.Point) (h : ¬ ∃ c : F, WeilSpec W n P Q c) :
    weilPairingVal W n P Q = 1 := by
  classical
  rw [weilPairingVal, dif_neg h]

/-! ## ★★★★★★★一意性 -/

variable [IsAlgClosed F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★★**Weil 対の値は一意である**(選択に依らない)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 119(`μ` の一意性)・第 177(`τ` の一意性、`g` の定数倍不変)・
第 128(生成元は定数倍)がここで合流する。 -/
theorem weilSpec_unique (h2 : IsUnit (2 : F)) {n : ℕ} (hn : 1 ≤ n)
    {P Q : W.Point} {c c' : F}
    (hc : WeilSpec W n P Q c) (hc' : WeilSpec W n P Q c') : c = c' := by
  obtain ⟨x, y, hP, x₀, y₀, hQ, fP, μ, g, τ, xn, yn, hns,
    hPe, hQe, hfP, hinj, hμF, hμx, hμy, hμP, hg, hτx, hτy, hdiv⟩ := hc
  obtain ⟨x', y', hP', x₀', y₀', hQ', fP', μ', g', τ', xn', yn', hns',
    hPe', hQe', hfP', hinj', hμF', hμx', hμy', hμP', hg', hτx', hτy', hdiv'⟩ := hc'
  rw [hPe] at hPe'
  rw [hQe] at hQe'
  rw [hμP] at hμP'
  obtain ⟨hxx, hyy⟩ : x = x' ∧ y = y' := by injection hPe' with a b; exact ⟨a, b⟩
  obtain ⟨hx0, hy0⟩ : x₀ = x₀' ∧ y₀ = y₀' := by injection hQe' with a b; exact ⟨a, b⟩
  obtain ⟨hxn, hyn⟩ : xn = xn' ∧ yn = yn' := by injection hμP' with a b; exact ⟨a, b⟩
  subst hxx; subst hyy; subst hx0; subst hy0; subst hxn; subst hyn
  have hμμ : μ = μ' :=
    coordinateRing_hom_ext _ _ (fun d => by rw [hμF d, hμF' d]) (by rw [hμx, hμx'])
      (by rw [hμy, hμy'])
  subst hμμ
  have hττ : τ = τ' := aut_ext W (by rw [hτx, hτx']) (by rw [hτy, hτy'])
  subst hττ
  have hassoc : Associated fP fP' :=
    Ideal.span_singleton_eq_span_singleton.1 (by rw [← hfP, hfP'])
  obtain ⟨u, hu⟩ := hassoc
  obtain ⟨a, ha0, hau⟩ := isUnit_coordinateRing u.isUnit
  have hfPne : fP ≠ 0 := fP_ne_zero W n fP hfP
  have hz : μ fP ≠ 0 := fun h0 => hfPne (hinj (by rw [h0, map_zero]))
  have hg'2 : g' ^ n = algebraMap F W.FunctionField a * μ fP := by
    rw [hg', ← hu, map_mul, hau, hμF a, mul_comm]
  have hkey := aut_div_eq_of_const_mul W (τ := τ) h2 hn ha0 hz hg hg'2
  rw [hdiv, hdiv'] at hkey
  exact (algebraMap F W.FunctionField).injective hkey.symm

/-- ★★★★★★★**特徴づけを満たす値は Weil 対の値そのもの**。 -/
theorem weilPairingVal_eq (h2 : IsUnit (2 : F)) {n : ℕ} (hn : 1 ≤ n)
    {P Q : W.Point} {c : F} (h : WeilSpec W n P Q c) :
    weilPairingVal W n P Q = c :=
  weilSpec_unique W h2 hn (weilPairingVal_spec W n P Q ⟨c, h⟩) h

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の値の定義)",
    sectionId := "genell-thm-3-8" }

def weilSpec_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の値が選択に依らないこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
