import ABC3.Found.GaloisRep.WeilBilinCore

/-!
# Galois (G5) 第 184 ブロック —— **★★★★★★★★双線型性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★5 性質の 2 つ目

    e_n(P₁, Q) · e_n(P₂, Q) = e_n(P₁+P₂, Q)

★第 183 の中核を `weilPairingVal` の言葉に翻訳する。

### ★★★★★同じ `μ`・同じ `τ` を使うのが要点

3 つの点 `P₁`・`P₂`・`P₁+P₂` について `WeilSpec` の witness を作るとき、
**同じ `μ` と同じ `τ`** を使う。★`μ` は `n` と曲線だけで決まり、`τ` は `Q` だけで
決まるので、これは自然に揃う。★★そのうえで第 183 の

    (τg₁/g₁)(τg₂/g₂) = τg₃/g₃

を当てると、`τg_i/g_i = algebraMap c_i` から `c₁ c₂ = c₃` が出る。

### ★★★これで `hμinj` を仮定に置かなくてよくなった

第 180 の `mulByN_injective` を内部で呼ぶので、単射性は仮定から消えた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `h_ne_zero_of_id` | ★イデアルの関係式から `h ≠ 0` |
| `weilPairingVal_mul` | ★★★★★★★★**双線型性(第 1 変数)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★イデアルの関係式から `h ≠ 0`。 -/
theorem h_ne_zero_of_id
    {x₁ y₁ x₂ y₂ x₃ y₃ : F} (h₁ : W.Nonsingular x₁ y₁) (h₂ : W.Nonsingular x₂ y₂)
    (h₃ : W.Nonsingular x₃ y₃) {hh : W.FunctionField}
    (hid : (CoordinateRing.XYIdeal' h₁ : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
      * (CoordinateRing.XYIdeal' h₂)
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ hh * (CoordinateRing.XYIdeal' h₃)) :
    hh ≠ 0 := by
  intro h0
  rw [h0, FractionalIdeal.spanSingleton_zero, zero_mul, ← Units.val_mul] at hid
  exact (Units.ne_zero _) hid

variable [IsAlgClosed F] [Infinite F] [W.IsElliptic]

/-- ★★★★★★★★**双線型性(第 1 変数)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★3 つの点について**同じ `μ`・同じ `τ`**で `WeilSpec` の witness を作り、
第 183 の中核を当てる。★★単射性は第 180 で内部的に得られるので仮定に無い。 -/
theorem weilPairingVal_mul (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    {x₁ y₁ x₂ y₂ x₃ y₃ : F} (hP₁ : W.Nonsingular x₁ y₁) (hP₂ : W.Nonsingular x₂ y₂)
    (hP₃ : W.Nonsingular x₃ y₃)
    (hsum : Point.some x₁ y₁ hP₁ + Point.some x₂ y₂ hP₂ = Point.some x₃ y₃ hP₃)
    (ht₁ : n • Point.some x₁ y₁ hP₁ = 0) (ht₂ : n • Point.some x₂ y₂ hP₂ = 0)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (hQt : n • Point.some x₀ y₀ hQ = 0) :
    weilPairingVal W n (Point.some x₁ y₁ hP₁) (Point.some x₀ y₀ hQ)
      * weilPairingVal W n (Point.some x₂ y₂ hP₂) (Point.some x₀ y₀ hQ)
      = weilPairingVal W n (Point.some x₃ y₃ hP₃) (Point.some x₀ y₀ hQ) := by
  have ht₃ : n • Point.some x₃ y₃ hP₃ = 0 := by
    rw [← hsum, nsmul_add, ht₁, ht₂, add_zero]
  have hinj := mulByN_injective W h2 n hn (hchar n hn le_rfl) μ hμF hns hμx hμy hμP
  obtain ⟨f₁, hf₁ne, hf₁⟩ := xyIdeal_pow_isPrincipal_integral hP₁ n ht₁
  obtain ⟨f₂, hf₂ne, hf₂⟩ := xyIdeal_pow_isPrincipal_integral hP₂ n ht₂
  obtain ⟨f₃, hf₃ne, hf₃⟩ := xyIdeal_pow_isPrincipal_integral hP₃ n ht₃
  obtain ⟨hh, hid⟩ := exists_h_of_add W hP₁ hP₂ hP₃ hsum
  have hhne := h_ne_zero_of_id W hP₁ hP₂ hP₃ hid
  obtain ⟨u, hu0, hrel⟩ := elem_relation_of_add W hP₁ hP₂ hP₃ hf₁ hf₂ hf₃ hid
  obtain ⟨g₁, hg₁⟩ := exists_nthRoot_mulByN W h2 hP₁ n hn hchar ht₁ μ hinj hμF hns hμP hμx hμy
    f₁ hf₁
  obtain ⟨g₂, hg₂⟩ := exists_nthRoot_mulByN W h2 hP₂ n hn hchar ht₂ μ hinj hμF hns hμP hμx hμy
    f₂ hf₂
  obtain ⟨g₃, hg₃⟩ := exists_nthRoot_mulByN W h2 hP₃ n hn hchar ht₃ μ hinj hμF hns hμP hμx hμy
    f₃ hf₃
  obtain ⟨τ, hτx, hτy⟩ := exists_translateAut_all W h4 hQ
  have hcomp := aut_comp_mulByN W τ hQ hτx hτy n hQt μ hμF hns hμx hμy hμP
  have hτz : ∀ r : W.CoordinateRing, τ (μ r) = μ r := fun r =>
    congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hcomp) r
  have hne : ∀ f : W.CoordinateRing, f ≠ 0 → μ f ≠ 0 :=
    fun f hf h0 => hf (hinj (by rw [h0, map_zero]))
  have hcore := aut_div_mul_of_relation W h2 hn τ hinj hμF hcomp hhne hu0 hrel hg₁ hg₂ hg₃
    (hne f₃ hf₃ne)
  obtain ⟨c₁, hc₁0, hd₁⟩ := exists_const_aut_div W h2 hn hg₁ (hne f₁ hf₁ne) (hτz f₁)
  obtain ⟨c₂, hc₂0, hd₂⟩ := exists_const_aut_div W h2 hn hg₂ (hne f₂ hf₂ne) (hτz f₂)
  obtain ⟨c₃, hc₃0, hd₃⟩ := exists_const_aut_div W h2 hn hg₃ (hne f₃ hf₃ne) (hτz f₃)
  have hval : ∀ {x y : F} (hPx : W.Nonsingular x y) (f : W.CoordinateRing)
      (hf : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {f})
      (g : W.FunctionField) (hg : g ^ n = μ f) {c : F}
      (hd : τ g / g = algebraMap F W.FunctionField c),
      weilPairingVal W n (Point.some x y hPx) (Point.some x₀ y₀ hQ) = c := by
    intro x y hPx f hf g hg c hd
    exact weilPairingVal_eq W h2 hn
      ⟨x, y, hPx, x₀, y₀, hQ, f, μ, g, τ, xn, yn, hns, rfl, rfl, hf, hinj, hμF,
        hμx, hμy, hμP, hg, hτx, hτy, hd⟩
  rw [hval hP₁ f₁ hf₁ g₁ hg₁ hd₁, hval hP₂ f₂ hf₂ g₂ hg₂ hd₂, hval hP₃ f₃ hf₃ g₃ hg₃ hd₃]
  rw [hd₁, hd₂, hd₃, ← map_mul] at hcore
  exact (algebraMap F W.FunctionField).injective hcore

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
