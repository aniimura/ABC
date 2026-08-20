import ABC3.Found.GaloisRep.Divisible

/-!
# Galois (G5) 第 187 ブロック —— **★★★★★★★★無限遠の極の位数は `−n` ちょうど**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★交代性を止めていた唯一の量

交代性は `∏_{i=0}^{n-1} τ_{iP}^*(f_P)` が定数であることを使う。★指数で書くと

    count_v(τ_T f_P) = n·[Q_v + T = P] − n·[Q_v + T = O]

の第 2 項——**平行移動で無限遠に落ちる素点での極の位数**——だけが未知だった。
★第 185 の `count_translate_gen` は `Q_v + T = O` のとき `v'` が存在しないので使えない。
★★Abel–Jacobi(第 166)も**この量を決められない**——類群には次数の情報が無いからである
(`(n + a)·T = 0` しか出ず、`a = −n` は出ない)。

## ★★★★★★★★値群の全射性で決まる

`w := count_v ∘ τ` は `W.FunctionField` 上の**全射**離散付値である(`τ` は自己同型)。
★一方、第 143-149(場合 B)の道具立てから、`count_v(τ x) = 2k` と置くと

| 材料 | 出どころ |
|---|---|
| `count(τ (p(x))) = deg p · 2k` | `count_algebraMap_poly`(第 145) |
| `2·count(τ y) = 3·count(τ x)`、すなわち `count(τ y) = 3k` | `count_genY_relation` |
| `2 deg p ≠ 2 deg q + 3` で 2 項は決して打ち消さない | `valuation_poly_ne_genY_term` |

より **`count_v(τ a) ∈ k·ℤ`**(`a ∈ F[W]∖0`)。★分数体に延ばしても同じ。
★★`w` は全射なので `k ∣ 1`、`count_v(τ x) < 0` から `k < 0`、よって **`k = −1`**。

### ★★★★★★そこから `−n` が出る

第 146-148 の場合 B の計算は、実は**等式**を与えている:

    2·count(τ f_P) = n·count(τ x)

(`f_P · f_{−P} = c(x − x_P)^n` に `count` を当て、超楕円対合で `f_{−P}` を `f_P` に潰す)。
★`count(τ x) = −2` を入れて **`count(τ f_P) = −n`**。

### ★★★場合 B に入ることは幾何で判定できる

`Q_v + T = O` なら `redPoint(autFF τ (generic)) = 0`、第 163 の `redPoint_eq_zero_iff` で
`1 < v(τ x)`——すなわち**場合 B**。★第 170 の `aut_transport_hyps` の中で使われていた
補題が、そのまま逆向きに効く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `count_eq_of_valuation_eq` | 付値が等しければ `count` も等しい |
| `count_div'` | `count(a/b) = count a − count b` |
| `exists_count_eq_one` | `count` は 1 を取る(単項化元) |
| `autHom` | 自己同型から作る `μ` |
| `dvd_count_case_B_all` | ★★★★★★**場合 B の値群は `k·ℤ`** |
| `count_autHom_genX` | ★★★★★★★★**`count(τ x) = −2`** |
| `two_mul_count_case_B` | ★★★★★★**場合 B の等式版** |
| `count_translate_infty` | ★★★★★★★★**極の位数は `−n`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★`count` の細工 -/

/-- 付値が等しければ `count` も等しい。 -/
theorem count_eq_of_valuation_eq (v : HeightOneSpectrum W.CoordinateRing)
    {a b : W.FunctionField} (ha : a ≠ 0) (hb : b ≠ 0)
    (h : v.valuation W.FunctionField a = v.valuation W.FunctionField b) :
    FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton W.CoordinateRing⁰ a)
      = FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ b) := by
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v ha,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hb, WithZero.exp_inj] at h
  omega

/-- `count(a/b) = count a − count b`。 -/
theorem count_div' (v : HeightOneSpectrum W.CoordinateRing) {a b : W.FunctionField}
    (ha : a ≠ 0) (hb : b ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (a / b))
      = FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ a)
        - FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ b) := by
  have hid : a = (a / b) * b := by field_simp
  have hdiv : a / b ≠ 0 := div_ne_zero ha hb
  have hmul := count_mul' W v hdiv hb
  rw [← hid] at hmul
  omega

/-- ★★★★★**`count` は 1 を取る**——単項化元があるから。 -/
theorem exists_count_eq_one (v : HeightOneSpectrum W.CoordinateRing) :
    ∃ z : W.FunctionField, z ≠ 0 ∧
      FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z) = 1 := by
  obtain ⟨π, hπ⟩ := v.intValuation_exists_uniformizer
  have hπ0 : π ≠ 0 := by
    intro h0
    rw [h0, (v.intValuation).map_zero] at hπ
    exact absurd hπ.symm (by simp)
  have hz : algebraMap W.CoordinateRing W.FunctionField π ≠ 0 := fun h0 =>
    hπ0 (IsFractionRing.injective W.CoordinateRing W.FunctionField (by rw [h0, map_zero]))
  refine ⟨algebraMap W.CoordinateRing W.FunctionField π, hz, ?_⟩
  have hval : v.valuation W.FunctionField (algebraMap W.CoordinateRing W.FunctionField π)
      = WithZero.exp (-1) := by
    rw [v.valuation_of_algebraMap]; exact hπ
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hz, WithZero.exp_inj] at hval
  omega

/-! ## ★自己同型から作る `μ` -/

/-- 自己同型 `τ` から作る `μ : F[W] →+* F(W)`。 -/
noncomputable def autHom (τ : W.FunctionField ≃ₐ[F] W.FunctionField) :
    W.CoordinateRing →+* W.FunctionField :=
  (τ.toAlgHom.toRingHom).comp (algebraMap W.CoordinateRing W.FunctionField)

theorem autHom_apply (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (a : W.CoordinateRing) :
    autHom W τ a = τ (algebraMap W.CoordinateRing W.FunctionField a) := rfl

theorem autHom_injective (τ : W.FunctionField ≃ₐ[F] W.FunctionField) :
    Function.Injective (autHom W τ) :=
  τ.injective.comp (IsFractionRing.injective W.CoordinateRing W.FunctionField)

theorem autHom_algebraMap (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (c : F) :
    autHom W τ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c := by
  rw [autHom_apply, ← IsScalarTower.algebraMap_apply, AlgEquiv.commutes]

/-! ## ★★★★★★場合 B の値群 -/

/-- ★★★★★★**場合 B の値群は `k·ℤ` に入る**(`count(μ x) = 2k`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★基底 `{1, y}` で分解し、2 項の付値が決して一致しないことから和の `count` は
どちらか一方に一致する。★★`count(p(x)) = deg p·2k`、`count(q(x)y) = (2 deg q + 3)k`。 -/
theorem dvd_count_case_B_all (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) {k : ℤ}
    (hk : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) = 2 * k)
    {a : W.CoordinateRing} (ha : a ≠ 0) :
    k ∣ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ a)) := by
  set w := v.valuation W.FunctionField with hw
  have hcy : 2 * FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genY W))) = 3 * (2 * k) := by
    rw [← hk]; exact count_genY_relation W h2 v μ hμF hx
  have hcy' : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genY W))) = 3 * k := by omega
  obtain ⟨p, q, hpq⟩ := CoordinateRing.exists_smul_basis_eq a
  have hdec : a = algebraMap (Polynomial F) W.CoordinateRing p
      + algebraMap (Polynomial F) W.CoordinateRing q * genY W := by
    rw [← hpq, Algebra.smul_def, Algebra.smul_def, mul_one]; rfl
  set A := μ (algebraMap (Polynomial F) W.CoordinateRing p) with hA
  set B := μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W) with hB
  have hsum : μ a = A + B := by rw [hdec, map_add, map_mul]
  have hcA : ∀ _ : p ≠ 0, FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ A) = (p.natDegree : ℤ) * (2 * k) := by
    intro hp; rw [hA, count_algebraMap_poly W v μ hμF hx hp, hk]
  have hyne : μ (genY W) ≠ 0 := fun h0 =>
    valuation_genY_ne_zero W h2 v μ hμF hx (by rw [h0, Valuation.map_zero])
  have hcB : ∀ _ : q ≠ 0, FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ B)
        = (q.natDegree : ℤ) * (2 * k) + 3 * k := by
    intro hq
    have hqne : μ (algebraMap (Polynomial F) W.CoordinateRing q) ≠ 0 := fun h0 =>
      valuation_algebraMap_poly_ne_zero W v μ hμF hx hq (by rw [h0, Valuation.map_zero])
    rw [hB, count_mul' W v hqne hyne, count_algebraMap_poly W v μ hμF hx hq, hk, hcy']
  by_cases hp : p = 0
  · by_cases hq : q = 0
    · exact absurd (by rw [hdec, hp, hq]; simp) ha
    · have hval : μ a = B := by rw [hsum, hA, hp, map_zero, map_zero, zero_add]
      rw [hval, hcB hq]; exact ⟨(q.natDegree : ℤ) * 2 + 3, by ring⟩
  · by_cases hq : q = 0
    · have hval : μ a = A := by rw [hsum, hB, hq, map_zero, map_zero, zero_mul, add_zero]
      rw [hval, hcA hp]; exact ⟨(p.natDegree : ℤ) * 2, by ring⟩
    · have hAne : A ≠ 0 := fun h0 =>
        valuation_algebraMap_poly_ne_zero W v μ hμF hx hp
          (by rw [hA] at h0; rw [h0, Valuation.map_zero])
      have hqne : μ (algebraMap (Polynomial F) W.CoordinateRing q) ≠ 0 := fun h0 =>
        valuation_algebraMap_poly_ne_zero W v μ hμF hx hq (by rw [h0, Valuation.map_zero])
      have hBne : B ≠ 0 := mul_ne_zero hqne hyne
      have hne : w A ≠ w B := valuation_poly_ne_genY_term W h2 v μ hμF hx hp hq
      have hmax : w (μ a) = max (w A) (w B) := by
        rw [hsum]; exact Valuation.map_add_of_distinct_val w hne
      have hane : μ a ≠ 0 := fun h0 => ha (hμinj (by rw [h0, map_zero]))
      rcases lt_or_gt_of_ne hne with hlt | hgt
      · rw [max_eq_right (le_of_lt hlt)] at hmax
        rw [count_eq_of_valuation_eq W v hane hBne hmax, hcB hq]
        exact ⟨(q.natDegree : ℤ) * 2 + 3, by ring⟩
      · rw [max_eq_left (le_of_lt hgt)] at hmax
        rw [count_eq_of_valuation_eq W v hane hAne hmax, hcA hp]
        exact ⟨(p.natDegree : ℤ) * 2, by ring⟩

/-- ★★★★★★★★**自己同型の場合、`count(τ x) = −2` ちょうど**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`count_v ∘ τ` は全射な離散付値なので、値群が `k·ℤ` に入ることから `k ∣ 1`。
★★`count(τ x) < 0` と合わせて `k = −1`。 -/
theorem count_autHom_genX (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    (hx : 1 < v.valuation W.FunctionField (autHom W τ (genX W))) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (autHom W τ (genX W))) = -2 := by
  set μ := autHom W τ with hμ
  have hμF := autHom_algebraMap W τ
  have hμinj := autHom_injective W τ
  obtain ⟨k, hk⟩ := even_count_genX W h2 v μ hμF hx
  have hall : ∀ z : W.FunctionField, z ≠ 0 → k ∣ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (τ z)) := by
    intro z hz
    obtain ⟨a, b, hb, hab⟩ := IsFractionRing.div_surjective (A := W.CoordinateRing) z
    have hb0 : algebraMap W.CoordinateRing W.FunctionField b ≠ 0 := by
      intro h0
      have hbz : b = 0 :=
        IsFractionRing.injective W.CoordinateRing W.FunctionField (by rw [h0, map_zero])
      exact absurd (hbz ▸ hb) (by simpa using (zero_notMem_nonZeroDivisors (M₀ := W.CoordinateRing)))
    have ha0 : algebraMap W.CoordinateRing W.FunctionField a ≠ 0 := by
      intro h0; rw [h0, zero_div] at hab; exact hz hab.symm
    have hbne : b ≠ 0 := fun h0 => hb0 (by rw [h0, map_zero])
    have hane : a ≠ 0 := fun h0 => ha0 (by rw [h0, map_zero])
    have hτ : τ z = μ a / μ b := by
      rw [← hab, map_div₀, hμ, autHom_apply, autHom_apply]
    rw [hτ, count_div' W v (fun h0 => hane (hμinj (by rw [h0, map_zero])))
      (fun h0 => hbne (hμinj (by rw [h0, map_zero])))]
    exact dvd_sub (dvd_count_case_B_all W h2 v μ hμinj hμF hx hk hane)
      (dvd_count_case_B_all W h2 v μ hμinj hμF hx hk hbne)
  obtain ⟨z₀, hz₀, hc₀⟩ := exists_count_eq_one W v
  have hk1 : k ∣ 1 := by
    have hd := hall (τ.symm z₀) (by simp [hz₀])
    rwa [AlgEquiv.apply_symm_apply, hc₀] at hd
  have hneg : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) < 0 :=
    count_genX_neg W v μ hx
  rw [hk] at hneg ⊢
  have hkm : k = -1 := by
    rcases Int.isUnit_iff.1 (isUnit_of_dvd_one hk1) with hone | hmone
    · omega
    · exact hmone
  omega

/-! ## ★★★★★★場合 B の等式版 -/

/-- ★★★★★★**場合 B の等式版**——`2·count(μ f_P) = n·count(μ x)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 149 の `dvd_count_case_B` は可除性だけを返していたが、中身は等式である。 -/
theorem two_mul_count_case_B (h2 : IsUnit (2 : F)) [W.IsElliptic]
    (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {x y : F} (h : W.Nonsingular x y) (n : ℕ) (fP fN : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (hfN : (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN}) :
    2 * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP))
      = (n : ℤ) * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
  obtain ⟨c, hc0, hprod⟩ := mu_fP_mul_mu_fNegP W h n fP fN hfP hfN μ hμF
  have hfP0 : fP ≠ 0 := fP_ne_zero W n fP hfP
  have hfN0 : fN ≠ 0 := fP_ne_zero W n fN hfN
  have hmP : μ fP ≠ 0 := fun h0 => hfP0 (hμinj (by rw [h0, map_zero]))
  have hmN : μ fN ≠ 0 := fun h0 => hfN0 (hμinj (by rw [h0, map_zero]))
  have hcK : algebraMap F W.FunctionField c ≠ 0 := fun h0 => hc0
    ((algebraMap F W.FunctionField).injective (by rw [h0, map_zero]))
  have hxc : μ (genX W) - algebraMap F W.FunctionField x ≠ 0 := by
    intro h0
    have hcs := count_genX_sub_const W v μ hx x
    rw [h0] at hcs
    simp only [FractionalIdeal.spanSingleton_zero, FractionalIdeal.count_zero] at hcs
    have hneg := count_genX_neg W v μ hx
    omega
  have hcount := congrArg (fun t => FractionalIdeal.count W.FunctionField v
    (FractionalIdeal.spanSingleton W.CoordinateRing⁰ t)) hprod
  rw [count_mul' W v hmP hmN, count_mul' W v hcK (pow_ne_zero n hxc), count_pow' W v,
    count_algebraMap_const W v hc0, count_genX_sub_const W v μ hx x, zero_add] at hcount
  obtain ⟨u, hu⟩ := hyperInv_fP_assoc W n fP fN hfP hfN
  have hNeq : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fN))
      = FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) := by
    rw [← hu, map_mul, count_mul' W v
      (fun h0 => hmN (by rw [← hu, map_mul, h0, zero_mul]))
      (fun h0 => by
        have hun : IsUnit (μ (u : W.CoordinateRing)) := u.isUnit.map μ
        rw [h0] at hun
        exact absurd hun not_isUnit_zero),
      count_unit_image_eq_zero W v μ hμF u.isUnit, add_zero,
      count_hyperInv W h2 v μ hμinj hμF hx hfP0]
  rw [hNeq] at hcount
  omega

/-! ## ★★★★★★★★極の位数 -/

/-- ★★★★★★★★**無限遠に落ちる素点での指数は `−n` ちょうど**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Q_v + T = O` なら `redPoint = 0`、第 163 で `1 < v(τ x)`——場合 B に入る。
★★そこで `count(τ x) = −2` と等式版を合わせる。 -/
theorem count_translate_infty [DecidableEq F] [W.IsElliptic] (h2 : IsUnit (2 : F))
    (v : HeightOneSpectrum W.CoordinateRing)
    {c y₀ : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hxτ : τ (coordX W) = translateX W x₀ y₀T) (hyτ : τ (coordY W) = translateY W x₀ y₀T)
    (hzero : Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQ = 0)
    {x y : F} (hP : W.Nonsingular x y) (n : ℕ) (hPt : n • (Point.some x y hP) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (τ (algebraMap W.CoordinateRing W.FunctionField fP))) = -(n : ℤ) := by
  obtain ⟨hns1, hpt⟩ := autFF_some W τ (nonsingular_coord W)
  have hred : redPoint W v h hv (Point.some (τ (coordX W)) (τ (coordY W)) hns1) = 0 := by
    rw [← hpt, ← ABC3.Found.GaloisRep.genericPoint,
      redPoint_aut_generic W v h hv h2 τ hQ hxτ hyτ, hzero]
  have hx : 1 < v.valuation W.FunctionField (autHom W τ (genX W)) :=
    (redPoint_eq_zero_iff W v h hv hns1).1 hred
  have hgx := count_autHom_genX W h2 v τ hx
  obtain ⟨fN, _, hfN⟩ := exists_fNeg W hP n hPt
  have hkey := two_mul_count_case_B W h2 v (autHom W τ) (autHom_injective W τ)
    (autHom_algebraMap W τ) hx hP n fP fN hfP hfN
  rw [hgx] at hkey
  show FractionalIdeal.count W.FunctionField v
    (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (autHom W τ fP)) = -(n : ℤ)
  omega

/-! ## ★出典の紐付け(`.src`) -/

def count_autHom_genX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——無限遠の上の count(x) = −2)",
    sectionId := "genell-thm-3-8" }

def count_translate_infty.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——極の位数が −n であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
