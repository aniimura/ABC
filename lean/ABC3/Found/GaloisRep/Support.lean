import ABC3.Found.GaloisRep.FiberSum

/-!
# Galois (G5) 第 151 ブロック —— **★★★★★★因子の台の同定**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★剰余体を経由せずに台が書けた

D2 の残りは「`count_v(μ f_P) ≠ 0` となる `v` はどれか」である。
★§9-461 では**点の還元(特殊化)写像**が要ると見積もった(8-20 ブロック)。
★★**剰余体を経由しない書き方がある**——**引き戻した素イデアル**

    P' := { a ∈ F[W] | w(μ a) < 1 }

を使う。これは `w(μ ·) ≤ 1`(場合 A)のとき**素イデアル**であり、

    count_v(μ f_P) ≠ 0   ⟺   I_P = P'

★★★これで台が**イデアルの等式**として書ける。剰余体も点の還元も要らない。

## ★★★★機構

第 143 の `w(μ f) = max(w(μa), w(μb))^n` を `count` に直すと

    count(μ f_P) = n · min(count(μ a), count(μ b))

場合 A では `count ≥ 0` なので

    count(μ f_P) ≠ 0 ⟺ count(μ a) > 0 かつ count(μ b) > 0
                     ⟺ a, b ∈ P'  ⟺  I_P ⊆ P'  ⟺  I_P = P'

★最後の同値は `I_P` が極大(第 138)で `P' ≠ ⊤` だから。

## ★★★D2 に残るもの

`P'` が定める点が `n·Q_v` であること——ここで初めて群法則が要る。
★台の同定そのものは本ブロックで済んだ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackPrime` | ★★★★引き戻した素イデアル `{a \| w(μ a) < 1}` |
| `pullbackPrime_isPrime` | ★★★素であること |
| `count_pos_iff` | ★★`0 < count ⟺ w < 1` |
| `count_eq_min_of_span_pair_pow` | ★★★★★`count(μ f) = n·min(count(μa), count(μb))` |
| `count_ne_zero_iff` | ★★★★★★**台の同定 `count ≠ 0 ⟺ I_P = P'`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★★★引き戻した素イデアル -/

section Pullback

variable {R K Γ₀ : Type} [CommRing R] [Field K] [LinearOrderedCommGroupWithZero Γ₀]

/-- ★★★★**引き戻した素イデアル** `{a | w(μ a) < 1}`。 -/
def pullbackPrime (w : Valuation K Γ₀) (μ : R →+* K) (hnn : ∀ r : R, w (μ r) ≤ 1) : Ideal R where
  carrier := {r | w (μ r) < 1}
  zero_mem' := by simp
  add_mem' := by
    intro x y hx hy
    simp only [Set.mem_setOf_eq] at *
    have hμ : μ (x + y) = μ x + μ y := map_add μ x y
    rw [hμ]
    exact lt_of_le_of_lt (Valuation.map_add w _ _) (max_lt hx hy)
  smul_mem' := by
    intro c x hx
    simp only [Set.mem_setOf_eq] at *
    have hμ : μ (c • x) = μ c * μ x := by rw [smul_eq_mul, map_mul]
    rw [hμ, Valuation.map_mul]
    calc w (μ c) * w (μ x) ≤ 1 * w (μ x) := mul_le_mul' (hnn c) le_rfl
      _ = w (μ x) := one_mul _
      _ < 1 := hx

theorem mem_pullbackPrime {w : Valuation K Γ₀} {μ : R →+* K} {hnn : ∀ r : R, w (μ r) ≤ 1}
    {r : R} : r ∈ pullbackPrime w μ hnn ↔ w (μ r) < 1 := Iff.rfl

/-- ★★★引き戻した素イデアルは実際に素である。 -/
theorem pullbackPrime_isPrime (w : Valuation K Γ₀) (μ : R →+* K) (hnn : ∀ r : R, w (μ r) ≤ 1) :
    (pullbackPrime w μ hnn).IsPrime := by
  constructor
  · intro htop
    have h1 : (1 : R) ∈ pullbackPrime w μ hnn := by rw [htop]; trivial
    rw [mem_pullbackPrime, map_one, Valuation.map_one] at h1
    exact absurd h1 (lt_irrefl 1)
  · intro a b hab
    rw [mem_pullbackPrime, map_mul, Valuation.map_mul] at hab
    by_contra hcon
    rw [not_or] at hcon
    obtain ⟨ha, hb⟩ := hcon
    rw [mem_pullbackPrime, not_lt] at ha hb
    have ha1 : w (μ a) = 1 := le_antisymm (hnn a) ha
    have hb1 : w (μ b) = 1 := le_antisymm (hnn b) hb
    rw [ha1, hb1, one_mul] at hab
    exact absurd hab (lt_irrefl 1)

end Pullback

/-! ## ★★`count` と付値 -/

section Count

variable {R K : Type} [CommRing R] [IsDedekindDomain R] [Field K] [Algebra R K]
  [IsFractionRing R K]

theorem count_pos_iff (v : HeightOneSpectrum R) (μ : R →+* K) {a : R} (ha : μ a ≠ 0) :
    0 < FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ a))
      ↔ v.valuation K (μ a) < 1 := by
  rw [valuation_eq_exp_neg_count (K := K) v ha]
  constructor
  · intro hpos
    have hlt : WithZero.exp (-FractionalIdeal.count K v
        (FractionalIdeal.spanSingleton R⁰ (μ a))) < WithZero.exp (0 : ℤ) :=
      WithZero.exp_lt_exp.2 (by omega)
    rwa [WithZero.exp_zero] at hlt
  · intro hlt
    rw [← WithZero.exp_zero] at hlt
    have := WithZero.exp_lt_exp.1 hlt
    omega

/-- ★★★★★`count(μ f) = n · min(count(μ a), count(μ b))`(場合 A)。 -/
theorem count_eq_min_of_span_pair_pow (v : HeightOneSpectrum R) (μ : R →+* K)
    (hμinj : Function.Injective μ) (hnn : ∀ r : R, v.valuation K (μ r) ≤ 1)
    {a b f : R} (ha : a ≠ 0) (hb : b ≠ 0) (hf0 : f ≠ 0) {n : ℕ}
    (hf : (Ideal.span ({a, b} : Set R)) ^ n = Ideal.span {f}) :
    FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ f))
      = n * min (FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ a)))
        (FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ b))) := by
  have hne : ∀ {r : R}, r ≠ 0 → μ r ≠ 0 := fun {r} hr h => hr (hμinj (by rw [h, map_zero]))
  set ca := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ a)) with hca
  set cb := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ b)) with hcb
  set cf := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ f)) with hcf
  have hkey := valuation_pow_of_span_pair_pow (v.valuation K) μ hnn hf
  rw [valuation_eq_exp_neg_count (K := K) v (hne ha),
    valuation_eq_exp_neg_count (K := K) v (hne hb),
    valuation_eq_exp_neg_count (K := K) v (hne hf0)] at hkey
  have hmax : max (WithZero.exp (-ca)) (WithZero.exp (-cb)) = WithZero.exp (-(min ca cb)) := by
    rcases le_total ca cb with hle | hle
    · rw [min_eq_left hle, max_eq_left (WithZero.exp_le_exp.2 (by omega))]
    · rw [min_eq_right hle, max_eq_right (WithZero.exp_le_exp.2 (by omega))]
  rw [hmax, ← WithZero.exp_nsmul, WithZero.exp_inj, nsmul_eq_mul, mul_neg] at hkey
  exact neg_inj.1 hkey

end Count

/-! ## ★★★★★★台の同定 -/

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

theorem count_nonneg_of_nonneg (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1)
    {a : W.CoordinateRing} (ha : μ a ≠ 0) :
    0 ≤ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ a)) := by
  have hle := hnn a
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v ha, ← WithZero.exp_zero,
    WithZero.exp_le_exp] at hle
  omega

/-- ★★★★★★**因子の台の同定**(場合 A)——`count_v(μ f_P) ≠ 0 ⟺ I_P = P'`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★剰余体も点の還元も経由しない。**イデアルの等式**として書ける。
★★`I_P` が極大(第 138)で `P' ≠ ⊤` だから、包含から等式が出る。 -/
theorem count_ne_zero_iff (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1)
    {x y : F} (h : W.Nonsingular x y) {n : ℕ} (hn : 1 ≤ n) (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) ≠ 0
      ↔ CoordinateRing.XYIdeal W x (Polynomial.C y)
        = pullbackPrime (v.valuation W.FunctionField) μ hnn := by
  have hne : ∀ {r : W.CoordinateRing}, r ≠ 0 → μ r ≠ 0 :=
    fun {r} hr hz => hr (hμinj (by rw [hz, map_zero]))
  have hane : CoordinateRing.XClass W x ≠ 0 := CoordinateRing.XClass_ne_zero x
  have hbne : CoordinateRing.YClass W (Polynomial.C y) ≠ 0 :=
    CoordinateRing.YClass_ne_zero (Polynomial.C y)
  have hspan : (Ideal.span ({CoordinateRing.XClass W x,
      CoordinateRing.YClass W (Polynomial.C y)} : Set W.CoordinateRing)) ^ n
      = Ideal.span {fP} := by rw [← hfP, CoordinateRing.XYIdeal]
  have hcnt := count_eq_min_of_span_pair_pow v μ hμinj hnn hane hbne
    (fP_ne_zero W n fP hfP) hspan
  have hca := count_nonneg_of_nonneg W v μ hnn (hne hane)
  have hcb := count_nonneg_of_nonneg W v μ hnn (hne hbne)
  rw [hcnt]
  constructor
  · intro hnz
    have hmin : 0 < min (FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (CoordinateRing.XClass W x))))
        (FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (μ (CoordinateRing.YClass W (Polynomial.C y))))) := by
      rcases lt_or_eq_of_le (le_min hca hcb) with hlt | heq
      · exact hlt
      · exfalso; apply hnz; rw [← heq]; ring
    refine (xyIdeal_isMaximal W h.1).eq_of_le
      (pullbackPrime_isPrime (v.valuation W.FunctionField) μ hnn).ne_top ?_
    rw [CoordinateRing.XYIdeal, Ideal.span_le]
    rintro z (rfl | rfl)
    · exact (count_pos_iff v μ (hne hane)).1 (lt_of_lt_of_le hmin (min_le_left _ _))
    · exact (count_pos_iff v μ (hne hbne)).1 (lt_of_lt_of_le hmin (min_le_right _ _))
  · intro heq
    have hma : CoordinateRing.XClass W x ∈ pullbackPrime (v.valuation W.FunctionField) μ hnn := by
      rw [← heq, CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert _ _)
    have hmb : CoordinateRing.YClass W (Polynomial.C y)
        ∈ pullbackPrime (v.valuation W.FunctionField) μ hnn := by
      rw [← heq, CoordinateRing.XYIdeal]
      exact Ideal.subset_span (Set.mem_insert_of_mem _ rfl)
    have h1 := (count_pos_iff v μ (hne hane)).2 hma
    have h2 := (count_pos_iff v μ (hne hbne)).2 hmb
    have hminpos : 0 < min (FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (CoordinateRing.XClass W x))))
        (FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰
          (μ (CoordinateRing.YClass W (Polynomial.C y))))) := lt_min h1 h2
    have hnpos : 0 < (n : ℤ) := by exact_mod_cast hn
    positivity

/-! ## ★出典の紐付け(`.src`) -/

def count_ne_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——div(f_P ∘ [n]) の台の同定)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
