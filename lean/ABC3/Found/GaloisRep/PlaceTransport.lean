import ABC3.Found.GaloisRep.AutInvariant

/-!
# Galois (G5) 第 169 ブロック —— **★★★★★★★素点の輸送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★自己同型は指数をそのまま移す

`e : K ≃+* K` が Dedekind 環 `A` の分数体の自己同型で、`w_v ∘ e` が `A` 上非負、
その中心が `v'` であるとする。★このとき

    count_v(e z) = count_{v'}(z)      (∀ z ≠ 0)

これが `e_v` の一定性の最後の道具である。

### ★★★★★★中心が `v'` の付値は `v'`-進付値の冪

一般に、`W'` を `K` 上の付値で

* `W'(a) ≤ 1`  (∀ a ∈ A)
* `W'(a) < 1 ⟺ a ∈ v'`  (∀ a ∈ A)

を満たすものとすると、一意化元 `π`(`v'.intValuation π = exp(-1)`)を取って

    W'(z) = W'(π) ^ ord_{v'}(z)

が成り立つ(`valuation_eq_zpow_of_center`)。★鍵は
**`ord_{v'}(u) = 0` なら `u = c/d`(`c, d ∉ v'`)** であり、
これは mathlib の `valuationSubringAtPrime_eq_valuationSubring` と
`IsLocalization (v'.asIdeal.primeCompl) (valuationSubringAtPrime K v')` から出る
(`exists_frac_of_valuation_le_one`)。

### ★★★★★`m = 1` は対称性から出る——全射性を示す必要がない

上の表示から `count_v(e z) = m · count_{v'}(z)`(`m = count_v(e π)`)。
★`e.symm` にも同じことを当てて `count_{v'}(e.symm w) = m' · count_v(w)`。
★★`w := e π` と置くと `1 = m'·m`、そして `π ∈ v'` から `m > 0` なので **`m = 1`**。

★★★通常の証明は「`W'` が全射」を経由するが、**`e` が同型であることを 2 回使えば
全射性は要らない**。

## ★★★これは (G7) でも効く

素点の輸送は半安定モデルの議論でもそのまま使う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_frac_of_valuation_le_one` | ★★★★★付値 ≤ 1 の元は `v` の外の分母で書ける |
| `count_algebraMap_uniformizer` | ★★一意化元の指数は 1 |
| `valuation_eq_zpow_of_center` | ★★★★★★★**中心が `v'` の付値は冪** |
| `count_comp_eq_mul` | ★★★★★★指数は `m` 倍になる |
| `count_comp_eq` | ★★★★★★★**素点の輸送**(`m = 1`) |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain nonZeroDivisors

variable {A K : Type} [CommRing A] [IsDedekindDomain A] [Field K] [Algebra A K]
  [IsFractionRing A K]

/-! ## ★★★★★局所化からの分数表示 -/

/-- ★★★★★**付値 ≤ 1 の元は `v` の外の分母で書ける**。

★mathlib の `valuationSubringAtPrime_eq_valuationSubring` と
`IsLocalization (v.asIdeal.primeCompl) (valuationSubringAtPrime K v)` から。 -/
theorem exists_frac_of_valuation_le_one (v : HeightOneSpectrum A) {u : K}
    (hu : v.valuation K u ≤ 1) :
    ∃ c d : A, d ∉ v.asIdeal ∧ u * algebraMap A K d = algebraMap A K c := by
  have hmem : u ∈ HeightOneSpectrum.valuationSubringAtPrime K v := by
    rw [HeightOneSpectrum.valuationSubringAtPrime_eq_valuationSubring]
    exact hu
  obtain ⟨⟨c, d⟩, hcd⟩ :=
    IsLocalization.surj (M := v.asIdeal.primeCompl)
      (S := HeightOneSpectrum.valuationSubringAtPrime K v) (⟨u, hmem⟩)
  refine ⟨c, d, d.2, ?_⟩
  have hq := congrArg (fun t : HeightOneSpectrum.valuationSubringAtPrime K v => (t : K)) hcd
  simpa [IsScalarTower.algebraMap_apply A (HeightOneSpectrum.valuationSubringAtPrime K v) K]
    using hq

/-- ★★一意化元の指数は 1。 -/
theorem count_algebraMap_uniformizer (v : HeightOneSpectrum A) {π : A}
    (hπ : v.intValuation π = WithZero.exp (-1)) :
    FractionalIdeal.count K v
      (FractionalIdeal.spanSingleton A⁰ (algebraMap A K π)) = 1 := by
  have hne : algebraMap A K π ≠ 0 := by
    intro h0
    rw [show π = 0 from IsFractionRing.injective A K (by rw [h0, map_zero]),
      Valuation.map_zero] at hπ
    exact absurd hπ.symm (by simp)
  have h1 := valuation_eq_exp_neg_count (K := K) v hne
  rw [HeightOneSpectrum.valuation_of_algebraMap, hπ] at h1
  have h2 := WithZero.exp_injective h1
  omega

/-! ## ★★★★★★★中心が `v'` の付値は `v'`-進付値の冪 -/

/-- ★★★★★★★**中心が `v` の付値は `v`-進付値の冪である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ord_v(u) = 0` なら `u = c/d`(`c, d ∉ v`)なので `W'(u) = 1`。
★★あとは一意化元の冪で書き分けるだけ。 -/
theorem valuation_eq_zpow_of_center {Γ₀ : Type} [LinearOrderedCommGroupWithZero Γ₀]
    (W' : Valuation K Γ₀) (v : HeightOneSpectrum A)
    (hle : ∀ a : A, W' (algebraMap A K a) ≤ 1)
    (hlt : ∀ a : A, W' (algebraMap A K a) < 1 ↔ a ∈ v.asIdeal)
    {π : A} (hπ : v.intValuation π = WithZero.exp (-1))
    {z : K} (hz : z ≠ 0) :
    W' z = (W' (algebraMap A K π))
      ^ (FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ z)) := by
  set k := FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ z) with hk
  set p := algebraMap A K π with hp
  have hp0 : p ≠ 0 := by
    intro h0
    have hz0 : v.intValuation π = 0 := by
      rw [show π = 0 from IsFractionRing.injective A K (by rw [← hp, h0, map_zero]),
        Valuation.map_zero]
    rw [hz0] at hπ
    exact absurd hπ.symm (by simp)
  have hwp : v.valuation K p = WithZero.exp (-1) := by
    rw [hp, HeightOneSpectrum.valuation_of_algebraMap, hπ]
  have hwz : v.valuation K z = WithZero.exp (-k) := valuation_eq_exp_neg_count v hz
  set u := z * p ^ (-k) with hu
  have hu0 : u ≠ 0 := mul_ne_zero hz (zpow_ne_zero _ hp0)
  have hwu : v.valuation K u = 1 := by
    rw [hu, Valuation.map_mul, map_zpow₀, hwz, hwp,
      show (WithZero.exp (-1 : ℤ)) ^ (-k) = WithZero.exp ((-k) * (-1)) from rfl,
      ← WithZero.exp_add]
    simp
  obtain ⟨c, d, hd, hcd⟩ := exists_frac_of_valuation_le_one v (le_of_eq hwu)
  have hwd : v.valuation K (algebraMap A K d) = 1 := by
    rcases lt_or_eq_of_le (HeightOneSpectrum.valuation_le_one (K := K) v d) with hlt' | heq
    · exact absurd ((HeightOneSpectrum.valuation_lt_one_iff_mem v d).1 hlt') hd
    · exact heq
  have hwc : v.valuation K (algebraMap A K c) = 1 := by
    rw [← hcd, Valuation.map_mul, hwu, hwd, one_mul]
  have hc : c ∉ v.asIdeal := by
    intro hmem
    exact absurd ((HeightOneSpectrum.valuation_lt_one_iff_mem v c).2 hmem) (by rw [hwc]; simp)
  have hWd : W' (algebraMap A K d) = 1 := by
    rcases lt_or_eq_of_le (hle d) with h' | h'
    · exact absurd ((hlt d).1 h') hd
    · exact h'
  have hWc : W' (algebraMap A K c) = 1 := by
    rcases lt_or_eq_of_le (hle c) with h' | h'
    · exact absurd ((hlt c).1 h') hc
    · exact h'
  have hWu : W' u = 1 := by
    have hq := congrArg W' hcd
    rw [Valuation.map_mul, hWd, mul_one, hWc] at hq
    exact hq
  have hz' : z = u * p ^ k := by
    rw [hu, mul_assoc, ← zpow_add₀ hp0]
    simp
  rw [hz', Valuation.map_mul, hWu, one_mul, map_zpow₀]

/-! ## ★★★★★★★素点の輸送 -/

/-- ★★★★★★指数は `m` 倍になる。 -/
theorem count_comp_eq_mul (σ : K →+* K) (hσ : Function.Injective σ)
    (v v' : HeightOneSpectrum A)
    (hle : ∀ a : A, v.valuation K (σ (algebraMap A K a)) ≤ 1)
    (hlt : ∀ a : A, v.valuation K (σ (algebraMap A K a)) < 1 ↔ a ∈ v'.asIdeal)
    {π : A} (hπ : v'.intValuation π = WithZero.exp (-1))
    {z : K} (hz : z ≠ 0) :
    FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ (σ z))
      = (FractionalIdeal.count K v
          (FractionalIdeal.spanSingleton A⁰ (σ (algebraMap A K π))))
        * FractionalIdeal.count K v' (FractionalIdeal.spanSingleton A⁰ z) := by
  have hkey := valuation_eq_zpow_of_center ((v.valuation K).comap σ) v' hle hlt hπ hz
  have hσz : σ z ≠ 0 := fun h0 => hz (hσ (by rw [h0, map_zero]))
  have hp0 : algebraMap A K π ≠ 0 := by
    intro h0
    have hq : v'.intValuation π = 0 := by
      rw [show π = 0 from IsFractionRing.injective A K (by rw [h0, map_zero]), Valuation.map_zero]
    rw [hq] at hπ
    exact absurd hπ.symm (by simp)
  have hσp : σ (algebraMap A K π) ≠ 0 := fun h0 => hp0 (hσ (by rw [h0, map_zero]))
  rw [show ((v.valuation K).comap σ) z = v.valuation K (σ z) from rfl,
    show ((v.valuation K).comap σ) (algebraMap A K π)
      = v.valuation K (σ (algebraMap A K π)) from rfl,
    valuation_eq_exp_neg_count v hσz, valuation_eq_exp_neg_count v hσp,
    show ∀ (a n : ℤ), (WithZero.exp a) ^ n = WithZero.exp (n * a) from fun _ _ => rfl] at hkey
  have hq := WithZero.exp_injective hkey
  linarith [hq]

/-- ★★★★★★★**素点の輸送**——自己同型は指数をそのまま移す。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`e` と `e.symm` の両方に `count_comp_eq_mul` を当てると `m·m' = 1`。
★★`π ∈ v'` から `m > 0` なので `m = 1`。**全射性を経由しない**。 -/
theorem count_comp_eq (e : K ≃+* K) (v v' : HeightOneSpectrum A)
    (hle : ∀ a : A, v.valuation K (e (algebraMap A K a)) ≤ 1)
    (hlt : ∀ a : A, v.valuation K (e (algebraMap A K a)) < 1 ↔ a ∈ v'.asIdeal)
    (hle' : ∀ a : A, v'.valuation K (e.symm (algebraMap A K a)) ≤ 1)
    (hlt' : ∀ a : A, v'.valuation K (e.symm (algebraMap A K a)) < 1 ↔ a ∈ v.asIdeal)
    {z : K} (hz : z ≠ 0) :
    FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ (e z))
      = FractionalIdeal.count K v' (FractionalIdeal.spanSingleton A⁰ z) := by
  obtain ⟨π, hπ⟩ := HeightOneSpectrum.intValuation_exists_uniformizer v'
  obtain ⟨ρ, hρ⟩ := HeightOneSpectrum.intValuation_exists_uniformizer v
  set m := FractionalIdeal.count K v
    (FractionalIdeal.spanSingleton A⁰ (e (algebraMap A K π))) with hm
  set m' := FractionalIdeal.count K v'
    (FractionalIdeal.spanSingleton A⁰ (e.symm (algebraMap A K ρ))) with hm'
  have key : ∀ {w : K}, w ≠ 0 →
      FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ (e w))
        = m * FractionalIdeal.count K v' (FractionalIdeal.spanSingleton A⁰ w) :=
    fun {w} hw => count_comp_eq_mul (e : K →+* K) e.injective v v' hle hlt hπ hw
  have key' : ∀ {w : K}, w ≠ 0 →
      FractionalIdeal.count K v' (FractionalIdeal.spanSingleton A⁰ (e.symm w))
        = m' * FractionalIdeal.count K v (FractionalIdeal.spanSingleton A⁰ w) :=
    fun {w} hw => count_comp_eq_mul (e.symm : K →+* K) e.symm.injective v' v hle' hlt' hρ hw
  have hp0 : algebraMap A K π ≠ 0 := by
    intro h0
    have hq : v'.intValuation π = 0 := by
      rw [show π = 0 from IsFractionRing.injective A K (by rw [h0, map_zero]), Valuation.map_zero]
    rw [hq] at hπ
    exact absurd hπ.symm (by simp)
  have hep : e (algebraMap A K π) ≠ 0 := by simpa using hp0
  have hcπ : FractionalIdeal.count K v' (FractionalIdeal.spanSingleton A⁰ (algebraMap A K π)) = 1 :=
    count_algebraMap_uniformizer v' hπ
  have hround := key' hep
  rw [RingEquiv.symm_apply_apply, hcπ, key hp0, hcπ, mul_one] at hround
  have hπlt : v'.intValuation π < 1 := by
    rw [hπ, show (1 : WithZero (Multiplicative ℤ)) = WithZero.exp 0 from rfl]
    exact WithZero.exp_lt_exp.2 (by norm_num)
  have hπmem : π ∈ v'.asIdeal := (HeightOneSpectrum.intValuation_lt_one_iff_mem v' π).1 hπlt
  have hmpos : 0 < m := by
    have hlt1 := (hlt π).2 hπmem
    rw [valuation_eq_exp_neg_count v hep, ← hm,
      show (1 : WithZero (Multiplicative ℤ)) = WithZero.exp 0 from rfl,
      WithZero.exp_lt_exp] at hlt1
    omega
  have hdvd : m ∣ 1 := ⟨m', by rw [hround]; ring⟩
  have hm1 : m = 1 := Int.eq_one_of_dvd_one (le_of_lt hmpos) hdvd
  rw [key hz, hm1, one_mul]

/-! ## ★出典の紐付け(`.src`) -/

def valuation_eq_zpow_of_center.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——中心が v の付値が v-進付値の冪であること)",
    sectionId := "genell-thm-3-8" }

def count_comp_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——自己同型による素点の輸送)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
