import ABC3.Found.HerbrandIndex
import ABC3.Found.PGC.UnitsGroupInvariants
import ABC3.Found.PGC.QuotientCardinality

/-!
# `p` 乗指数——`[𝒪_K^× : (𝒪_K^×)^p] = p^{[K:ℚ_p]} · #μ_p(K)`

[pGC] Proposition 1.2 への**経路 C**(`ResearchPaper/pgc-goal.md`)の (C-d) 側の
計数部分。`Found/PGC/UnitsPowPrimeToP.lean`(第 1016)が `p ∤ n` の場合
(`[𝒪_K^× : (𝒪_K^×)^n] = gcd(n, q−1)`)を片付けたので、残った `n = p` をここで
片付ける。**この二つで `n` のすべての素因数が尽きる**わけではないが、経路 C が
`q`(→ `p ∤ n` 側)と `[K:ℚ_p]`(→ `p` 側)を分けて読むのに必要な二本は揃う。

## なぜ `p ∤ n` と同じ手が使えないか

`p ∤ n` のときは Hensel が主単数群の中で `n` 乗写像を全射にした
(`UnitsPowPrimeToP.lean::pow_surjective_principalUnits`)。`n = p` では
還元 `X^p − 1 = (X−1)^p` が `1` を**重根**に持つので Hensel が効かず、
実際に主単数群の `p` 乗指数は `1` ではなく `p^{[K:ℚ_p]}` になる。
そこで別の道具——**Herbrand 商の有限指数不変性**(`Found/HerbrandIndex.lean`)
——を使う。

## 筋

1. `smallPrincipalUnits K`(`‖u−1‖ ≤ 1/4` の層、`Found/PGC/PrincipalUnitsLog.lean`)は
   `K^×` の部分群として定義されているが、`‖u−1‖ ≤ 1/4` なら `‖u‖ = 1` なので
   実は `𝒪_K^×` の中に住んでいる。これを `𝒪_K^×` の部分群として引き戻したものが
   `smallPrincipalUnitsO K`(`unitsToCarrier` による `comap`)で、
   `smallPrincipalUnitsOEquiv` が `≃* smallPrincipalUnits K` を与える。
2. **有限指数**:素元 `ϖ` を取り `‖ϖ‖^m ≤ 1/4` となる `m` を選ぶと
   `principalUnits K ϖ m ⊆ smallPrincipalUnitsO K`。前者は
   `Found/PGC/QuotientCardinality.lean::finite_principalUnitsQuotient` により
   有限指数だから、後者も有限指数(`Subgroup.index_dvd_of_le`)。
3. **`ℤ_p^d` 側の二つの量**:`smallPrincipalUnitsEquivPi`(在庫)で
   `≃* Multiplicative (Fin d → ℤ_[p])`。`ℤ_p` は整域で `p ≠ 0` だから
   `p` 捩れは自明(`#= 1`)、`p` 乗指数は `p^d`
   (`index_powRange_smallPrincipalUnits`、在庫)。
4. **Herbrand で移す**:`index_pow_mul_card_torsion` を
   `A = 𝒪_K^×`・`B = smallPrincipalUnitsO K` に当てると

   ```
   [𝒪_K^× : (𝒪_K^×)^p] · 1  =  p^{[K:ℚ_p]} · #μ_p(K)
   ```

   ——左辺の `#B[p]` が `1`、右辺の `[B : B^p]` が `p^d` である。

## ★退化の自己検査——どこを変えると偽になるか

この定理には仮定らしい仮定が無い(`K` が `p` 進局所体であることだけ)。
落とせるのは**結論の各因子**なので、そこを検査する。

* **`#μ_p(K)` を落として `= p^{[K:ℚ_p]}` にすると偽**。`K = ℚ_2` を取る。
  `d = [K:ℚ_2] = 1`、`μ_2(ℤ_2^×) = {±1}` なので右辺は `2 · 2 = 4`。
  実際 `ℤ_2^× ≅ {±1} × (1 + 4ℤ_2) ≅ ℤ/2 × ℤ_2` だから平方の指数は `2 · 2 = 4` で
  合う。`#μ_2` を落とした `2^1 = 2` は誤り。
  (`1 + 4ℤ_2` はまさに本ファイルの `smallPrincipalUnits`——`‖u−1‖ ≤ 1/4`——である。)
* **`p^{[K:ℚ_p]}` を落として `= #μ_p(K)` にすると偽**。`p` を奇素数、`K = ℚ_p`
  とすると `μ_p(ℤ_p^×) = {1}` で右辺 `1`、左辺は
  `ℤ_p^× ≅ ℤ/(p−1) × ℤ_p` から `gcd(p, p−1) · p = p`。
* **`p` を `p ∤ n` の `n` に取り替えると偽**。その場合の正しい値は
  `gcd(n, q−1)`(`UnitsPowPrimeToP.lean::index_powRange_units`)であって
  `n^{[K:ℚ_p]} · #μ_n(K)` ではない。`n` が `p` とも `q−1` とも素なら
  左辺は `1`、後者は `n^{[K:ℚ_p]} > 1`。**`p` 乗と `p` 以外の乗は本質的に別の式**で、
  それが経路 C で `q` と `[K:ℚ_p]` を分離して読める理由でもある。
* **`n = 1` では自明**。両辺とも `1`。中身を持つのは `n = p`(と `p ∤ n` 側)だけ。
* **有限指数(`index_smallPrincipalUnitsO_ne_zero`)を落とすと橋が落ちる。**
  `Found/HerbrandIndex.lean` の主定理は `B.index ≠ 0` を落とすと**偽になる**
  (`A = Multiplicative ℤ`・`B = ⊥`・`n = 2` で `2 ≠ 1`)。ここで
  `smallPrincipalUnits` を「もっと小さい層」に取り替えても指数は有限のままなので
  結論は変わらないが、たとえば `B` を「`1` に収束する列が張る部分群」のような
  有限指数でないものに取り替えると、`p^d` の側だけが動いて等式は成り立たない。

## 逸脱の記録

`UnitsPowPrimeToP.lean` と同じく、原典 [pGC] が Serre の局所類体論(相互律)に
投げている論拠を経由しない(`ResearchPaper/pgc-goal.md` の「経路 C」で記録済み)。
本ファイルはさらに `Found/HerbrandIndex.lean`(mathlib にも ABC3 にも無かった
一般の群論)を新設して使う。消費する結論は「`[K:ℚ_p]` が群論的に決まる」だけなので
影響は無い。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found Subgroup
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## `smallPrincipalUnits` を `𝒪_K^×` の中に置き直す -/

/-- `𝒪_K^× → K^×`。 -/
noncomputable def unitsToCarrier (K : PAdicLocalField p) : (𝒪[K.carrier])ˣ →* (K.carrier)ˣ :=
  Units.map (algebraMap 𝒪[K.carrier] K.carrier : 𝒪[K.carrier] →+* K.carrier).toMonoidHom

@[simp] theorem unitsToCarrier_val (K : PAdicLocalField p) (u : (𝒪[K.carrier])ˣ) :
    ((unitsToCarrier K u : (K.carrier)ˣ) : K.carrier) = ((u : 𝒪[K.carrier]) : K.carrier) := rfl

theorem unitsToCarrier_injective (K : PAdicLocalField p) :
    Function.Injective (unitsToCarrier K) := by
  intro u v h
  exact Units.ext (Subtype.ext (congrArg (fun z : (K.carrier)ˣ => (z : K.carrier)) h))

/-- **主単数群(小さい層)を `𝒪_K^×` の部分群として見たもの**。 -/
noncomputable def smallPrincipalUnitsO (K : PAdicLocalField p) : Subgroup (𝒪[K.carrier])ˣ :=
  Subgroup.comap (unitsToCarrier K) (smallPrincipalUnits K)

/-- `‖w−1‖ ≤ 1/4` なら `‖w‖ = 1`、つまり `w` は `𝒪_K` の単数から来る。 -/
theorem exists_unitsToCarrier (K : PAdicLocalField p) (w : (K.carrier)ˣ)
    (hw : w ∈ smallPrincipalUnits K) : ∃ u : (𝒪[K.carrier])ˣ, unitsToCarrier K u = w := by
  have hw' : ‖(w : K.carrier) - 1‖ ≤ 1 / 4 := hw
  have hnorm : ‖(w : K.carrier)‖ = 1 :=
    norm_eq_one_of_norm_sub_one_lt K (lt_of_le_of_lt hw' (by norm_num))
  have hv : Valued.v (w : K.carrier) ≤ 1 := by
    have h : Valued.v (w : K.carrier) = (‖(w : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [h]
    exact_mod_cast le_of_eq hnorm
  have hy : IsUnit (⟨(w : K.carrier), hv⟩ : 𝒪[K.carrier]) :=
    Valued.integer.isUnit_iff_norm_eq_one.mpr hnorm
  refine ⟨hy.unit, Units.ext ?_⟩
  rw [unitsToCarrier_val, IsUnit.unit_spec]

/-- **`smallPrincipalUnitsO K ≃* smallPrincipalUnits K`**——`𝒪_K ↪ K` が単射で、
`‖u−1‖ ≤ 1/4` の単数は必ず `𝒪_K^×` から来るから。 -/
noncomputable def smallPrincipalUnitsOEquiv (K : PAdicLocalField p) :
    smallPrincipalUnitsO K ≃* smallPrincipalUnits K :=
  MulEquiv.ofBijective
    (MonoidHom.codRestrict ((unitsToCarrier K).comp (smallPrincipalUnitsO K).subtype)
      (smallPrincipalUnits K) (fun u => u.2))
    (by
      constructor
      · intro u v h
        exact Subtype.ext (unitsToCarrier_injective K (congrArg Subtype.val h))
      · rintro ⟨w, hw⟩
        obtain ⟨u, hu⟩ := exists_unitsToCarrier K w hw
        have hmem : u ∈ smallPrincipalUnitsO K := by
          show unitsToCarrier K u ∈ smallPrincipalUnits K
          rw [hu]
          exact hw
        exact ⟨⟨u, hmem⟩, Subtype.ext hu⟩)

/-! ## 有限指数——`1 + ϖ^m 𝒪_K` を挟む -/

/-- `𝒪_K` の元のノルムは `1` 以下。 -/
theorem norm_le_one_carrierIntegers (K : PAdicLocalField p) (c : 𝒪[K.carrier]) :
    ‖(c : K.carrier)‖ ≤ 1 := by
  have h : Valued.v ((c : K.carrier)) = (‖(c : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
  have h2 : Valued.v ((c : K.carrier)) ≤ 1 := c.2
  rw [h] at h2
  exact_mod_cast h2

/-- `‖ϖ‖^m ≤ 1/4` なら `1 + ϖ^m 𝒪_K ⊆ smallPrincipalUnitsO K`。 -/
theorem principalUnits_le_smallPrincipalUnitsO (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    {m : ℕ} (hm : ‖ϖ‖ ^ m ≤ 1 / 4) :
    principalUnits K ϖ m ≤ smallPrincipalUnitsO K := by
  intro v hv
  obtain ⟨c, hc⟩ := (mem_principalUnits_iff K ϖ m v).mp hv
  show ‖((unitsToCarrier K v : (K.carrier)ˣ) : K.carrier) - 1‖ ≤ 1 / 4
  rw [unitsToCarrier_val]
  have hcv : ((v : 𝒪[K.carrier]) : K.carrier) - 1
      = (c : K.carrier) * (ϖ : K.carrier) ^ m := by
    rw [hc]; push_cast; ring
  rw [hcv, norm_mul, norm_pow]
  calc ‖(c : K.carrier)‖ * ‖(ϖ : K.carrier)‖ ^ m
      ≤ 1 * ‖(ϖ : K.carrier)‖ ^ m :=
        mul_le_mul_of_nonneg_right (norm_le_one_carrierIntegers K c)
          (pow_nonneg (norm_nonneg _) m)
    _ ≤ 1 / 4 := by rw [one_mul]; exact hm

/-- **★★★★`smallPrincipalUnitsO K` は `𝒪_K^×` の有限指数部分群**。

`‖ϖ‖ < 1` なので `‖ϖ‖^n < 1/4` となる `n` があり、`m := n+1 ≥ 1` に対して
`principalUnits K ϖ m` が中に入る。後者は `𝒪_K/ϖ^m` が有限であることから
有限指数(`finite_principalUnitsQuotient`)なので、指数は約数の関係で伝わる。
**Herbrand の橋(`index_pow_mul_card_torsion`)の唯一の仮定がこれ**である。 -/
theorem index_smallPrincipalUnitsO_ne_zero (K : PAdicLocalField p) :
    (smallPrincipalUnitsO K).index ≠ 0 := by
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hlt : ‖ϖ‖ < 1 := Valued.integer.norm_irreducible_lt_one hϖ
  obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one (show (0:ℝ) < 1 / 4 by norm_num) hlt
  have hnn : (0:ℝ) ≤ ‖ϖ‖ ^ n := pow_nonneg (norm_nonneg _) n
  have hpow : ‖ϖ‖ ^ (n + 1) ≤ 1 / 4 := by
    have hrw : ‖ϖ‖ ^ (n + 1) = ‖ϖ‖ ^ n * ‖ϖ‖ := by ring
    rw [hrw]
    nlinarith [norm_nonneg ϖ]
  haveI := finite_principalUnitsQuotient K hϖ.maximalIdeal_eq (n + 1) (by omega)
  have hidx : (principalUnits K ϖ (n + 1)).index ≠ 0 := Subgroup.index_ne_zero_of_finite
  intro hzero
  exact hidx (Nat.eq_zero_of_zero_dvd
    (hzero ▸ Subgroup.index_dvd_of_le (principalUnits_le_smallPrincipalUnitsO K hpow)))

/-! ## `ℤ_p^d` 側の二つの量 -/

/-- `ℤ_p^d` に `p` 捩れは無い——`ℤ_p` は整域で `p ≠ 0`。 -/
theorem card_torsion_pi_padicInt (p d : ℕ) [Fact p.Prime] :
    Nat.card { x : Multiplicative (Fin d → ℤ_[p]) // x ^ p = 1 } = 1 := by
  have hp0 : (p : ℤ_[p]) ≠ 0 := Nat.cast_ne_zero.mpr (Fact.out : p.Prime).ne_zero
  have key : ∀ x : Multiplicative (Fin d → ℤ_[p]), x ^ p = 1 → x = 1 := by
    intro x hx
    have h : (p : ℕ) • (Multiplicative.toAdd x) = 0 := hx
    have h2 : Multiplicative.toAdd x = 0 := by
      funext i
      have hi : (p : ℕ) • (Multiplicative.toAdd x i) = 0 := congrFun h i
      rw [nsmul_eq_mul] at hi
      exact (mul_eq_zero.mp hi).resolve_left hp0
    exact Multiplicative.toAdd.injective h2
  haveI : Subsingleton { x : Multiplicative (Fin d → ℤ_[p]) // x ^ p = 1 } :=
    ⟨fun a b => Subtype.ext ((key a a.2).trans (key b b.2).symm)⟩
  exact Nat.card_eq_one_iff_unique.mpr ⟨inferInstance, ⟨⟨1, one_pow p⟩⟩⟩

/-- 主単数群(小さい層)に `p` 捩れは無い。 -/
theorem card_torsion_smallPrincipalUnitsO (K : PAdicLocalField p) :
    Nat.card { u : smallPrincipalUnitsO K // u ^ p = 1 } = 1 := by
  rw [card_torsion_of_mulEquiv (smallPrincipalUnitsOEquiv K) p,
    card_torsion_of_mulEquiv (smallPrincipalUnitsEquivPi K) p]
  exact card_torsion_pi_padicInt p _

/-- 主単数群(小さい層)の `p` 乗指数は `p^{[K:ℚ_p]}`
(`index_powRange_smallPrincipalUnits` の言い換え)。 -/
theorem index_powRange_smallPrincipalUnitsO (K : PAdicLocalField p) :
    ((powMonoidHom p : smallPrincipalUnitsO K →* smallPrincipalUnitsO K).range).index
      = p ^ Module.finrank ℚ_[p] K.carrier := by
  rw [index_powRange_of_mulEquiv (smallPrincipalUnitsOEquiv K) p]
  exact index_powRange_smallPrincipalUnits K

/-! ## 本体 -/

/-- **★★★★★★★★`[𝒪_K^× : (𝒪_K^×)^p] = p^{[K:ℚ_p]} · #μ_p(K)`**。

Herbrand 商の有限指数不変性(`Found/HerbrandIndex.lean::index_pow_mul_card_torsion`)を
`A = 𝒪_K^×`・`B = smallPrincipalUnitsO K` に当てる。`B ≅ ℤ_p^{[K:ℚ_p]}` 側では
`p` 捩れが自明・`p` 乗指数が `p^{[K:ℚ_p]}` なので、そのまま `A` 側へ移る。

`p ∤ n` の場合(`index_powRange_units`、`= gcd(n, q−1)`)と合わせて、
[pGC] 経路 C が `q` と `[K:ℚ_p]` を分けて読むための二本が揃った。 -/
theorem index_powRange_units_p (K : PAdicLocalField p) :
    ((powMonoidHom p : (𝒪[K.carrier])ˣ →* (𝒪[K.carrier])ˣ).range).index
      = p ^ Module.finrank ℚ_[p] K.carrier
        * Nat.card { u : (𝒪[K.carrier])ˣ // u ^ p = 1 } := by
  have h := index_pow_mul_card_torsion (smallPrincipalUnitsO K)
    (index_smallPrincipalUnitsO_ne_zero K) p
  rw [card_torsion_smallPrincipalUnitsO K, mul_one, index_powRange_smallPrincipalUnitsO K] at h
  exact h

end ABC3.Found.PGC
