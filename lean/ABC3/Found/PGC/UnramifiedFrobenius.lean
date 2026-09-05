import ABC3.Found.PGC.UnramifiedRootsOfUnity

/-!
# 不分岐拡大の Frobenius は `p` と素な 1 の冪根に `q` 乗で作用する

`Gal(K^ur/K) ≅ Ẑ` の生成元(Frobenius)が、`p` と素な位数の 1 の冪根に
どう作用するかを決める。`q` を群論的に読む(LCFT 非経由)道の核心の算術。

## 筋

`σ` を `exists_frobenius`(`UnramifiedExtension.lean`)が与える Frobenius、
`ζ : adjoinIntegers K x` を `ζ^m = 1`(`p ∤ m`)とする。

1. `σ ζ` と `ζ^q` はどちらも `m` 乗すると `1`。
2. 剰余体では `residue (σ ζ) = residueAlgEquiv K x σ (residue ζ)
   = (residue ζ)^q = residue (ζ^q)`——これが Frobenius の定義そのもの。
3. したがって `u := (σ ζ)·(ζ^q)⁻¹` は `u^m = 1` かつ `u ≡ 1 (mod 𝔪)`。
4. 局所環で `m` が単元なら主単数群に `m` 捩れは無い
   (`eq_one_of_pow_eq_one_of_sub_one_mem_maximalIdeal`)ので `u = 1`。

`p ∤ m` から `(m : adjoinIntegers K x)` が単元であることは、
`‖(m : K)‖ = 1`(`norm_natCast_eq_one_of_not_dvd`)をスペクトルノルムで
`K(x)` へ延長して得る(`isUnit_natCast_adjoinIntegers`)。

## ★★★測定と、そこから引いた結論の訂正(2026-09-05)

以前この節には次のように書いていた——

> `adjoinIntegers` と `𝒪[(adjoinField K x).carrier]` の境界は越えられない。
> `(integersEquivAdjoinIntegers K x w).1 = w.1` は定義上ほとんど自明だが、
> 両辺の型が `↥K⟮x⟯` と `(adjoinField K x).carrier` で異なるため
> 既定の heartbeats で `(kernel) deterministic timeout`(212 秒)、
> `maxHeartbeats 1000000` で `(deterministic) timeout at whnf`(126 秒)。

**この測定値は本物**である(`𝒪[·]` が `Valuation.integer` +
スペクトルノルム由来の `Valued` インスタンスで、`whnf` が展開しきれない)。
しかしそこから引いた結論——「`exists_frobenius` を `𝒪[·]` 側で述べ直すしか
ない」——は**寄せる側が逆だった**。

`exists_frobenius` も `algEquivIntegers` も `residueAlgEquiv` も
`residueGalHom` も、**最初から `adjoinIntegers` 側で完結している**。
橋を渡る必要は無い。**`𝒪[(adjoinField K x).carrier]` を一度も書かなければ
壁は発生しない。** 以前ここにあった `integersEquivOf` / `residueEquivOf` は
`algEquivIntegers` / `residueAlgEquiv` を `𝒪[·]` 側で作り直したもので、
二重化はそこで生じていた。両者とも削除した。

教訓は `tools/lean-idioms.md` #69 にも反映してある——
**寄せる先は `adjoinIntegers` 側**(`exists_frobenius` がいる側)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 局所環の主単数群に「`m` が単元」の捩れは無い

付値もノルムも使わない汎用の補題。`μ_m(R) → μ_m(R/𝔪)` の単射性。 -/

/-- **局所環 `R` で `(m : R)` が単元なら、主単数の `m` 捩れは自明**。

`u = 1 + x`(`x ∈ 𝔪`)とすると `0 = u^m - 1 = x·S`、`S = ∑_{i<m} u^i`。
`S - m = ∑_{i<m}(u^i - 1) ∈ 𝔪` で `m` が単元だから `S ∉ 𝔪`、
局所環なので `S` は単元。よって `x = 0`。整域性は使わない。 -/
theorem eq_one_of_pow_eq_one_of_sub_one_mem_maximalIdeal
    {R : Type*} [CommRing R] [IsLocalRing R] {m : ℕ}
    (hm : IsUnit ((m : ℕ) : R)) {u : R} (hu : u ^ m = 1)
    (h1 : u - 1 ∈ IsLocalRing.maximalIdeal R) : u = 1 := by
  set S : R := ∑ i ∈ Finset.range m, u ^ i with hS
  have hfac : (u - 1) * S = 0 := by
    rw [hS, mul_comm, geom_sum_mul, hu, sub_self]
  have hui : ∀ i : ℕ, u ^ i - 1 ∈ IsLocalRing.maximalIdeal R := by
    intro i
    rw [← geom_sum_mul u i]
    exact Ideal.mul_mem_left _ _ h1
  have hSm : S - ((m : ℕ) : R) ∈ IsLocalRing.maximalIdeal R := by
    have hexp : S - ((m : ℕ) : R) = ∑ i ∈ Finset.range m, (u ^ i - 1) := by
      rw [Finset.sum_sub_distrib, hS]; simp
    rw [hexp]
    exact Ideal.sum_mem _ (fun i _ => hui i)
  have hSunit : IsUnit S := by
    rw [← IsLocalRing.residue_ne_zero_iff_isUnit]
    intro h0
    rw [IsLocalRing.residue_eq_zero_iff] at h0
    rw [← IsLocalRing.residue_ne_zero_iff_isUnit, ne_eq,
      IsLocalRing.residue_eq_zero_iff] at hm
    exact hm (by simpa using Ideal.sub_mem _ h0 hSm)
  obtain ⟨v, hv⟩ := hSunit
  rw [← hv] at hfac
  exact sub_eq_zero.mp ((Units.mul_left_eq_zero v).mp hfac)

/-! ## `p ∤ m` なら `m` は `adjoinIntegers K x` の単元 -/

/-- `algebraMap K → K(x)` はノルムを保つ——スペクトルノルムが基点上の
ノルムを延長すること(`spectralNorm_extends`)から。 -/
theorem norm_algebraMap_adjoin (K : PAdicLocalField p) (x : K.closure) (y : K.carrier) :
    ‖algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) y‖ = ‖y‖ := by
  rw [NormedAlgebra.norm_eq_spectralNorm K.carrier
    (algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) y)]
  exact spectralNorm_extends y

/-- ノルムがちょうど `1` の整数は単元——逆元のノルムも `1` なので
`adjoinIntegers K x` に入る。 -/
theorem isUnit_of_norm_eq_one (K : PAdicLocalField p) (x : K.closure) (y : adjoinIntegers K x)
    (h : ‖(y : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1) : IsUnit y := by
  set L := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hL
  have hne : (y : L) ≠ 0 := by
    intro h0; rw [h0, norm_zero] at h; exact zero_ne_one h
  have hmem : ‖(y : L)⁻¹‖ ≤ 1 := by rw [norm_inv, h]; simp
  refine ⟨⟨y, ⟨(y : L)⁻¹, hmem⟩, ?_, ?_⟩, rfl⟩
  · apply Subtype.ext; show (y : L) * (y : L)⁻¹ = 1; exact mul_inv_cancel₀ hne
  · apply Subtype.ext; show (y : L)⁻¹ * (y : L) = 1; exact inv_mul_cancel₀ hne

/-- **`p ∤ m` なら `(m : adjoinIntegers K x)` は単元**——
`‖(m : K)‖ = 1`(`norm_natCast_eq_one_of_not_dvd`)を
`algebraMap` で `K(x)` へ運ぶだけ。剰余体の標数を経由しない。 -/
theorem isUnit_natCast_adjoinIntegers (K : PAdicLocalField p) (x : K.closure) {m : ℕ}
    (hm : ¬ (p ∣ m)) : IsUnit ((m : ℕ) : adjoinIntegers K x) := by
  refine isUnit_of_norm_eq_one K x _ ?_
  have h1 : (((m : ℕ) : adjoinIntegers K x) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ((m : ℕ) : K.carrier) := by
    push_cast; ring
  rw [h1, norm_algebraMap_adjoin K x, norm_natCast_eq_one_of_not_dvd K hm]

/-! ## 本題——Frobenius は `p` と素な 1 の冪根に `q` 乗で作用する -/

/-- **★★★剰余体で `q` 乗として作用する `σ` は、`p` と素な位数の 1 の冪根に
そのまま `q` 乗として作用する**。

`Gal(K^ur/K) ≅ Ẑ` の生成元が `μ^{(p')}` に与える作用が、剰余体側の
Frobenius `z ↦ z^q` と**完全に一致する**ということ。持ち上げの一意性は
「主単数群に `p` と素な捩れが無い」ことから来る。 -/
theorem algEquivIntegers_eq_pow_of_pow_eq_one (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hσ : ∀ z : IsLocalRing.ResidueField (adjoinIntegers K x),
      residueAlgEquiv K x σ z = z ^ (Nat.card 𝓀[K.carrier]))
    {m : ℕ} (hm : ¬ (p ∣ m)) {ζ : adjoinIntegers K x} (hζ : ζ ^ m = 1) :
    algEquivIntegers K x σ ζ = ζ ^ (Nat.card 𝓀[K.carrier]) := by
  have hm0 : m ≠ 0 := by rintro rfl; exact hm (dvd_zero p)
  set q := Nat.card 𝓀[K.carrier] with hq
  set a : adjoinIntegers K x := algEquivIntegers K x σ ζ with ha
  set b : adjoinIntegers K x := ζ ^ q with hb
  have hbm : b ^ m = 1 := by rw [hb, ← pow_mul, mul_comm, pow_mul, hζ, one_pow]
  have ham : a ^ m = 1 := by rw [ha, ← map_pow, hζ, map_one]
  obtain ⟨B, hB⟩ := IsUnit.of_pow_eq_one hbm hm0
  have hBm : B ^ m = 1 := Units.ext (by rw [Units.val_pow_eq_pow_val, hB, hbm, Units.val_one])
  set u : adjoinIntegers K x := a * (B⁻¹ : (adjoinIntegers K x)ˣ) with hu
  have hum : u ^ m = 1 := by
    rw [hu, mul_pow, ham, one_mul, ← Units.val_pow_eq_pow_val, inv_pow, hBm, inv_one,
      Units.val_one]
  have hres : IsLocalRing.residue (adjoinIntegers K x) a
      = IsLocalRing.residue (adjoinIntegers K x) b := by
    rw [ha, ← residueAlgEquiv_apply K x σ ζ, hσ, hb, map_pow]
  have hu1 : u - 1 ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
    refine (IsLocalRing.residue_eq_zero_iff _).mp ?_
    rw [map_sub, map_one, hu, map_mul, hres, ← hB, ← map_mul, ← Units.val_mul, mul_inv_cancel,
      Units.val_one, map_one, sub_self]
  have h1 := eq_one_of_pow_eq_one_of_sub_one_mem_maximalIdeal
    (isUnit_natCast_adjoinIntegers K x hm) hum hu1
  have h2 : (a * (B⁻¹ : (adjoinIntegers K x)ˣ)) * (B : adjoinIntegers K x)
      = 1 * (B : adjoinIntegers K x) := by rw [← hu, h1]
  rw [mul_assoc, ← Units.val_mul, inv_mul_cancel, Units.val_one, mul_one, one_mul, hB] at h2
  exact h2

/-! ### 体の側での言い換え

`ζ^m = 1` なら自動的に `‖ζ‖ = 1` なので、`ζ` が整数環に入ることを
仮定に書く必要は無い。`μ_m` は体の中にいるので、こちらの形の方が使いやすい。 -/

/-- 1 の冪根のノルムは `1`——`‖ζ‖^m = 1` と `‖ζ‖ ≥ 0` から。

★第 3 の場合(`‖ζ‖ = -1`)を潰すのに `linarith` / `norm_num at h` を使うと
文脈全体(`↥K⟮x⟯` の巨大な項)を前処理して止まる。項で書く
(`tools/lean-idioms.md` #70)。 -/
theorem norm_eq_one_of_pow_eq_one (K : PAdicLocalField p) (x : K.closure)
    {m : ℕ} (hm0 : m ≠ 0)
    {ζ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)} (hζ : ζ ^ m = 1) :
    ‖ζ‖ = 1 := by
  have h1 : ‖ζ‖ ^ m = 1 := by rw [← norm_pow, hζ, norm_one]
  rcases (pow_eq_one_iff_cases.mp h1) with h | h | h
  · exact absurd h hm0
  · exact h
  · exact absurd (h.1 ▸ norm_nonneg ζ) (by norm_num)

/-- **★★★体の側での本題**——`σ ζ = ζ^q`(`ζ ∈ K(x)`、`ζ^m = 1`、`p ∤ m`)。 -/
theorem algEquiv_eq_pow_of_pow_eq_one (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hσ : ∀ z : IsLocalRing.ResidueField (adjoinIntegers K x),
      residueAlgEquiv K x σ z = z ^ (Nat.card 𝓀[K.carrier]))
    {m : ℕ} (hm : ¬ (p ∣ m))
    {ζ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)} (hζ : ζ ^ m = 1) :
    σ ζ = ζ ^ (Nat.card 𝓀[K.carrier]) := by
  have hm0 : m ≠ 0 := by rintro rfl; exact hm (dvd_zero p)
  have hmem : ζ ∈ adjoinIntegers K x := le_of_eq (norm_eq_one_of_pow_eq_one K x hm0 hζ)
  have hzm : (⟨ζ, hmem⟩ : adjoinIntegers K x) ^ m = 1 := by
    apply Subtype.ext
    show ζ ^ m = 1
    exact hζ
  have hmain := algEquivIntegers_eq_pow_of_pow_eq_one K x σ hσ hm hzm
  have hcoe : ((algEquivIntegers K x σ ⟨ζ, hmem⟩ : adjoinIntegers K x) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = σ ζ := rfl
  have hcoe2 : (((⟨ζ, hmem⟩ : adjoinIntegers K x) ^ (Nat.card 𝓀[K.carrier])
      : adjoinIntegers K x) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = ζ ^ (Nat.card 𝓀[K.carrier]) := rfl
  rw [← hcoe, ← hcoe2]
  exact congrArg Subtype.val hmain

/-- **★★★★不分岐拡大の Frobenius は `p` と素な 1 の冪根に `q` 乗で作用する**。

`exists_frobenius`(位数ちょうど `n` の生成元)に
`algEquivIntegers_eq_pow_of_pow_eq_one` / `algEquiv_eq_pow_of_pow_eq_one`
を当てただけ。 -/
theorem exists_frobenius_pow_of_pow_eq_one (K : PAdicLocalField p) (n : ℕ) (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ ∃ σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
          orderOf σ = n
            ∧ (∀ {m : ℕ}, ¬ (p ∣ m) → ∀ ζ : adjoinIntegers K x, ζ ^ m = 1 →
                algEquivIntegers K x σ ζ = ζ ^ (Nat.card 𝓀[K.carrier]))
            ∧ ∀ {m : ℕ}, ¬ (p ∣ m) →
                ∀ ζ : IntermediateField.adjoin K.carrier ({x} : Set K.closure), ζ ^ m = 1 →
                  σ ζ = ζ ^ (Nat.card 𝓀[K.carrier]) := by
  obtain ⟨x, hrank, hu, σ, hσ, hord⟩ := exists_frobenius K n hn
  exact ⟨x, hrank, hu, σ, hord,
    fun hm ζ hζ => algEquivIntegers_eq_pow_of_pow_eq_one K x σ hσ hm hζ,
    fun hm ζ hζ => algEquiv_eq_pow_of_pow_eq_one K x σ hσ hm hζ⟩

end ABC3.Found.PGC
