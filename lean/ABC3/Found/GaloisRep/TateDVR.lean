import ABC3.Found.GaloisRep.TateInjective

/-!
# Galois (G6) 第 260 ブロック —— **★★★★★★付値と `x` 座標の相異性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★付値と極大イデアルの冪の橋

第 259 までは `I` 進の言葉だけで進んだが、`(u, q/u)` の対を作る段になると
**「`v(u) ≤ v(q)` なら `q/u` が環に入る」**が要る。離散付値環では

    x ∈ 𝔪^n  ↔  n ≤ addVal x                    (`mem_maximalIdeal_pow_iff`)
    v(u) < v(q)  ⟹  ∃ w ∈ 𝔪, u·w = q            (`exists_complement`)

★mathlib の `IsDiscreteValuationRing.addVal_le_iff_dvd`(`v(a) ≤ v(b) ↔ a ∣ b`)と
`irreducible_iff_uniformizer`(`𝔪 = span {ϖ}`)を合わせるだけで出る。

## ★★★★★★`x` 座標が相異なる条件は `u ≠ v` だけ

第 257 の群法則は「3 点の `x` 座標が相異なること」を要求した。
第 259 の `tate_inj_X` を 3 点

    (u, vw),  (v, uw),  (w, uv)        (`u v w = q`)

に当てると、`X₁ = X₂` は

    (u = v かつ vw = uw)   または   (u = uw かつ vw = v)

を意味する。★★後者は `u(1 − w) = 0` すなわち `w = 1` を強いるが、
`w ∈ I` なので `1 − w` は単元、`w = 1` なら `0` が単元になって矛盾する。
★したがって **`X₁ = X₂ ⟺ u = v`**——非退化条件は `u, v, w` が相異なることだけになった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_maximalIdeal_pow_iff` | ★★★★★付値と `𝔪^n` の橋 |
| `exists_complement` | ★★★★★相方 `q/u` が環に入る |
| `tateXpair_ne_of_ne` | ★★★★★★**`x` 座標の相異性は `u ≠ v` から** |
-/

namespace ABC3.Found.GaloisRep

/-! ## ★★★★★付値と極大イデアルの冪 -/

section DVR

open IsDiscreteValuationRing IsLocalRing

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]

/-- ★★★★★**付値と極大イデアルの冪の橋**。 -/
theorem mem_maximalIdeal_pow_iff {ϖ : R} (hϖ : Irreducible ϖ) (n : ℕ) (x : R) :
    x ∈ (maximalIdeal R) ^ n ↔ (n : ℕ∞) ≤ addVal R x := by
  rw [(irreducible_iff_uniformizer ϖ).1 hϖ, Ideal.span_singleton_pow,
    Ideal.mem_span_singleton, ← addVal_le_iff_dvd, addVal_pow, addVal_uniformizer hϖ]
  simp

theorem mem_maximalIdeal_pow_iff' (n : ℕ) (x : R) :
    x ∈ (maximalIdeal R) ^ n ↔ (n : ℕ∞) ≤ addVal R x := by
  obtain ⟨ϖ, hϖ⟩ := exists_irreducible R
  exact mem_maximalIdeal_pow_iff hϖ n x

theorem mem_maximalIdeal_iff_one_le (x : R) : x ∈ maximalIdeal R ↔ 1 ≤ addVal R x := by
  have h := mem_maximalIdeal_pow_iff' 1 x
  rw [pow_one] at h
  simpa using h

theorem exists_mul_eq_of_addVal_le {a q : R} (h : addVal R a ≤ addVal R q) :
    ∃ w : R, a * w = q := by
  obtain ⟨c, hc⟩ := addVal_le_iff_dvd.1 h
  exact ⟨c, hc.symm⟩

/-- ★★★★★**相方の存在**——`v(u) < v(q)` なら `u·w = q` となる `w` が極大イデアルに在る。 -/
theorem exists_complement {u q : R} (hlt : addVal R u < addVal R q) :
    ∃ w : R, u * w = q ∧ w ∈ maximalIdeal R := by
  obtain ⟨w, hw⟩ := exists_mul_eq_of_addVal_le hlt.le
  refine ⟨w, hw, ?_⟩
  rw [mem_maximalIdeal_iff_one_le]
  rcases eq_or_ne (addVal R w) 0 with h0 | h0
  · exfalso
    have hmul : addVal R u + addVal R w = addVal R q := by
      rw [← hw]
      exact (AddValuation.map_mul (addVal R) u w).symm
    rw [h0, add_zero] at hmul
    exact absurd hmul (ne_of_lt hlt)
  · exact Order.one_le_iff_ne_zero.2 h0

end DVR

/-! ## ★★★★★★`x` 座標の相異性 -/

section Distinct

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★**3 点の `x` 座標が相異なる条件は `u ≠ v` だけ**。

★`u = u·w` の枝は `u(1 − w) = 0` を強いるが、`w ∈ I` なので `1 − w` は単元であり、
`w = 1` なら `0` が単元になって矛盾する。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXpair_ne_of_ne [IsAdicComplete I R] [IsDomain R] (u v w q : R) (hq : q ∈ I)
    (hu : u ∈ I) (hv : v ∈ I) (hw : w ∈ I) (huvw : u * v * w = q)
    (hu0 : u ≠ 0) (huv : u ≠ v) :
    tateXpair u (v * w) q hq ≠ tateXpair v (u * w) q hq := by
  intro heq
  have hvw : v * w ∈ I := Ideal.mul_mem_right _ _ hv
  have huw : u * w ∈ I := Ideal.mul_mem_right _ _ hu
  have h1 : u * (v * w) = q := by rw [← huvw]; ring
  have h2 : v * (u * w) = q := by rw [← huvw]; ring
  have hw1 : w ≠ 1 := by
    intro h
    have hun := isUnit_one_sub (I := I) hw
    rw [h, sub_self] at hun
    exact not_isUnit_zero hun
  rcases tate_inj_X u (v * w) v (u * w) q hq hu hvw hv huw h1 h2 heq with ⟨h, _⟩ | ⟨h, _⟩
  · exact huv h
  · have hz : u * (1 - w) = 0 := by linear_combination h
    rcases mul_eq_zero.1 hz with h' | h'
    · exact hu0 h'
    · exact hw1 (by linear_combination -h')

end Distinct

/-! ## ★出典の紐付け(`.src`) -/

def mem_maximalIdeal_pow_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——付値と極大イデアルの冪の橋)",
    sectionId := "genell-def-3-3" }

def tateXpair_ne_of_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——x 座標の相異性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
