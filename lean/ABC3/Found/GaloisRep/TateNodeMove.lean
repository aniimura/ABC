import ABC3.Found.GaloisRep.TateSplitRoot

/-!
# Galois (G6) 第 306 ブロック —— **★★★★★★★★節点は `R` に持ち上がる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `c₄` が単元で `Δ ∈ 𝔪` なら、平行移動 `(r, t)` が在って
> **`a₃, a₄, a₆` がすべて `𝔪` に入る**(`exists_node_translation`)

★★★これは「還元の**節点を原点に移す**」ことに他ならない。
★G6 第 2 段(Tate 標準形)の骨格である。

## ★★★★★★節点の座標は `c₄` で割れば出る

    x₀ = (18b₆ − b₂b₄)/c₄                    (`nodeX` は分子)

★★これは**すべての標数で使える**——`c₄` が単元なのは乗法還元の仮定そのもの。
★★★mathlib の分裂条件に現れる 2 次式 `c₄T² + a₁c₄T − (54b₆ − 3b₂b₄ + a₂c₄)` は、
まさに `3x₀ + a₂` を `c₄` で書き直したものである。

## ★★★★★★★`y₀` は標数で 2 通り——しかし場合分けは 2 つで済む

`y₀` は `∂F/∂y = 0`(すなわち `2y + a₁x + a₃ = 0`)から出したいが、
**剰余標数 2 では `2` が割れない**。そこで `∂F/∂x = 0` を使う。

| 場合 | `y₀` の分子 | 分母 |
|---|---|---|
| `2` が単元 | `−(a₁x₀ + a₃c₄)`(`nodeYA`) | `2c₄` |
| `a₁` が単元 | `3x₀² + 2a₂c₄x₀ + a₄c₄²`(`nodeYB`) | `a₁c₄²` |

★★★★**この 2 つで尽きる**:`2` と `a₁` がともに `𝔪` に入れば
`c₄ = a₁·a₁³ + 2·(…)` も `𝔪` に入り、単元でなくなる(`isUnit_two_or_isUnit_a1`)。

## ★★★★★★分母を払えば `ring` で verify できる

証明の要は 3 本の多項式恒等式である(分母は払ってある):

    6x₀² + 4c₄a₂x₀ + 2c₄²a₄ + a₁²c₄x₀ + a₁a₃c₄² = −72Δ            (`node_id1`)
    c₄y₀² + 2c₄a₁x₀y₀ + 2c₄²a₃y₀ − 4x₀³ − ⋯ − 4c₄³a₆ = −4c₆Δ       (`node_id2`)
    2y₀ᴮ + a₁²c₄x₀ + a₁a₃c₄² = −72Δ                                (`node_id3`)

★★★どれも `b`・`c`・`Δ` の定義を開けば **`ring` が出す**。
★★★★★**分母を払ってから `ring` に渡す**——第 274 以来の型を、ここでも使った。

## ★★★剰余標数 2 の `a₆` は判別式から出る

`a₃, a₄ ∈ 𝔪` と `2 ∈ 𝔪` のとき、剰余体では `Δ = −a₁⁶a₆` になる。
★これも恒等式 `Δ + a₁⁶a₆ = 2P + a₃Q + a₄S`(`delta_add_a1_pow_six`)として
**`ring` で verify できる**形に書いた。★★`a₁` が単元なので `a₆ ∈ 𝔪` が出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `node_id1`・`node_id2`・`node_id3` | ★★★★★節点の 3 恒等式 |
| `delta_add_a1_pow_six` | ★★★★剰余標数 2 用の恒等式 |
| `isUnit_two_or_isUnit_a1` | ★★★★★`2` か `a₁` は単元 |
| `exists_node_translation` | ★★★★★★★★**節点を原点に移せる** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

section Node

variable {R : Type} [CommRing R]

/-! ## ★★★★★節点の座標(分母を払った形) -/

/-- ★★節点の `x` 座標の分子(分母は `c₄`)。 -/
def nodeX (E : WeierstrassCurve R) : R := 18 * E.b₆ - E.b₂ * E.b₄

/-- ★★節点の `y` 座標の分子(`2` が単元の場合、分母は `2c₄`)。 -/
def nodeYA (E : WeierstrassCurve R) : R := -(E.a₁ * nodeX E + E.a₃ * E.c₄)

/-- ★★節点の `y` 座標の分子(`a₁` が単元の場合、分母は `a₁c₄²`)。 -/
def nodeYB (E : WeierstrassCurve R) : R :=
  3 * (nodeX E) ^ 2 + 2 * E.a₂ * E.c₄ * nodeX E + E.a₄ * E.c₄ ^ 2

/-- ★★★★★**第 1 恒等式**——`∂F/∂x` の側。 -/
theorem node_id1 (E : WeierstrassCurve R) :
    6 * (nodeX E) ^ 2 + 4 * E.c₄ * E.a₂ * nodeX E + 2 * E.c₄ ^ 2 * E.a₄
      + E.a₁ ^ 2 * E.c₄ * nodeX E + E.a₁ * E.a₃ * E.c₄ ^ 2 = -72 * E.Δ := by
  simp only [nodeX, WeierstrassCurve.c₄, WeierstrassCurve.Δ, WeierstrassCurve.b₂,
    WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- ★★★★★**第 2 恒等式**——曲線の式の側。 -/
theorem node_id2 (E : WeierstrassCurve R) :
    E.c₄ * (nodeYA E) ^ 2 + 2 * E.c₄ * E.a₁ * nodeX E * nodeYA E
      + 2 * E.c₄ ^ 2 * E.a₃ * nodeYA E - 4 * (nodeX E) ^ 3
      - 4 * E.c₄ * E.a₂ * (nodeX E) ^ 2 - 4 * E.c₄ ^ 2 * E.a₄ * nodeX E
      - 4 * E.c₄ ^ 3 * E.a₆ = -4 * E.c₆ * E.Δ := by
  simp only [nodeYA, nodeX, WeierstrassCurve.c₄, WeierstrassCurve.c₆, WeierstrassCurve.Δ,
    WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- ★★★★★**第 3 恒等式**——`a₁` が単元の場合の `∂F/∂y` の側。 -/
theorem node_id3 (E : WeierstrassCurve R) :
    2 * nodeYB E + E.a₁ ^ 2 * E.c₄ * nodeX E + E.a₁ * E.a₃ * E.c₄ ^ 2 = -72 * E.Δ := by
  simp only [nodeYB, nodeX, WeierstrassCurve.c₄, WeierstrassCurve.Δ, WeierstrassCurve.b₂,
    WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

/-- ★★★★**剰余標数 2 用の恒等式**——`Δ + a₁⁶a₆` は `(2, a₃, a₄)` に入る。 -/
theorem delta_add_a1_pow_six (E : WeierstrassCurve R) :
    E.Δ + E.a₁ ^ 6 * E.a₆
      = 2 * (-6 * E.a₁ ^ 4 * E.a₂ * E.a₆ - 24 * E.a₁ ^ 2 * E.a₂ ^ 2 * E.a₆
          - 32 * E.a₂ ^ 3 * E.a₆ - 216 * E.a₆ ^ 2)
        + E.a₃ * (E.a₁ ^ 5 * E.a₄ - E.a₁ ^ 4 * E.a₂ * E.a₃ + 8 * E.a₁ ^ 3 * E.a₂ * E.a₄
          + E.a₁ ^ 3 * E.a₃ ^ 2 + 36 * E.a₁ ^ 3 * E.a₆ - 8 * E.a₁ ^ 2 * E.a₂ ^ 2 * E.a₃
          - 30 * E.a₁ ^ 2 * E.a₃ * E.a₄ + 16 * E.a₁ * E.a₂ ^ 2 * E.a₄
          + 36 * E.a₁ * E.a₂ * E.a₃ ^ 2 + 144 * E.a₁ * E.a₂ * E.a₆ - 96 * E.a₁ * E.a₄ ^ 2
          - 16 * E.a₂ ^ 3 * E.a₃ + 72 * E.a₂ * E.a₃ * E.a₄ - 27 * E.a₃ ^ 3
          - 216 * E.a₃ * E.a₆)
        + E.a₄ * (E.a₁ ^ 4 * E.a₄ + 8 * E.a₁ ^ 2 * E.a₂ * E.a₄ + 72 * E.a₁ ^ 2 * E.a₆
          + 16 * E.a₂ ^ 2 * E.a₄ + 288 * E.a₂ * E.a₆ - 64 * E.a₄ ^ 2) := by
  simp only [WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈]
  ring

/-! ## ★平行移動の係数 -/

theorem trans_a₃ (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₃ = E.a₃ + r * E.a₁ + 2 * t := by
  rw [WeierstrassCurve.variableChange_a₃]
  simp only [inv_one, Units.val_one, one_pow, one_mul]

theorem trans_a₄ (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₄ = E.a₄ + 2 * r * E.a₂ - t * E.a₁ + 3 * r ^ 2 := by
  rw [WeierstrassCurve.variableChange_a₄]
  simp only [inv_one, Units.val_one, one_pow, one_mul]
  ring

theorem trans_a₆ (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₆
      = E.a₆ + r * E.a₄ + r ^ 2 * E.a₂ + r ^ 3 - t * E.a₃ - t ^ 2 - r * t * E.a₁ := by
  rw [WeierstrassCurve.variableChange_a₆]
  simp only [inv_one, Units.val_one, one_pow, one_mul]

theorem trans_a₁ (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₁ = E.a₁ := by
  rw [WeierstrassCurve.variableChange_a₁]
  simp only [inv_one, Units.val_one, one_mul]
  ring

theorem trans_Delta (E : WeierstrassCurve R) (r t : R) :
    ((⟨1, r, 0, t⟩ : VariableChange R) • E).Δ = E.Δ := by
  rw [WeierstrassCurve.variableChange_Δ]
  simp only [inv_one, Units.val_one, one_pow, one_mul]

end Node

/-! ## ★★★★★`2` か `a₁` は単元 -/

section Local

variable {R : Type} [CommRing R] [IsLocalRing R]

/-- ★★★★★**`c₄` が単元なら `2` か `a₁` が単元**。 -/
theorem isUnit_two_or_isUnit_a1 (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄) :
    IsUnit (2 : R) ∨ IsUnit E.a₁ := by
  by_cases h2 : IsUnit (2 : R)
  · exact Or.inl h2
  right
  by_contra h1
  have hm2 : (2 : R) ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal _).2 h2
  have hm1 : E.a₁ ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal _).2 h1
  have heq : E.c₄ = E.a₁ * E.a₁ ^ 3
      + 2 * (4 * E.a₁ ^ 2 * E.a₂ + 8 * E.a₂ ^ 2 - 24 * E.a₄ - 12 * E.a₁ * E.a₃) := by
    simp only [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
    ring
  have : E.c₄ ∈ IsLocalRing.maximalIdeal R := by
    rw [heq]
    exact Ideal.add_mem _ (Ideal.mul_mem_right _ _ hm1) (Ideal.mul_mem_right _ _ hm2)
  exact (IsLocalRing.notMem_maximalIdeal.2 hc4) this

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**節点を原点に移せる**——`a₃, a₄, a₆` がすべて `𝔪` に入る平行移動がある。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_node_translation (E : WeierstrassCurve R) (hc4 : IsUnit E.c₄)
    (hΔ : E.Δ ∈ IsLocalRing.maximalIdeal R) :
    ∃ r t : R, ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₃ ∈ IsLocalRing.maximalIdeal R
      ∧ ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₄ ∈ IsLocalRing.maximalIdeal R
      ∧ ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  set cinv : R := ((hc4.unit⁻¹ : Rˣ) : R) with hcinv
  have hc : cinv * E.c₄ = 1 := by
    have h2 : ((hc4.unit⁻¹ : Rˣ) : R) * ((hc4.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [hc4.unit_spec] at h2
    exact h2
  set r : R := nodeX E * cinv with hrdef
  have hr : E.c₄ * r = nodeX E := by
    rw [hrdef]
    linear_combination (nodeX E) * hc
  by_cases h2 : IsUnit (2 : R)
  · -- ★`2` が単元の場合
    set ti : R := ((h2.unit⁻¹ : Rˣ) : R) with hti
    have hi : ti * 2 = 1 := by
      have h3 : ((h2.unit⁻¹ : Rˣ) : R) * ((h2.unit : Rˣ) : R) = 1 := by
        rw [← Units.val_mul, inv_mul_cancel]
        rfl
      rw [h2.unit_spec] at h3
      exact h3
    set t : R := nodeYA E * cinv * ti with htdef
    have ht : 2 * E.c₄ * t = nodeYA E := by
      rw [htdef]
      linear_combination (nodeYA E * 2 * ti) * hc + (nodeYA E) * hi
    refine ⟨r, t, ?_, ?_, ?_⟩
    · rw [trans_a₃]
      refine (Ideal.unit_mul_mem_iff_mem _ hc4).1 ?_
      have hzero : E.c₄ * (E.a₃ + r * E.a₁ + 2 * t) = 0 := by
        simp only [nodeYA] at ht
        linear_combination E.a₁ * hr + ht
      rw [hzero]
      exact Ideal.zero_mem _
    · rw [trans_a₄]
      refine (Ideal.unit_mul_mem_iff_mem _ (h2.mul (hc4.pow 2))).1 ?_
      have hval : 2 * E.c₄ ^ 2 * (E.a₄ + 2 * r * E.a₂ - t * E.a₁ + 3 * r ^ 2) = -72 * E.Δ := by
        simp only [nodeYA] at ht
        linear_combination node_id1 E
          + (4 * E.c₄ * E.a₂ + 6 * (E.c₄ * r + nodeX E)) * hr - E.a₁ * E.c₄ * ht
      rw [hval]
      exact Ideal.mul_mem_left _ _ hΔ
    · rw [trans_a₆]
      refine (Ideal.unit_mul_mem_iff_mem _ ((h2.pow 2).mul (hc4.pow 3))).1 ?_
      have hval : 2 ^ 2 * E.c₄ ^ 3
          * (E.a₆ + r * E.a₄ + r ^ 2 * E.a₂ + r ^ 3 - t * E.a₃ - t ^ 2 - r * t * E.a₁)
          = 4 * E.c₆ * E.Δ := by
        linear_combination -node_id2 E
          + (4 * E.c₄ ^ 2 * E.a₄ + 4 * E.c₄ * E.a₂ * (E.c₄ * r + nodeX E)
            + 4 * ((E.c₄ * r) ^ 2 + (E.c₄ * r) * nodeX E + (nodeX E) ^ 2)
            - 2 * E.a₁ * E.c₄ * (2 * E.c₄ * t)) * hr
          + (-2 * E.c₄ ^ 2 * E.a₃ - E.c₄ * (2 * E.c₄ * t + nodeYA E)
            - 2 * E.a₁ * E.c₄ * nodeX E) * ht
      rw [hval]
      exact Ideal.mul_mem_left _ _ hΔ
  · -- ★`a₁` が単元の場合(剰余標数 2)
    have hm2 : (2 : R) ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal _).2 h2
    have ha1 : IsUnit E.a₁ := (isUnit_two_or_isUnit_a1 E hc4).resolve_left h2
    set ainv : R := ((ha1.unit⁻¹ : Rˣ) : R) with hainv
    have ha : ainv * E.a₁ = 1 := by
      have h3 : ((ha1.unit⁻¹ : Rˣ) : R) * ((ha1.unit : Rˣ) : R) = 1 := by
        rw [← Units.val_mul, inv_mul_cancel]
        rfl
      rw [ha1.unit_spec] at h3
      exact h3
    set t : R := nodeYB E * cinv ^ 2 * ainv with htdef
    have ht : E.a₁ * E.c₄ ^ 2 * t = nodeYB E := by
      rw [htdef]
      linear_combination (nodeYB E * ainv * E.a₁ * (cinv * E.c₄ + 1)) * hc + (nodeYB E) * ha
    have h3F : ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₃ ∈ IsLocalRing.maximalIdeal R := by
      rw [trans_a₃]
      refine (Ideal.unit_mul_mem_iff_mem _ (ha1.mul (hc4.pow 2))).1 ?_
      have hval : E.a₁ * E.c₄ ^ 2 * (E.a₃ + r * E.a₁ + 2 * t) = -72 * E.Δ := by
        linear_combination node_id3 E + (E.a₁ ^ 2 * E.c₄) * hr + 2 * ht
      rw [hval]
      exact Ideal.mul_mem_left _ _ hΔ
    have h4F : ((⟨1, r, 0, t⟩ : VariableChange R) • E).a₄ ∈ IsLocalRing.maximalIdeal R := by
      rw [trans_a₄]
      refine (Ideal.unit_mul_mem_iff_mem _ (hc4.pow 2)).1 ?_
      have hzero : E.c₄ ^ 2 * (E.a₄ + 2 * r * E.a₂ - t * E.a₁ + 3 * r ^ 2) = 0 := by
        simp only [nodeYB] at ht
        linear_combination (2 * E.a₂ * E.c₄ + 3 * (E.c₄ * r + nodeX E)) * hr - ht
      rw [hzero]
      exact Ideal.zero_mem _
    refine ⟨r, t, h3F, h4F, ?_⟩
    set F := (⟨1, r, 0, t⟩ : VariableChange R) • E with hF
    have ha1F : F.a₁ = E.a₁ := trans_a₁ E _ _
    have hΔF : F.Δ = E.Δ := trans_Delta E _ _
    have hkey := delta_add_a1_pow_six F
    have hmem : F.a₁ ^ 6 * F.a₆ ∈ IsLocalRing.maximalIdeal R := by
      have hrw : F.a₁ ^ 6 * F.a₆ = (F.Δ + F.a₁ ^ 6 * F.a₆) - F.Δ := by ring
      rw [hrw, hkey, hΔF]
      refine Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _ ?_ ?_) ?_) hΔ
      · exact Ideal.mul_mem_right _ _ hm2
      · exact Ideal.mul_mem_right _ _ h3F
      · exact Ideal.mul_mem_right _ _ h4F
    rw [ha1F] at hmem
    exact (Ideal.unit_mul_mem_iff_mem _ (ha1.pow 6)).1 hmem

end Local

/-! ## ★出典の紐付け(`.src`) -/

def node_id1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(節点の座標——第 1 恒等式)",
    sectionId := "genell-def-3-3" }

def exists_node_translation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(節点を原点に移す平行移動)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
