import ABC3.Found.GaloisRep.RemovableGlue

/-!
# Galois (G6) 第 248 ブロック —— **★★★★★★★★℘ の共線性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★葉 (c) 段 2 —— `℘` の加法定理(共線性の形)

第 246 の `collDet`、第 247 の貼り合わせ、第 241 の Liouville を組み合わせて

> `z₁ + z₂ + z₃ ∈ L` なら 3 点 `(℘(zᵢ), ℘'(zᵢ))` は一直線上にある

が出る(`weierstrass_collinear`)。これが**葉 (c) 段 2 の到達点**である。

    ℘(z₁)(℘'(z₂) − ℘'(z₃)) + ℘(z₂)(℘'(z₃) − ℘'(z₁)) + ℘(z₃)(℘'(z₁) − ℘'(z₂)) = 0

## ★★★組み立て

| 補題 | 役割 |
|---|---|
| `analyticAt_collDet` | 良い点(`z ∉ L` かつ `z + w ∉ L`)で解析的 |
| `exists_local_collDet` | ★★★★★★**全点で局所延長を持つ** |
| `collDet_eq_zero` | ★★★★★★★★`collDet = 0` |
| `weierstrass_collinear` | ★★★★★★★★対称な 3 点の形 |

★`exists_local_collDet` の第 2 種の悪い点(`z₀ + w ∈ L`)は、**`collDet_neg` で
第 1 種に帰着する**——`z₀ + w = l` なら `−z₀ − w = −l ∈ L` なので、
第 1 種の延長 `g` から `v ↦ −g(−v − w)` を作ればよい。★新しい極の解析は要らない。

## ★★共線性で十分である

葉 (c) が要求するのは「一意化が**準同型**」であり、それは

> `u₁ u₂ u₃ = 1` なら 3 点 `P(uᵢ)` は一直線上にある

という形である。★**`℘(z+w)` の明示式は要らない**——共線性そのものが目的の形であり、
段 3(形式的べき級数への移送)でもこの多項式の形のまま扱える。

## ★次(葉 (c) 段 3)

`X(u) = (2πi)^{-2}(℘(z) − …)` の対応で、上の共線性を Tate 級数の言葉に翻訳し、
第 221–240 の普遍環による降下を 3 変数 `(u₁, u₂, q)` に拡張する。
-/

namespace ABC3.Found.GaloisRep

open Complex Real Filter Topology PeriodPair

/-! ## ★良い点での解析性 -/

theorem analyticAt_collDet (L : PeriodPair) (w z₀ : ℂ)
    (h1 : z₀ ∉ L.lattice) (h2 : z₀ + w ∉ L.lattice) :
    AnalyticAt ℂ (collDet L w) z₀ := by
  have hP : AnalyticAt ℂ L.weierstrassP z₀ := L.analyticOnNhd_weierstrassP z₀ h1
  have hP' : AnalyticAt ℂ L.derivWeierstrassP z₀ := L.analyticOnNhd_derivWeierstrassP z₀ h1
  have hQ : AnalyticAt ℂ (fun z : ℂ => L.weierstrassP (z + w)) z₀ := analyticAt_shiftP L w z₀ h2
  have hQ' : AnalyticAt ℂ (fun z : ℂ => L.derivWeierstrassP (z + w)) z₀ := by
    have h3 := analyticAt_deriv_shiftP L w z₀ h2
    rwa [deriv_shiftP] at h3
  exact ((analyticAt_const.mul (hP.sub hQ)).sub (analyticAt_const.mul (hP'.add hQ'))).add
    ((hP.mul hQ').add (hP'.mul hQ))

/-! ## ★★★★★★全点で局所延長を持つ -/

/-- ★★★★★★**全点で局所延長を持つ**。

★第 2 種の悪い点(`z₀ + w ∈ L`)は `collDet_neg` で第 1 種に帰着する——
`z₀ + w = l` なら `−z₀ − w = −l ∈ L` なので、第 1 種の延長 `g` から
`v ↦ −g(−v − w)` を作ればよい。**新しい極の解析は要らない**。 -/
theorem exists_local_collDet (L : PeriodPair) (w : ℂ) (hw : w ∉ L.lattice) (z₀ : ℂ) :
    ∃ g : ℂ → ℂ, AnalyticAt ℂ g z₀ ∧ ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → collDet L w z = g z := by
  by_cases h1 : z₀ ∈ L.lattice
  · obtain ⟨g, hg, hgeq⟩ := exists_analytic_collDet L w hw ⟨z₀, h1⟩
    exact ⟨g, hg, Filter.Eventually.of_forall hgeq⟩
  by_cases h2 : z₀ + w ∈ L.lattice
  · have hl : -z₀ - w ∈ L.lattice := by
      have h3 : -(z₀ + w) ∈ L.lattice := neg_mem h2
      rwa [show -(z₀ + w) = -z₀ - w by ring] at h3
    obtain ⟨g, hg, hgeq⟩ := exists_analytic_collDet L w hw ⟨-z₀ - w, hl⟩
    refine ⟨fun v => -g (-v - w), ?_, ?_⟩
    · have hmap : AnalyticAt ℂ (fun v : ℂ => -v - w) z₀ :=
        (analyticAt_id.neg).sub analyticAt_const
      exact (AnalyticAt.comp (f := fun v : ℂ => -v - w) (x := z₀) hg hmap).neg
    · refine Filter.Eventually.of_forall ?_
      intro v hv
      have hvne : -v - w ≠ (⟨-z₀ - w, hl⟩ : L.lattice) := by
        intro h
        exact hv (by simpa using (by linear_combination -h : v = z₀))
      have h4 := hgeq (-v - w) hvne
      have hkey : collDet L w v = -collDet L w (-v - w) := by
        have h5 := collDet_neg L w (-v - w)
        rwa [show -(-v - w) - w = v by ring] at h5
      rw [hkey, h4]
  · exact ⟨collDet L w, analyticAt_collDet L w z₀ h1 h2,
      Filter.Eventually.of_forall (fun z _ => rfl)⟩

/-! ## ★★★★★★★★共線性 -/

/-- ★★★★★★★★**共線性**——`P(z)`, `P(w)`, `−P(z+w)` は一直線上にある。

第 247 の `eq_zero_of_periodic_antisymm` に、周期性(`collDet_add_lattice`)と
反対称性(`collDet_neg`)を流し込むだけ。 -/
theorem collDet_eq_zero (L : PeriodPair) (w : ℂ) (hw : w ∉ L.lattice) (z : ℂ)
    (h1 : z ∉ L.lattice) (h2 : z + w ∉ L.lattice) : collDet L w z = 0 := by
  refine eq_zero_of_periodic_antisymm L (collDet L w) w (exists_local_collDet L w hw) ?_ ?_ z
    (analyticAt_collDet L w z h1 h2)
  · intro l hl u
    exact collDet_add_lattice L w ⟨l, hl⟩ u
  · intro u
    exact collDet_neg L w u

/-- ★★★★★★★★**`℘` の共線性(対称な形)**——`z₁ + z₂ + z₃ ∈ L` なら 3 点
`(℘(zᵢ), ℘'(zᵢ))` は一直線上にある。

★これが**葉 (c) 段 2 の到達点**である。一意化が準同型であることは、
Tate 側でこの形に対応する多項式恒等式そのものである。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem weierstrass_collinear (L : PeriodPair) (z₁ z₂ z₃ : ℂ)
    (hs : z₁ + z₂ + z₃ ∈ L.lattice)
    (h₁ : z₁ ∉ L.lattice) (h₂ : z₂ ∉ L.lattice) (h₃ : z₃ ∉ L.lattice) :
    L.weierstrassP z₁ * (L.derivWeierstrassP z₂ - L.derivWeierstrassP z₃)
      + L.weierstrassP z₂ * (L.derivWeierstrassP z₃ - L.derivWeierstrassP z₁)
      + L.weierstrassP z₃ * (L.derivWeierstrassP z₁ - L.derivWeierstrassP z₂) = 0 := by
  have h12 : z₁ + z₂ ∉ L.lattice := by
    intro h
    apply h₃
    have h4 : (z₁ + z₂ + z₃) - (z₁ + z₂) ∈ L.lattice := sub_mem hs h
    rwa [show (z₁ + z₂ + z₃) - (z₁ + z₂) = z₃ by ring] at h4
  have hkey := collDet_eq_zero L z₂ h₂ z₁ h₁ h12
  have hQ : L.weierstrassP (z₁ + z₂) = L.weierstrassP z₃ := by
    have h5 := L.weierstrassP_add_coe (-z₃) ⟨z₁ + z₂ + z₃, hs⟩
    rw [show -z₃ + ((⟨z₁ + z₂ + z₃, hs⟩ : L.lattice) : ℂ) = z₁ + z₂ from by
      show -z₃ + (z₁ + z₂ + z₃) = z₁ + z₂; ring] at h5
    rw [h5, L.weierstrassP_neg]
  have hQ' : L.derivWeierstrassP (z₁ + z₂) = -L.derivWeierstrassP z₃ := by
    have h5 := L.derivWeierstrassP_add_coe (-z₃) ⟨z₁ + z₂ + z₃, hs⟩
    rw [show -z₃ + ((⟨z₁ + z₂ + z₃, hs⟩ : L.lattice) : ℂ) = z₁ + z₂ from by
      show -z₃ + (z₁ + z₂ + z₃) = z₁ + z₂; ring] at h5
    rw [h5, L.derivWeierstrassP_neg]
  simp only [collDet, hQ, hQ'] at hkey
  linear_combination hkey

/-! ## ★出典の紐付け(`.src`) -/

def collDet_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の行列式は 0)",
    sectionId := "genell-def-3-3" }

def weierstrass_collinear.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘ の共線性、対称な形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
