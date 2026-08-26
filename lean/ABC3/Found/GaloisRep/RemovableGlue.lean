import ABC3.Found.GaloisRep.CollinearDet

/-!
# Galois (G6) 第 247 ブロック —— **★★★★★★★局所延長の貼り合わせ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★穴あき近傍の極限で一様に貼る

第 246 で `collDet` は各格子点のまわりで解析的に延びると分かった。しかし
**その延長どうしを貼り合わせて整関数を作る**必要がある。悪い点は無限個(格子の 2 つの
平行移動)なので、`Function.update` で 1 点ずつ潰す形は使えない。

★★★そこで**場所によらない一様な定義**を採る:

    glue f z := limUnder (𝓝[≠] z) f

つまり「`z` を除いた近傍での極限」。これなら 1 本の式で全点を定義できる。

- **良い点**(`f` がそこで解析的)では `glue f z = f z`。
- **悪い点**(`f` が穴あき近傍で解析的な `g` に一致する)では `glue f z = g z`。

★どちらも同じ補題 `glue_eq` で出る——場合分けが要らないのがこの定義の利点である。

さらに `glue f` は `z₀` の**穴のない**近傍で `g` に一致する(`glue_eventuallyEq`)ので、
`AnalyticAt.congr` で `glue f` の解析性がそのまま出る。

## ★★★★★アフィン変換は `glue` を通り抜ける

`f (a z + b) = ε f(z)` なら `glue f (a z + b) = ε · glue f z`(`glue_affine`)。
`(a, b, ε) = (1, l, 1)` が**周期性**、`(-1, -c, -1)` が**反対称性**である。

★逆写像 `v ↦ (v - b)/a` で局所延長を運ぶだけなので、フィルターの押し出しを
書く必要がない。`Tendsto.eventually` 一発で済む。

## ★★★★★★★周期的かつ反対称なら 0

    eq_zero_of_periodic_antisymm : (局所延長) → (周期性) → (反対称性) → f z₀ = 0

★第 241 の Liouville で `glue f` が定数と分かり、反対称性から
`glue f 0 = -glue f 0`、よって **0**。**値を一点で計算する必要がない**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `glue` | ★★★★★★穴あき近傍の極限による一様な貼り合わせ |
| `glue_eventuallyEq` | ★★★★★★`glue f` は近傍で局所延長に一致する |
| `glue_eq`・`analyticAt_glue` | ★★点での一致と解析性 |
| `glue_affine` | ★★★★★アフィン変換は通り抜ける |
| `eq_zero_of_periodic_antisymm` | ★★★★★★★**周期的かつ反対称なら 0** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real Filter Topology PeriodPair

/-! ## ★★★★★★穴あき近傍の極限 -/

/-- ★★★★★★**局所延長を貼り合わせる**——各点で穴あき近傍の極限を取る。

場所によらない 1 本の式なので、悪い点が無限個あっても定義できる。 -/
noncomputable def glue (f : ℂ → ℂ) : ℂ → ℂ := fun z => limUnder (𝓝[≠] z) f

/-- ★★★★★★`f` が `z₀` の**穴あき**近傍で解析的な `g` に一致するなら、
`glue f` は `z₀` の**穴のない**近傍で `g` に一致する。

★穴あきの点 `z ≠ z₀` では `f` が `z` の近傍全体で `g` に一致するので `glue f z = g z`、
`z = z₀` では `g` の連続性から極限が `g z₀` になる。 -/
theorem glue_eventuallyEq (f g : ℂ → ℂ) (z₀ : ℂ) (hg : AnalyticAt ℂ g z₀)
    (heq : ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → f z = g z) :
    ∀ᶠ z in 𝓝 z₀, glue f z = g z := by
  have hV : {z : ℂ | (z ≠ z₀ → f z = g z) ∧ AnalyticAt ℂ g z} ∈ 𝓝 z₀ :=
    heq.and hg.eventually_analyticAt
  obtain ⟨V, hVsub, hVopen, hz₀V⟩ := mem_nhds_iff.1 hV
  refine Filter.eventually_of_mem (hVopen.mem_nhds hz₀V) ?_
  intro z hzV
  have hcg : ContinuousAt g z := ((hVsub hzV).2).continuousAt
  have hmem : V ∩ {z₀}ᶜ ∈ 𝓝[≠] z := by
    refine Filter.inter_mem (nhdsWithin_le_nhds (hVopen.mem_nhds hzV)) ?_
    by_cases hz : z = z₀
    · subst hz; exact self_mem_nhdsWithin
    · exact nhdsWithin_le_nhds (isOpen_compl_singleton.mem_nhds hz)
  have hfg : f =ᶠ[𝓝[≠] z] g := by
    filter_upwards [hmem] with u hu
    exact (hVsub hu.1).1 hu.2
  exact ((hcg.continuousWithinAt.tendsto).congr' hfg.symm).limUnder_eq

theorem glue_eq (f g : ℂ → ℂ) (z₀ : ℂ) (hg : AnalyticAt ℂ g z₀)
    (heq : ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → f z = g z) : glue f z₀ = g z₀ :=
  (glue_eventuallyEq f g z₀ hg heq).self_of_nhds

theorem analyticAt_glue (f g : ℂ → ℂ) (z₀ : ℂ) (hg : AnalyticAt ℂ g z₀)
    (heq : ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → f z = g z) : AnalyticAt ℂ (glue f) z₀ := by
  have h : (glue f) =ᶠ[𝓝 z₀] g := glue_eventuallyEq f g z₀ hg heq
  exact hg.congr h.symm

/-! ## ★★★★★アフィン変換は通り抜ける -/

/-- ★★★★★**アフィン変換は `glue` を通り抜ける**——`(a,b,ε) = (1,l,1)` が周期性、
`(-1,-c,-1)` が反対称性。

★逆写像 `v ↦ (v - b)/a` で局所延長を運ぶだけなので、フィルターの押し出しは要らない。 -/
theorem glue_affine (f : ℂ → ℂ) (a b ε : ℂ) (ha : a ≠ 0)
    (hf : ∀ z : ℂ, f (a * z + b) = ε * f z)
    (hloc : ∀ z₀ : ℂ, ∃ g, AnalyticAt ℂ g z₀ ∧ ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → f z = g z) (z : ℂ) :
    glue f (a * z + b) = ε * glue f z := by
  obtain ⟨g, hg, heq⟩ := hloc z
  have hinv : ∀ v : ℂ, a * ((v - b) / a) + b = v := by
    intro v; field_simp; ring
  have hz' : (a * z + b - b) / a = z := by field_simp; ring
  have hgz : AnalyticAt ℂ g ((a * z + b - b) / a) := by rw [hz']; exact hg
  have hmap : AnalyticAt ℂ (fun v : ℂ => (v - b) / a) (a * z + b) :=
    (analyticAt_id.sub analyticAt_const).div analyticAt_const ha
  have hg' : AnalyticAt ℂ (fun v : ℂ => ε * g ((v - b) / a)) (a * z + b) :=
    analyticAt_const.mul
      (AnalyticAt.comp (f := fun v : ℂ => (v - b) / a) (x := a * z + b) hgz hmap)
  have htend : Filter.Tendsto (fun v : ℂ => (v - b) / a) (𝓝 (a * z + b)) (𝓝 z) := by
    have h0 : Filter.Tendsto (fun v : ℂ => (v - b) / a) (𝓝 (a * z + b))
        (𝓝 ((a * z + b - b) / a)) :=
      ((continuous_id.sub continuous_const).div_const a).tendsto _
    rwa [hz'] at h0
  have heq' : ∀ᶠ v in 𝓝 (a * z + b), v ≠ a * z + b → f v = ε * g ((v - b) / a) := by
    filter_upwards [htend.eventually heq] with v hv hvne
    have hne : (v - b) / a ≠ z := by
      intro h
      exact hvne (by rw [← hinv v, h])
    have hfv : f v = ε * f ((v - b) / a) := by rw [← hf ((v - b) / a), hinv v]
    rw [hfv, hv hne]
  rw [glue_eq f (fun v => ε * g ((v - b) / a)) (a * z + b) hg' heq',
    glue_eq f g z hg heq, hz']

/-! ## ★★★★★★★周期的かつ反対称なら 0 -/

set_option maxHeartbeats 400000 in
/-- ★★★★★★★**周期的かつ反対称なら 0**。

局所延長を貼り合わせて整関数 `glue f` を作り、第 241 の Liouville で定数、
反対称性から `glue f 0 = -glue f 0`、よって **0**。
★値を一点で計算する必要がない。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem eq_zero_of_periodic_antisymm (L : PeriodPair) (f : ℂ → ℂ) (c : ℂ)
    (hloc : ∀ z₀ : ℂ, ∃ g, AnalyticAt ℂ g z₀ ∧ ∀ᶠ z in 𝓝 z₀, z ≠ z₀ → f z = g z)
    (hper : ∀ l ∈ L.lattice, ∀ z : ℂ, f (z + l) = f z)
    (hneg : ∀ z : ℂ, f (-z - c) = -f z)
    (z₀ : ℂ) (hz₀ : AnalyticAt ℂ f z₀) : f z₀ = 0 := by
  have hdiff : Differentiable ℂ (glue f) := by
    intro x
    obtain ⟨g, hg, heq⟩ := hloc x
    exact (analyticAt_glue f g x hg heq).differentiableAt
  have hper' : ∀ l ∈ L.lattice, ∀ z : ℂ, glue f (l + z) = glue f z := by
    intro l hl z
    have h1 : ∀ u : ℂ, f (1 * u + l) = 1 * f u := by
      intro u; rw [one_mul, one_mul]; exact hper l hl u
    have h2 := glue_affine f 1 l 1 one_ne_zero h1 hloc z
    rw [one_mul, one_mul] at h2
    rw [add_comm l z]
    exact h2
  have h2 : ∀ u : ℂ, f ((-1) * u + (-c)) = (-1) * f u := by
    intro u
    rw [show (-1 : ℂ) * u + -c = -u - c by ring, neg_one_mul]
    exact hneg u
  have hanti : ∀ z : ℂ, glue f (-z - c) = -glue f z := by
    intro z
    have h3 := glue_affine f (-1) (-c) (-1) (by norm_num) h2 hloc z
    rw [show (-1 : ℂ) * z + -c = -z - c by ring, neg_one_mul] at h3
    exact h3
  have hconst := eq_of_periodic_differentiable L (glue f) hdiff hper'
  have h0 : glue f 0 = -glue f 0 := by
    have h4 := hanti 0
    rw [show -(0 : ℂ) - c = -c by ring] at h4
    rw [← h4]
    exact hconst 0 (-c)
  have hz : glue f 0 = 0 := by
    have h5 : (2 : ℂ) * glue f 0 = 0 := by linear_combination h0
    simpa using h5
  have hgz : glue f z₀ = f z₀ :=
    glue_eq f f z₀ hz₀ (Filter.Eventually.of_forall (fun z _ => rfl))
  rw [← hgz, hconst z₀ 0, hz]

/-! ## ★出典の紐付け(`.src`) -/

def glue.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——局所延長の貼り合わせ)",
    sectionId := "genell-def-3-3" }

def eq_zero_of_periodic_antisymm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——周期的かつ反対称なら 0)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
