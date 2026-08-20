import ABC3.Found.GaloisRep.RootUnit

/-!
# Galois (G5) 第 166 ブロック —— **★★★★★★因子の類は点の和である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★D2 が 1 つの等式に凝縮された

第 165 で `[J] = ∏ᶠ v, [v]^{e_v}` と `[v] = toClass(Q_v)` まで来た。
★本ブロックで有限台の `Finset` に落とし、`Point.toClass` が `AddMonoidHom` であることを
使って**点の和**にまとめる:

    [J] = toClass( Σ_{v ∈ s} e_v · Q_v )

★★`toClass_eq_zero`(mathlib)から

    J が単項  ⟸  Σ_{v ∈ s} e_v · Q_v = 0

となり、**D2 はこの 1 つの等式だけになった**。

## ★★★★★Abel–Jacobi も同じ機構で出る

`z ≠ 0` なら `(z)` は単項なので `[(z)] = 1`、したがって

    Σ_{v ∈ s} count_v(z) · Q_v = 0

★これは**アフィン曲線版の Abel–Jacobi**である。★★`(μ f_P)` に当てると
`count_v = n·e_v` から `n · (Σ e_v · Q_v) = 0`、すなわち
**`S := Σ e_v · Q_v` は `n` 等分点**であることまでは無料で出る。

★★★残るのは `S = 0`、すなわち `e_v` がファイバー上で一定であることだけ。

## ★★`Additive` の扱い

`Point.toClass` は `W.Point →+ Additive (ClassGroup F[W])` である。★`finprod` と
`finsum` の橋渡しは有限台の `Finset` に落としてから `toMul_sum`・`toMul_zsmul` で行う
(`finprod`/`finsum` のままだと `Additive` の変換補題が mathlib に無い)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `classGroup_rootUnit_prod` | ★★★有限台の積として書き直す |
| `classGroup_rootUnit_eq_toClass` | ★★★★★★**類は点の和の類** |
| `exists_spanSingleton_rootUnit` | ★★★★★★**和が 0 なら単項** |
| `finsum_count_smul_eq_zero` | ★★★★★★**Abel–Jacobi(アフィン版)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {R : Type} [CommRing R] [IsDedekindDomain R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★有限台の上での積として書き直す。 -/
theorem classGroup_rootUnit_prod (I : FractionalIdeal R⁰ K) (n : ℕ)
    (s : Finset (HeightOneSpectrum R))
    (hs : ∀ v : HeightOneSpectrum R, FractionalIdeal.count K v I ≠ 0 → v ∈ s) :
    ClassGroup.mk K (rootUnit I n)
      = ∏ v ∈ s, (ClassGroup.mk K (primeUnit K v)) ^ (FractionalIdeal.count K v I / (n : ℤ)) := by
  rw [classGroup_rootUnit]
  refine finprod_eq_prod_of_mulSupport_subset _ ?_
  intro v hv
  simp only [Function.mem_mulSupport] at hv
  refine Finset.mem_coe.2 (hs v ?_)
  intro h0
  exact hv (by rw [h0, Int.zero_ediv, zpow_zero])

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]
  (Q : HeightOneSpectrum W.CoordinateRing → W.Point)
  (hQ : ∀ v : HeightOneSpectrum W.CoordinateRing,
    ClassGroup.mk W.FunctionField (primeUnit W.FunctionField v)
      = Additive.toMul (Point.toClass (Q v)))

include hQ in
/-- ★★★★★★**因子の類は点の和の類である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Point.toClass` が `AddMonoidHom` であることと `toMul_sum`・`toMul_zsmul` だけ。 -/
theorem classGroup_rootUnit_eq_toClass
    (I : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) (n : ℕ)
    (s : Finset (HeightOneSpectrum W.CoordinateRing))
    (hs : ∀ v : HeightOneSpectrum W.CoordinateRing,
      FractionalIdeal.count W.FunctionField v I ≠ 0 → v ∈ s) :
    ClassGroup.mk W.FunctionField (rootUnit I n)
      = Additive.toMul (Point.toClass
        (∑ v ∈ s, (FractionalIdeal.count W.FunctionField v I / (n : ℤ)) • Q v)) := by
  rw [classGroup_rootUnit_prod I n s hs, map_sum, toMul_sum]
  refine Finset.prod_congr rfl fun v _ => ?_
  rw [hQ v, AddMonoidHom.map_zsmul, toMul_zsmul]

include hQ in
/-- ★★★★★★**点の和が 0 なら `rootUnit` は単項である**——D2 の帰着先。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `toClass_eq_zero` と第 140 の `isPrincipal_of_classGroup_eq_one` を繋ぐだけ。 -/
theorem exists_spanSingleton_rootUnit
    (I : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) (n : ℕ)
    (s : Finset (HeightOneSpectrum W.CoordinateRing))
    (hs : ∀ v : HeightOneSpectrum W.CoordinateRing,
      FractionalIdeal.count W.FunctionField v I ≠ 0 → v ∈ s)
    (hsum : (∑ v ∈ s, (FractionalIdeal.count W.FunctionField v I / (n : ℤ)) • Q v) = 0) :
    ∃ g : W.FunctionField, ((rootUnit I n : (FractionalIdeal W.CoordinateRing⁰ W.FunctionField)ˣ)
      : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g :=
  isPrincipal_of_classGroup_eq_one _ (by
    rw [classGroup_rootUnit_eq_toClass W Q hQ I n s hs, hsum, map_zero]; rfl)

include hQ in
/-- ★★★★★★**Abel–Jacobi(アフィン版)**——単項イデアルの因子は点の和として 0。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`(μ f_P)` に当てると `count_v = n·e_v` から `n · (Σ e_v · Q_v) = 0` が無料で出る。 -/
theorem finsum_count_smul_eq_zero {z : W.FunctionField} (hz : z ≠ 0)
    (s : Finset (HeightOneSpectrum W.CoordinateRing))
    (hs : ∀ v : HeightOneSpectrum W.CoordinateRing,
      FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z) ≠ 0 → v ∈ s) :
    (∑ v ∈ s, (FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z)) • Q v) = 0 := by
  set I := FractionalIdeal.spanSingleton W.CoordinateRing⁰ z with hI
  have hIne : I ≠ 0 := by
    rw [hI, ne_eq, FractionalIdeal.spanSingleton_eq_zero_iff]; exact hz
  have hcoe : ((rootUnit I 1 : (FractionalIdeal W.CoordinateRing⁰ W.FunctionField)ˣ)
      : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) = I := by
    have hp := rootUnit_pow (K := W.FunctionField) hIne (n := 1) (fun v => one_dvd _)
    rwa [pow_one] at hp
  have hone : ClassGroup.mk W.FunctionField (rootUnit I 1) = 1 := by
    rw [ClassGroup.mk_eq_one_iff, hcoe, hI]
    exact (FractionalIdeal.isPrincipal_iff _).2 ⟨z, rfl⟩
  rw [classGroup_rootUnit_eq_toClass W Q hQ I 1 s hs] at hone
  have h0 : Point.toClass (∑ v ∈ s,
      (FractionalIdeal.count W.FunctionField v I / ((1 : ℕ) : ℤ)) • Q v) = 0 := hone
  rw [Point.toClass_eq_zero] at h0
  simpa using h0

/-! ## ★出典の紐付け(`.src`) -/

def classGroup_rootUnit_eq_toClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の類が点の和の類であること)",
    sectionId := "genell-thm-3-8" }

def finsum_count_smul_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——Abel–Jacobi、アフィン版)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
