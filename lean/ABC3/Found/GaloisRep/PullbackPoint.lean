import ABC3.Found.GaloisRep.RedKernel
import ABC3.Found.GaloisRep.Support

/-!
# Galois (G5) 第 164 ブロック —— **★★★★★★★引き戻した素点は `n·Q_v` である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★D2 の台の同定が完成した

第 151 で `count_v(μ f_P) ≠ 0 ⟺ I_P = P'`(`P'` は引き戻した素イデアル)まで来ていた。
本ブロックで **`P' = XYIdeal(n·Q_v)`** が出るので、台は

    { v : n·Q_v = P }   =   [n]⁻¹(P)

すなわち**ファイバーそのもの**になる。

### ★★★★★仕組みは 3 行である

1. 生成点の還元は `v` に対応する点そのもの(`redPoint_generic`)——
   `genX − c ∈ v.asIdeal` を `valuation_lt_one_iff_mem` で読むだけ。
2. 第 163 の `redPoint_nsmul` で `red(n·生成点) = n·Q_v`。
3. `n·生成点 = (μ(genX), μ(genY))` なので、`redConst(μ genX) = x(n·Q_v)`。
   ★これは `w(μ(genX − x(n·Q_v))) < 1`、すなわち **`XClass ∈ P'`** に他ならない。
   ★★`XYIdeal(n·Q_v) ⊆ P'` と極大性(第 138)から等式が出る。

★★★**剰余体を経由しない**——`redConst` が `F` に値を取ることが効いている。

## ★★★★`XYIdeal` は点を決める

`XYIdeal(x,y) = XYIdeal(x',y')` なら差が定数になり、極大イデアルに定数は 0 しか入らない。
★これで台を**点の等式**として書ける(`count_ne_zero_iff_nsmul`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `xyIdeal_inj` | ★★★★`XYIdeal` は点を決める |
| `redConst_coordX` / `redConst_coordY` | ★★★生成点の座標の還元 |
| `redPoint_generic` | ★★★★★**生成点の還元は `Q_v`** |
| `pullbackPrime_eq_xyIdeal` | ★★★★★★★**`P' = XYIdeal(n·Q_v)`** |
| `exists_pullbackPrime_eq_xyIdeal` | ★★★★★★★同上(存在形、`hnn` だけで済む) |
| `count_ne_zero_iff_nsmul` | ★★★★★★★**台は `[n]⁻¹(P)` そのもの** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-! ## ★★★★`XYIdeal` は点を決める -/

/-- ★★★★**`XYIdeal` は点を決める**。

★差 `x' − x` は定数で、極大イデアルに 0 でない定数は入らない。 -/
theorem xyIdeal_inj {x y x' y' : F} (h' : W.Equation x' y')
    (heq : CoordinateRing.XYIdeal W x (Polynomial.C y)
      = CoordinateRing.XYIdeal W x' (Polynomial.C y')) : x = x' ∧ y = y' := by
  have hne : CoordinateRing.XYIdeal W x' (Polynomial.C y') ≠ ⊤ :=
    (xyIdeal_isMaximal W h').ne_top
  have hconst : ∀ d : F, algebraMap F W.CoordinateRing d
      ∈ CoordinateRing.XYIdeal W x' (Polynomial.C y') → d = 0 := by
    intro d hd
    by_contra hd0
    exact hne (Ideal.eq_top_of_isUnit_mem _ hd
      ((IsUnit.mk0 d hd0).map (algebraMap F W.CoordinateRing)))
  constructor
  · have hx : CoordinateRing.XClass W x
        ∈ CoordinateRing.XYIdeal W x' (Polynomial.C y') := by
      rw [← heq, CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert _ _)
    have hx' : CoordinateRing.XClass W x'
        ∈ CoordinateRing.XYIdeal W x' (Polynomial.C y') := by
      rw [CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert _ _)
    have hsub : algebraMap F W.CoordinateRing (x' - x)
        = CoordinateRing.XClass W x - CoordinateRing.XClass W x' := by
      rw [XClass_eq, XClass_eq, map_sub]; ring
    exact (sub_eq_zero.1 (hconst _ (hsub ▸ Ideal.sub_mem _ hx hx'))).symm
  · have hy : CoordinateRing.YClass W (Polynomial.C y)
        ∈ CoordinateRing.XYIdeal W x' (Polynomial.C y') := by
      rw [← heq, CoordinateRing.XYIdeal]
      exact Ideal.subset_span (Set.mem_insert_of_mem _ rfl)
    have hy' : CoordinateRing.YClass W (Polynomial.C y')
        ∈ CoordinateRing.XYIdeal W x' (Polynomial.C y') := by
      rw [CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert_of_mem _ rfl)
    have hsub : algebraMap F W.CoordinateRing (y' - y)
        = CoordinateRing.YClass W (Polynomial.C y)
          - CoordinateRing.YClass W (Polynomial.C y') := by
      rw [YClass_eq, YClass_eq, map_sub]; ring
    exact (sub_eq_zero.1 (hconst _ (hsub ▸ Ideal.sub_mem _ hy hy'))).symm

variable [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

/-! ## ★★★生成点の還元 -/

/-- ★★★生成点の `x` 座標の還元は `c`。 -/
theorem redConst_coordX : redConst W v h hv (coordX W) = c := by
  refine redConst_eq W v h hv
    (show v.valuation W.FunctionField (coordX W) ≤ 1 from
      HeightOneSpectrum.valuation_le_one v _) ?_
  have hsub : coordX W - algebraMap F W.FunctionField c
      = algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.XClass W c) := by
    rw [XClass_eq, map_sub, coordX,
      IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField c]
  rw [hsub]
  exact (HeightOneSpectrum.valuation_lt_one_iff_mem v _).2
    (by rw [hv, CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert _ _))

/-- ★★★生成点の `y` 座標の還元は `y₀`。 -/
theorem redConst_coordY : redConst W v h hv (coordY W) = y₀ := by
  refine redConst_eq W v h hv
    (show v.valuation W.FunctionField (coordY W) ≤ 1 from
      HeightOneSpectrum.valuation_le_one v _) ?_
  have hsub : coordY W - algebraMap F W.FunctionField y₀
      = algebraMap W.CoordinateRing W.FunctionField
        (CoordinateRing.YClass W (Polynomial.C y₀)) := by
    rw [YClass_eq, map_sub, coordY,
      IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField y₀]
  rw [hsub]
  exact (HeightOneSpectrum.valuation_lt_one_iff_mem v _).2
    (by rw [hv, CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert_of_mem _ rfl))

variable [W.IsElliptic]

/-- ★★★★★**生成点の還元は `v` に対応する点そのもの**。 -/
theorem redPoint_generic :
    redPoint W v h hv (ABC3.Found.GaloisRep.genericPoint W)
      = Point.some c y₀ (equation_iff_nonsingular.mp h) := by
  rw [ABC3.Found.GaloisRep.genericPoint, redPoint_some W v h hv (nonsingular_coord W)
    (HeightOneSpectrum.valuation_le_one v _) (HeightOneSpectrum.valuation_le_one v _)]
  simp only [redConst_coordX W v h hv, redConst_coordY W v h hv]

variable [DecidableEq F]

/-! ## ★★★★★★★引き戻した素点の同定 -/

include hv in
/-- ★★★★★★★**引き戻した素イデアルは `n·Q_v` の素イデアルである**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`redConst(μ genX) = x(n·Q_v)` は `w(μ(XClass)) < 1`、つまり `XClass ∈ P'` そのもの。
★★包含 `XYIdeal(n·Q_v) ⊆ P'` と極大性(第 138)から等式になる。 -/
theorem pullbackPrime_eq_xyIdeal (h2 : IsUnit (2 : F)) (n : ℕ)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1)
    {xR yR : F} (hR : W.Nonsingular xR yR)
    (hQ : n • Point.some c y₀ (equation_iff_nonsingular.mp h) = Point.some xR yR hR) :
    pullbackPrime (v.valuation W.FunctionField) μ hnn
      = CoordinateRing.XYIdeal W xR (Polynomial.C yR) := by
  have hred : redPoint W v h hv (Point.some xn yn hns) = Point.some xR yR hR := by
    rw [← hμP, redPoint_nsmul W v h hv h2, redPoint_generic W v h hv, hQ]
  have hxle : v.valuation W.FunctionField xn ≤ 1 := by
    by_contra hc
    rw [not_le] at hc
    rw [(redPoint_eq_zero_iff W v h hv hns).2 hc] at hred
    exact absurd hred.symm (by simp)
  have hyle : v.valuation W.FunctionField yn ≤ 1 := val_y_le_of_val_x_le W v hns.1 hxle
  rw [redPoint_some W v h hv hns hxle hyle] at hred
  have hxR : redConst W v h hv xn = xR := by injection hred
  have hyR : redConst W v h hv yn = yR := by injection hred
  refine ((xyIdeal_isMaximal W hR.1).eq_of_le
    (pullbackPrime_isPrime (v.valuation W.FunctionField) μ hnn).ne_top ?_).symm
  rw [CoordinateRing.XYIdeal, Ideal.span_le]
  rintro z (rfl | rfl)
  · show v.valuation W.FunctionField (μ (CoordinateRing.XClass W xR)) < 1
    rw [XClass_eq, map_sub, hμx, hμF, ← hxR]
    exact redConst_spec W v h hv hxle
  · show v.valuation W.FunctionField
      (μ (CoordinateRing.YClass W (Polynomial.C yR))) < 1
    rw [YClass_eq, map_sub, hμy, hμF, ← hyR]
    exact redConst_spec W v h hv hyle

include hv in
/-- ★★★★★★★**引き戻した素イデアルは `n·Q_v` の素イデアルである**(存在形)。

★`hnn` があれば `w(μ genX) ≤ 1` なので `n·Q_v` は自動的に無限遠点ではない。 -/
theorem exists_pullbackPrime_eq_xyIdeal (h2 : IsUnit (2 : F)) (n : ℕ)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1) :
    ∃ (x' y' : F) (h' : W.Nonsingular x' y'),
      Point.some x' y' h' = n • Point.some c y₀ (equation_iff_nonsingular.mp h) ∧
      pullbackPrime (v.valuation W.FunctionField) μ hnn
        = CoordinateRing.XYIdeal W x' (Polynomial.C y') := by
  have hxle : v.valuation W.FunctionField xn ≤ 1 := by rw [← hμx]; exact hnn _
  have hyle : v.valuation W.FunctionField yn ≤ 1 := by rw [← hμy]; exact hnn _
  have hred : redPoint W v h hv (Point.some xn yn hns)
      = Point.some (redConst W v h hv xn) (redConst W v h hv yn)
        (equation_iff_nonsingular.mp
          (equation_redConst W v h hv hns.1 hxle hyle)) :=
    redPoint_some W v h hv hns hxle hyle
  have hQ : Point.some (redConst W v h hv xn) (redConst W v h hv yn)
      (equation_iff_nonsingular.mp (equation_redConst W v h hv hns.1 hxle hyle))
      = n • Point.some c y₀ (equation_iff_nonsingular.mp h) := by
    rw [← hred, ← hμP, redPoint_nsmul W v h hv h2, redPoint_generic W v h hv]
  exact ⟨_, _, _, hQ,
    pullbackPrime_eq_xyIdeal W v h hv h2 n μ hμF hns hμx hμy hμP hnn _ hQ.symm⟩

include hv in
/-- ★★★★★★★**因子の台は `[n]⁻¹(P)` そのものである**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 151(`count ≠ 0 ⟺ I_P = P'`)と本ブロック(`P' = XYIdeal(n·Q_v)`)、
そして `xyIdeal_inj` を繋ぐと、台が**点の等式**として書ける。 -/
theorem count_ne_zero_iff_nsmul (h2 : IsUnit (2 : F)) {n : ℕ} (hn : 1 ≤ n)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hnn : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1)
    {x y : F} (hP : W.Nonsingular x y) (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    {xR yR : F} (hR : W.Nonsingular xR yR)
    (hQ : n • Point.some c y₀ (equation_iff_nonsingular.mp h) = Point.some xR yR hR) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) ≠ 0
      ↔ Point.some x y hP = Point.some xR yR hR := by
  rw [count_ne_zero_iff W v μ hμinj hnn hP hn fP hfP,
    pullbackPrime_eq_xyIdeal W v h hv h2 n μ hμF hns hμx hμy hμP hnn hR hQ]
  constructor
  · intro heq
    obtain ⟨rfl, rfl⟩ := xyIdeal_inj W hR.1 heq
    rfl
  · intro heq
    have hx : x = xR := by injection heq
    have hy : y = yR := by injection heq
    rw [hx, hy]

/-! ## ★出典の紐付け(`.src`) -/

def pullbackPrime_eq_xyIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——引き戻した素点が n·Q_v であること)",
    sectionId := "genell-thm-3-8" }

def count_ne_zero_iff_nsmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の台が [n]⁻¹(P) であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
