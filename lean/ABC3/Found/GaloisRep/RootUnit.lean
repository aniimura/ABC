import ABC3.Found.GaloisRep.PullbackPoint
import ABC3.Found.GaloisRep.DivisorTools

/-!
# Galois (G5) 第 165 ブロック —— **★★★★★★`n` 乗根イデアルを明示する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★存在では足りない——類を計算するには明示が要る

第 142 の `exists_pow_of_dvd_count` は `∃ J, J^n = I` しか返さない。
★D2(`J` が単項であること)を言うには **`J` の類を計算**しなければならず、
そのためには `J` を**明示的に**書き下す必要がある:

    J  :=  ∏ᶠ v,  v ^ (count_v(I) / n)

★★しかも `ClassGroup.mk` は**単元**にしか当たらないので、
分数イデアルの単元群 `(FractionalIdeal R⁰ K)ˣ` の中で組む(`rootUnit`)。
Dedekind 環では 0 でない分数イデアルは可逆(`CommGroupWithZero`)なので
`Units.mk0` で作れる。

## ★★★★★類は素点の類の積になる

`ClassGroup.mk` は `MonoidHom` なので `MonoidHom.map_finprod` がそのまま効く:

    [J] = ∏ᶠ v, [v] ^ (count_v(I) / n)

★そして **`[v] = toClass(Q_v)`**(`classGroup_primeUnit`)——
mathlib の `XYIdeal'` と `Point.toClass_some` に繋がる。

★★これで D2 は「`Σᶠ v, e_v · Q_v = 0`」に帰着する。

## ★★★残るのはファイバーの総和だけ

`e_v ≠ 0` となるのは `n·Q_v = P` のとき(第 164)。★`e_v` がファイバー上で一定
(平行移動 `τ_T` で不変)であれば、第 150 の `sum_fiber_eq` から総和は 0 になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `primeUnit` | ★★素点に対応する分数イデアルの単元 |
| `rootUnit` | ★★★★★**明示的な `n` 乗根イデアル** |
| `coe_rootUnit` | ★★分数イデアルとしての表示 |
| `rootUnit_pow` | ★★★★★**本当に `n` 乗根である** |
| `classGroup_rootUnit` | ★★★★★★**類は素点の類の積** |
| `primeUnit_eq_xyIdeal'` | ★★★★mathlib の `XYIdeal'` と一致する |
| `classGroup_primeUnit` | ★★★★★**`[v] = toClass(Q_v)`** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain nonZeroDivisors

/-! ## ★★素点に対応する単元 -/

variable {R : Type} [CommRing R] [IsDedekindDomain R]
  (K : Type) [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★素点に対応する分数イデアルの単元。

★Dedekind 環では 0 でない分数イデアルは可逆なので `Units.mk0` で作れる。 -/
noncomputable def primeUnit (v : HeightOneSpectrum R) : (FractionalIdeal R⁰ K)ˣ :=
  Units.mk0 ((v.asIdeal : Ideal R) : FractionalIdeal R⁰ K)
    (by rw [ne_eq, FractionalIdeal.coeIdeal_eq_zero]; exact v.ne_bot)

@[simp] theorem primeUnit_coe (v : HeightOneSpectrum R) :
    ((primeUnit K v : (FractionalIdeal R⁰ K)ˣ) : FractionalIdeal R⁰ K)
      = ((v.asIdeal : Ideal R) : FractionalIdeal R⁰ K) := rfl

variable {K}

/-! ## ★★★★★明示的な `n` 乗根 -/

/-- ★★★★★**明示的な `n` 乗根イデアル**——各素点の指数を `n` で割ったもの。

★第 142 の存在証明の中身をそのまま単元群の中で組み直したもの。 -/
noncomputable def rootUnit (I : FractionalIdeal R⁰ K) (n : ℕ) : (FractionalIdeal R⁰ K)ˣ :=
  ∏ᶠ v : HeightOneSpectrum R, (primeUnit K v) ^ (FractionalIdeal.count K v I / (n : ℤ))

theorem hasFiniteMulSupport_primeUnit (I : FractionalIdeal R⁰ K) (n : ℕ) :
    Function.HasFiniteMulSupport
      (fun v : HeightOneSpectrum R =>
        (primeUnit K v) ^ (FractionalIdeal.count K v I / (n : ℤ))) := by
  refine Set.Finite.subset (FractionalIdeal.finite_factors (K := K) I) ?_
  intro v hv
  simp only [Function.mem_mulSupport] at hv
  intro h0
  exact hv (by rw [h0, Int.zero_ediv, zpow_zero])

theorem hasFiniteMulSupport_coe (I : FractionalIdeal R⁰ K) (n : ℕ) :
    Function.HasFiniteMulSupport
      (fun v : HeightOneSpectrum R =>
        ((v.asIdeal : Ideal R) : FractionalIdeal R⁰ K)
          ^ (FractionalIdeal.count K v I / (n : ℤ))) := by
  refine Set.Finite.subset (FractionalIdeal.finite_factors (K := K) I) ?_
  intro v hv
  simp only [Function.mem_mulSupport] at hv
  intro h0
  exact hv (by rw [h0, Int.zero_ediv, zpow_zero])

/-- ★★分数イデアルとしての表示。 -/
theorem coe_rootUnit (I : FractionalIdeal R⁰ K) (n : ℕ) :
    ((rootUnit I n : (FractionalIdeal R⁰ K)ˣ) : FractionalIdeal R⁰ K)
      = ∏ᶠ v : HeightOneSpectrum R,
        (((v.asIdeal : Ideal R) : FractionalIdeal R⁰ K)
          ^ (FractionalIdeal.count K v I / (n : ℤ))) := by
  rw [rootUnit, ← Units.coeHom_apply,
    MonoidHom.map_finprod _ (hasFiniteMulSupport_primeUnit I n)]
  exact finprod_congr fun v => by
    rw [Units.coeHom_apply, Units.val_zpow_eq_zpow_val, primeUnit_coe]

/-- ★★★★★**`rootUnit` は本当に `n` 乗根である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `finprod_heightOneSpectrum_factorization'` に `Int.ediv_mul_cancel` を当てるだけ。 -/
theorem rootUnit_pow {I : FractionalIdeal R⁰ K} (hI : I ≠ 0) {n : ℕ}
    (h : ∀ v : HeightOneSpectrum R, (n : ℤ) ∣ FractionalIdeal.count K v I) :
    ((rootUnit I n : (FractionalIdeal R⁰ K)ˣ) : FractionalIdeal R⁰ K) ^ n = I := by
  rw [coe_rootUnit, finprod_pow (hasFiniteMulSupport_coe I n)]
  conv_rhs => rw [← FractionalIdeal.finprod_heightOneSpectrum_factorization' K hI]
  refine finprod_congr fun v => ?_
  rw [← zpow_natCast (((v.asIdeal : Ideal R) : FractionalIdeal R⁰ K)
    ^ (FractionalIdeal.count K v I / (n : ℤ))) n, ← zpow_mul]
  congr 1
  exact Int.ediv_mul_cancel (h v)

/-- ★★★★★★**類は素点の類の積になる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ClassGroup.mk` は `MonoidHom` なので `MonoidHom.map_finprod` がそのまま効く。 -/
theorem classGroup_rootUnit (I : FractionalIdeal R⁰ K) (n : ℕ) :
    ClassGroup.mk K (rootUnit I n)
      = ∏ᶠ v : HeightOneSpectrum R,
        (ClassGroup.mk K (primeUnit K v)) ^ (FractionalIdeal.count K v I / (n : ℤ)) := by
  rw [rootUnit, MonoidHom.map_finprod (ClassGroup.mk K) (hasFiniteMulSupport_primeUnit I n)]
  exact finprod_congr fun v => map_zpow (ClassGroup.mk K) _ _

/-! ## ★★★★★曲線の座標環での読み替え -/

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★**素点の単元は mathlib の `XYIdeal'` そのものである**。 -/
theorem primeUnit_eq_xyIdeal' (v : HeightOneSpectrum W.CoordinateRing) {c y₀ : F}
    (hns : W.Nonsingular c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀)) :
    primeUnit W.FunctionField v = CoordinateRing.XYIdeal' hns :=
  Units.ext (by rw [primeUnit_coe, CoordinateRing.XYIdeal'_eq, hv])

/-- ★★★★★**素点の類は対応する点の類である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これで因子の類が **`Point.toClass` の言葉**になる。 -/
theorem classGroup_primeUnit [DecidableEq F] (v : HeightOneSpectrum W.CoordinateRing) {c y₀ : F}
    (hns : W.Nonsingular c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀)) :
    ClassGroup.mk W.FunctionField (primeUnit W.FunctionField v)
      = Additive.toMul (Point.toClass (Point.some c y₀ hns)) := by
  rw [primeUnit_eq_xyIdeal' W v hns hv]
  exact (Point.toClass_some hns).symm

/-! ## ★出典の紐付け(`.src`) -/

def rootUnit_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——n 乗根イデアルを明示すること)",
    sectionId := "genell-thm-3-8" }

def classGroup_rootUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の類が素点の類の積になること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
