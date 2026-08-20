import ABC3.Found.GaloisRep.WeilExist

/-!
# Galois (G5) 第 180 ブロック —— **★★★★★★★★`[n]` の引き戻しは単射である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★長らく仮定に置いていた `hμinj` が定理になった

`Skeleton` の `dvd_count_pullback` の `.needs` に

> ★★`μ` の単射性を仮定に置いた。第 118 の `pointHom_injective_of_transcendental`
> から出るが、`x([n]·)` の超越性を別途示す必要がある(3-8 ブロック)

と記録していた項目である。★**超越性を経由せずに済んだ。**

### ★★★★★★機構——第 171 の還元を使う

1. `μ` の核は素イデアル。0 でなければ第 138 の `prime_eq_xyIdeal` で
   `XYIdeal(R)` の形になる。
2. すると `XClass ∈ ker μ` から **`x_n = μ(genX)` が定数**になる(`y_n` も同様)。
3. ところが第 171 の `redPoint_mulByN` は、**どの素点 `v` でも**
   `redPoint_v(n·生成点) = n·Q_v` を与える。★定数点の還元はそれ自身(第 172)なので

       n·Q_v = R        (すべての `v` について)

4. `v` は `F` のすべての点を走る(第 172 の `placeOf`)。★★したがって
   任意の 2 点の差が `n` 等分点になり、`E(F)` が `E[n]`(有限)に単射する。
5. `F` は代数閉なので**各 `x` に対して点がある**(第 122 `exists_point_of_isAlgClosed`)。
   ★`F` は無限なので矛盾。

★★★**D2 のために積んだ第 171・172 が、別の未解決事項をそのまま解いた。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `yOf` / `somePt` | ★各 `x` に対する曲線上の点(選択) |
| `somePt_injective` | ★★`x` が違えば点も違う |
| `exists_ne_nsmul` | ★★★★★**`n` 倍写像は定数ではない** |
| `mulByN_injective` | ★★★★★★★★**`[n]` の引き戻しは単射** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [IsAlgClosed F] [DecidableEq F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-! ## ★各 `x` に対する点 -/

/-- ★与えられた `x` に対する `y` 座標(選択)。 -/
noncomputable def yOf (x : F) : F := Classical.choose (exists_point_of_isAlgClosed W x)

theorem equation_yOf (x : F) : W.Equation x (yOf W x) :=
  Classical.choose_spec (exists_point_of_isAlgClosed W x)

theorem nonsingular_yOf (x : F) : W.Nonsingular x (yOf W x) :=
  equation_iff_nonsingular.mp (equation_yOf W x)

/-- ★与えられた `x` に対する曲線上の点。 -/
noncomputable def somePt (x : F) : W.Point :=
  Point.some x (yOf W x) (nonsingular_yOf W x)

theorem somePt_eq (x : F) : somePt W x = Point.some x (yOf W x) (nonsingular_yOf W x) := rfl

theorem somePt_ne_zero (x : F) : somePt W x ≠ 0 := by rw [somePt_eq]; simp

theorem somePt_injective : Function.Injective (somePt W) := by
  intro a b hab
  rw [somePt_eq, somePt_eq] at hab
  injection hab

/-- ★★★★★**`n` 倍写像は定数ではない**。

★定数なら任意の 2 点の差が `n` 等分点になり、`E(F)` が有限集合に単射してしまう。 -/
theorem exists_ne_nsmul [Infinite F] (n : ℕ) (hn : 1 ≤ n) (hchar : (n : F) ≠ 0) :
    ∃ a b : F, n • somePt W a ≠ n • somePt W b := by
  by_contra hcon
  push_neg at hcon
  have hinj : Function.Injective (fun x : F => somePt W x - somePt W 0) :=
    fun a b hab => somePt_injective W (sub_left_injective hab)
  have hsub : ∀ x : F, (somePt W x - somePt W 0) ∈ {P : W.Point | n • P = 0} := by
    intro x
    show n • (somePt W x - somePt W 0) = 0
    rw [smul_sub, hcon x 0, sub_self]
  exact (Set.infinite_of_injective_forall_mem hinj hsub) (finite_torsion W n hn hchar)

/-! ## ★★★★★★★★単射性 -/

variable [Infinite F] [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★★★★**`[n]` の引き戻しは単射である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★核が 0 でなければ `x_n` が定数になり、第 171 の `redPoint_mulByN` から
**どの点でも `n` 倍が同じ点**になってしまう。★★`E(F)` は無限、`E[n]` は有限。 -/
theorem mulByN_injective (h2 : IsUnit (2 : F)) (n : ℕ) (hn : 1 ≤ n) (hchar : (n : F) ≠ 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns) :
    Function.Injective μ := by
  rw [RingHom.injective_iff_ker_eq_bot]
  by_contra hker
  haveI : (RingHom.ker μ).IsPrime := RingHom.ker_isPrime μ
  obtain ⟨c, y₀, heq, hkerEq⟩ := prime_eq_xyIdeal h2 W (RingHom.ker μ) hker
  have hxc : xn = algebraMap F W.FunctionField c := by
    have hmem : CoordinateRing.XClass W c ∈ RingHom.ker μ := by
      rw [hkerEq, CoordinateRing.XYIdeal]; exact Ideal.subset_span (Set.mem_insert _ _)
    rw [RingHom.mem_ker, XClass_eq, map_sub, hμx, hμF] at hmem
    exact sub_eq_zero.1 hmem
  have hyc : yn = algebraMap F W.FunctionField y₀ := by
    have hmem : CoordinateRing.YClass W (Polynomial.C y₀) ∈ RingHom.ker μ := by
      rw [hkerEq, CoordinateRing.XYIdeal]
      exact Ideal.subset_span (Set.mem_insert_of_mem _ rfl)
    rw [RingHom.mem_ker, YClass_eq, map_sub, hμy, hμF] at hmem
    exact sub_eq_zero.1 hmem
  have hcast : Point.some xn yn hns
      = toFF W (Point.some c y₀ (equation_iff_nonsingular.mp heq)) := by
    rw [toFF_some]
    exact point_some_congr hxc hyc hns _
  have hall : ∀ a : F, n • somePt W a = Point.some c y₀ (equation_iff_nonsingular.mp heq) := by
    intro a
    have hvid : (placeOf W (somePt W a) (somePt_ne_zero W a)).asIdeal
        = CoordinateRing.XYIdeal W a (Polynomial.C (yOf W a)) :=
      placeOf_some W (nonsingular_yOf W a) _
    have hred := redPoint_mulByN W (placeOf W (somePt W a) (somePt_ne_zero W a))
      (equation_yOf W a) hvid h2 n hns hμP
    rw [hcast, redPoint_toFF W _ (equation_yOf W a) hvid] at hred
    exact hred.symm
  obtain ⟨a, b, hab⟩ := exists_ne_nsmul W n hn hchar
  exact hab ((hall a).trans (hall b).symm)

/-! ## ★出典の紐付け(`.src`) -/

def mulByN_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——[n] の引き戻しが単射であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
