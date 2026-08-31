/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightIdealNorm
import ABC3.Meta.Claim

/-!
# ★★★★★★★★整数表現に直しても素朴高さは変わらない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは何か

`§9-853`（段 C2c-2）は**整数座標**についての公式

    `H(x) · N(span{x_i}) = ∏_{v 無限素点} (sup_i v(x_i))^{mult v}`

を与えた。★しかし消費側（`projPointCoord`、`§9-C2b`）が出すのは
**正規化座標 `x_j/x_i ∈ F`** であり、整数とは限らない。

★★本ファイルはその橋である:

> 任意の（全部 0 でない）組 `y : ι → K` に対し、**整数表現 `x : ι → 𝓞_K` が取れて**
> 同じ公式が `H(y)` について成り立つ。

★★★機構は 2 本だけ:

| 道具 | 役割 |
|---|---|
| `Height.mulHeight_smul_eq_mulHeight`（mathlib） | ★**素朴高さはスカラー倍で変わらない**（射影不変） |
| `IsLocalization.exist_integer_multiples`（mathlib） | ★★共通分母が取れる |

## ★測定の記録

★どちらも mathlib にあった。**段 C2c-2 の周りは在庫が厚い**
——`§9-853` でも `absNorm_mul_finprod_finitePlace_eq_one` がそのまま効いた。

## ★残っている段（明示）

★★段 C2c-1（超平面因子の引き戻しが座標の生成するイデアルであること）と
段 C2c-3（与えられた計量と Fubini–Study 型の差が一様に有界であること）が残る。
★C2c-3 は Prop 1.4 (iii)（`§9-806`）が本体を持っているので、
実質は **C2c-1（幾何の側）** だけである。
-/

namespace ABC3.Found.GenEll

open NumberField Ideal Height

/-! ## ★★★★★スカラー倍で整数に直す -/

/-- ★★★★★**整数表現に直しても公式は `H(y)` で成り立つ**。

★機構は `Height.mulHeight_smul_eq_mulHeight`（素朴高さは射影不変）だけである。 -/
theorem mulHeight_mul_absNorm_of_scaled (K : Type) [Field K] [NumberField K]
    {ι : Type} [Finite ι] (y : ι → K) (c : K) (hc : c ≠ 0) (x : ι → 𝓞 K)
    (hcx : ∀ i, (x i : K) = c * y i) (hx : x ≠ 0) :
    Height.mulHeight y * ((Ideal.span (Set.range x)).absNorm : ℝ)
      = ∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult := by
  have h1 : (fun i => (x i : K)) = c • y := by
    funext i; rw [hcx i]; rfl
  have h2 := mulHeight_mul_absNorm K x hx
  rw [h1, Height.mulHeight_smul_eq_mulHeight y hc] at h2
  exact h2

/-! ## ★★★★★★★★整数表現の存在 -/

/-- ★★★★★★★★**どんな組にも整数表現があり、そこで公式が成り立つ**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

`y : ι → K` が全部 0 でなければ、★**整数表現 `x : ι → 𝓞_K`** が取れて

    `H(y) · N(span{x_i}) = ∏_{v 無限素点} (sup_i v(x_i))^{mult v}`

★★これで `projPointCoord`（正規化座標、`§9-C2b`）が出す組にも `§9-853` が使える。
★★★機構は `IsLocalization.exist_integer_multiples`（共通分母）＋
`Height.mulHeight_smul_eq_mulHeight`（射影不変）だけである。 -/
theorem exists_integral_repr (K : Type) [Field K] [NumberField K]
    {ι : Type} [Fintype ι] (y : ι → K) (hy : y ≠ 0) :
    ∃ x : ι → 𝓞 K, x ≠ 0 ∧
      Height.mulHeight y * ((Ideal.span (Set.range x)).absNorm : ℝ)
        = ∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult := by
  obtain ⟨b, hb⟩ := IsLocalization.exist_integer_multiples
    (nonZeroDivisors (𝓞 K)) (Finset.univ : Finset ι) y
  have hc0 : ((b.1 : 𝓞 K) : K) ≠ 0 := by
    intro h
    exact nonZeroDivisors.coe_ne_zero b (by exact_mod_cast h)
  have hval : ∀ i, ((((hb i (Finset.mem_univ i)).choose : 𝓞 K)) : K)
      = ((b.1 : 𝓞 K) : K) * y i := by
    intro i
    have := (hb i (Finset.mem_univ i)).choose_spec
    simpa [Algebra.smul_def] using this
  have hxne : (fun i => (hb i (Finset.mem_univ i)).choose) ≠ (0 : ι → 𝓞 K) := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hy
    refine Function.ne_iff.mpr ⟨i, ?_⟩
    intro h0
    rw [show ((0 : ι → 𝓞 K) i) = 0 from rfl] at h0
    have h := hval i
    rw [h0] at h
    simp only [map_zero] at h
    exact (mul_ne_zero hc0 (by simpa using hi)) h.symm
  exact ⟨_, hxne,
    mulHeight_mul_absNorm_of_scaled K y ((b.1 : 𝓞 K) : K) hc0 _ hval hxne⟩

/-! ## ★出典の紐付け(`.src`) -/

def mulHeight_mul_absNorm_of_scaled.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(整数表現に直しても公式は H(y) で成り立つ)",
    sectionId := "genell-prop-1-4" }

def exists_integral_repr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(どんな組にも整数表現があり、そこで公式が成り立つ)",
    sectionId := "genell-prop-1-4" }

def exists_integral_repr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Height.mulHeight_smul_eq_mulHeight(素朴高さは射影不変)"
      (.inMathlib "Height.mulHeight_smul_eq_mulHeight") 2,
    .citation "[mathlib]" "IsLocalization.exist_integer_multiples(共通分母)"
      (.inMathlib "IsLocalization.exist_integer_multiples") 2,
    .citation "[ABC3]" "mulHeight_mul_absNorm(段 C2c-2、§9-853)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mulHeight_mul_absNorm") 2,
    .implicitStep
      ("★どちらの道具も mathlib にあった。**段 C2c-2 の周りは在庫が厚い**" ++
       "——§9-853 でも absNorm_mul_finprod_finitePlace_eq_one がそのまま効いた") 3,
    .implicitStep
      ("★★残るのは段 C2c-1(超平面因子の引き戻しが座標の生成するイデアルであること)と " ++
       "段 C2c-3(計量の差が一様に有界)である。" ++
       "★C2c-3 は Prop 1.4 (iii)(§9-806)が本体を持っているので、実質は **C2c-1(幾何の側)** だけ") 5 ]

end ABC3.Found.GenEll
