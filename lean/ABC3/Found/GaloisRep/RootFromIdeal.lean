import ABC3.Found.GaloisRep.PointSpectrum

/-!
# Galois (G5) 第 139 ブロック —— **★★★★★★分数イデアルから `n` 乗根へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★目標が 2 節点に還元される

層 3 の目標は `∃ g, g^n = μ(f_P)` である。★本ブロックは**その最後の一歩**を出す:

    J^n = (μ f_P)   かつ   J = (g)     ⟹   ∃ g', g'^n = μ f_P

★★残るのは次の 2 つだけになる:

| 節点 | 内容 | 見積もり |
|---|---|---|
| D1 | `∃ J, J^n = (μ f_P)`(因子が `n` で割れる) | 20-50 |
| D2 | その `J` が単項であること(Abel–Jacobi) | 10-25 |

## ★★★★機構——単元は定数だから吸収できる

`(g^n) = (μ f_P)` なら両者は `F[W]` の単元倍で一致する。
★第 128 ブロック(`isUnit_coordinateRing`)により**単元は定数**であり、
★★`F` が代数閉だからその `n` 乗根 `b` が取れて、`g·b` に取り替えれば消える。

## ★★★★★★見通しが 1 つ変わった——**不分岐性は要らない**

当初は「`[n]` が不分岐だから `ord_Q(μ f) = ord_{[n]Q}(f)`」を示す予定だった。
★しかし一般の付値論で `ord_Q(μ f) = e_Q · ord_{[n]Q}(f)` であり、
`div(f_P) = n(P) − n(O)` から `ord_{[n]Q}(f_P) ∈ {n, −n, 0}` なので、
★★**`e_Q` が何であっても `n` で割れる**。不分岐性の証明は不要である。

★★★D2(類が自明であること)でも `e_Q` は消える——平行移動不変性から
`e_Q` は各ファイバー上で一定であり、`Σ_{Q ∈ [n]⁻¹(P)} Q = n²Q₀ = nP = 0` および
`Σ_{T ∈ E[n]} T = 0` により、共通因子 `e` を括り出しても類は 0 のままである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_nthRoot_of_spanSingleton_eq` | ★★★★★`(g^n) = (u)` ⟹ `∃ g', g'^n = u` |
| `exists_nthRoot_of_fractionalIdeal` | ★★★★★★**`J^n = (u)` かつ `J = (g)` ⟹ `∃ g', g'^n = u`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine nonZeroDivisors

variable {F : Type} [Field F]

/-- ★★★★★**`(g^n) = (u)` なら `u` は `n` 乗根を持つ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★単項分数イデアルが等しいことから `u = z · g^n`(`z` は `F[W]` の単元)。
★★第 128 で単元は定数 `c` であり、`F` が代数閉だから `c = b^n`。`g·b` が求めるもの。 -/
theorem exists_nthRoot_of_spanSingleton_eq [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    [W.IsElliptic] {u g : W.FunctionField} {n : ℕ} (hn : 1 ≤ n)
    (hg : FractionalIdeal.spanSingleton W.CoordinateRing⁰ (g ^ n)
        = FractionalIdeal.spanSingleton W.CoordinateRing⁰ u) :
    ∃ g' : W.FunctionField, g' ^ n = u := by
  obtain ⟨z, hz⟩ := FractionalIdeal.spanSingleton_eq_spanSingleton.1 hg
  obtain ⟨c, hc0, hc⟩ := isUnit_coordinateRing (u := (z : W.CoordinateRing)) z.isUnit
  obtain ⟨b, hb⟩ := IsAlgClosed.exists_pow_nat_eq (k := F) c (n := n) hn
  refine ⟨g * algebraMap F W.FunctionField b, ?_⟩
  rw [mul_pow, ← map_pow, hb, ← hz, Units.smul_def, Algebra.smul_def, hc,
    ← IsScalarTower.algebraMap_apply, mul_comm]

/-- ★★★★★★**`J^n = (u)` かつ `J` が単項なら `u` は `n` 乗根を持つ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが層 3(`g_P` の存在)の**最後の一歩**である。 -/
theorem exists_nthRoot_of_fractionalIdeal [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    [W.IsElliptic] {u g : W.FunctionField} {n : ℕ} (hn : 1 ≤ n)
    {J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField}
    (hJ : J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ u)
    (hgJ : J = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g) :
    ∃ g' : W.FunctionField, g' ^ n = u := by
  refine exists_nthRoot_of_spanSingleton_eq W (g := g) hn ?_
  rw [← FractionalIdeal.spanSingleton_pow, ← hgJ, hJ]

/-! ## ★出典の紐付け(`.src`) -/

def exists_nthRoot_of_fractionalIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——分数イデアルから n 乗根を取り出す最後の一歩)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
