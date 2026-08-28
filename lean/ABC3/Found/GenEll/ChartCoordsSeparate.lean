/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayGenerated
import ABC3.Found.GenEll.ChartPoints
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★正規化座標が点を分ける —— Northcott 路の穴が塞がった（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6-7。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★これは何か

`§9-849` で「段 E1d は Prop 1.4 (iv) の必須路ではない」と測り、
「座標が点を分ける」（`ext_of_chart`）を出したが、そこは
**`A⁰_{x_i}` 全体での一致**を仮定にしていた。
★`§9-850` で「`A⁰_{x_i}` は定数と正規化座標で生成される」（`ext_of_projCoord`）を示した。

★★本ファイルはその 2 本を繋ぐ:

> `X_{s_i}` の 2 つの `Spec F`-点は、**定数と正規化座標 `x_j/x_i` での値が一致すれば等しい**。

★★★これが `northcott_of_projModel`（`Found/GenEll/NorthcottCoord.lean`）が要求する
`hinj`（正規化座標が単射）の**そのままの形**である。

## ★これで Northcott 路の幾何側は塞がった

| 段 | 出典 |
|---|---|
| `ample` ⟹ 有限被覆 | `§9-817` |
| 分母を払う ⟹ 大域切断 | `§9-822`〜`§9-831` |
| チャートの座標環は有限生成 | `§9-832` |
| `A⁰_{x_i} → Γ(X, X_{s_i})` が全射 | `§9-847` |
| 座標が点を分ける | `§9-849` ＋ `§9-850` ＋ 本ファイル |

★★**残っているのは高さの側**である——C2c 前半（超平面因子の高さ ＝ 素朴高さ）と
`northcott_of_projModel` の `hcmp`（`mulHeight ≤ exp(ht + const)`）を繋ぐ段。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-- ★★★★★★★★★★**正規化座標が点を分ける** —— Northcott 路の幾何側の終点。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

`X_{s_i}` がアフィンで `A⁰_{x_i} → Γ(X, X_{s_i})` が全射（`§9-847`）なら、
★`X_{s_i}` の 2 つの `Spec F`-点は
**定数と正規化座標 `x_j/x_i` での値が一致すれば等しい**。

★★これが `northcott_of_projModel` の `hinj` の形である
——**スキームの射 `ψ : X ⟶ ℙᴺ_R` を作らずに済む**（`§9-849` の測定）。

★★★機構は `ext_of_projCoord`（`§9-850`——生成）と `ext_of_chart`（`§9-849`——全射性）の合成だけである。 -/
theorem ext_of_chart_coords (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i)))
    (hsurj : Function.Surjective (globalAwayHom M hM φ s i))
    {F : Type} [CommRing F]
    (p q : Spec (CommRingCat.of F) ⟶ (nonVanishing M (s i)).toScheme)
    (hc : ∀ c : R, p.appTop.hom (chartRingHom M hM φ s i (awayConst N R i c))
        = q.appTop.hom (chartRingHom M hM φ s i (awayConst N R i c)))
    (hx : ∀ j, p.appTop.hom (chartRingHom M hM φ s i (projCoord N R i j))
        = q.appTop.hom (chartRingHom M hM φ s i (projCoord N R i j))) :
    p = q := by
  have h : ((p.appTop).hom).comp (chartRingHom M hM φ s i)
      = ((q.appTop).hom).comp (chartRingHom M hM φ s i) :=
    ext_of_projCoord i _ _ hc hx
  exact ext_of_chart M hM φ s i haff hsurj p q (fun z => RingHom.congr_fun h z)

/-! ## ★出典の紐付け(`.src`) -/

def ext_of_chart_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(正規化座標が点を分ける)",
    sectionId := "genell-prop-1-4" }

def ext_of_chart_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ext_of_projCoord(A⁰_{x_i} は定数と座標で生成される、§9-850)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_projCoord") 2,
    .citation "[ABC3]" "ext_of_chart(座標が点を分ける、§9-849)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_chart") 2,
    .citation "[ABC3]" "northcott_of_projModel(消費側)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_projModel") 2,
    .implicitStep
      ("★これで Northcott 路の**幾何側**は塞がった" ++
       "(§9-817 の有限被覆 → §9-822〜§9-831 の大域化 → §9-832 の有限生成 →" ++
       " §9-847 の全射性 → 本ファイルの点分離)。" ++
       "★★残っているのは**高さの側**である——C2c 前半(超平面因子の高さ ＝ 素朴高さ)と " ++
       "northcott_of_projModel の hcmp(mulHeight ≤ exp(ht + const))を繋ぐ段") 6 ]

end ABC3.Found.GenEll
