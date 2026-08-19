import ABC3.Found.GaloisRep.TorsionIdeal

/-!
# スケルトン —— **関数体への引き戻し(Weil 対の葉)**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★(G5) の**葉**はここである

第 113 ブロックで `f_P`(`div(f_P) = n(P) − n(O)`)は `Found` に入った。
★Weil 対の残りの構成が要求するのは、**関数体への 2 つの引き戻し**である:

| 節点 | 内容 | mathlib |
|---|---|---|
| `exists_translateAut` | 平行移動 `τ_Q : f ↦ f∘(·+Q)` | ★**0 件**(2026-08-20 実測) |
| `exists_mulByNPullback` | 乗法 `[n]^* : f ↦ f∘[n]` | ★**0 件**(同上) |

★★どちらも**層番号最小の葉**である——`Found` と mathlib だけに依存する。

## ★★★★どちらも「式」で固定してある

空虚な posit を避けるため、両方とも**値を明示式で固定**した:

- `τ_Q` は加法公式そのもの——`addX`・`addY`(mathlib にある)で座標関数の像を指定。
- `[n]^*` は**分点多項式**——`x([n]P) = Φ_n(x)/Ψ_n(x)²`。
  `WeierstrassCurve.Φ` と `WeierstrassCurve.ΨSq` は mathlib にある。

★★★したがって「存在するとだけ言って中身が空」という退化は起きない。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

/-! ## ★座標関数 -/

/-- ★座標関数 `x`。 -/
noncomputable def coordX (W : WeierstrassCurve.Affine F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.mk W (Polynomial.C Polynomial.X))

/-- ★座標関数 `y`。 -/
noncomputable def coordY (W : WeierstrassCurve.Affine F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.mk W Polynomial.X)

/-- ★`F[X]` の元を関数体に送る。 -/
noncomputable def polyToFF (W : WeierstrassCurve.Affine F) (p : Polynomial F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.mk W (Polynomial.C p))

/-- ★点 `(x₀, y₀)` を通る割線の傾き(関数体の元として)。 -/
noncomputable def slopeFF (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) : W.FunctionField :=
  (coordY W - algebraMap F W.FunctionField y₀) / (coordX W - algebraMap F W.FunctionField x₀)

/-! ## ★★★★★葉 1 —— 平行移動 -/

/-- ★★★★★**平行移動 `τ_Q` の関数体への引き戻し**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★Weil 対 `e_n(P,Q) = g_P(·+Q)/g_P(·)` の「`·+Q`」がこれである。
★★座標関数の像を**加法公式で固定**してあるので、空虚な存在主張ではない。 -/
theorem exists_translateAut (W : WeierstrassCurve.Affine F) (x₀ y₀ : F)
    (hQ : W.Nonsingular x₀ y₀) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField,
      τ (coordX W)
        = (W.map (algebraMap F W.FunctionField)).toAffine.addX (coordX W)
            (algebraMap F W.FunctionField x₀) (slopeFF W x₀ y₀)
      ∧ τ (coordY W)
        = (W.map (algebraMap F W.FunctionField)).toAffine.addY (coordX W)
            (algebraMap F W.FunctionField x₀) (coordY W) (slopeFF W x₀ y₀) := by
  sorry

/-! ## ★★★★★葉 2 —— 乗法 `[n]` -/

/-- ★★★★★**乗法 `[n]` の関数体への引き戻し**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★Weil 対の構成の `g_P^n = f_P ∘ [n]` の「`∘ [n]`」がこれである。
★★`x([n]P) = Φ_n(x)/Ψ_n(x)²` で固定してある——mathlib の分点多項式を使う。 -/
theorem exists_mulByNPullback (W : WeierstrassCurve.Affine F) (n : ℤ) :
    ∃ μ : W.FunctionField →ₐ[F] W.FunctionField,
      μ (coordX W) = polyToFF W (W.Φ n) / polyToFF W (W.ΨSq n) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def coordX.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(関数体の座標関数 x)",
    sectionId := "genell-thm-3-8" }

def coordY.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(関数体の座標関数 y)",
    sectionId := "genell-thm-3-8" }

def polyToFF.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(F[X] の元を関数体に送る写像)",
    sectionId := "genell-thm-3-8" }

def slopeFF.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(割線の傾きの関数体版)",
    sectionId := "genell-thm-3-8" }

def exists_translateAut.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_translateAut.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3(平行移動が関数体の F 自己同型を誘導すること)"
      (.absent "mathlib に平行移動の関数体への引き戻しは 0 件(2026-08-20、EllipticCurve/ 配下を `translat|mulByN|scalarMul` で全文検索して 0 件)") 19,
    .implicitStep
      "★平行移動は**座標環の**自己同型ではない(無限遠点を動かす)。関数体の側でしか定義できず、その事実自体が段である(5-15 ブロック)" 19,
    .implicitStep
      "★★加法公式が関数体の元として well-defined であること——分母 `x − x₀` が 0 でないこと(座標環が整域であることは mathlib にある: `instIsDomainCoordinateRing`)(5-10 ブロック)" 19,
    .implicitStep
      "★★★τ が全単射であること(逆は `τ_{−Q}`)。群作用 `τ_{Q+Q'} = τ_Q ∘ τ_{Q'}` は下流の双線型性で消費される(5-15 ブロック)" 19 ]

def exists_mulByNPullback.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_mulByNPullback.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.4 / Exercise 3.7(乗法 [n] の関数体への作用)"
      (.absent "mathlib に [n] の関数体への引き戻しは 0 件(2026-08-20、同上の検索)") 19,
    .implicitStep
      "★分点多項式 `Φ_n` / `ΨSq_n` は mathlib にある(`WeierstrassCurve.Φ`・`WeierstrassCurve.ΨSq`、2026-08-20 実測)。**足場はある**(0 ブロック)" 19,
    .implicitStep
      "★★`y([n]P)` の側の式(`ω_n/ψ_n³` 型)——mathlib の在庫は未測定(5-15 ブロック)" 19,
    .implicitStep
      "★★★引き戻しが**環準同型である**こと、すなわち Weierstrass 多項式を殺すこと。これが本体の計算である(15-40 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
