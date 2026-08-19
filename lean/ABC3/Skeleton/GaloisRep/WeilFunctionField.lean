import ABC3.Found.GaloisRep.Translate

/-!
# スケルトン —— **関数体への引き戻し(Weil 対の葉)**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★葉は 2 件に絞れた

Weil 対の残りの構成が要求するのは、**関数体への 2 つの引き戻し**である。
★2026-08-20 の作業で、平行移動の側は**本体が `Found` に入った**:

| 段 | 場所 | 状態 |
|---|---|---|
| 生成点が曲線の非特異点であること | `Found/GaloisRep/GenericPoint.lean`(第 114) | ✅ |
| 平行移動した座標も非特異点であること | `Found/GaloisRep/Translate.lean`(第 115) | ✅ |
| 環準同型 `τ_Q : F[W] →+* F(W)` | 同上 | ✅ |
| **`τ_Q` が単射であること** | 本ファイル | ★**葉** |
| **乗法 `[n]` の引き戻し** | 本ファイル | ★**葉** |

★★平行移動の側は「分数体へ延ばすための単射性」だけになった
——`nonsingular_add` が効いて、環準同型そのものは mathlib から出た。

## ★★★★どちらも「式」で固定してある

空虚な posit を避けるため、値を明示式で固定している:

- `τ_Q` は加法公式そのもの——`Found` の `translateX`・`translateY`。
- `[n]^*` は**分点多項式**——`x([n]P) = Φ_n(x)/Ψ_n(x)²`。
  `WeierstrassCurve.Φ` と `WeierstrassCurve.ΨSq` は mathlib にある。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F]

/-- ★`F[X]` の元を関数体に送る。 -/
noncomputable def polyToFF (W : WeierstrassCurve.Affine F) (p : Polynomial F) : W.FunctionField :=
  algebraMap W.CoordinateRing W.FunctionField
    (CoordinateRing.mk W (Polynomial.C p))

/-! ## ★★★★★葉 1 —— 平行移動を分数体へ延ばす -/

/-- ★★★★★**平行移動の環準同型は単射である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが言えれば `IsFractionRing.lift` で `F(W) →ₐ[F] F(W)` に延び、
逆写像は `τ_{−Q}` で与えられるので**自己同型**になる。
★★環準同型そのものは第 115 ブロックで `Found` に入っている。 -/
theorem translateHom_injective (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    Function.Injective (translateHom W hQ) := by
  sorry

/-- ★★★★★★**平行移動 `τ_Q` は関数体の `F` 自己同型を誘導する**。 -/
theorem exists_translateAut (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField,
      τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀ := by
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

def polyToFF.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(F[X] の元を関数体に送る写像)",
    sectionId := "genell-thm-3-8" }

def translateHom_injective.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の環準同型が単射であること)",
    sectionId := "genell-thm-3-8" }

def translateHom_injective.needs : List ProofObligation :=
  [ .implicitStep
      "★★★`F[W]` は `F[X]` 上**階数 2 の自由加群**(mathlib の `CoordinateRing.basis`)なので整拡大である。整拡大では (0) の上にある素イデアルは (0) しかないので核は 0 になる——**mathlib での経路は未特定**(2026-08-20)(5-15 ブロック)" 19,
    .implicitStep
      "★別解: `translateX` が `F` 上超越的であることを直接示す。`genX_ne_const`(第 115、**Found に済**)の一般化になる(10-25 ブロック)" 19 ]

def exists_translateAut.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_translateAut.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3(平行移動が関数体の F 自己同型を誘導すること)"
      (.absent "mathlib に平行移動の関数体への引き戻しは 0 件(2026-08-20、EllipticCurve/ 配下を `translat|mulByN|scalarMul` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——平行移動の環準同型が単射であること)" 19,
    .implicitStep
      "★★★★**環準同型そのものは Found に入った**(第 115 ブロック `translateHom`)。生成点が曲線の非特異点であること(第 114)から mathlib の `nonsingular_add` がそのまま効き、`AdjoinRoot.lift` で構成できた。当初「平行移動は座標環の自己同型ではないので段が要る」と見積もっていた部分は**これで解消**(0 ブロック)" 19,
    .implicitStep
      "★★単射性から `IsFractionRing.lift` で `F(W) →ₐ[F] F(W)` へ延ばす(3-8 ブロック)" 19,
    .implicitStep
      "★★★全単射性——逆は `τ_{−Q}`。合成が恒等であることは生成元での計算に帰着する(5-15 ブロック)" 19 ]

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
      "★★★引き戻しが**環準同型である**こと、すなわち Weierstrass 多項式を殺すこと。★★★★平行移動の側と同じ道が使えるはずである——`[n]`(生成点)が曲線の点であることを言えば `AdjoinRoot.lift` に流せる。mathlib の `Point` の `nsmul` を生成点に適用する形になる(10-30 ブロック、当初 15-40 から下方修正)" 19 ]

end ABC3.Skeleton.GaloisRep
