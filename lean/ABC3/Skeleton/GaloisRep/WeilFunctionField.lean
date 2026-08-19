import ABC3.Found.GaloisRep.TranslateAut

/-!
# スケルトン —— **関数体への引き戻し(Weil 対の葉)**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★葉 1(平行移動)は**ほぼ閉じた**

2026-08-20 の作業で、平行移動の側は次まで `Found` に入った:

| 段 | 場所 | 状態 |
|---|---|---|
| 生成点が曲線の非特異点であること | `Found/GaloisRep/GenericPoint.lean`(第 114) | ✅ |
| 平行移動した座標も非特異点であること | `Found/GaloisRep/Translate.lean`(第 115) | ✅ |
| 環準同型 `τ_Q : F[W] →+* F(W)` | 同上 | ✅ |
| 単射性を超越性に帰着 | `Found/GaloisRep/Transcendence.lean`(第 116) | ✅ |
| **`translateX` の超越性**(`Q` が 2 等分点でない) | `Found/GaloisRep/TranslateAut.lean`(第 117) | ✅ |
| **関数体の自己準同型 `translateFieldHom`** | 同上 | ✅ |
| `Q` が **2 等分点**のときの単射性 | 本ファイル | ★**葉** |
| **全単射性**(自己同型であること) | 本ファイル | ★**葉** |
| **乗法 `[n]` の引き戻し** | 本ファイル | ★**葉** |

★超越性の証明は **`−Q` での 1 点評価**だけで済んだ——
`u = x − x₀` は `−Q` で消えるが `A = (x−x₀)²·X(τ_Q)` は消えない、という一点で決まる。
★★2 等分点では `−Q = Q` なので `A` も消え、この議論が使えない。

## ★★★★式で固定してある

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

/-! ## ★★★★★葉 1a —— 2 等分点の場合 -/

/-- ★★★★★**`Q` が 2 等分点のときの単射性**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 117 ブロックの `−Q` 評価は `−Q = Q` になるので使えない。
★★別解: 2 等分点 `Q` を `Q = Q' + Q''`(どちらも 2 等分点でない)と分解し、
`τ_Q = τ_{Q'} ∘ τ_{Q''}` として単射性を合成で得る。 -/
theorem translateHom_injective_twoTorsion (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ = y₀) :
    Function.Injective (translateHom W hQ) := by
  sorry

/-! ## ★★★★★葉 1b —— 全単射性 -/

/-- ★★★★★★**平行移動 `τ_Q` は関数体の `F` 自己同型である**。

★自己準同型 `translateFieldHom` は第 117 ブロックで `Found` に入っている。
★★逆は `τ_{−Q}` であり、合成が恒等であることは**生成元での計算**に帰着する。 -/
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

def translateHom_injective_twoTorsion.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——2 等分点での平行移動の単射性)",
    sectionId := "genell-thm-3-8" }

def translateHom_injective_twoTorsion.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★2 等分点でない場合は **Found に済**(第 117 `translateHom_injective`)。`−Q` での 1 点評価で決まった(0 ブロック)" 19,
    .implicitStep
      "★★2 等分点では `−Q = Q` なので `ev_{−Q}(A) = 0` になり、第 117 の議論が使えない。分解 `Q = Q' + Q''`(どちらも 2 等分点でない)を取って合成で回す(5-15 ブロック)" 19,
    .implicitStep
      "★分解の存在——体が無限なら 2 等分点でない点は十分ある。`E[2]` は高々 4 点(第 65-72 の `E[n] ≃ (ℤ/n)²` から)(3-8 ブロック)" 19 ]

def exists_translateAut.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_translateAut.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3(平行移動が関数体の F 自己同型を誘導すること)"
      (.absent "mathlib に平行移動の関数体への引き戻しは 0 件(2026-08-20、EllipticCurve/ 配下を `translat|mulByN|scalarMul` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——2 等分点での平行移動の単射性)" 19,
    .implicitStep
      "★★★★★★**環準同型・単射性・分数体への延長はすべて Found に入った**(第 114-117)。生成点が曲線の非特異点であることから mathlib の `nonsingular_add` が効き、超越性は `−Q` での 1 点評価で決まった。当初の見積もり(20-60 ブロック)は**大きく下方修正**される(0 ブロック)" 19,
    .implicitStep
      "★★全単射性——逆は `τ_{−Q}`。合成 `τ_{−Q} ∘ τ_Q` が恒等であることは、生成元 `x`・`y` での計算(mathlib の群法則の結合則)に帰着する(5-15 ブロック)" 19,
    .implicitStep
      "★`F` 代数射であること——`translateHom` が `F` 代数射なので延長も `F` を固定する(2-5 ブロック)" 19 ]

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
      "★★★★**平行移動と同じ道が使える**——`[n]`(生成点)が曲線の点であること(第 114 の `genericPoint` に mathlib の `nsmul` を当てる)を言えば `AdjoinRoot.lift` に流せる。単射性も同様に 1 点評価で出るはずである(10-25 ブロック、当初 15-40 から下方修正)" 19,
    .implicitStep
      "★★`y([n]P)` の側の式(`ω_n/ψ_n³` 型)——mathlib の在庫は未測定(5-15 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
