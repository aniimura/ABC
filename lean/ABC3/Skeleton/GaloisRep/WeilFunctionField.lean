import ABC3.Found.GaloisRep.TranslateComp
import ABC3.Found.GaloisRep.PointHom

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
| 合成則 `τ_{Q₁} ∘ τ_{Q₂} = τ_{Q₁+Q₂}` | `Found/GaloisRep/TranslateComp.lean`(第 121) | ✅ |
| 2 等分点の単射性(**分解があれば**) | 同上 | ✅ |
| **2 等分点の分解の存在** | 本ファイル | ★**葉** |
| **全単射性**(2 等分点以外) | `Found/GaloisRep/TranslateEquiv.lean`(第 120) | ✅ |
| 2 等分点での自己同型 | 本ファイル | ★**葉** |
| `[n]` の環準同型(`pointHom` の一般形から) | `Found/GaloisRep/PointHom.lean`(第 118) | ✅ |
| **生成点が捻れ点でないこと** | 本ファイル | ★**葉** |
| **`x([n]P) = Φ_n/ΨSq_n`** | 本ファイル | ★**葉** |

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

/-! ## ★★★★★葉 1b —— 2 等分点での自己同型 -/

/-- ★★★★★★**平行移動 `τ_Q` は関数体の `F` 自己同型である**(`Q` が 2 等分点のとき)。

★★★★★★2 等分点でない場合は **`Found` に入った**(第 120 ブロック
`exists_translateAut_of_not_twoTorsion`)——逆は `τ_{−Q}` で、合成が恒等であることは
`functionField_algHom_ext` により生成元での計算に落ち、`(G + Z) + Q = G` という
mathlib の群法則で片付いた。
★2 等分点では単射性がまだ無いので、ここだけが残っている。 -/
theorem exists_translateAut_twoTorsion (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ = y₀) :
    ∃ τ : W.FunctionField ≃ₐ[F] W.FunctionField,
      τ (coordX W) = translateX W x₀ y₀ ∧ τ (coordY W) = translateY W x₀ y₀ := by
  sorry

/-! ## ★★★★★葉 2a —— 生成点は捩れ点でない -/

/-- ★★★★★★**生成点は捩れ点ではない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 118 ブロックの `exists_mulByNHom` はこれを仮定として受けている。
★★証明の道: `n` 等分点の `x` 座標は **`F` の代数閉包の元**である
(曲線が `F` 上定義されているから)。★★★我々は `E[n] ≃ (ℤ/n)²`(第 65-72、**Found に済**)を
持つので、`E[n](F̄)` と `E[n](F(W)‾)` がどちらも `n²` 個であり、
`Point.map` の単射性から**両者は一致する**。
★★★★すると `coordX` が `F` 上代数的になり、第 116 の `coordX_transcendental` に矛盾する。 -/
theorem genericPoint_not_torsion (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (n : ℕ) (hn : 1 ≤ n) : n • genericPoint W ≠ 0 := by
  sorry

/-! ## ★★★★★葉 2b —— 乗法 `[n]` -/

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
      "★★★★★★**分解があれば単射性は出る**——第 121 ブロック `translateHom_injective_of_decomp` が **Found に入った**。合成則 `τ_{Q₁} ∘ τ_{Q₂} = τ_{Q₁+Q₂}` と「体の自己準同型は単射」から出る(0 ブロック)" 19,
    .implicitStep
      "★★**残るのは分解の存在だけ**—— 2 等分点 `Q₃` に対し `Q₃ = Q₁ + Q₂`(どちらも 2 等分点でなく、どちらもアフィン点)を取る。`E[2]` は高々 4 点(第 65-72 の `E[n] ≃ (ℤ/n)²` から)なので、体が無限なら取れる(3-8 ブロック)" 19 ]

def exists_translateAut_twoTorsion.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_translateAut_twoTorsion.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3(平行移動が関数体の F 自己同型を誘導すること)"
      (.absent "mathlib に平行移動の関数体への引き戻しは 0 件(2026-08-20、EllipticCurve/ 配下を `translat|mulByN|scalarMul` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——2 等分点での平行移動の単射性)" 19,
    .implicitStep
      "★★★★★★**環準同型・単射性・分数体への延長はすべて Found に入った**(第 114-117)。生成点が曲線の非特異点であることから mathlib の `nonsingular_add` が効き、超越性は `−Q` での 1 点評価で決まった。当初の見積もり(20-60 ブロック)は**大きく下方修正**される(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**2 等分点以外の全単射性は Found に入った**(第 120 ブロック `translateAlgEquiv`)。`functionField_algHom_ext`(第 119)で生成元の計算に落ち、`(G + Z) + Q = G` という mathlib の群法則で片付いた(0 ブロック)" 19,
    .implicitStep
      "★残るのは 2 等分点の場合だけであり、それは単射性(上の辺)が出れば同じ道で閉じる(3-8 ブロック)" 19 ]

def genericPoint_not_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——生成点が捩れ点でないこと)",
    sectionId := "genell-thm-3-8" }

def genericPoint_not_torsion.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★`E[n] ≃+ (ℤ/n)²`(第 65-72 ブロック、**Found に済**)がそのまま使える。代数閉体上で `#E[n] = n²` であることが両側で効く(0 ブロック)" 19,
    .implicitStep
      "★★`Point.map` が単射であること(mathlib の `Point.map_injective`)と、有限集合の間の単射が同数なら全単射であること(5-10 ブロック)" 19,
    .implicitStep
      "★★★`coordX` が `F` 上超越的であること(第 116 `coordX_transcendental`、**Found に済**)が矛盾を出す(0 ブロック)" 19,
    .implicitStep
      "★代数閉包の取り扱い——`F(W)` の代数閉包の中で `F̄` を見る。mathlib の `AlgebraicClosure` と `IsAlgClosed` の橋が要る(5-15 ブロック)" 19 ]

def exists_mulByNPullback.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の関数体への引き戻し)",
    sectionId := "genell-thm-3-8" }

def exists_mulByNPullback.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.4 / Exercise 3.7(乗法 [n] の関数体への作用)"
      (.absent "mathlib に [n] の関数体への引き戻しは 0 件(2026-08-20、同上の検索)") 19,
    .implicitStep
      "★分点多項式 `Φ_n` / `ΨSq_n` は mathlib にある(`WeierstrassCurve.Φ`・`WeierstrassCurve.ΨSq`、2026-08-20 実測)。**足場はある**(0 ブロック)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——生成点が捻れ点でないこと)" 19,
    .implicitStep
      "★★★★★★**環準同型は Found に入った**(第 118 ブロック `exists_mulByNHom`)。`pointHom`(点 ⇒ 環準同型)として一般化したところ、平行移動と**完全に同じ道**が `[n]` にも効いた。単射性も `pointHom_injective_of_transcendental` に帰着済み(0 ブロック)" 19,
    .implicitStep
      "★★`x([n]P)` が `Φ_n/ΨSq_n` と一致すること——mathlib の分点多項式と群法則を結ぶ段。mathlib にはこのリンクが無い(2026-08-20 実測)(10-25 ブロック)" 19,
    .implicitStep
      "★`x([n]P)` の超越性——`n` 等分点での 1 点評価(`ΨSq_n` が消えて `Φ_n` が消えない)で出る見込み。第 117 と同じ型(5-15 ブロック)" 19,
    .implicitStep
      "★★`y([n]P)` の側の式(`ω_n/ψ_n³` 型)——mathlib の在庫は未測定(5-15 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
