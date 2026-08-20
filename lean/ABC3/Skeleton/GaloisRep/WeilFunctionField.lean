import ABC3.Found.GaloisRep.GenericNotTorsion
import ABC3.Found.GaloisRep.PointHom

/-!
# スケルトン —— **関数体への引き戻し(Weil 対の葉)**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★葉 1(平行移動)は**閉じた**

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
| 2 等分点の分解の存在 | `Found/GaloisRep/TwoTorsionDecomp.lean`(第 122) | ✅ |
| 2 等分点でない点の存在(代数閉・標数≠2) | `Found/GaloisRep/NotTwoTorsionPoint.lean`(第 123) | ✅ |
| **平行移動の単射性(仮定なし)** | 同上 | ✅ |
| **全単射性**(2 等分点以外) | `Found/GaloisRep/TranslateEquiv.lean`(第 120) | ✅ |
| **2 等分点での自己同型** | `Found/GaloisRep/TranslateAutAll.lean`(第 124) | ✅ |
| `[n]` の環準同型(`pointHom` の一般形から) | `Found/GaloisRep/PointHom.lean`(第 118) | ✅ |
| **生成点が捻れ点でないこと** | `Found/GaloisRep/GenericNotTorsion.lean`(第 125) | ✅ |
| `[n]` の引き戻し(群法則で固定) | `Found/GaloisRep/GenericNotTorsion.lean`(第 125) | ✅ |
| `x([n]P) = Φ_n/ΨSq_n`(式の側) | 本ファイル | ★臨界路外 |

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

/-! ## ★臨界路外 —— `x([n]P) = Φ_n/ΨSq_n`(式の側) -/

/-- ★**`x([n]P) = Φ_n/ΨSq_n`**——古典的な式。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★★**2026-08-20: この節点は臨界路から外れた。**
Weil 対の構成が実際に消費するのは `f_P ∈ F[W]` への作用だけであり、
それは第 118・125 ブロックで**群法則によって固定された形で Found に入った**
(`exists_mulByNHom_charZero`)。
★式の側は真であるが、現在の道では消費されない——将来の層が要求したときのために残す。 -/
theorem exists_mulByNPullback (W : WeierstrassCurve.Affine F) (n : ℤ) :
    ∃ μ : W.FunctionField →ₐ[F] W.FunctionField,
      μ (coordX W) = polyToFF W (W.Φ n) / polyToFF W (W.ΨSq n) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def polyToFF.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(F[X] の元を関数体に送る写像)",
    sectionId := "genell-thm-3-8" }

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
      "★★★★★★★★**2026-08-20: 本節点は臨界路に戻った**。(G5) の非退化性が `deg[n] = n²` を要求し(第 196・197)、それが `[F(x) : F(x∘[n])] ≤ n²`、すなわち `x` が `F(x_n)` 上でモニック多項式 `Φ_n(X) − x_n·ΨSq_n(X)`(第 198 で次数ちょうど `n²`)の根であることに帰着したためである。★消費されるのは **`x` の側だけ**——`y([n]P)` の式は要らない" 19,
    .implicitStep
      "★★★★★★**2026-08-20 の実測(足場)**: `x(2P) = Φ₂/ΨSq₂` は第 199 で証明できた(`ΨSq₂(x) = (2y+a₁x+a₃)²` を出してから分母を払う)。★`y` を消す Kummer の 2 公式(和・積)は第 200・201 で取れた。★★`n = 3` の梯子恒等式 `Φ₃·X = Φ₂²X² − b₄Φ₂X·ΨSq₂ − …` も `ring` で通った(第 202、2.6 秒)。★★★**mathlib の `preΨ₄` の定義が群法則と噛み合っている**ことは確認済みである" 19,
    .implicitStep
      "★★★★★★★★**2026-08-20 の実測(壁の所在)**: `n = 4` の段(`n → 2n`)を `b₂ b₄ b₆ b₈` を自由変数にして `ring` に掛けると **124 秒かけて失敗**した——`b₂b₆ − b₄² = 4b₈` の関係が要る。★`a₁ … a₆` に展開すると既定の heartbeats を超える。★★つまり**一般の `n` を多項式展開で押すのは実際的でない**。★★★構造的には EDS 恒等式 `W(m+n)W(m−n)W(r)² = W(m+r)W(m−r)W(n)² − W(n+r)W(n−r)W(m)²` が要るが、これは **mathlib 自身の TODO** である(`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean` の `TODO: prove that normEDS satisfies IsEllDivSequence`、2026-08-20 実測)。★★★★見積もりは上振れする(30-80 ブロック)。上流に入れるのが本筋である" 19,
    .implicitStep
      "★`x([n]P)` の超越性——`n` 等分点での 1 点評価(`ΨSq_n` が消えて `Φ_n` が消えない)で出る見込み。第 117 と同じ型(5-15 ブロック)" 19,
    .implicitStep
      "★★`y([n]P)` の側の式(`ω_n/ψ_n³` 型)は **Weil 対には要らない**ことが分かった(上記)。なお mathlib も `ω_n` を持たない(docstring に `TODO: the bivariate polynomials ωₙ`、2026-08-20 実測)" 19 ]

end ABC3.Skeleton.GaloisRep
