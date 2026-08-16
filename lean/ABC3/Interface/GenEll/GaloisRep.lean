import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.InfinitePlace.Basic

/-!
# [GenEll] §3 —— l 捩れ点への Galois 表現の `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.19。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★なぜ `Interface` なのか —— statement を書くだけで語彙が足りない

`Theorem 3.8` の statement に現れる語のうち、**mathlib に無いもの**を目視で数えた:

| statement の語 | 由来 | mathlib(2026-08-16 実測) |
|---|---|---|
| `ht^Falt([E_L])`(Faltings 高さ) | Arakelov 理論(§1 の枠組み) | ★**0 件** |
| `Gal(ℚ̄/L) → GL₂(ℤ_l)` | l 進 Tate 加群への Galois 表現 | ★**0 件**(`EllipticCurve/*.lean` の `torsion` は 2-torsion polynomial のみ) |
| compactly bounded subset | Example 1.3, (ii) | —(§1 待ち) |
| Galois-finite | Example 1.3, (i) | —(§1 待ち) |
| local heights / potentially multiplicative reduction | Remark 3.3.1 / §3 | Tate 曲線が無い |

★**したがって `Theorem 3.8` のスケルトンは、それ自体が `Interface` に支えられる。**
本ファイルはその `Interface` である。

## ★これは「先送り」であって「解決」ではない

`Interface` は**何を作らねばならないかを型で固定する**装置であって、
作る手間を減らすものではない。`node tools/check.mjs` の
「Interface 実装待ち」に本構造体が並ぶことが、その事実の機械的な表示である。

★**公開プロジェクトとの関係**: l 捩れ点への Galois 表現は
**FLT プロジェクト**(Buzzard、EPSRC 2029-09 まで)のブループリント
§2.5「Galois representations and elliptic curves」・§3「Reducibility of p-torsion of the Frey curve」
と重なる(2026-08-16 取得)。★**独立に作ると重複投資になる。**
ただし「ブループリントに章がある」は「Lean で完成している」ではない——
使う前に clone して `sorry` を数えること(`ResearchPaper/lean-ecosystem.json` の規律)。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **`Theorem 3.8` の statement を書くのに要る語彙**を受ける `Interface`。

★意図的に「述語の束」として取っている——`M_ell(ℚ̄)` や `Gal(ℚ̄/L)` の**中身**を型に出すと、
Arakelov 理論と Tate 加群を抱え込む。**境界は薄く取る。**

★各フィールドは原文 p.19 の statement に**実際に現れる語**に 1 対 1 で対応させた。
現れない語は入れていない(仮説の強化を避けるため)。 -/
structure TorsionGaloisRepData where
  /-- `M_ell(ℚ̄)` —— `ℚ̄` 上の楕円曲線の同型類。 -/
  EllClass : Type
  /-- 数体 `L ⊆ ℚ̄` 上の楕円曲線 `E_L`。 -/
  Curve : Type
  /-- `E_L ↦ [E_L] ∈ M_ell(ℚ̄)`。 -/
  cls : Curve → EllClass
  /-- `d ≝ [L : ℚ]`。 -/
  degOfDefinition : Curve → ℕ
  /-- Faltings 高さ `ht^Falt`。★これが Arakelov 理論を要求する箇所である。 -/
  faltingsHeight : EllClass → ℝ
  /-- 原文「`E_L` has at least one prime of potentially multiplicative reduction」。 -/
  HasPotMultRed : Curve → Prop
  /-- 原文「`l` is prime to the local heights of `E_L` at all of its primes of
  potentially multiplicative reduction」。 -/
  PrimeToLocalHeights : Curve → ℕ → Prop
  /-- 原文「compactly bounded subset [cf. Example 1.3, (ii)]」。 -/
  CompactlyBounded : Set EllClass → Prop
  /-- 原文「Galois-finite [cf. Example 1.3, (i)] subset」。 -/
  GaloisFinite : Set EllClass → Prop
  /-- ★**結論の述語**——「`Gal(ℚ̄/L) → GL₂(ℤ_l)` の像が `SL₂(ℤ_l)` を含む」。

  ★ここが l 進 Tate 加群への Galois 表現を要求する箇所である。
  ★`Found/GenEll/Lemma31.lean` の `lemma_3_1_iii` は
  **この述語の `mod l` 版を出すための群論**であり、原文もそこだけを群論に頼る。 -/
  ImageContainsSL2 : Curve → ℕ → Prop

/-- Track B は何を作らねばならないか。 -/
def TorsionGaloisRepData.waiting : WaitingFor :=
  { what := "楕円曲線の l 進 Tate 加群への Galois 表現 Gal(ℚ̄/L) → GL₂(ℤ_l)、Faltings 高さ ht^Falt、局所高さ・潜在的乗法還元(Tate 曲線)、および compactly bounded / Galois-finite(Example 1.3)"
    trackB := "Found/GenEll — ★mathlib は楕円曲線の群構造(`Affine.Point` の `AddCommGroup`)と分点多項式(`DivisionPolynomial`)までは持つが、`E[n] ≅ (ℤ/n)²` と Galois 作用は無い(2026-08-16 実測)。★FLT プロジェクト(EPSRC 2029 まで)のブループリント §2.5・§3 と重なるので、独立に作る前に突き合わせること" }

/-! ## ★出典の紐付け(`.src`) -/

def TorsionGaloisRepData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8",
    sectionId := "genell-thm-3-8" }

end ABC3.Interface.GenEll
