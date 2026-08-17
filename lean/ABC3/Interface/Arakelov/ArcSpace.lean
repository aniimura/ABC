import ABC3.Interface.Arakelov.LineBundle
import Mathlib.Analysis.Complex.Basic

/-!
# Arakelov 理論のスケルトン(2/3)—— **アルキメデス側: `X^arc` と計量**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★2026-08-17 の撤回を反映する

`Interface/GenEll/ArithLineBundle.lean` は層 C を
**「`X^arc`(スキームの解析化)・GAGA」**と書いていた。
★★★**これは誤りだった。**同日、`Found/GenEll/ProjClosed.lean` /
`ArcModel.lean` で `X^arc` の位相とコンパクト性を
**商位相と多項式の連続性だけ**で構成した——正則構造も連接層も GAGA も現れない。

★したがって層 C の残りは以下の 3 件に縮む。

| # | 受けるもの | 現状(2026-08-17 実測) |
|---|---|---|
| C1 | `X^arc` の位相・`ι_X`・コンパクト性 | ★★**`ArcModel` を与えれば構成済み**(`Found/GenEll/ArcModel.lean`) |
| C2 | ℤ-固有 ⇒ 射影埋め込み(`ArcModel` の構成) | ★**未取得**。`Proj` の点の関手が要る |
| C3 | 可逆層の解析化 `L^arc` と hermitian 計量 | ★**未取得**。(B1) に従属 |

★★★**C2 が層 C の律速である**——`X^arc` そのものではなく、
**`X` が ℙⁿ に入ることを言う段**が残っている。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★C1 —— `X^arc` -/

/-- **(C1)** `X^arc`(複素点のなす位相空間)と複素共役 `ι_X`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Found/GenEll/ArcModel.lean` が**射影モデルを与えられたときに**
これを構成する(`ArcModel.topology` / `ArcModel.compactSpace`)。
★★★本 obligation はそれを `X` について一様に受ける形である。 -/
structure ArcSpaceData where
  /-- `X^arc` の台。 -/
  Arc : Scheme.{0} → Type
  /-- その位相。 -/
  topology : (X : Scheme.{0}) → TopologicalSpace (Arc X)
  /-- ★★**コンパクト**——原文が `Proposition 1.6` の証明で使う性質。 -/
  compact : (X : Scheme.{0}) → @CompactSpace (Arc X) (topology X)
  /-- 複素共役 `ι_X`。 -/
  conj : (X : Scheme.{0}) → Arc X → Arc X
  /-- ★`ι_X` は対合。 -/
  conj_involutive : ∀ (X : Scheme.{0}) (p : Arc X), conj X (conj X p) = p
  /-- ★`ι_X` は連続。 -/
  conj_continuous : ∀ (X : Scheme.{0}), @Continuous _ _ (topology X) (topology X) (conj X)

def ArcSpaceData.waiting : WaitingFor :=
  { what := "(C1) X^arc(複素点のなす位相空間)とその複素共役 ι_X、およびコンパクト性"
    trackB := "Found/Arakelov — ★★**射影モデル(ArcModel)を与えれば構成済み**(`Found/GenEll/ArcModel.lean` の `topology` / `compactSpace`、`ProjClosed.lean` の `isCompact_projZeroSet`。すべて sorry 0、2026-08-17)。★★★2026-08-16 に書いた『複素解析空間と GAGA が要る』という判定は**撤回した**——要るのは商位相と多項式の連続性だけである。★残るのは (C2)" }

/-! ## ★★★C2 —— ℤ-固有から射影モデルへ(層 C の律速) -/

/-- **(C2)** `X` が ℤ 上固有であることから**射影埋め込み**を作る段。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが層 C の律速である。**
`Found/GenEll/ArcModel.lean` の `ArcModel` は埋め込みを**データとして受けている**——
原文はそれを `X` の ℤ-固有性から**導いている**。

★具体的には `ℙⁿ` の**点の関手** `Hom(Spec ℂ, ℙⁿ) ≅ ℙ(ℂ^{n+1})` が要る。 -/
structure ProjectiveModelData where
  /-- 台となる `X^arc`。 -/
  toArcSpaceData : ArcSpaceData
  /-- `X` が ℤ 上固有であること(原文の仮定)。 -/
  ProperOverZ : Scheme.{0} → Prop
  /-- ★★★**ℤ-固有なら `X^arc` はコンパクトな射影モデルを持つ**。

  ★これを埋めるには `ℙⁿ` の点の関手が要る(`ResearchPaper/genell-goal.md` §9-18)。 -/
  hasProjectiveModel : ∀ X : Scheme.{0}, ProperOverZ X →
    ∃ (n : ℕ) (emb : toArcSpaceData.Arc X → (Fin (n + 1) → ℂ)),
      @Continuous (toArcSpaceData.Arc X) (Fin (n + 1) → ℂ) (toArcSpaceData.topology X)
        inferInstance emb ∧ Function.Injective emb

def ProjectiveModelData.waiting : WaitingFor :=
  { what := "(C2) X が ℤ 上固有であることから射影埋め込み X^arc ↪ ℙⁿ(ℂ) を得る段"
    trackB := "Found/Arakelov — ★★★**これが層 C の律速**。要るのは `ℙⁿ` の**点の関手** `Hom(Spec ℂ, Proj 𝒜) ≅ ℙ(ℂ^{n+1})`。★mathlib は `ProjectiveSpectrum/Basic.lean` に `Proj.awayι`(開埋め込み)と `basicOpenIsoSpec`、`ProjectiveSpectrum/Functor.lean` に次数付き環準同型からの関手性を持つが、**点の関手は無い**(2026-08-17 実測)。★構成 5 段は `ResearchPaper/genell-goal.md` §9-18 に記録済み" }

/-! ## ★★C3 —— 可逆層の解析化と hermitian 計量 -/

/-- **(C3)** 可逆層 `L` の解析化 `L^arc` と、その上の **`ι_X` と両立する計量**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★原文の `| − |_L` である。★★★**`ι_X` との両立が要である**——
それが無いと高さが `X(ℚ̄)` の上で well-defined にならない
(理由は `Found/GenEll/ArchBaseChange.lean` の測定: `v.embedding` は
`w.embedding|_F` **またはその共役**である)。 -/
structure HermitianMetricData where
  /-- 台となる `Pic`。 -/
  toPicardData : PicardData
  /-- 台となる `X^arc`。 -/
  toArcSpaceData : ArcSpaceData
  /-- 可逆層 `L` の点 `p ∈ X^arc` における計量の対数(Green 関数)。 -/
  logMetric : (X : Scheme.{0}) → toPicardData.Pic X → toArcSpaceData.Arc X → ℝ
  /-- ★連続であること(原文の「continuous function `|s|_L`」)。 -/
  logMetric_continuous : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X),
    @Continuous (toArcSpaceData.Arc X) ℝ (toArcSpaceData.topology X) inferInstance (logMetric X L)
  /-- ★★**`ι_X` と両立する**。 -/
  logMetric_conj : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (p : toArcSpaceData.Arc X),
    logMetric X L (toArcSpaceData.conj X p) = logMetric X L p
  /-- ★テンソル積で足し算になる(高さの加法性の源)。 -/
  logMetric_mul : ∀ (X : Scheme.{0}) (L M : toPicardData.Pic X) (p : toArcSpaceData.Arc X),
    logMetric X
      (@HMul.hMul _ _ _
        (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M) p
      = logMetric X L p + logMetric X M p

def HermitianMetricData.waiting : WaitingFor :=
  { what := "(C3) 可逆層の解析化 L^arc と、その上の ι_X-両立な hermitian 計量"
    trackB := "Found/Arakelov — ★(B1) の可逆層に従属する。★★**ι_X との両立**は我々が既に定式化し(`Found/GenEll/ArchConj.lean` の `IsConjInvariant`)、**それが要る理由も測った**——`v.embedding` は `w.embedding|_F` **またはその共役**なので、両立を課さないと高さが X(ℚ̄) 上で well-defined にならない(`ArchBaseChange.lean`、2026-08-17)。★連続性と加法性は因子表示では既に扱えている" }

/-! ## ★出典の紐付け(`.src`) -/

def ArcSpaceData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——X^arc の位相と ι_X のみ)",
    sectionId := "genell-def-1-1-i" }

def ProjectiveModelData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——ℤ-固有から射影モデルへ)",
    sectionId := "genell-def-1-1-i" }

def HermitianMetricData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——解析化と hermitian 計量)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
