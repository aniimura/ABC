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
| C1 | `X^arc` の位相と `ι_X` | ★★**`ArcModel` を与えれば構成済み**(`Found/GenEll/ArcModel.lean`) |
| C2 | **ℤ-固有 ⇒ `X^arc` コンパクト** | ★**未取得**。★★★固有は射影を含意しないので Chow の補題が要る |
| C3 | 可逆層の解析化 `L^arc` と hermitian 計量 | ★**未取得**。(B1) に従属 |

★★★**C2 が層 C の律速である**——`X^arc` そのものではなく、
**固有性からコンパクト性を出す段**が残っている。

## ★★★2026-08-17: 自分の定式化の誤りを 2 つ直した

1. C1 に `compact : ∀ X, CompactSpace (Arc X)` と書いていた。
   ★★**`X = 𝔸¹` で偽である。**コンパクト性は固有性の帰結であって、
   `X^arc` の定義の一部ではない。C2 へ移した。
2. C2 を「ℤ-固有 ⇒ **射影埋め込み**」と書いていた。
   ★★★**固有は射影を含意しない。**原文は
   > Let X be a normal scheme which is proper and flat over Spec(Z)

   と書いているだけである(物理 p.3、実測)。
   ★正しい結論は「コンパクト」であり、射影の場合は**その一つの道**にすぎない
   (一般には Chow の補題を経る)。
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
  /-- 複素共役 `ι_X`。 -/
  conj : (X : Scheme.{0}) → Arc X → Arc X
  /-- ★`ι_X` は対合。 -/
  conj_involutive : ∀ (X : Scheme.{0}) (p : Arc X), conj X (conj X p) = p
  /-- ★`ι_X` は連続。 -/
  conj_continuous : ∀ (X : Scheme.{0}), @Continuous _ _ (topology X) (topology X) (conj X)
  /-- ★★★**台は複素点の集合そのもの**——これを言わないと `Arc` は何でもよくなる。

  ★`Interface/` は `Found/` を import できないので、`complexPoints` を参照せず
  ここに書き下す(中身は同じ)。 -/
  equivComplexPoints : (X : Scheme.{0}) → Arc X ≃ (Spec (CommRingCat.of ℂ) ⟶ X)
  /-- ★アフィン `Spec A` の点における、切断 `a : A` の**評価** `a(p) ∈ ℂ`。 -/
  evalAffine : (A : CommRingCat.{0}) → Arc (Spec A) → A → ℂ
  /-- ★★★**位相を固定する(その 1)**——アフィンでは**各点収束の位相**である。

  ★★★これが無いと**離散位相で埋まってしまう**
  (`conj_continuous` は離散位相なら自明に成り立つ)。
  ★原文の `X^arc` は複素多様体の位相であって、離散ではない。 -/
  topology_affine : ∀ A : CommRingCat.{0},
    topology (Spec A) = TopologicalSpace.induced (evalAffine A) Pi.topologicalSpace
  /-- 射に沿った複素点の移送。 -/
  map : {X Y : Scheme.{0}} → (X ⟶ Y) → Arc X → Arc Y
  /-- ★★★**位相を固定する(その 2)**——開埋め込みは位相を誘導する。

  ★★アフィンだけを縛っても、貼り合わせた `X`(例えば `ℙⁿ`)の位相は決まらない。
  ★開被覆に沿って誘導であることを課すと、そこも決まる。
  ★★★これが無いと、射影的な `X` で**密着位相**が通ってしまう。 -/
  topology_openImmersion : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y), IsOpenImmersion f →
    topology X = TopologicalSpace.induced (map f) (topology Y)

def ArcSpaceData.waiting : WaitingFor :=
  { what := "(C1) X^arc(複素点のなす位相空間)とその複素共役 ι_X"
    trackB := "Found/Arakelov — ★★**射影モデル(ArcModel)を与えれば構成済み**(`Found/GenEll/ArcModel.lean` の `topology` / `compactSpace`、`ProjClosed.lean` の `isCompact_projZeroSet`。すべて sorry 0、2026-08-17)。★★★2026-08-16 に書いた『複素解析空間と GAGA が要る』という判定は**撤回した**——要るのは商位相と多項式の連続性だけである。★残るのは (C2)" }

/-! ## ★★★C2 —— ℤ-固有からコンパクト性へ(層 C の律速) -/

/-- **(C2)** `X` が ℤ 上固有・平坦であることから **`X^arc` のコンパクト性**を出す段。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが層 C の律速である。**
`Found/GenEll/ArcModel.lean` の `ArcModel` は埋め込みを**データとして受けている**——
原文はコンパクト性を `X` の ℤ-固有性から**導いている**。

★★★**「射影埋め込みを作る」と書いてはならない**——原文の仮定は
`proper and flat over Spec(Z)` であって、**固有は射影を含意しない**。
★射影の場合は道の一つであり(`projectiveCase`、`ArcModel` で証明済)、
一般へ渡すには Chow の補題が要る。 -/
structure ProjectiveModelData where
  /-- 台となる `X^arc`。 -/
  toArcSpaceData : ArcSpaceData
  /-- `X` が ℤ 上固有・平坦であること(原文の仮定、物理 p.3)。 -/
  ProperFlatOverZ : Scheme.{0} → Prop
  /-- ★★★**原文が実際に使う結論**——ℤ-固有なら `X^arc` はコンパクト。

  ★★★**「射影的」ではない。**原文は
  > Let X be a normal scheme which is proper and flat over Spec(Z)

  と書いており、**固有は射影を含意しない**。 -/
  properFlat_compact : ∀ X : Scheme.{0}, ProperFlatOverZ X →
    @CompactSpace (toArcSpaceData.Arc X) (toArcSpaceData.topology X)
  /-- ★★**射影的な場合**——`ℙⁿ` への埋め込みからコンパクト性が出る道。

  ★`Found/GenEll/ArcModel.lean` は**この形の入力からコンパクト性を証明済み**である。
  ★★一般の固有 `X` へは Chow の補題(固有 → 射影的な変更を持つ)を経る。 -/
  projectiveCase : ∀ X : Scheme.{0}, ∀ n : ℕ,
    ∀ emb : toArcSpaceData.Arc X → (Fin (n + 1) → ℂ),
      @Continuous (toArcSpaceData.Arc X) (Fin (n + 1) → ℂ) (toArcSpaceData.topology X)
        inferInstance emb → Function.Injective emb →
      IsClosed (Set.range emb) → Bornology.IsBounded (Set.range emb) →
      @CompactSpace (toArcSpaceData.Arc X) (toArcSpaceData.topology X)

def ProjectiveModelData.waiting : WaitingFor :=
  { what := "(C2) X が ℤ 上固有・平坦であることから X^arc のコンパクト性を得る段。★固有は射影を含意しないので、射影の場合(ArcModel、実装済)から一般へ渡すには Chow の補題が要る"
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
  /-- ★★★**可逆層 `L` の上の計量の全体**。

  ★★★**2026-08-17 の設計修正**: 以前は計量を `L` ごとに **1 つ固定**していたので、
  「すべて 0」で埋まった(実測)。★**1 つの直線束には計量が多数ある**——
  それを型に出さないと `APic`(層と計量の**対**)の意味が消える。 -/
  Metric : (X : Scheme.{0}) → toPicardData.Pic X → Type
  /-- ★計量は存在する。 -/
  metric_nonempty : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X), Nonempty (Metric X L)
  /-- 計量の対数(Green 関数)。 -/
  logMetric : (X : Scheme.{0}) → (L : toPicardData.Pic X) → Metric X L →
    toArcSpaceData.Arc X → ℝ
  /-- ★連続であること(原文の「continuous function `|s|_L`」)。 -/
  logMetric_continuous : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (m : Metric X L),
    @Continuous (toArcSpaceData.Arc X) ℝ (toArcSpaceData.topology X) inferInstance
      (logMetric X L m)
  /-- ★★**定数倍による作用**——計量は正の定数倍で動かせる。 -/
  scale : (X : Scheme.{0}) → (L : toPicardData.Pic X) → ℝ → Metric X L → Metric X L
  /-- ★★★**定数倍は Green 関数を平行移動する**。

  ★★★**これが退化を殺す。**`logMetric ≡ 0` だと `0 = 0 + c` を全 `c` で要求し、矛盾する。 -/
  logMetric_scale : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (c : ℝ) (m : Metric X L)
    (p : toArcSpaceData.Arc X),
    logMetric X L (scale X L c m) p = logMetric X L m p + c
  /-- ★★**`ι_X` と両立する計量**であること。 -/
  IsConjCompatible : (X : Scheme.{0}) → (L : toPicardData.Pic X) → Metric X L → Prop
  /-- ★その意味——共役で値が変わらない。 -/
  isConjCompatible_iff : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (m : Metric X L),
    IsConjCompatible X L m ↔
      ∀ p : toArcSpaceData.Arc X,
        logMetric X L m (toArcSpaceData.conj X p) = logMetric X L m p
  /-- ★テンソル積で計量は掛かる(Green 関数は足される)——高さの加法性の源。 -/
  tensorMetric : (X : Scheme.{0}) → (L M : toPicardData.Pic X) → Metric X L → Metric X M →
    Metric X (@HMul.hMul _ _ _
      (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
  /-- ★★Green 関数の加法性。 -/
  logMetric_tensor : ∀ (X : Scheme.{0}) (L M : toPicardData.Pic X)
    (mL : Metric X L) (mM : Metric X M) (p : toArcSpaceData.Arc X),
    logMetric X _ (tensorMetric X L M mL mM) p
      = logMetric X L mL p + logMetric X M mM p

def HermitianMetricData.waiting : WaitingFor :=
  { what := "(C3) 可逆層の解析化 L^arc と、その上の ι_X-両立な hermitian 計量"
    trackB := "Found/Arakelov — ★(B1) の可逆層に従属する。★★**ι_X との両立**は我々が既に定式化し(`Found/GenEll/ArchConj.lean` の `IsConjInvariant`)、**それが要る理由も測った**——`v.embedding` は `w.embedding|_F` **またはその共役**なので、両立を課さないと高さが X(ℚ̄) 上で well-defined にならない(`ArchBaseChange.lean`、2026-08-17)。★連続性と加法性は因子表示では既に扱えている" }

/-! ## ★出典の紐付け(`.src`) -/

def ArcSpaceData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——X^arc の位相と ι_X のみ)",
    sectionId := "genell-def-1-1-i" }

def ProjectiveModelData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——ℤ-固有からコンパクト性へ)",
    sectionId := "genell-def-1-1-i" }

def HermitianMetricData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 C——解析化と hermitian 計量)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
