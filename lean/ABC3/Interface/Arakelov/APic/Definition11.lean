/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Interface.Arakelov.ArcSpace

/-!
# APic —— `[GenEll] Definition 1.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★D1 —— `APic(X)` -/

/-- **(D1)** **算術直線束の群** `APic(X)` —— 可逆層と計量の対。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の `L̄ = (L, |−|_L)` そのものである。 -/
structure APicData where
  /-- 台となる計量(その中に `Pic` も入っている)。 -/
  toHermitianMetricData : HermitianMetricData
  /-- `APic(X)` の台。

  ★★★★2026-08-19 の修正——**宇宙が `Pic` と合っていなかった**。
  それまで `Type`(= `Type 0`)と書いてあったが、(B1) の `Pic : Scheme.{0} → Type 1` と
  `forgetMetric` / `forgetMetric_mk` で結ばれる以上、`APic` も `Type 1` でなければならない
  ——`Metric X L` は常に空でないので `L ↦ ofMetric X L m` は `Pic X ↪ APic X` を与え、
  `Type 1 ↪ Type 0` を要求してしまう。★これは型の誤記であって数学の変更ではない。 -/
  APic : Scheme.{0} → Type 1
  /-- テンソル積による群構造。 -/
  group : (X : Scheme.{0}) → CommGroup (APic X)
  /-- 下にある可逆層を忘れる写像。 -/
  forgetMetric : (X : Scheme.{0}) → APic X → toHermitianMetricData.toPicardData.Pic X
  /-- ★忘れる写像は群準同型。 -/
  forgetMetric_mul : ∀ (X : Scheme.{0}) (L M : APic X),
    forgetMetric X
      (@HMul.hMul _ _ _
        (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _
          (@instHMul _ (toHermitianMetricData.toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (forgetMetric X L) (forgetMetric X M)
  /-- ★★**層と計量の対から算術直線束を作る**(`mk` は構造体コンストラクタ名なので `ofMetric`)。 -/
  ofMetric : (X : Scheme.{0}) → (L : toHermitianMetricData.toPicardData.Pic X) →
    toHermitianMetricData.Metric X L → APic X
  /-- ★★★**対を忘れると元の層に戻る**。

  ★★★**これが退化を殺す。**`APic := PUnit` だと `forgetMetric` は定数だが、
  この式はすべての `L` について `forgetMetric (mk X L m) = L` を要求する。
  ★`Pic` は (B1) の `equivPicRing` で非自明が強制されているので、矛盾する。 -/
  forgetMetric_mk : ∀ (X : Scheme.{0}) (L : toHermitianMetricData.toPicardData.Pic X)
    (m : toHermitianMetricData.Metric X L), forgetMetric X (ofMetric X L m) = L
  /-- ★★算術直線束はすべて対から来る(`APic` に余計な元が無い)。 -/
  mk_surjective : ∀ (X : Scheme.{0}) (x : APic X),
    ∃ (L : toHermitianMetricData.toPicardData.Pic X)
      (m : toHermitianMetricData.Metric X L), ofMetric X L m = x
  /-- 射に沿った引き戻し。 -/
  pullback : {X Y : Scheme.{0}} → (X ⟶ Y) → APic Y → APic X
  /-- ★引き戻しは群準同型。 -/
  pullback_mul : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : APic Y),
    pullback f
      (@HMul.hMul _ _ _
        (@instHMul _ (group Y).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _
          (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (pullback f L) (pullback f M)

def APicData.waiting : WaitingFor :=
  { what := "(D1) 算術直線束の群 APic(X) —— 可逆層と ι_X-両立な計量の対、そのテンソル積と引き戻し"
    trackB := "Found/Arakelov — ★(B1) の可逆層と (C3) の計量に従属する。★★因子表示では既に対応物を構成済み(`Found/GenEll/ArchPoint.lean` の `ArithCartier` = 因子 + Green 関数、テンソル積 `ArithCartier.tensor`、引き戻し `pullbackADiv` とその積保存 `pullbackADiv_tensor_unconditional`。すべて sorry 0、2026-08-17)" }

/-! ## ★出典の紐付け(`.src`) -/

def APicData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 D——APic(X) の群構造のみ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
