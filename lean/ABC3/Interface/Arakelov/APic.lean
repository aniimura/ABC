import ABC3.Interface.Arakelov.ArcSpace

/-!
# Arakelov 理論のスケルトン(3/3)—— **合流: `APic(X)` と高さ**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3–p.6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★本ファイルが受ける 3 件

| # | 受けるもの | 現状(2026-08-17 実測) |
|---|---|---|
| D1 | `APic(X)`(算術直線束の群) | ★(B1)+(C3) に従属 |
| D2 | `APic(Spec 𝓞_F) ≅ ADiv(F)/APrc(F)` と `deg_F` | ★★`ADiv` / `deg_F` / `APrc` は**実装済**(`Found/GenEll/ArithDiv.lean`)。橋が無い |
| D3 | 高さ `ht_L̄ : X(ℚ̄) → ℝ` と `Proposition 1.4` | ★★★**因子表示・`U_X(ℚ̄)` の範囲では構成済**(2026-08-17)。全域化が残る |

★★★**D3 は「まったく無い」ではない。**
`Found/GenEll/UPoint.lean` の `htU` が原文の高さ関数そのものであり、
`Prop 1.4, (i)(ii)(iii)` も因子表示では取れている。
★残るのは **`X(ℚ̄)` 全体への拡張**であり、それは (B2) が埋まれば出る
——可逆層はいつでも引き戻せるので「`x` が `D` を通らない」が消えるからである。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★D1 —— `APic(X)` -/

/-- **(D1)** **算術直線束の群** `APic(X)` —— 可逆層と計量の対。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の `L̄ = (L, |−|_L)` そのものである。 -/
structure APicData where
  /-- 台となる `Pic`。 -/
  toPicardData : PicardData
  /-- 台となる計量。 -/
  toHermitianMetricData : HermitianMetricData
  /-- `APic(X)` の台。 -/
  APic : Scheme.{0} → Type
  /-- テンソル積による群構造。 -/
  group : (X : Scheme.{0}) → CommGroup (APic X)
  /-- 下にある可逆層を忘れる写像。 -/
  forgetMetric : (X : Scheme.{0}) → APic X → toPicardData.Pic X
  /-- ★忘れる写像は群準同型。 -/
  forgetMetric_mul : ∀ (X : Scheme.{0}) (L M : APic X),
    forgetMetric X
      (@HMul.hMul _ _ _
        (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _
          (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (forgetMetric X L) (forgetMetric X M)
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

/-! ## ★★D2 —— `Spec 𝓞_F` 上での `APic` と `deg_F` -/

/-- **(D2)** `APic(Spec 𝓞_F) ≅ ADiv(F)/APrc(F)` と、正規化次数 `deg_F`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★**原文自身がこの同型を使っている。**
`ADiv(F)` / `deg_F` / `APrc(F)` は `Found/GenEll/ArithDiv.lean` に
**実装済み(sorry 0)**なので、残るのは `APic` 側と橋だけである。 -/
structure APicSpecData where
  /-- 台となる `APic`。 -/
  toAPicData : APicData
  /-- 正規化次数 `deg_F : APic(Spec 𝓞_F) → ℝ`。 -/
  degF : (F : Type) → [Field F] → [NumberField F] →
    toAPicData.APic (Spec (CommRingCat.of (𝓞 F))) → ℝ
  /-- ★`deg_F` はテンソル積を和に移す。 -/
  degF_mul : ∀ (F : Type) [Field F] [NumberField F]
    (L M : toAPicData.APic (Spec (CommRingCat.of (𝓞 F)))),
    degF F (@HMul.hMul _ _ _
        (@instHMul _ (toAPicData.group _).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = degF F L + degF F M
  /-- ★★**底変換で値が変わらない**(原文 p.4 の `deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`)。 -/
  degF_baseChange : ∀ (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] (φ : Spec (CommRingCat.of (𝓞 K)) ⟶ Spec (CommRingCat.of (𝓞 F)))
    (L : toAPicData.APic (Spec (CommRingCat.of (𝓞 F)))),
    degF K (toAPicData.pullback φ L) = degF F L

def APicSpecData.waiting : WaitingFor :=
  { what := "(D2) APic(Spec 𝓞_F) ≅ ADiv(F)/APrc(F) と正規化次数 deg_F、およびその底変換不変性"
    trackB := "Found/Arakelov — ★★`ADiv(F)` / `deg_F` / `ord_v` / 主因子 / `APrc(F)` は `Found/GenEll/ArithDiv.lean` に**実装済み(sorry 0)**。★★★**底変換不変性も因子表示では証明済み**(`Found/GenEll/HeightInvariant.lean` の `htArith_baseChange_natural`、2026-08-17)。★残るのは (D1) の `APic` 側と、その間の橋である" }

/-! ## ★★★D3 —— 高さ -/

/-- **(D3)** 高さ `ht_L̄ : X(ℚ̄) → ℝ` と `Proposition 1.4`。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**因子表示・`U_X(ℚ̄)` の範囲では構成済みである**
(`Found/GenEll/UPoint.lean` の `htU`)。
★残るのは **`X(ℚ̄)` 全体への拡張**であり、それは (B2) が埋まれば出る。 -/
structure ArakelovHeightData where
  /-- 台となる `APic`。 -/
  toAPicSpecData : APicSpecData
  /-- `X(ℚ̄)` の点。 -/
  AlgPoint : Scheme.{0} → Type
  /-- ★★**高さ** `ht_L̄`。 -/
  height : (X : Scheme.{0}) → toAPicSpecData.toAPicData.APic X → AlgPoint X → ℝ
  /-- ★★**`Proposition 1.4, (i)`** —— 加法性(`X(ℚ̄)` **全体**で)。 -/
  height_mul : ∀ (X : Scheme.{0}) (L M : toAPicSpecData.toAPicData.APic X) (x : AlgPoint X),
    height X (@HMul.hMul _ _ _
        (@instHMul _ (toAPicSpecData.toAPicData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
        L M) x
      = height X L x + height X M x
  /-- ★★**`Proposition 1.4, (ii)`** —— 有効な算術直線束の高さは下に一様有界。 -/
  IsEffective : (X : Scheme.{0}) → toAPicSpecData.toAPicData.APic X → Prop
  height_bddBelow : ∀ (X : Scheme.{0}) (L : toAPicSpecData.toAPicData.APic X),
    IsEffective X L → ∃ C : ℝ, ∀ x : AlgPoint X, -C ≤ height X L x

def ArakelovHeightData.waiting : WaitingFor :=
  { what := "(D3) 高さ ht_L̄ : X(Q̄) → R と Proposition 1.4 —— X(Q̄) 全体で"
    trackB := "Found/Arakelov — ★★★**因子表示・U_X(Q̄) の範囲では構成済み**(`Found/GenEll/UPoint.lean` の `htU`、`Prop 1.4 (i)(ii)(iii)` は `HeightAdditive` / `HeightNonneg` / `HeightClass`。すべて sorry 0、2026-08-17)。★残るのは **X(Q̄) 全体への拡張**で、(B2)(Cartier 因子 → 可逆層)が埋まれば出る——可逆層はいつでも引き戻せるので『x が D を通らない』が消えるからである。★★`Prop 1.4, (iv)`(Northcott)は別筋で、数体上は `Found/GenEll/NorthcottRat.lean` に取得済み" }

/-! ## ★出典の紐付け(`.src`) -/

def APicData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 D——APic(X) の群構造のみ)",
    sectionId := "genell-def-1-1-i" }

def APicSpecData.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "Definition 1.2, (i)(層 D——Spec 𝓞_F 上の APic と deg_F)",
    sectionId := "genell-def-1-2-i" }

def ArakelovHeightData.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2, (i)(層 D——X(ℚ̄) 全体での高さ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
