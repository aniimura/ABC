/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Interface.Arakelov.ArcSpace
import ABC3.Interface.Arakelov.APic.Definition11

/-!
# APic —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

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
  /-- ★★**底変換で値が変わらない**(原文 p.4 の `deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`)。

  ★★★★2026-08-19 の修正——**射を代数構造に結んだ**。

  それまで `φ` は**任意の**射 `Spec 𝓞_K ⟶ Spec 𝓞_F` であり、
  同時に宣言されている `[Algebra F K]` は**一切使われていなかった**
  ——これは書き落としである。★原文 p.4 が言っているのは
  「`F ⊆ K` に沿った制限で次数が変わらない」であって、任意の射ではない。

  ★★任意の `φ` 版も**真ではある**(どの環準同型 `𝓞_F → 𝓞_K` も体の埋め込みに延び、
  各埋め込みがちょうど `[K:F]` 回現れる)が、
  **原文に無い強い主張**であり、証明には分数体への延長の機構が要る。
  ★★★下流(D3)が使うのは代数写像に沿った場合だけなので、そちらに揃える。 -/
  degF_baseChange : ∀ (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] [IsScalarTower (𝓞 F) (𝓞 K) K]
    (L : toAPicData.APic (Spec (CommRingCat.of (𝓞 F)))),
    degF K (toAPicData.pullback
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L) = degF F L
  /-- ★★★**計量を定数倍すると次数が動く**——アルキメデス側の寄与そのもの。

  原文 (GenEll p.4):
  > — where xF : Spec(OF ) → X is any morphism that gives rise to x.

  ★★★**これが「`deg_F ≡ 0`」の退化を殺す。**`deg_F` は計量に依存しなければならない。
  ★係数 `1/[F:ℚ]` は原文の**正規化**である(`Found/GenEll/ArithDiv.lean` の
  `degNormalized` と同じ規約)。 -/
  degF_scale : ∀ (F : Type) [Field F] [NumberField F]
    (L : toAPicData.toHermitianMetricData.toPicardData.Pic
      (Spec (CommRingCat.of (𝓞 F))))
    (m : toAPicData.toHermitianMetricData.Metric _ L) (c : ℝ),
    degF F (toAPicData.ofMetric _ L
        (toAPicData.toHermitianMetricData.scale _ L c m))
      = degF F (toAPicData.ofMetric _ L m) - c

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
  /-- `X(ℚ̄)` の点。

  ★★2026-08-19 の修正——宇宙を `Type 1` にした。`IsPointOf` が数体 `F : Type` に
  わたって量化するので、`AlgPoint` が `Type 0` だと代表元を保持できない
  ((D1) の `APic` と同じ誤記)。 -/
  AlgPoint : Scheme.{0} → Type 1
  /-- `x ∈ X(ℚ̄)` が射 `x_F : Spec 𝓞_F → X` で表されること。 -/
  IsPointOf : (X : Scheme.{0}) → (F : Type) → [Field F] → [NumberField F] →
    (Spec (CommRingCat.of (𝓞 F)) ⟶ X) → AlgPoint X → Prop
  /-- ★★★★★**どの射も点を定める**。

  原文 (GenEll p.5):
  > as the height function associated to the arithmetic line bundle M.

  ★★★★**2026-08-19 に見つけた穴。**それまで `IsPointOf` が**空でもよかった**ので、
  `AlgPoint := PUnit`、`IsPointOf := fun _ _ _ _ _ => False`、`height := 0` が
  すべての欄を満たしてしまった(`height_eq_degF` が空虚に成り立つ)。

  ★★★塞ぎ方: **`Spec 𝓞_F ⟶ X` はすべて `X(ℚ̄)` の点を定める**ことを課す。
  これで `height_eq_degF` が効き、高さは `deg_F` の引き戻しに**強制される**。 -/
  isPointOf_exists : ∀ (X : Scheme.{0}) (F : Type) [Field F] [NumberField F]
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X), ∃ x : AlgPoint X, IsPointOf X F xF x
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
  /-- ★★★★★2026-08-20 の追加(§9-405)——**自明束はこの類に入る**。

  ★これが無いと `IsEffective := fun _ _ => False` で
  `height_bddBelow` が**空虚に成立**してしまう。
  ★★逆に `IsEffective := True` は `height_bddBelow` 自身が落とす
  (任意の算術直線束の高さは下に有界でない)。
  ★★★この 2 本で `IsEffective` は上下から挟まれる。

  ★★★★**未塗りの穴**:原文の前提は「有効」ではなく
  「`L_Q` のある正のテンソル冪が大域切断で生成される」である。
  ★★★★★**「有効」では偽になる**——反例は膨らみ上げ:
  `L = O(E)`(`E` は例外因子)は有効だが、
  `L|_E = O_{P^1}(-1)` なので `E` 上で高さは `-∞` へ行く。
  ★§9-405 に記録。 -/
  isEffective_one : ∀ (X : Scheme.{0}),
    IsEffective X (@One.one _
      (toAPicSpecData.toAPicData.group X).toDivInvMonoid.toMonoid.toOne)
  height_bddBelow : ∀ (X : Scheme.{0}) (L : toAPicSpecData.toAPicData.APic X),
    IsEffective X L → ∃ C : ℝ, ∀ x : AlgPoint X, -C ≤ height X L x
  /-- ★★★**高さは `deg_F` の引き戻しである**——原文 p.4 の定義そのもの。

  原文 (GenEll p.5):
  > as the height function associated to the arithmetic line bundle M.

  ★★★**これが「`height ≡ 0`」の退化を殺す。**
  `deg_F` は (D2) の `degF_scale` で非自明が強制されているので、
  高さもそれに従って動かねばならない。 -/
  height_eq_degF : ∀ (X : Scheme.{0}) (L : toAPicSpecData.toAPicData.APic X)
    (F : Type) [Field F] [NumberField F]
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) (x : AlgPoint X),
    IsPointOf X F xF x →
    height X L x = toAPicSpecData.degF F (toAPicSpecData.toAPicData.pullback xF L)

def ArakelovHeightData.waiting : WaitingFor :=
  { what := "(D3) 高さ ht_L̄ : X(Q̄) → R と Proposition 1.4 —— X(Q̄) 全体で"
    trackB := "Found/Arakelov — ★★★**因子表示・U_X(Q̄) の範囲では構成済み**(`Found/GenEll/UPoint.lean` の `htU`、`Prop 1.4 (i)(ii)(iii)` は `HeightAdditive` / `HeightNonneg` / `HeightClass`。すべて sorry 0、2026-08-17)。★残るのは **X(Q̄) 全体への拡張**で、(B2)(Cartier 因子 → 可逆層)が埋まれば出る——可逆層はいつでも引き戻せるので『x が D を通らない』が消えるからである。★★`Prop 1.4, (iv)`(Northcott)は別筋で、数体上は `Found/GenEll/NorthcottRat.lean` に取得済み" }

def APicSpecData.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "Definition 1.2, (i)(層 D——Spec 𝓞_F 上の APic と deg_F)",
    sectionId := "genell-def-1-2-i" }

def ArakelovHeightData.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2, (i)(層 D——X(ℚ̄) 全体での高さ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
