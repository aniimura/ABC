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
