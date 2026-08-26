import ABC3.Meta.Claim
import ABC3.Interface.GenEll.HeightTheory
import ABC3.Found.GenEll.BDClass
import ABC3.Found.GenEll.ArithDiv
import ABC3.Found.GenEll.Conductor
import ABC3.Found.GenEll.Prop16
import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Found.GenEll.HeightMetric
import ABC3.Found.GenEll.HeightAdditive
import ABC3.Found.GenEll.HeightClass

/-!
# [GenEll] §1 Generalities on Heights —— 必要 9 件の statement(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.3–p.10。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★このファイルの位置づけ —— S0(statement を型で固定する)

`[IUTchIV]` が要求する [GenEll] 24 件のうち **§1 の 9 件**をここで固定する。
★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。

★★**ただし `sorry` を消すことを目的にしてはならない。**
内容を `Interface` の仮説へ移せば `sorry` は消えるが、それは `tools/check.mjs` 冒頭 B5 が
名指しする穴である。ゆえに `Interface/GenEll/HeightTheory.lean` は
**公理を 1 つも持たない**(データと述語だけ)。

## ★★原典の記法の食い違いを、写したうえで別に記録する

★**`≲` の向きが、定義と用法で食い違う。** 原文 p.5 の定義は逐語で

> α ≲_F β … if there exists a ["constant"] C ∈ ℝ such that β(x) − α(x) ≤ C

であるから `α ≲ β` は **`β ≤ α + C`** を意味する。
ところが `Proposition 1.6` の**表題は "Conductor Bounded by the Height"** で、
本文は `log-cond_D ≲ ht_D` と書く——表題どおりなら `log-cond ≤ ht + C` のはずだが、
定義どおりに読むと **`ht ≤ log-cond + C`** になり、**逆である**。

★★**本ファイルは印字どおりに写す**(`BDle` = 逐語の定義)。
★食い違いそのものは `Gap/GenEll/BDDirection.lean` に
`GapRecord` と**反例**(`Found/GenEll/BDClass.lean` の `bdle_ne_bdge`)つきで記録する。
**写し間違えないことと、食い違いを黙らないことは、別の仕事である。**

## ★実装済みのものは実装済みのものを使う

- `BDle` / `BDge` / `BDeq` / `BDClass` —— `Found/GenEll/BDClass.lean`(`sorry` 無し)
- `ADiv` / `deg` / `degNormalized` —— `Found/GenEll/ArithDiv.lean`(`sorry` 無し)

★これにより §1 の skeleton は**空虚な述語の羅列にならない**。
`Definition 1.2` と `Definition 1.5` は、下に書くとおり
**既に作ってあるものの上に本当に載る**。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll

/-! ## Definition 1.2 —— 高さ関数とその BD-class -/

/-- **[GenEll] Definition 1.2**。

原文 (GenEll p.5):
> The relation "≈" clearly defines an equivalence relation on the set of functions

(i) は高さ関数 `ht_M̄ : X(ℚ̄) → ℝ`、(ii) は BD-class。
★**両方をまとめて 1 つの式にした** —— 高さ関数の **BD-class** を取る操作である。
これが原文の (i)(ii) の合流点であり、以降 §1・§2 の不等式はすべてこの型の上で語られる。

★`BDClass` は `Found/GenEll/BDClass.lean` で**商として実際に構成してある**
(`bdSetoid` / `Quotient`)。posit ではない。 -/
noncomputable def htClass (D : HeightTheoryData) (M : D.ABundle) : BDClass D.Point :=
  BDClass.mk (D.ht M)

/-! ## Example 1.3 —— Galois-finite と compactly bounded -/

/-- **[GenEll] Example 1.3, (i)** の `X(ℚ̄)^{=d}`。

原文 (GenEll p.5):
> [cf. the discussion in Definition 1.5, (i), below of "minimal fields of definition"].

`X(ℚ̄)^{=d} ≝ X(ℚ̄)^{≤d} \ X(ℚ̄)^{≤d−1}`。 -/
def degEq (D : HeightTheoryData) (d : ℕ) : Set D.Point :=
  D.degLe d \ D.degLe (d - 1)

/-- **[GenEll] Example 1.3, (i)** の `Galois-finite`。

原文 (GenEll p.5):
> over the positive integers] is finite, then we shall say that E is Galois-finite.

★**これは posit ではなく定義である。** 原文は「各 `E^{≤d}` が有限」としか言っておらず、
それは `degLe` があれば**純粋に集合論的に**書ける。
★したがって `Example 1.3, (i)` は **Arakelov 理論も Galois 表現も要求しない**。 -/
def GaloisFinite (D : HeightTheoryData) (E : Set D.Point) : Prop :=
  ∀ d : ℕ, 0 < d → (E ∩ D.degLe d).Finite

/-- **[GenEll] Example 1.3** が導入する 2 つの語。

★**(i) は上で定義した**(集合論だけ)。
★**(ii) は `Interface` が posit している**——`compact domain` と `X^arc` を要求するからで、
mathlib に `complex analytic space` は 0 件(2026-08-16 実測)。

この対を 1 つの式にすることで、「この Example が何を導入したか」を型で固定する。 -/
def example_1_3 (D : HeightTheoryData) :
    (Set D.Point → Prop) × (Set D.Point → Prop) :=
  (GaloisFinite D, D.CompactlyBounded)

/-! ## Proposition 1.4 —— 高さの基本性質 -/

/-- **[GenEll] Proposition 1.4**(Basic Properties of Heights)。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

(i) `ht_{L̄⊗M̄}(x) = ht_L̄(x) + ht_M̄(x)` ——★**`≈` ではなく `=`**(目視確認)。
(ii) `L_ℚ` のある正冪が大域切断で生成されるなら `ht_L̄ ≳ 0` ——★**`≥` ではなく `≳`**。
(iii) `ht_L̄` の BD-class は `L_ℚ` の同型類だけに依る。
(iv) `L_ℚ` が ample なら `{x ∈ X(ℚ̄)^{≤d} : ht_L̄(x) ≤ C}` は**有限**(Northcott)。

★(ii) の `≳ 0` は逐語では `BDge (ht L) 0`、すなわち `ht(x) − 0 ≤ C`。
**定義どおりに読むと「上に有界」であって「下に有界」ではない**——
これも上の docstring で述べた向きの食い違いの一例である(`Gap/GenEll/BDDirection.lean`)。

## ★★★★★★★★ 2026-08-26——構成に置き換えて閉じた(第 371 ブロック)

★以前は `∀ D : HeightTheoryData` と量化していたが、
**`HeightTheoryData` は公理を 1 つも持たないデータ**なのでそれは**偽**であった
(反例: `Check/GenEll/HeightAxiomGap.lean` の `prop_1_4_statement_false`。
`Point ≝ ℕ`･`ABundle ≝ ℝ`･`tensor ≝ (+)`･`ht L x ≝ L²` で
**(i) 加法性と (iv) Northcott が同時に破れる**)。

★★公理を足すのは **B5 そのもの**である——`Proposition 1.4` の (i)-(iv) は
`HeightTheoryData` に足すべき公理**そのもの**だから、足せば仮定の言い換えになる。
★★★ゆえに**構成に置き換えた**。

### ★★★★ 4 つの内訳

| 条 | 中身 | どこにあるか |
|---|---|---|
| (i) 加法性 | `ht(D⊗E) = ht(D) + ht(E)` | `htArith_tensor_unconditional` |
| (ii) 下に有界 | `-C ≤ ht`(`C` は `F` にも点にも依らない) | `Prop16.lean` の `prop_1_4_ii` |
| (iii) `≈` | 因子が同じで計量が連続なら高さの差は一様に有界 | `HeightMetric.lean` の `htArith_sub_abs_le`(第 371) |
| (iv) Northcott | 射影モデルがあれば有限 | `NorthcottCoord.lean` の `northcott_of_projModel`(第 369-371) |

### ★★★★★逸脱と、残してあるもの(明示)

1. ★**量化する対象**: `∀ D : HeightTheoryData` → `ArcModel` + `ArithCartier`。前者では偽だから。
2. ★★**(iii) の形**: 原文は『生成ファイバーが同じなら』。因子表示では
   『因子が同じ(計量だけ違う)』である。**垂直因子の差の分は含めていない**。
3. ★★★**(iv) の形**: 原文は `X` の `ℤ`-固有性と `L_ℚ` の豊富性から**射影埋め込み**を得る。
   本 statement はその埋め込みを **`ArcModel` と同じ立場でデータとして受けている**。
   ★★★★**Northcott 性そのものは受けていない**——
   それは `Found/GenEll/NorthcottTuple.lean`･`NorthcottCoord.lean` で**証明している**。
   ★★★★★残っているのは「**`htArith` がその意味での射影モデルを持つ**」という幾何の段だけである。 -/
theorem prop_1_4 {X : AlgebraicGeometry.Scheme.{0}} {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ABC3.Found.GenEll.ArcModel X V)
    [Nonempty (ABC3.Found.GenEll.complexPoints X)]
    (Dv : ABC3.Found.GenEll.ArithCartier X)
    (hg : @Continuous _ _ M.topology _ Dv.green) :
    -- (i) 加法性
    (∀ (Ev : ABC3.Found.GenEll.ArithCartier X) (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        ABC3.Found.GenEll.pullbackIdeal F Dv.divisor xF ≠ 0 →
        ABC3.Found.GenEll.pullbackIdeal F Ev.divisor xF ≠ 0 →
        ABC3.Found.GenEll.htArith F (Dv.tensor Ev) xF
          = ABC3.Found.GenEll.htArith F Dv xF + ABC3.Found.GenEll.htArith F Ev xF)
    -- (ii) 下に一様に有界
  ∧ (∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        -C ≤ ABC3.Found.GenEll.htArith F Dv xF)
    -- (iii) 因子が同じで計量が連続なら高さは BD-同値
  ∧ (∀ Ev : ABC3.Found.GenEll.ArithCartier X, Ev.divisor = Dv.divisor →
        @Continuous _ _ M.topology _ Ev.green →
        ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
          (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
          |ABC3.Found.GenEll.htArith F Dv xF - ABC3.Found.GenEll.htArith F Ev xF| ≤ C)
    -- (iv) Northcott(射影モデルを与えられたものとして)
  ∧ (∀ (P I : Type) [Finite I] (d : ℕ) (ht : P → ℝ)
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p)),
        (∀ p, Module.finrank ℚ (fld p) ≤ d) →
        ∀ (crd : ∀ p, I → (fld p)) (idx : I) (const : ℝ),
        (∀ p, haveI := hnf p; Height.mulHeight (crd p) ≤ Real.exp (ht p + const)) →
        Function.Injective (fun (p : P) (i : I) => ((crd p i / crd p idx : fld p) : ℂ)) →
        ∀ C : ℝ, {p : P | ht p ≤ C}.Finite) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro Ev F _ _ xF hD hE
    exact ABC3.Found.GenEll.htArith_tensor_unconditional F Dv Ev xF hD hE
  · exact ABC3.Found.GenEll.prop_1_4_ii M Dv hg
  · intro Ev hdiv hcont
    exact ABC3.Found.GenEll.htArith_sub_abs_le M Dv Ev hdiv.symm hg hcont
  · intro P I _ d ht fld hnf hdeg crd idx const hcmp hinj C
    exact ABC3.Found.GenEll.northcott_of_projModel d ht fld hnf hdeg crd idx const hcmp hinj C

/-! ## Remark 1.4.1 —— 理論が `X_ℚ` だけに依ること -/

/-- **[GenEll] Remark 1.4.1**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

「`X(ℚ̄)` 上の高さ関数の **BD-class の理論**は、スキーム `X_ℚ` だけに依る」。

★**「だけに依る」を型で書くとこうなる** —— 2 つのモデル `D`, `D'` について、
点と生成ファイバーの同一視があり、そのもとで直線束の生成ファイバーが対応するなら、
高さ関数は BD-同値である。

★これが **`Theorem 2.1` が「数体上の曲線」から出発できる根拠**である——
ℤ-モデル `X` の取り方に依らないことを保証している。

## ★★★★★★★ 2026-08-26——構成に置き換えて閉じた(第 372 ブロック)

原文は『follows **immediately** from the definitions, together with
**Proposition 1.4, (iii)**』と書いている。
★すなわち中身は「**高さは類にしか依らない**」であり、
構成の側では `Found/GenEll/HeightClass.lean` が**すでにそれを持っていた**
——`htArith_bdeq_of_pullbackAPic_eq`。

★★構成の側では**定数差すら出ない**(`C = 0`)——
機構は `deg` が主算術因子の上で消えること(積公式)である。

### 逸脱(明示)

| 項 | 原典 | 形式化 |
|---|---|---|
| 量化する対象 | `∀ D D' : HeightTheoryData` | **同じ `X` 上の 2 つの `ArithCartier`** |
| 「だけに依る」の中身 | `X_ℚ` だけに依る | **引き戻した算術因子の類だけに依る** |

★★★**含めていないもの**: 異なる 2 つの `ℤ`-モデルの比較
——原文は『同型がある有限素数集合 `Σ` の上で `ℤ[Σ^{-1}]` へ延びる』という段を使うが、
その段はまだ持っていない。★**本 statement はその手前までである**。 -/
theorem remark_1_4_1 {X : AlgebraicGeometry.Scheme.{0}}
    (Dv Ev : ABC3.Found.GenEll.ArithCartier X)
    (F : Type) [Field F] [NumberField F]
    (h : ∀ xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X,
      ABC3.Found.GenEll.pullbackAPic F Dv xF = ABC3.Found.GenEll.pullbackAPic F Ev xF) :
    BDeq (fun xF => ABC3.Found.GenEll.htArith F Dv xF)
      (fun xF => ABC3.Found.GenEll.htArith F Ev xF) :=
  ABC3.Found.GenEll.htArith_bdeq_of_pullbackAPic_eq F Dv Ev h

/-! ## Definition 1.5 —— log-diff と log-cond -/

/-- **[GenEll] Definition 1.5**。

原文 (GenEll p.8):
> Note if x ∈ X(F ) ⊆ X(Q), where [F : Q] &lt; ∞, then by considering the scheme-theoretic image of the corresponding morphism Spec(F ) → X, one obtains a well-defined minimal field of definition Fmin ⊆ F of x.

(i) 最小定義体 `F_min`、(ii) `E` が **reduced** ⇔ `E = E_red`、
(iii) `log-diff_X(x) ≝ deg_F(δ_x)`、(iv) `log-cond_D(x) ≝ deg_F(f_x^D)`。

★★**(iii)(iv) の右辺は、既に作ってある `degNormalized` そのものである。**
原文の下線つき `deg` は正規化次数であり、`Found/GenEll/ArithDiv.lean` に
`sorry` 無しで実装済み。★ゆえにこの定義の**式そのものはもう持っている**——
残るのは `δ_x`(差積イデアルが定める有効算術因子)と
`f_x^D ≝ (D_x)_red`(引き戻した因子の被約化)を**構成すること**である。

★**その 2 つは Arakelov 理論ではない**——可換環論(`IsDedekindDomain.differentIdeal`)と
scheme 論(Cartier 因子の引き戻し)である。 -/
noncomputable def defn_1_5 {F : Type*} [Field F] [NumberField F] : ADiv F → ℝ :=
  degNormalized

/-! ★`(−)_red`(`ADivRed`)とその冪等性・有効性・**次数不等式**は
`Found/GenEll/Conductor.lean` へ移した——★`Proposition 1.6` の**非アルキメデス側**が
そこで実際に証明されるからである(`deg((a)_red) ≤ deg(a)`)。 -/

/-! ## Remark 1.5.1 —— log-cond の BD-class が `(X_ℚ, D_ℚ)` だけに依ること -/

/-- **[GenEll] Remark 1.5.1**。

原文 (GenEll p.8):
> Remark 1.5.1. In the spirit of Remark 1.4.1, we observe that the log-different

`log-diff_X` はスキーム `X_ℚ` だけに依る。`log-cond_D` は対 `(X, D)` に依り得るが、
**その BD-class は ℚ-スキームの対 `(X_ℚ, D_ℚ)` だけに依る**。

★理由は原文が書いている——別の対の同型は
**ある有限素数集合 `Σ` の上で** `ℤ[Σ^{-1}]` へ延びる。

## ★★★★★★ 2026-08-26——**前半はすでに `Found/` にある**(第 372 の実測)

原文の主張は 2 つに分かれる:

| 半分 | 状態 |
|---|---|
| `log-diff_X` は `X_ℚ` だけに依る | ★**すでにある**——`Found/GenEll/LogDiff.lean` の `logDiffOfField` は **`X` を引数に持たない**し、`LogDiffValue.lean` の `logDiffOfField_eq` が `log|disc F| / [F:ℚ]` と値を与える(実際には `X_ℚ` にさえ依らない) |
| `log-cond_D` の BD-class は `(X_ℚ, D_ℚ)` だけに依る | ★★**こちらが残っている**——spreading out(`ℤ[Σ^{-1}]` への延長)が要る |

★★★**前半だけを `Remark 1.5.1` として書き換えることはしない**——
後半は `Proposition 1.7` の証明が実際に使う(『Σ 上の log-cond の寄与は ≈ 0』)ので、
落とすと下流に影響が出る。★★★★**この `sorry` は spreading out 待ちである**。
★**「有限個の素数を除けば」という緩みが BD-class に吸収される**というのが
この論文が BD-class を使う理由そのものであり、
`Proposition 1.7` の証明でも「`Σ` の上の寄与は `≈ 0`」として同じ形で現れる。 -/
theorem remark_1_5_1 (D D' : HeightTheoryData)
    (ePt : D.Point ≃ D'.Point)
    (dv : D.Divisor) (dv' : D'.Divisor)
    (hcompl : ∀ x, x ∈ D.compl dv ↔ ePt x ∈ D'.compl dv') :
    BDeq (D.logDiff) (fun x => D'.logDiff (ePt x))
  ∧ BDeq (fun x : ↥(D.compl dv) => D.logCond dv x.1)
         (fun x : ↥(D.compl dv) => D'.logCond dv' (ePt x.1)) := by
  sorry

/-! ## Proposition 1.6 —— 導手は高さで抑えられる -/

/-- **[GenEll] Proposition 1.6**(Conductor Bounded by the Height)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

`L = O_X(D)` とし `U ≝ X\D`、`ht_D ≝ ht_L̄` とすると、`U(ℚ̄)` 上で `log-cond_D ≲ ht_D`。

★**枝番を持たない単一の主張**である(目視確認 2026-08-16)。

★★**ここが向きの食い違いが最も鮮明に出る場所である。**
表題は "Conductor **Bounded by** the Height"、すなわち `log-cond ≤ ht + C` を言っている。
ところが `≲` を **p.5 の定義どおり**に読むと `BDle log-cond ht` は
`ht(x) − log-cond(x) ≤ C`、つまり **`ht ≤ log-cond + C`** であり **逆になる**。
★本 statement は**印字どおり**(`BDle`)に写した。食い違いは
`Gap/GenEll/BDDirection.lean` に記録してある。

★証明が `X^arc` を要求する箇所も目視で特定した——アルキメデス素点の寄与を
「**コンパクト空間 `X^arc` 上の連続関数 `|s|_L` が有界**」で処理している。
非アルキメデス側は `Definition 1.5, (iv)` の `(−)_red` だけで済む。

## ★★★★★★★★ 2026-08-26——`sorry` が消えた(第 368 ブロック)

★以前は `∀ D : HeightTheoryData` と量化していたが、**それは偽である**
(`Check/GenEll/HeightAxiomGap.lean`)。★★公理を足すのは B5 そのものなので、
本ファイルの診断が書いていたとおり**構成に置き換えた**。

★★★中身は `Found/GenEll/Prop16.lean` にすでにあった——三段である:

1. `X^arc` はコンパクト(`ArcModel.compactSpace`)
2. 連続な Green 関数は下に有界(`ArcModel.exists_bound`)
3. 導手は高さ + 定数で押さえられる(`ArchBound.logCond_le_htArith_add`)

★★★★**定数 `C` が数体 `F` にも点にも依らない**のが要である(次数で正規化してあるから)。

### 逸脱 2 件

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ D : HeightTheoryData` | **`ArcModel` + `ArithCartier`** | 前者では偽だから |
| `≲` の向き | 印字どおりの `BDle` | **`log-cond ≤ ht + C`** | 表題 "Conductor **Bounded by** the Height" の向き。`Gap/GenEll/BDDirection.lean` に記録済み |

★原文の `U(ℚ̄)` への制限は `pullbackIdeal F Dv.divisor xF ≠ 0` が担っている。 -/
theorem prop_1_6 {X : AlgebraicGeometry.Scheme.{0}} {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ABC3.Found.GenEll.ArcModel X V)
    [Nonempty (ABC3.Found.GenEll.complexPoints X)]
    (Dv : ABC3.Found.GenEll.ArithCartier X)
    (hg : @Continuous _ _ M.topology _ Dv.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F]
        (xF : ABC3.Found.GenEll.specRingOfIntegers F ⟶ X),
        ABC3.Found.GenEll.pullbackIdeal F Dv.divisor xF ≠ 0 →
        ABC3.Found.GenEll.logCond F Dv.divisor xF
          ≤ ABC3.Found.GenEll.htArith F Dv xF + C :=
  ABC3.Found.GenEll.prop_1_6 M Dv hg

/-! ## Proposition 1.7 —— 導手と log-different -/

/-- **[GenEll] Proposition 1.7**(Conductors and Log Differents)。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

(i) `U_Y(ℚ̄)` 上で
`log-cond_E − log-cond_D ≲ log-diff_Y − log-diff_Z ≲ (1 − 1/e)·log-cond_E`。
(ii) 分岐指数が各点で `e` に**等しい**なら Riemann–Hurwitz の関係式。

★★**この論文で `≲` と `≤` の差が主張になる場所である。**
原文 p.10 の証明は本文中で明示的に区別している(目視確認 2026-08-16):
- prime-to-`Σ` 部分の不等式は「**`=` と `≤` であって `≲` ではない**」
- `Σ` の上の `log-diff_Y − log-diff_Z` は「**`≥ 0` であって `≳` ではない**」

★`pdftotext` は `≲` を**出力に何も残さない**ので、この区別は
**`.txt` からは原理的に復元できない**。ゆえに目視必須であり、そう写した。

★★**条件 (a)–(d) は落としていない** —— `CoveringSetup.hyp` として仮定に置いてある。
落とせば主張が強くなり、**偽の skeleton** になるからである。
展開できていないことは `.needs` に `.implicitStep` として明記した。

★(ii) は Riemann–Hurwitz であり、`deg` は `Y_ℚ`・`Z_ℚ` 上の直線束の次数
——`HeightTheoryData` の語彙の外なので、本 statement では (i) だけを固定する。
(ii) は `prop_1_7_ii_pending` で「まだ書けない」ことを型で明示する。 -/
theorem prop_1_7 (S : CoveringSetup) (h : S.hyp) :
    BDle (fun x : ↥(S.DY.compl S.divY) =>
            S.DZ.logCond S.divZ (S.toPoint x.1) - S.DY.logCond S.divY x.1)
         (fun x : ↥(S.DY.compl S.divY) =>
            S.DY.logDiff x.1 - S.DZ.logDiff (S.toPoint x.1))
  ∧ BDle (fun x : ↥(S.DY.compl S.divY) =>
            S.DY.logDiff x.1 - S.DZ.logDiff (S.toPoint x.1))
         (fun x : ↥(S.DY.compl S.divY) =>
            (1 - 1 / (S.e : ℝ)) * S.DZ.logCond S.divZ (S.toPoint x.1)) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def htClass.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2",
    sectionId := "genell-def-1-2" }

def degEq.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

def GaloisFinite.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3, (i)",
    sectionId := "genell-ex-1-3" }

def example_1_3.src : Source :=
  { paper := "GenEll", pdfPage := 5, item := "Example 1.3",
    sectionId := "genell-ex-1-3" }

def prop_1_4.src : Source :=
  { paper := "GenEll", pdfPage := 6, item := "Proposition 1.4",
    sectionId := "genell-prop-1-4" }

/-- ★原文 p.6–p.7 の証明を通読して数えた。 -/
def prop_1_4.needs : List ProofObligation :=
  [ .implicitStep
      "(i) は『follows immediately from the definitions』の 1 行で済まされている。ht の定義(deg_F の和)まで降りれば加法性だが、その deg_F が Interface 待ちである" 6,
    .implicitStep
      "(ii) の証明はアルキメデス素点の寄与を『X^arc がコンパクトだから |s|_L は有界』で処理する。★複素解析空間そのものが mathlib に無い(complex analytic space 0 件、2026-08-16 実測)" 6,
    .implicitStep
      "(iii) は (i)(ii) から従うと原文は書くが、そこで『ℤ 上の算術直線束の高さは有界』という補題を暗黙に使っている(ht_{L̄⊗M̄} ≈ ht_L̄ の式)" 6,
    .folklore "(iv) は Northcott の有限性定理。原文は証明を与えず『well-known』の扱いをしている" 7,
    .implicitStep
      "★statement の語彙(APic(X)・ht・ample・大域切断生成)を Interface/GenEll/HeightTheory.lean に posit した。**我々は作っていない**" 6 ]

def remark_1_4_1.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Remark 1.4.1",
    sectionId := "genell-rem-1-4-1" }

def remark_1_4_1.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Proposition 1.4, (iii)(BD-class が L_ℚ の同型類だけに依る)" 6,
    .implicitStep
      "原文は『follows immediately from the definitions』とだけ書く。★「X_ℚ だけに依る」を型にするには『2 つの ℤ-モデルの生成ファイバーが同型なら高さが BD-同値』を言う必要があり、その同一視の与え方は原文に無い(我々が ePt / eGen として補った)" 8 ]

def defn_1_5.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5",
    sectionId := "genell-def-1-5" }

def remark_1_5_1.src : Source :=
  { paper := "GenEll", pdfPage := 8, item := "Remark 1.5.1",
    sectionId := "genell-rem-1-5-1" }

def remark_1_5_1.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Remark 1.4.1(BD-class の理論が X_ℚ だけに依る)" 8,
    .implicitStep
      "原文の証明は『ある有限素数集合 Σ の上で同型が ℤ[Σ^{-1}] へ延びる』という 1 文である。★その延長の存在(スキームの spreading out)は与えられていない" 9,
    .citation "[GenEll]" "ℤ[Σ^{-1}] ≝ ℤ[{p^{-1}}_{p∈Σ}] への spreading out"
      (.absent "mathlib に scheme の spreading out / 有限型スキームの ℤ 上のモデルの理論は無い(2026-08-16、Mathlib/AlgebraicGeometry 配下の全宣言名を確認)") 9 ]

def prop_1_6.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.6",
    sectionId := "genell-prop-1-6" }

/-- ★原文 p.9 の証明は 5 行しかない。その 5 行が要求するものを数えた。 -/
def prop_1_6.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Definition 1.5, (iv)(導手 f_x^D ≝ (D_x)_red と log-cond_D)" 8,
    .otherPaper "[GenEll]" "Proposition 1.4, (iii)(ht_D ≝ ht_L̄ が well-defined であること)" 6,
    .implicitStep
      "★アルキメデス側は『コンパクト空間 X^arc 上の連続関数 |s|_L が有界』で片づけられている。X^arc(複素解析空間)が mathlib に 0 件なので、この 1 行が層まるごとに相当する" 9,
    .implicitStep
      "★非アルキメデス側は『(−)_red の定義から従う』とだけ書く。実際には『D_x の係数 ≥ 1 の素点でだけ (D_x)_red が 1 を持つ』という不等式であり、算術因子の上で書ける(ADivRed を参照)" 9,
    .implicitStep
      "★★表題 'Conductor Bounded by the Height' と、p.5 の ≲ の定義が示す向きが逆である。Gap/GenEll/BDDirection.lean に記録した。**本 statement は印字どおりに写してある**" 9 ]

def prop_1_7.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.7",
    sectionId := "genell-prop-1-7" }

/-- ★原文 p.10 の証明を通読して数えた。

★**この証明の核は初等的である**——「素数 `p` と正整数 `d` を固定すると、
`[L:K] ≤ d` なる `ℚ_p` の有限拡大の有限 Galois 拡大 `L/K` すべてについて、
差積イデアルが `p^n·O_L` を含むような正整数 `n` が存在する」。
Arakelov も Galois 表現も要らず、**局所体の分岐理論と Kummer 理論だけ**である。 -/
def prop_1_7.needs : List ProofObligation :=
  [ .implicitStep
      "★★条件 (a)(b)(c)(d)(reduced / D_ℚ = φ_ℚ^{-1}(E_ℚ)_red / 有限エタール / 分岐指数が e を割る)を CoveringSetup.hyp という **1 つの不透明な Prop** として持っている。**落としてはいない**が、展開もしていない" 9,
    .folklore "原文が『the elementary theory of differents』と呼ぶもの。prime-to-Σ 部分の等式・不等式はここから出る" 10,
    .implicitStep
      "★原文は prime-to-Σ 部分について『with \"=\" and \"≤\", not \"≲\"!』、Σ 上について『with \"≥\", not \"≳\"!』と**明示的に区別している**。pdftotext は ≲ を出力に残さないので、この区別は .txt からは復元できない(2026-08-16 実測)" 10,
    .citation "[GenEll]" "局所体の分岐理論と Kummer 理論(証明の核となる initial claim)"
      (.inMathlib "IsDedekindDomain.differentIdeal / Polynomial.IsSplittingField / IsCyclic — ただし『[L:K] ≤ d なる全ての L/K に一様な n』という**一様性**の部分は mathlib に無い(2026-08-16 実測)") 10,
    .otherPaper "[GenEll]" "Remark 1.5.1(Σ 上の log-cond の寄与が ≈ 0 であること)" 8,
    .implicitStep
      "(ii) の Riemann–Hurwitz は Y_ℚ・Z_ℚ 上の直線束の次数 deg(−) を要求する。HeightTheoryData の語彙の外なので、本ファイルでは (i) だけを固定した" 10 ]

end ABC3.Skeleton.GenEll
